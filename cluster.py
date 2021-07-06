import re
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

'''
Parse arguments.
Input file is required.
'''
def parseArgs():
    parser = argparse.ArgumentParser(description='Cluster cell types using mcmicro marker expression data.')
    parser.add_argument('-i', '--input', help="Input CSV of mcmicro marker expression data for cells", type=str, required=True)
    parser.add_argument('-o', '--output', help='The directory to which output files will be saved', type=str, required=False)
    parser.add_argument('-m', '--markers', help='A text file with a marker on each line to specify which markers to use for clustering', type=str, required=False)
    parser.add_argument('-v', '--verbose', help='Flag to print out progress of script', action="store_true", required=False)
    parser.add_argument('-k', '--neighbors', help='the number of nearest neighbors to use when clustering. The default is 30.', default=30, type=int, required=False)
    args = parser.parse_args()
    return args


'''
Get input data file name
'''
def getDataName(path):
    fileName = path.split('/')[-1]
    return fileName


'''
Get markers to use for clustering from a text file where each marker is on a line and corresponds exactly to the column name in the input data file.
Returns a list of markers to use for clustering.
'''
def get_markers(markers_file):
    markers = [] # list of markers in file

    # read markers from file
    f = open(markers_file, 'r')
    for line in f:
        markers.append(line.strip())

    return markers


'''
Clean data in input file.
NOTE: Currently we are doing this with pandas however, using csv might be faster, or numpy.

Exclude the following data from clustering:
    - X_centroid, …, Extent, Orientation - morphological features
    - Any of the following DNA stains
        - DNA0, DNA1, DNA2, …
        - Hoechst0, Hoechst1, ....
        - DAPI0, DAPI1, …
    - AF* (e.g., AF488, AF555, etc.) - autofluorescence
    - A* (e.g., A488, A555, etc.) - secondary antibody staining only

To include any of these markers in the clustering, provide their exact names in a file passed in with the '-m' flag
'''
def clean(input_file):

    # constants

    # a default list of features to exclude from clustering
    FEATURES_TO_REMOVE = ['X_centroid', 'Y_centroid', # morphological features
                        'column_centroid', 'row_centroid', 
                        'Area', 'MajorAxisLength', 
                        'MinorAxisLength', 'Eccentricity', 
                        'Solidity', 'Extent', 'Orientation', 
                        'DNA.*', 'Hoechst.*', 'DAP.*', # DNA stain
                        'AF.*', # autofluorescence
                        'A\d{3}.*'] # secondary antibody staining only (iy has to have 3 digist after)

    if args.verbose:
        print('Cleaning data...')

    # load csv
    data = pd.read_csv(input_file)

    # if markers provided, keep only those features and the Cell IDs. It is important that the CellID column is first.
    if args.markers:
        if CELL_ID not in markers: # add cell ID to list of columns to keep
            markers.insert(0, CELL_ID)
        elif markers.index(CELL_ID) != 0: # if cell ID column is included but not first, move it to the front
            markers.insert(0, markers.pop(markers.index(CELL_ID)))
        data = data[markers]
    else:
        # find any columns in the input csv that should be excluded from clustering be default
        # NOTE: may want to replace this with regex, it might be faster.
        col_to_remove = []
        cols = data.columns
        for feature in FEATURES_TO_REMOVE:
            r = re.compile(feature)
            col_to_remove.extend(list(filter(r.match, cols)))
        
        # drop all columns that should be excluded
        data = data.drop(columns=col_to_remove, axis=1)

    # save cleaned data to csv
    data.to_csv(f'{output}/{CLEAN_DATA_FILE}', index=False)

    if args.verbose:
        print(f'Done. Cleaned data is in {output}/clean_data.csv.')


'''
Write CELLS_FILE from leidenCluster() adata
'''
def writeCells(adata):
    cells = pd.DataFrame(adata.obs[CELL_ID].astype(int)) # extract cell IDs to dataframe
    cells[CLUSTER] = adata.obs[LEIDEN] # extract and add cluster assignments to cells dataframe
    cells.to_csv(f'{output}/{cells_file}', index=False)


'''
Write CLUSTERS_FILE from leidenCluster() adata
'''
def writeClusters(adata):
    clusters = pd.DataFrame(columns=adata.var_names, index=adata.obs[LEIDEN].cat.categories)   
    clusters.index.name = CLUSTER # name indices as cluster column
    for cluster in adata.obs.leiden.cat.categories: # this assumes that LEIDEN = 'leiden' if the name is changed, replace it for 'leiden' in this line
        clusters.loc[cluster] = adata[adata.obs[LEIDEN].isin([cluster]),:].X.mean(0)
    clusters.to_csv(f'{output}/{clusters_file}')


'''
Cluster data using the Leiden algorithm via scanpy
'''
def leidenCluster():

    sc.settings.verbosity = 3 # print out information
    adata_init = sc.read(CLEAN_DATA_FILE, cache=False) # load in clean data

    # move CellID info into .obs
    # this assumes that 'CELL_ID' is the first column in the csv
    adata_init.obs[CELL_ID] = adata_init.X[:,0]
    adata = ad.AnnData(np.delete(adata_init.X, 0, 1), obs=adata_init.obs, var=adata_init.var.drop([CELL_ID]))

    # compute neighbors and cluster
    sc.pp.neighbors(adata, n_neighbors=args.neighbors, n_pcs=0) # compute neighbors, using the number of neighbors provided in the command line. Default is 30.
    sc.tl.leiden(adata, key_added = LEIDEN) # run leidan clustering. default resolution in 1.0

    # write cell/cluster information to 'CELLS_FILE'
    writeCells(adata)

    # write cluster mean feature expression to 'CLUSTERS_FILE'
    writeClusters(adata)


'''
Main.
'''
if __name__ == '__main__':
    args = parseArgs() # parse arguments

    # get user-defined output dir (strip last slash if present) or set to current
    if args.output is None:
        output = '.'
    elif args.output[-1] == '/':
        output = args.output[:-1]
    else:
        output = args.output

    # get list of markers if provided
    if args.markers is not None:
        markers = get_markers(args.markers)

    # constants
    CLEAN_DATA_FILE = 'clean_data.csv' # name of output cleaned data CSV file
    CELL_ID = 'CellID' # column name holding cell IDs
    CLUSTER = 'Cluster' # column name holding cluster number
    LEIDEN = 'leiden' # obs name for cluster assignment
    
    # output file names
    data_prefix = getDataName(args.input) # get the name of the input data file to add as a prefix to the output file names
    clusters_file = f'{data_prefix}_clusters.csv' # name of output CSV file that contains the mean expression of each feaute, for each cluster
    cells_file = f'{data_prefix}_cells.csv' # name of output CSV file that contains each cell ID and it's cluster assignation
    
    # clean input data file
    clean(args.input)

    # cluster using scanpy implementation of Leiden algorithm
    leidenCluster()
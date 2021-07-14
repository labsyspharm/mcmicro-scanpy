# mcmicro-scanpy
An MCMICRO module implementation of scanpy for clustering cells using the Leiden algorithm.

## Paramenter Reference
```
usage: cluster.py [-h] -i INPUT [-o OUTPUT] [-m MARKERS] [-k NEIGHBORS] [-c]

Cluster cell types using mcmicro marker expression data.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input CSV of mcmicro marker expression data for cells
  -o OUTPUT, --output OUTPUT
                        The directory to which output files will be saved
  -m MARKERS, --markers MARKERS
                        A text file with a marker on each line to specify which markers to use for clustering
  -k NEIGHBORS, --neighbors NEIGHBORS
                        the number of nearest neighbors to use when clustering. The default is 20.
  -c, --method          Include a column with the method name in the output files.
```
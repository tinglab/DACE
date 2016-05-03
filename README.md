# DACE
DACE: A Scalable DP-means Algorithm for Clustering Extremely Large Sequence Data

### INSTALLATION
##### Requirements
 cmake, MPI, Boost, g++

##### Compile
```
cd DACE
cmake . && make
```
Please find the executable file in DACE/bin/dace.


### USAGE

```
mpiexec -np [Host number] bin/dace [Options]             

   Note: Host number must be greater than 1: one for master, the others for slaves.

Options
   -i  input file in fasta format [required]
   -o  output prefix [required]
   -p  clustering threshold, default 0.97 
   -c  thread number, default 2
   -b  block size of DP-means, default 2000
   -l  maximum LSH iteration number, default 100
   -x  convergence threshold for LSH iteration, default 0.01
   -d  iteration number of DP-mean algorithm, default 5
   -s  suffix array filter threshold, default 3000
   -q  0/1, sort by length/weight in Big-Kmer Mapping, default 1
   -h  print this help

Example:
   $ mpiexec -np 2 bin/dace -i db.fa -o db_output -c 8
```
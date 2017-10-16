Readme for parallel io

gauss_col.c is the updated version for collective communication
gauss_io.c is the updated version for parallel io 


The preprocessing file could be compiled by

mpicc preprocessing.c -o preprocessing

And could be run by using

'mpiexec -n 1 ./preprocessing ./ge_data/sizeYYxYY'

This step will generate the zzz.mat and zzz.vec files with xxx.pre name.

The gaussian elimination could be run in the way we have been doing until now.
GE will generate the file in the binary format so to convert the file into
ASCII use the postprocessing file. It could be compiled by

mpicc postprocessing.c -o postprocessing

And the said file could be run by

'mpiexec -n 1 ./postprocessing ./ge_data/xxx.sol'  


The preprocessing and post processing steps have been written in a way that
they could use only one process at most.


We also provide a check.c file, which implement a check on binary solution
and original solution files. One can easily add the automatic check in 
the bash script and automatically check the solution for parallel io.



compbio_algos
=============

Various comp bio algorithms I've implemented.

Python 

Dead End Elimination is in dee.py.  It can use two matrices, pair.dat and pair_small.dat. 
The file pair.dat is a plain text, tab-delimited
file, containing the energies as a matrix of the following form:

E1    E12    E13  ...  E1N
0     E2     E23  ...  E2N
0     0      E3   ...  E3N
...  ...     ...  ...
0     0       0   ...  EN

The file pair_small.dat is a smaller version of the same file in the same format.


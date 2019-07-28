# SparsePR: Matlab Software for Sparse Phase Retrieval

This repository contains Matlab code for solving the 
sparse phase retrieval problem. Details of the method, 
theoretical guarantees and representative numerical 
results can be found in 

Robust Sparse Phase Retrieval Made Easy   
Mark Iwen, Aditya Viswanthan and Yang Wang    
arXiv    
2014 

This software was developed at the [Department of 
Mathematics][msumath], [Michigan State University][msu] 
and is released under the MIT license. 

The software was developed and tested using Matlab 
R2014a and uses [TFOCS][tfocs], a Matlab software 
package for the efficient construction and solution of 
convex optimization problems. A copy of the TFOCS 
package is included under the third party software 
directory at third/. A selection of scripts also uses 
the CVX optimization software, which can be downloaded 
[here][cvx].


## Directory Structure and Contents

The software package is organized under the following 
directory structure:

 - demos/    
   This folder contains a representative implementation 
   of the algorithm discussed in the paper. It also 
   contains implementations of related algorithms, such 
   as Basis Pursuit Compressive Sensing reconstruction, 
   PhaseLift and Compressive Phase Retrieval via Lifting
   (CPRL). If you have just downloaded the software, 
   execute sparsePR.m from this folder in Matlab.

 - src/
   This folder contains auxiliary functions necessary to 
   implement the algorithm. Examples include generating 
   test functions and measurements. 

 - tests/
   This folder contains convergence tests; for example, 
   scripts to generate noise vs error plots, runtime 
   plots and so on. These typically take a long time to 
   run and are provided for those looking to recreate 
   the results in the paper.

 - third/
   This software contains third party software used by 
   SparsePR. In particular, it contains an installation 
   of TFOCS with minor modifications (to the 
   tfocs_initialize.m routine). 

 - contrib/
   This folder will contain contributed code and 
   additional examples, modifications and applications.


## Instructions

Extract the contents of the zip file and execute, in 
Matlab, scripts from the demos/ and tests/ folder. 
Scripts which use CVX have _cvx.m in their file name. 
These require a working installation of CVX.


## Contact

Bug reports, comments and suggestions are welcome 
at the [SparsePR Bitbucket repository][bitbucket].


[msu]: http://www.msu.edu/
[msumath]: http://math.msu.edu/
[tfocs]: http://cvxr.com/tfocs/
[cvx]: http://cvxr.com/cvx/
[mark]: http://www.math.msu.edu/~markiwen/
[yang]: http://www.math.msu.edu/~ywang/
[bitbucket]: https://bitbucket.org/charms/sparsepr/

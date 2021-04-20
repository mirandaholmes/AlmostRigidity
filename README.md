# AlmostRigidity

This is Matlab code to perform the Almost Rigidity test introduced in 

 *Holmes‐Cerfon, M., Theran, L., & Gortler, S. J. (2019). Almost‐Rigidity of Frameworks. Communications on Pure and Applied Mathematics.*

The main script that performs the test is **driver_almostrigidity.m**. This script requires CVX, which you can download at [http://cvxr.com/cvx/](http://cvxr.com/cvx/).
    
The script runs two examples. You can choose which one by commenting out one of the following lines: 

    framework = load_n10;                % n=10 hypostatic framework
    framework = load_exampleTensegrity;  % example of a 2d tensegrity
    
To run it on your own example, create a file similar to **load_.m**, which loads the following data into a 'framework' struct:   
* framework.x:         vector of framework vertices (n x dim - dimensional vector)
* framework.n:         # of vertices
* framework.dim:       dimension of ambient space (2 or 3)  
* framework.a:         adjacency matrix (n x n)  
* framework.lengths:   distance matrix, same size as a. Set to NAN if lengths are determined automatically from given lengths  
* framework.ap:        periodic adjacency matrix (empty if not periodic)  
* framework.lattice:   lattice vectors for periodic frameworks (empty if not periodic)  
* framework.pfix:      vector of particle indices to fix ([i,j] if dim=2 or [i,j,k] if dim=3). Empty to fix COM & rotation constraints instead.  
* framework.types:     Bar types for a tensegrity (0=bar, 1=strut, -1=cable.)  
                       Order in types, must match order of bonds in `[rr,cc] = find(triu(a))`.  
                       Empty if a bar framework.   
                       
Other things you might want to change in driver_almostrigidity.m are:   
* tols = maximum singular value for constructing almost-flex-space and almost-stress space  
* plotting options -- see options in code.   

I am always interested in learning about applications of this code. 
If you use it, please don't hesitate to write with comments / questions / bugs / etc (mch328@nyu.edu).
If you use this code to help with a publication, please cite the paper it is based on: 

  *Holmes‐Cerfon, M., Theran, L., & Gortler, S. J. (2019). Almost‐Rigidity of Frameworks. Communications on Pure and Applied Mathematics.*
    

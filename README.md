# Sparse_matrices
Algorithms for solving equations involving sparse matrices 

In this project, different algorithms were developed to store sparse matrices, to compute their determinant and permanent and to use their sparsity in efficiently solving sparse system of linear algebraic equations using different methods such as:       
	1. Gauss-Seidel method                                                                      
	2. Cramer’s rule                                                                                    
	3. LU Decomposition method (Direct Solution and Iterative Refinement Method)                                            
Throughout the project the Compressed Row Storage (CRS) scheme was used to store the sparse matrices.

# The Algorithms in this project are:
1.Storing the given sparse matrix by Compressed Row Storage                                                                 
2.Efficient computation of the determinant and permanent of square sparse matrix                                          
3.Solving 'n x n' sparse linear system of equations using Gauss-Seidel method                                           
4.Solving 'n' variable sparse linear system of equations using Cramer’s rule                                        
5.1.Solving 'n x n' sparse linear system of equations by LU Decomposition (Direct Method)                                              
5.2.Solving 'n x n' sparse linear system of equations by LUIR (LU Decomposition Iterative Refinement)                        

# Following citations by Professor R.C. Mittal were used in developing these algorithms:
Application of the Cramer rule in the solution of sparse systems of linear algebraic equations:
http://www.sciencedirect.com/science/article/pii/S0377042700005720

LU-decomposition and numerical structure for solving large sparse nonsymmetric linear systems:
http://www.sciencedirect.com/science/article/pii/S0898122101002796?np=y&npKey=2791864f79d1d9af933ed32ea9414080e9458b6d966ec9d5a3621be04519e20e

Efficient computation of the permanent of a sparse matrix:
http://www.tandfonline.com/doi/abs/10.1080/00207160108805061?needAccess=true#aHR0cDovL3d3dy50YW5kZm9ubGluZS5jb20vZG9pL3BkZi8xMC4xMDgwLzAwMjA3MTYwMTA4ODA1MDYxP25lZWRBY2Nlc3M9dHJ1ZUBAQDA=

Efficient solution of a sparse non-symmetric system of linear equations:
http://www.tandfonline.com/doi/abs/10.1080/00207160210940#aHR0cDovL3d3dy50YW5kZm9ubGluZS5jb20vZG9pL3BkZi8xMC4xMDgwLzAwMjA3MTYwMjEwOTQwP25lZWRBY2Nlc3M9dHJ1ZUBAQDA=

# All the programs in this project are in Python 3
The Algorithms in the project are naive and may have a good scope for improvement in terms of efficiency.
Any suggestions and improvements in the algorithms are encouraged

# For Python users:
There are some excellent modules like "scipy.sparse" available in Python for the most efficient computations involving sparse matrices

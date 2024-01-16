# Computing Eigenvalues and Eigenvectors using Normalized Power Iteration together with Deflation
### Assignment: 
Implementing the normalized power iteration algorithm together with deflation to find two eigenvalues and corresponding eigenvector for the dominant eigenvalue of a given matrix A, where A is a square matrix. The program must read A from an input file and output the dominant eigenvalue, its corresponding eigenvector and the next eigenvalue as a text file.


### Compile the code

> g++ -std=c++11 assignment2.cpp -o run


### Run the program

**Usage:** ./run path_to_matrix_A   tolerance   output_file

**Example:** ./run A.txt 1e-6 result.txt

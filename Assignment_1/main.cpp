#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>

const double epsilon = 1e-10;

//read matrix from txt file
template<class T>
T** readMatrix(const char* file_name, int& n){ // n will be size of matrix
    T** retval = NULL; // retval is defined to terminate function when file cannot open
    std::fstream matrix_file(file_name);
    if (!matrix_file.is_open()) {
        std::cout << "Error opening file.\n";
        return retval;
    }
    n = 0;
    std::string row;

    //taking lines of matrix one by one and increasing n until line is empty. So we can find size of matrix
    while (std::getline(matrix_file, row)) {
        if (!row.empty()) n++;
    }
    
    //allocating memory for matrix
    T** matrix = new T*[n]; 
    for(int i = 0;i < n; i++){
        matrix[i] = new T[n];
    }

    //this part is to get back 'getline' function to the beginning of the txt file
    matrix_file.clear();
    matrix_file.seekg(0, std::ios::beg);

    //reading txt file and assigning elements to the allocated matrix
    int i = 0;
    while (std::getline(matrix_file, row)) {
        std::stringstream ss(row);
        for ( int j = 0; j < n; j++){   
            ss >> matrix[i][j];
        }
        i++;
    }

    matrix_file.close();
    return matrix;
}

//read vector from txt file
template<class T>
T* readVector(const char* file_name, int n){

    T* retval = NULL; // retval is defined to terminate function when file cannot open
    std::ifstream vector_file(file_name);
    if (!vector_file.is_open()) {
        std::cout << "Error opening file.\n";
        return retval;
    }
    int i=0;
    std::string row;
    T* b = new T[n];
    double value;   
    while (vector_file >> value && i < n) // Read each value from the file and assign it to the array
    {
        b[i] = value;
        i++;
    }
    vector_file.close();
    return b;
}

//print matrix for debugging issues
template<class T>
void printMatrix(T** matrix, int n){ 

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            std::cout << matrix[i][j] << " ";
        }
        std::cout << std::endl;
        
    }
    
}

// frees memory allocated by matrix
template<class T> 
void freeMatrix(T** matrix, int n){ 
    for (int i = 0; i < n; i++)
    {
        delete[] matrix[i];
    }
    delete[] matrix;
}

//calculates determinant of 2 dimensional matrix
template<class T>
double findDeterminant2D(T** submatrix){ 
    
    return submatrix[0][0] * submatrix[1][1] - submatrix[0][1] * submatrix[1][0]; 
}

//this function calculates determinant of unknown size matrix by
//calling function recursively
template<class T>
double findDeterminant(T** matrix,int n){ 
    double determinant= 0;
    
    if(n==1) return matrix[0][0];
    
    else if(n==2) return findDeterminant2D(matrix);
    
    else{
    for(int i=0; i<n ; i++){
        T** submatrix = new T*[n-1];
        for(int m=0; m<n-1; m++){
            //submatrix does not have ith and jth column and rows of upper matrix
            submatrix[m] = new T[n-1]; 
        }
            for (int j = 1; j < n; j++) {
                    for (int k = 0; k < i; k++) {
                        submatrix[j-1][k] = matrix[j][k];
                    }
                    for (int k = i+1; k < n; k++) {
                        submatrix[j-1][k-1] = matrix[j][k];
                    }
                }
        
        double c = findDeterminant(submatrix,n-1); //c is cofactor that multiplies submatrix
        // if else is to determine sign of cofactor * 
        if(i%2==0) determinant += c * matrix[0][i];
        else determinant -= c * matrix[0][i];  
        freeMatrix(submatrix,n-1);
    }
    return determinant;
    }
}

//Implements Gaussian Elimination with Partial Pivoting
template<class T>
T** GaussElimination(T** matrix, T* vector, int n){
    
    for(int i= 0; i<n; i++){
        int maxindex= i;
        for(int j=i+1; j<n; j++){
            if(fabs(matrix[j][i]) > fabs(matrix[maxindex][i])) maxindex = j; //check which pivot is bigger
       }
        if(maxindex != i){
            //swap rows by higher pivot coming up
            std::swap(matrix[i],matrix[maxindex]);
            std::swap(vector[i],vector[maxindex]);
        } 
            

      
        for (int j = i+1; j < n; j++) // this is for calculating factor to multiply rows below pivot
        {   

            double f = matrix[j][i] / matrix[i][i]; // multiplier to make zero elements under pivot
            vector[j] -= f * vector[i];
            for (int k = i; k < n; k++) // for iterating over row below pivot
            {
                matrix[j][k] -= f * matrix[i][k];
                
            }
        }

    }   
        
    return matrix;

}

//Implement Backward Substitution over Upper Triangle Matrix
template<class T>
T* BackwardSubstution(T** matrix, T* vector, int n){
    T* x = new T[n];
    x[n-1] = vector[n-1] / matrix[n-1][n-1];  //bottom element of x is determined.
    for (int i = n-2; i >= 0; i--) 
    {
    // known x's are summed and substracted from corresponding b value then divided by coefficient of unknown x value
        double sum = 0;
        for(int j = n-1; j>i; j--){
            sum += x[j] * matrix[i][j];
            
        }
        x[i] = (vector[i] - sum) / matrix[i][i];
      
    }
    return x;
}

//find condition number for 2x2 matrix
template<class T>
void findConditionNumber(T** matrix){
    
    //allocate memory for inverse matrix
    T** inverse = new T*[2];
    inverse[0] = new T[2];
    inverse[1] = new T[2];
    
    for(int i=0; i < 2; i++){
        for(int j=0; j<2; j++){
            inverse[i][j] = matrix[i][j];
        }      
    }

    //find determinant of matrix
    double determinant = findDeterminant2D<T>(matrix);
    
    //find inverse matrix
    std::swap(inverse[0][0],inverse[1][1]);
    inverse[0][0] /= determinant;
    inverse[0][1] /= (-1 * determinant);
    inverse[1][0] /= (-1 * determinant);
    inverse[1][1] /= determinant;

    //calculate norm-1 of both matrix and its inverse
    double norm1 = fmax(fabs(matrix[0][0]) + fabs(matrix[1][0]) , fabs(matrix[0][1])+fabs(matrix[1][1])); // calculate the norm-1 A
    double norm1_inverse = fmax(fabs(inverse[0][0])+fabs(inverse[1][0]), fabs(inverse[0][1])+fabs(inverse[1][1])); // calculate the norm-1 of inverse A
    
    //calculate norm-infinity of both matrix and its inverse
    double norminf = fmax(fabs(matrix[0][0])+fabs(matrix[0][1]), fabs(matrix[1][0])+fabs(matrix[1][1]));
    double norminf_inverse = fmax(fabs(inverse[0][0])+fabs(inverse[0][1]), fabs(inverse[1][0])+fabs(inverse[1][1]));
    
    double condition_1 = norm1 * norm1_inverse; // calculate the condition number of 1
    double condition_inf = norminf * norminf_inverse; // calculate the condition number of infinity
    std::cout << "Condition number of 1: " << condition_1 << std::endl;
    std::cout << "Condition number of infinity: " << condition_inf << std::endl;

}


int main(int argc, char* argv[]){
    int n;
    if(argc != 3){
        std::cout << "ERROR: " << argc << " arguments entered. 3 needed"<<std::endl;
        return 1;
    }

    double** matrix = readMatrix<double>(argv[1],n); //read the A matrix from txt file
    
    if(n==2){    //check if the dimension of the matrix is 2x2
    findConditionNumber<double>(matrix); //if dimension is 2x2, find condition numbers
    }


    double determinant = findDeterminant(matrix,n); // calculate determinant of matrix

    std::cout << "Determinant " << determinant << std::endl;
    if(fabs(determinant) < epsilon){ // check if the matrix is singular
        std::cout << "The Matrix is singular. \n Exiting...";
        return 1;
    }
    
    double * b = readVector<double>(argv[2],n); // read the b vector from txt file
    
    double ** eliminated = GaussElimination<double>(matrix,b,n); // implement gaussian elimination with partial pivoting
    double * x = BackwardSubstution<double>(eliminated,b,n); // Find x vector using backward substitution
    
    std::cout << "\nX values:" << std::endl;
    for(int i = 0; i<n; i++) std::cout << x[i] << std::endl; 

    // free allocated memories
    freeMatrix<double>(eliminated,n);
    delete[] x;
    delete[] b;

 }


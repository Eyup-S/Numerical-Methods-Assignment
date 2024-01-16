// Class: Vector
// Description: Represents a mathematical vector and its operations.
//
// Member Functions:
// - random(): Generates random values for each element in the vector.
// - print(): Prints the elements of the vector.
// - getValue(int index) const: Returns the value at the given index.
// - setValue(double value, int index): Sets the value at the given index.
// - findLength(): Returns the length (magnitude) of the vector.
// - infinityNorm(): Returns the infinity norm of the vector.
// - setVector(Vector vec): Sets the current vector equal to another vector.
// - operator/(double length): Returns a new vector, which is the current vector divided by the given length.
// - vecCross(Vector) : Returns cross product of two vectors
//
// Class: Matrix
// Description: Represents a mathematical square matrix and its operations.
//
// Member Functions:
// - returnSize(): Returns the size of the square matrix.
// - print(): Prints the elements of the matrix.
// - initializeEigenvector(): Initializes the eigenvector with random values.
// - printEigenvectors(): Prints the eigenvectors of the matrix.
// - operator*(Vector& vec): Multiplies the matrix by a given vector.
// - operator-(Matrix): Implement substraction between two matrices.
// - getEigenvalue(): Returns the maximum eigenvalue of the matrix.
// - setMaxEigenvalue(double eigen): Sets the maximum eigenvalue of the matrix.
// - getEigenvector(): Returns the eigenvector of the matrix.
// - setEigenvector(double value, int index): Sets the value of the eigenvector at the given index.
// - setValue(double value, int index1, int index2): Sets the value of the matrix at the given indices.
// - getValue(int index1, int index2): Returns the value of the matrix at the given indices.
// - cross(Vector& vec): Returns the cross product of the matrix and a given vector.
//
// Function: powerIteration(Matrix& matrix, double tolerance)
// Description: Computes the dominant eigenvalue and eigenvector of the given matrix using the power iteration method.
//
// Function: Deflation(Matrix& matrix)
// Description: Applies the deflation method to the given matrix and returns the deflated matrix.

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <math.h>
class Matrix;

class Vector{
private:
    int n;
    double * vector;
public:
    //constructor with size argument
    Vector(int a){
        n = a;
        vector = new double[n];
    }
    // ~Vector(){
    //     delete[] vector;
    // }
    //generate random vector using rand function from stdlib
    void random(){
        for(int i = 0; i<n; i++){
            //random numbers are normalized so that they are between [0,1]
            vector[i] = ((double)rand()) / RAND_MAX;
        }
        
    }
    //get size of vector
    int size(){
        return n;
    }

    //print the vector
    void print(){

        for(int i = 0; i<n;i++){
            std::cout << vector[i] << std::endl;
        }
    }
    //get the value of element of vector at 'index'
    double getValue(int index) const{
        return vector[index];
    }
    //set the value of element of vector at 'index'
    void setValue(double value, int index){
        vector[index] = value;
    }
    //find infinity norm of the vector
    double infinityNorm(){
        double max = 0; //max element is held here
        for(int i = 0; i<n;i++){
            if(fabs(vector[i]) > max ){
                max = fabs(vector[i]); // if any element is bigge than max, max is updated
            }
        }
        return max;
    }

    //equate two vectors
    void setVector(Vector vec){
        for(int i=0; i<n; i++){
            vector[i] = vec.getValue(i);
        }

    }
    //overloading / operator to divide each element of the vector with given double value
    Vector operator/(double length){
        Vector normalized(n);
        for(int i = 0; i<n; i++){
            normalized.setValue(vector[i] / length,i);
        }
        return normalized;
    }
    

    Matrix vecCross(Vector);

};

class Matrix{
private:
    int n = 0; // assumed the matrix will be square always
    double ** matrix; //double array to hold valuess
    Vector * eigenvector; 
    double maxEigenvalue; // dominant eigenvalue
public:
    //constructor for reading matrix from file
    Matrix(const char* path){
        //read the matrix from the txt file
        std::fstream file(path);
        if(!file.is_open()){
            std::cerr << "file could not open" << std::endl;
            return;
        }
        //determine size of matrix
        std::string row;
        while (std::getline(file, row)){
            if (!row.empty()) n++;
        }
        
        //get back to the beginning of txt file
        file.clear();
        file.seekg(0, std::ios::beg);
        
        //allocate memory for matrix
        matrix = new double*[n]; 
        for(int i = 0 ; i < n ; i++){
            matrix[i] = new double[n];
        }
        int i = 0;

        //assign values read from the file to the matrix
        while (std::getline(file, row)) {
        std::stringstream ss(row);
            for ( int j = 0; j < n; j++){   
                ss >> matrix[i][j];
            }
            i++;
        
        }
        
        //close txt file
        file.close();
        //create eigenvector Vector with the size of equal to matrix
        eigenvector =  new Vector(n);

    
    }
    //constructor with size argument
    Matrix(int size){
        n = size;
        matrix = new double*[size]; 
        for(int i = 0 ; i < size ; i++){
            matrix[i] = new double[size];
        }
        for(int i=0;i<size;i++){
            for(int j=0;j<size;j++){
                matrix[i][j] = 0;
            }
        }
        eigenvector =  new Vector(size);

    }
    
    //return size of matrix
    int returnSize(){
        return n;
    }

    //print the matrix
    void print(){
        for(int i = 0; i<n;i++){
            for(int j = 0; j<n; j++){

                std::cout << matrix[i][j] << " ";
            }
        std::cout << std::endl;
        }
    }
    //initialze eigenvector with random values
    void initializeEigenvector(){
        eigenvector->random();
    }

    //print Eigenvectors
    void printEigenvectors(){
        eigenvector->print();
    }

    //operator overloading for multiplying a matrix with a vector
    Vector operator*(Vector& vec ){ 
        Vector result(n);
        for(int i = 0; i<n; i++){
            double sum = 0;
            for(int j = 0; j < n; j++){
                sum += matrix[i][j] * vec.getValue(j);
                
            }
            result.setValue(sum,i);
        }
        result.print();
        return result; 
    }

    //get Eigenvalue of matrix
    double getEigenvalue(){
        return maxEigenvalue;
    }

    //set dominant eigenvalue of matrix
    void setMaxEigenvalue(double eigen){
        maxEigenvalue = eigen;       
    }

    //get eigenvector of matrix
    Vector* getEigenvector(){
        return eigenvector;
    }

    //set the value of Eigenvector of the matrix at given 'index'
    void setEigenvector(double value, int index){
            eigenvector->setValue(value,index);
    }

    //set value of matrix at given 'indexes'
    void setValue(double value, int index1,int index2){
        matrix[index1][index2] = value;
    }

    //get the value of matrix at given 'indexes'
    double getValue(int index1,int index2){
        return matrix[index1][index2];
    }

    //Implement cross product between a matrix and a vector and return a vector
    Vector cross(Vector& vec){
        Vector result(n);
        for(int i = 0; i<n; i++){
            double sum = 0;
            for(int j = 0; j < n; j++){
                //sum product of matrix row and vector column
                sum += matrix[i][j] * vec.getValue(j);
                
            }
            result.setValue(sum,i);
        }
        return result; 
    }

    Matrix operator * (double scalar){
        Matrix result(n);
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                result.setValue(scalar * matrix[i][j],i,j); 
            }
        }
        return result;

    }
    Matrix operator -(Matrix mat){
        Matrix result(n);
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                result.setValue(matrix[i][j] - mat.getValue(i,j),i,j); 
            }
        }
        return result;
        
    }    

};

Matrix Vector::vecCross(Vector vec){
        Matrix mat(vec.size());
        for(int i =0;i<n;i++){
            for(int j=0;j<n;j++){
                mat.setValue(vector[i] * vec.getValue(j) ,i,j);
            }
        }
        return mat;
    }
//Implement normalized power iteration method to find dominant eigenvalue and its eigenvector
void powerIteration(Matrix& matrix,double tolerance){
    //initialize random eigenvector of the matrix
    matrix.initializeEigenvector();
    //Create vectors x_k and y_k
    Vector x_k = *matrix.getEigenvector();
    Vector y_k(matrix.returnSize());    
    double difference = 1000; // 1000 is given arbitrarily
    double oldEigen= 0; //hold previous eigenvalue
    double newEigen; //hold new computed eigenvalue
    //power iteration is implemented 
    while(difference > tolerance){ //check convergence
        double oldInfinity = y_k.infinityNorm();  // hold old infinity norm
        y_k.setVector(matrix.cross(x_k)); //cross product A * x_k
        x_k = y_k / y_k.infinityNorm();
        newEigen = y_k.infinityNorm();
        difference = fabs(newEigen -  oldEigen); //find difference for convergence
        oldEigen = newEigen;
    }
    //set the dominant eigenvalue found above
    matrix.setMaxEigenvalue(newEigen);
    //set the Eigenvector found above
    for(int i = 0 ; i<matrix.returnSize(); i++){
        matrix.setEigenvector(x_k.getValue(i),i);
    }
    std::cout << "Eigenvalue: " << newEigen << std::endl;

} 

// Implement Deflation method to find next eigenvalue and eigenvector
Matrix Deflation(Matrix& matrix,double tolerance) {
	Vector first_row(matrix.returnSize());//form vector for first row of given matrix
    //assign values to first_row vector
    for(int i = 0;i<matrix.returnSize(); i++){
        first_row.setValue(matrix.getValue(0,i),i );
    }
    //take eigenvector of the matrix
    Vector *eigenvec = matrix.getEigenvector();
    double Vi= 0; 
    //check if the element of the eigenvector is zero
    for(int i = 0;i < matrix.returnSize();i++){
        if(eigenvec->getValue(i) !=0){
            Vi= eigenvec->getValue(i);
            break;
        } 
    }
    //find x vector
    first_row = first_row / matrix.getEigenvalue() / Vi;
    //find the matrix that is result of v * x
    Matrix mat = eigenvec->vecCross(first_row);
    
    mat = matrix - mat *  matrix.getEigenvalue();
    Matrix B(matrix.returnSize() - 1);
    //take B matrix from the transformed matrix
    for(int i =0; i<matrix.returnSize() - 1; i++){
        for(int j =0; j<matrix.returnSize() - 1;j++){
            B.setValue(mat.getValue(i+1,j+1),i,j);
        }
    }
    return B;

}


int main(int argc, char* argv[]){
    if(argc !=4){
        std::cout << "ERROR!\nUsage: ./run <path_to_matrix_A> <tolerance> <output_file>" << std::endl;
        return -1;
    }
    double tolerance = std::strtod(argv[2],nullptr); //char to double conversion
    Matrix matrix(argv[1]); //create matrix object
    std::ofstream output(argv[3]); //open output file
    if(!output.is_open()){
        std::cout <<"ERROR: Unable to open the output file" << std::endl;
        return -1;
    }
    if(matrix.returnSize() == 1){ //check the size of given matrix
        output << matrix.getValue(0,0);
        output.close();
        return 0;
    }
    else{
    powerIteration(matrix,tolerance); //implement power iteration method for created matrix

    //write to txt file the first eigenvalue and its corresponding eigenvector
    output << "Eigenvalue#1: " <<matrix.getEigenvalue() << std::endl;
    for(int i =0; i<matrix.returnSize(); i++){
        output << matrix.getEigenvector()->getValue(i) << std::endl;
    }
    matrix.getEigenvector()->print();
    Matrix B = Deflation(matrix,tolerance); //create Matrix object B for deflation method
    powerIteration(B, tolerance); //find dominant eigenvalue of deflated matrix
    output << "Eigenvalue#2: " << B.getEigenvalue() << std::endl; //write the second eigenvalue to txt file
    output.close();
    }
    
}
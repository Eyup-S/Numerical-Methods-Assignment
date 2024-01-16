#include <iostream>
#include <cstdlib>
#include <math.h>

//sign function returns -1 if the given value is negative and vice versa
int sign(double value){
    if(value < 0) {
        return -1;
    }
    else
    {
        return 1;
    }
}
//defining given polynomial
double f(double x,int n, char* coef[]){
    double y = 0;
    for(int i = 0; i<n;i++){
        y += std::strtod(coef[i+1],nullptr) * pow(x,n-i-1); //multiplies input value with corresponding order
    }
    return y;
}

//implement bisection method
double bisection(double interval1,double interval2,double tolerance,int n, char* coeff[],int &count){
    double mid = 0; //define mean value
    while(interval2-interval1 > tolerance){ //check the convergence
        mid = interval1 + (interval2-interval1)/2; //calculate mid value
        if (sign(f(interval1,n,coeff)) == sign(f(mid,n,coeff)) ){ //check the signs
            interval1 = mid; //if signs are same slide towards interval2
        }
        else {
            interval2 = mid; //vice versa
        }
        count++;
    }
    return mid;

}

//Implement secant method
double secant(double x_0,double x_1,double tolerance,int n, char* coeff[],int &count){
    double h = 1000; //1000 is arbitrarily chosen
    double x_km1 = x_0; //x_k-1
    double x_k = x_1; //x_k value
    while(fabs(h) > tolerance){
        h = f(x_k,n,coeff) * ( (x_k - x_km1) / (f(x_k,n,coeff) - f(x_km1,n,coeff) ) ); //calculating difference
        x_km1 = x_k;
        x_k = x_k - h; // calculating x_k+1 value
        count++;
    }
    return x_k;

}
//Implement bisection method for first two iterations and then implement secant method
double hybrid(double interval1,double interval2,double tolerance,int n, char* coeff[],int &count){
    double mid = 0; //define mean value
    while(count <2){ // check count for only two iterations
    mid = interval1 + (interval2-interval1)/2; //calculate mid value
    if (sign(f(interval1,n,coeff)) == sign(f(mid,n,coeff)) ){ //check the signs
            interval1 = mid; //if signs are same slide towards interval2
        }
        else {
            interval2 = mid; //vice versa
        }
        count++;
    }
    double h = 1000; // 1000 is chosen arbitrarily
    double x_km1 = interval1; //x_k-1
    double x_k = interval2; //x_k
    while(fabs(h) > tolerance){
        h = f(x_k,n,coeff) * ( (x_k - x_km1) / (f(x_k,n,coeff) - f(x_km1,n,coeff) ) ); //calculate difference
        x_km1 = x_k; //assigning next x_k-1 value 
        x_k = x_k - h; // calculate x_k+1 value
        count++;
    }
    return x_k;

}

int main(int argc, char* argv[]){
    if(argc < 5){
        std::cout << "ERROR: Not enough parameters!" << std::endl;
        std::cout << "Usage: ./run <coefficients of polynomial with space> <interval1> <interval2> <tolerance>" << std::endl;
        return 1;
    }

    // Extracting input values from command-line arguments
    double x_0 = std::strtod(argv[argc - 3],nullptr);
    double x_1 = std::strtod(argv[argc - 2],nullptr);
    double tolerance = std::strtod(argv[argc - 1],nullptr);
    int order = argc -4; //maximum order of polynomial

    int bisection_count = 0;
    int secant_count = 0;
    int hybrid_count = 0;

    // Bisection method
    std::cout << "---------Bisection--------" << std::endl;
    double bisection_value = bisection(x_0, x_1, tolerance, order, argv,bisection_count);
    std::cout << "value: " << bisection_value << std::endl;
    std::cout << "count: " << bisection_count << std::endl;

    //Secant method
    std::cout << "----------Secant----------" << std::endl;
    double secant_value = secant(x_0, x_1, tolerance, order, argv,secant_count);
    std::cout << "value: " << secant_value << std::endl;
    std::cout << "count: " << secant_count << std::endl;

    //Hybrid Method
    std::cout << "----------Hybrid----------" << std::endl;
    double hybrid_value = hybrid(x_0, x_1, tolerance, order, argv,hybrid_count);
    std::cout << "value: " << hybrid_value << std::endl;
    std::cout << "hybrid count: " << hybrid_count << std::endl;

}
# Secant and Bisection Methods

### Assignment:
Implementing secant and bisection algorithms in order to solve f(x)=0 for any given polynomial f.
The program should take the coefficients of the function, initial guesses and the tolerance value as command line arguments and return the resulting values of x as well as the numbers of iterations for each method.
You should implement both methods separately first. Then you should use a hybrid method where you start with bisection method for the first two iterations and then continue with the secant method for the rest of the iterations. Your program should print out the number of iterations required for each of the 3 methods

### Compile the code

> g++ -std=c++11 assignment3.cpp -o run


### Run the program

**Usage:** ./run coefficients_of_polynomial_with_space   interval1   interval2   tolerance

**Example:** ./run 2 2 -7 1 -7 1.5 1.8 0.001

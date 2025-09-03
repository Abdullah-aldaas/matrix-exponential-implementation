# Matrix Exponential Implementation – Bachelor Thesis Codes (MATLAB)
This repository contains the MATLAB codes used in my Bachelor thesis (Computing the Matrix Exponential via Scaling and Squaring with Subdiagonal Padé Approximation: Theory, Stability, and Numerical Experiments) to implement sexpm/sexpmv and generate test matrices with specific properties.

These codes were used in the numerical experiments of my bachelor thesis to test the stability and efficiency of the sexpm/sexpmv method.

## How to run?
```matlab
% Example for computing e^A:
A = rand(3,3);
sexpm(A)
```
```matlab
% Example for computing e^Av:
A = rand(3,3);
v = rand(3,1);
sexpmv(A,v)
```
```matlab
% Example for generate test matrices:
seed = 0;
n = 6;
A = generate_random(seed,n)
```

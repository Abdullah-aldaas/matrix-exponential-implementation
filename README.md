# Random Matrix Generator – Bachelor Thesis Code (MATLAB)
This repository contains the MATLAB code used in my Bachelor thesis (Computing the Matrix Exponential via Scaling and Squaring with Subdiagonal Padé Approximation: Theory, Stability, and Numerical Experiments) to generate test matrices with specific properties.

## Matrix Properties:
- **Essentially nonnegative**
- **Diagonal dominant**
- All eigenvalues **λ(A) ≤ 0**
- One eigenvalue is **close to 0**
- Matrix norm approximately **‖A‖ ≈ 10^(seed/5)**

These matrices were used in my thesis to test the stability and efficiency of the sexpm/sexpmv method.

## How to run?
```matlab
%Example:
seed = 0;
n = 6;
A = generate_mat(seed,n);

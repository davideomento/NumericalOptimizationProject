# NumericalOptimizationProject
Project for the exam Numerical Optimization for Large Scale Problems in collaboration with Geard Koci.

# Numerical Optimization Project

This MATLAB project explores several variants of Newton's methods for solving unconstrained optimization problems. In particular, it implements and compares:

- Modified Newton method with backtracking  
- Truncated Newton method using CG and PCG  
- Empirical analysis on well-known test functions from the literature

## 📁 Contents

- `mod_newton_backtrack.m`: Modified Newton method with regularized Hessian and backtracking  
- `tru_newton_backtrack_cg.m`: Truncated Newton method using Conjugate Gradient (CG)  
- `tru_newton_backtrack_pcg.m`: Truncated Newton method using Preconditioned Conjugate Gradient (PCG)  
- Implementations and tests on:  
  - Rosenbrock function  
  - Chained Rosenbrock (n=1000)  
  - Problems No. 76 and No. 82 from the reference book  



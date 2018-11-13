# Comp_math

It is a set of programs for Computational Mathematics on seven topics (see below). A program's number matches a topic's number.

For each program use `Python 3.x`. One should install `matplotlib` and `numpy` to use all programs.

The description of program's result you can find in the first their strings or below. 

1. Errors. The program `1_taylor_exp.py` calculates Taylor's series for exponent.

2. Linear systems. The program `2_linear.py` compares Fast Gradient Method and Method of Minimal Residual for multiple random linear system. Their number is 1000.

3. The method of least squares. The program `3_OLS.py` finds the best degree of aproximal polinomial for exponent and uses it for noisy exponent.

4. Non-linear systems. The program `4_nonlinear.py` gives three results:
  *  Show one disadvantage of Newton's Method: it does not converge to a solution if starting point is choosed badly.
  *  Show how many iterations needs for to find all roots of a random polynomial of constant degree
  *  Show how dichotomy works for different degrees of polinomial with same number of iteration

5. Interpolation (*one should install also* `sympy`). The program `5_interpolation.py` creates interpolate Newton's polinomial and cubic spline for three functions(exponent, polinomial, sinus).

6. Numerical Integration. The program `6_integration.py` compares three method of computional integration (trapezium, rectangle and Simpson) for three functions.

7. Differential equation (*it has not been written yet*)

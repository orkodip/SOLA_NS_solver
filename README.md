# SOLA_NS_solver
This repository contain the code (written in C++) to compute the fluid flow inside a square lid driven cavity, demonstrating the performance of the incompressible Navier-Stokes solver: SOLA. Finite Difference Method is used to discretize the governing equations in a staggered grid. For more details, please refer to the works of Hirt et al. (1975) and Brandt et al. (1980).

# Reference
Hirt, C.W., Nichols, B.D. and Romero, N.C., 1975. SOLA: A Numerical Solution Algorithm for Transient Fluid Flows; LA-5852, Los Alamos Scientific Lab. N. Mex.: Los Alamos, NM, USA.

Brandt, A., Dendy Jr, J.E. and Ruppel, H., 1980. The multigrid method for semi-implicit hydrodynamics codes. Journal of Computational Physics, 34(3), pp.348-370.

# Software requirements
This solver needs:

- gcc

# How to install the required packages (on a Linux system)

To install gcc (compiler for C++ code)

```bash
sudo apt install build-essential
```

# How to compile and run the code

To compile the code

```bash
g++ SOLA.cpp -o output
```
To run this code

```bash
./output
```

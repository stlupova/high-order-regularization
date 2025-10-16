# high-order-regularization
Code for https://arxiv.org/abs/2510.13639

For questions/comments contact Svetlana Tlupova (svetlana.tlupova@farmingdale.edu) 

Reference:
J. T. Beale and S. Tlupova. High order regularization of nearly singular surface integrals. arXiv; Cornell University Library, 2025. https://arxiv.org/abs/

## Testing instructions

### Test: Stokes flow around a translating body (extending integral values to the whole grid) (Fig. 7 and 8 in referenced paper):

1.	Dependencies:
    * The FFTW library (https://www.fftw.org/)
    * The OpenMP library
3.	Code subdirectory: `highreg_stokes_extend_to_grid`
4.	In main.cpp, set the following:
    * Line 72: grid size (h=3/N): N = 3* (32, 64, 128, 256) (choose one)
    * Lines 76-78 (optional): order of regularization and kappa_0 for delta=kappa * h^q
5.  Compile using `make`. The makefile is set up for a Mac OS.
6.  Run using `./main.out`

### Test: Surfaces close to each other in Stokes flow (Fig. 9 and 10 in referenced paper):

1.	Dependencies: the OpenMP library
2.	Code subdirectory: `highreg_stokes_2surf`
3.	In main.cpp, set the following:
    * Line 68: grid size (h=1/N): N = 32, 64, 128, 256 (choose one)
    * Lines 147-148 (optional): order of regularization and kappa_0 for delta=kappa * h^q
5.  In KITC.h, set the following treecode parameters:
    * degree of interpolating polynominal P (recommended: P=6 for h=1/32 and h=1/64, P=8 for h=1/128, P=10 for h=1/256, etc., add 2 for h/2)
    * leaf size N0 (recommended: N0=1000 for h=1/32, N0=2000 for h=1/64, N0=4000 for h=1/128, N0=8000 for h=1/256, etc., double for h/2)
    * MAC parameter theta in sq_theta=theta^2 (recommended: theta=0.6)
5.  Compile using `make`. The makefile is set up for a Mac OS.
6.  Run using `./main.out`
   

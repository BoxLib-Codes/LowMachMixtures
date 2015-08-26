Low Mach Number Fluctuating Hydrodynamics of Multispecies Liquid Mixtures
Authors: Andy J. Nonaka, CCSE, Lawrence Berkeley Labs and
         Aleksandar Donev, Courant Institute, NYU
Last updated: August 21st 2015

These codes use the BoxLib framework to develop a low-Mach number code for simulating diffusively-mixing multispecies liquid and gas mixtures in the presence of gravity. The algorithm is specifically tailored to low Reynolds number flows, including steady Stokes flow, and does not use a projection or splitting algorithm.

Details can be found in the paper (see doc/LowMachMultispecies.pdf):

1. "Low Mach Number Fluctuating Hydrodynamics of Multispecies Liquid Mixtures", A. Donev and A. J. Nonaka and A. K. Bhattacharjee and A. L. Garcia and J. B. Bell, Physics of Fluids, 27(3):037103, 2015 [ArXiv:1412.6503].

The algorithmic details are given in the paper (see doc/LowMachImplicit.pdf):

2. "Low Mach Number Fluctuating Hydrodynamics of Binary Liquid Mixtures", A. J. Nonaka and Y. Sun and J. B. Bell and A. Donev, to appear in CAMCOS, 2015 [ArXiv:1410.2300].

and the Stokes solver is described here (see doc/StokesPreconditioners.pdf):

3. "Efficient Variable-Coefficient Finite-Volume Stokes Solvers", M. Cai and A. J. Nonaka and J. B. Bell and B. E. Griffith and A. Donev, Commun. Comput. Phys. (CiCP), 16(5):1263-1297, 2014 [ArXiv:1308.4605].

To compile the code, you need to also download the Fortran BoxLib library from
https://github.com/BoxLib-Codes/BoxLib
You will also need the LAPACK Fortran 90 interface. We package a copy of this in
src_multiSpec/LAPACK95/
and you can edit the Makefile there and build the Fortran 90 interface library. It is possible to run the code without LAPACK support since we also provide iterative procedures instead of dense linear algebra.

The executable codes are in the directory exec. Edit the makefile in
exec/test/GNUMakefile
and then do "make" to build the code.

Example input files are provided in exec/test:

- inputs_thermodiffusion_2d is the thermodiffusion example described in Section III.C.3 in Ref. [1]
- inputs_centrifuge_2d is the ultra-centrifuge barodiffusion example described in Section III.C.3 in Ref. [1]
- inputs_mmi_3d is the mixed-mode instability described in Section IV.B in Ref. [1]
- inputs_dlc_shadow_3d is the DLC instability example described in Section IV.C in Ref. [1].

These tests were originally performed with a previous (binary only) version of the code and have been roughly ported to the multispecies code:

-inputs_liddriven_2d is similar to the lid-driven cavity test from 4.B in Ref. [2]
-inputs_square_2d is similar to the sharp interface square test from 4.C in Ref. [2]
-inputs_khlm_3d/2d is similar to the Kelvin-Helmholtz problem from 4.D in Ref. [2]

Caveats: 

1) The code used to perform spectral analysis of the results is not included in this distribution. The code writes files using the BoxLib file formats and these can be read and visualized using the tools described in the BoxLib manual. In particular the code can also write projections (averages) along a selected axes of the full 2D or 3D grid to reduce I/O requirements and simpify visualization and analysis.

2) BDS advection with non-periodic boundary conditions is not fully implemented, or, if implemented, may not actually be third-order accurate. Correct treatment of boundary conditions in BDS is an open research problem.

3) Quadratic BDS advection (advection_type=3 and 4) is still under development.


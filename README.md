# xSDK Examples v0.3.0

The example codes provided here demonstrate the use of of various xSDK libraries in tandem to solve problems of 
interest.  Each of the library folders has one or more examples codes that are built of that library 
and utilize code integrations with other xSDK libraries.  Running these example codes and
examining the output is a good way to better understand how these libraries can work together. The
code samples are a good place to start for new projects.  More details about the examples can be found 
in the README.md files in the library subfolders.  For more information on the xSDK see <https://xsdk.info/>.

## Example Summary

These examples were tested and verified against xsdk@0.7.0. 

| Example                                               | Libraries                | Description                                       | GPUs           |
|:------------------------------------------------------|:-------------------------|:--------------------------------------------------|:---------------|
| hypre/ij_laplacian.c                                  | HYPRE+SuperLU_Dist       | 2D Laplacian problem                              |                |
| libensemble/test_persistent_aposmm_tao.py             | libEnsemble+PETSc        | 2D constrained optimization problem               |                |
| mfem/hypre-superlu/convdiff.cpp                       | MFEM+HYPRE+SuperLU_Dist  | 2D steady state convective diffusion              |                |
| mfem/ginkgo/mfem_ex1_gko.cpp                          | MFEM+Ginkgo              | 2D Poisson problem with Ginko solver              | ![cuda]        |
| mfem/petsc/obstacle.cpp                               | MFEM+PETSc               | Membrane obstacle problem (min energy functional) |                |
| mfem/strumpack/diffusion-eigen.cpp                    | MFEM+STRUMPACK+HYPRE     | Diffusion eigenvalue problem                      |                |
| mfem/sundials/transient-heat.cpp                      | MFEM+SUNDIALS            | 2D Transient nonlinear heat conduction            |                |
| mfem/hypre/magnetic-diffusion.cpp                     | MFEM+HYPRE               | Steady state magnetic diffusion problem           | ![cuda]        |
| mfem/sundials/advection.cpp                           | MFEM+SUNDIALS            | 2D Time-dependent advection                       | ![cuda]        |
| petsc/ex19.c                                          | PETSc+HYPRE+SuperLU_Dist | 2D nonlinear driven cavity problem                | ![cuda]        |
| plasma/ex1solve.c                                     | PLASMA+SLATE+BLASPP      | Linear system direct solution                     | ![cuda]        |
| sundials/ark_brusselator1D_FEM_sludist.cpp            | SUNDIALS+SuperLU_Dist    | Chemical kinetics brusselator problem             |                |
| sundials/cv_petsc_ex7.c                               | SUNDIALS+PETSc           | 2D nonlinear PDE solution                         |                |
| sundials/cvRoberts_blockdiag_magma.cpp                | SUNDIALS+MAGMA           | Solves a group of chemical kinetics ODEs          | ![cuda] ![hip] | 
| trilinos/SimpleSolve_WithParameters.cpp               | Trilinos+SuperLU_Dist    | Small linear system direct solution               |                |
| strumpack/sparse.cpp                                  | STRUMPACK+ButterflyPACK  | 3D Poisson problem with STRUMPACK preconditioner  |                |

These examples are currently in the repo but will not be enabled in the xsdk-examples spack package until we release a new version of the xSDK.
They can still be built using CMake directly.

| Example                                               | Libraries                | Description                                       | GPUs           |
|:------------------------------------------------------|:-------------------------|:--------------------------------------------------|:---------------|
| amrex/sundials/amrex_sundials_advection_diffusion.cpp | AMReX+SUNDIALS           | 2D Advection-diffusion problem                    | ![cuda] ![hip] |
| mfem/hypre/magnetic-diffusion.cpp                     | MFEM+HYPRE               | Steady state magnetic diffusion problem           | ![hip]         |

## Installing the Examples

The examples can be installed along with the xSDK utilizing the Spack package.

```
spack install xsdk-examples
```

To install with CUDA support,

```
spack install xsdk-examples+cuda cuda_arch=<arch> ^xsdk+cuda cuda_arch=<arch> ^sundials+magma
```

Since `xsdk-examples` depends on the `xsdk` Spack package, Spack will also install `xsdk`. In many cases, it may be easier to install the `xsdk` package (separately) following https://xsdk.info/download/ prior to the `xsdk-examples` package. 

Alternatively the examples can be built and installed with CMake directly:

```
git clone https://github.com/xsdk-project/xsdk-examples
cmake -DCMAKE_PREFIX_PATH=/path/to/libraries -DENABLE_CUDA=<TRUE|FALSE> -DENABLE_HIP=<TRUE|FALSE> -S xsdk-examples/ -B xsdk-examples/builddir
cd xsdk-examples/builddir
make
make install
```

Note, that to build with HIP support CMake must be used directly.

## Running and Testing

xsdk-examples is setup to use `ctest`. Each example in the repository is tested with at least a set of default options. If CMake is used to build xsdk-examples, the tests can be run from the build directory (`builddir` above):
```
ctest .
```
or
```
make test
```
Details on how to run each example code manually (and with different options) can be found in each example folder's README.md file.


[cuda]: https://img.shields.io/badge/-cuda-brightgreen?style=flat "CUDA"
[hip]: https://img.shields.io/badge/-hip-red?style=flat "HIP"

# Tasmanian example

Usage:
```
    ./tasmanian_magma
```

Constructs a surrogate model for a simple reference model using random
samples.
The surrogate uses a spars grid basis, but the samples are not aligned
to the grid and the basis coefficients are computed using
a least-squares approximation.
A linear algebra engine is needed, either BLAS or the GPU native
cuBLAS/rocBLAS or the out-of-core methods implemented in MAGMA.
MAGMA is chosen by default, but the other two can be used as fallback.

#### MAGMA integration

The example is chosen to be small for reference purposes;
however, in practice, the size of the random data can easily exceed
the memory capacity of the GPU.
MAGMA offers out-of-core methods for QR factorization, which can be
used to solve the least-squares problem while leveraging
the computational capacity of the GPU but without the memory limitations.
In comparison, using CPU BLAS will be slow and cuBLAS/rocBLAS will
soon hit the memory limit of the GPU device.

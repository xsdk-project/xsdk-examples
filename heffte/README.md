# heFFTe example

Usage:
```
    mpirun -np 2 heffte_example
```

Computes a 3D FFT transform on data of size 4x4x4 using the GPU backend,
the data is distributed across two MPI ranks.
Two transforms are computed, forward and backward, and the result of the
forward transform is scaled by the full factor `1/4^3 = 1/64`.
Thus, the two transforms are mathematical inverses, which is verified
with a simple check in the end.

#### MAGMA integration

If MAGMA integration is enabled within heFFTe, then the MAGMA linear
algebra methods are used in place of the internal heFFTe kernels.
Most notably, the scaling BLAS method `xscal()` is implemented
much more efficiently in MAGMA compared to the reference kernels
in heFFTe.

MAGMA is enabled at compile time and is automatically used without
and additional runtime switches or options.

/* ex1solve.cpp */

#include <atomic>

#include <cstdio>
#include <cstdlib>

#include <plasma.h>
#include <blas.hh>
#include <lapack.hh>
#include <lapack/flops.hh>
#include <slate/slate.hh>

#include <unistd.h>
#include <sys/time.h>

static double
wtime() {
struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + 1e-6 * tv.tv_usec;
}


static inline blas::Op
t_p2bpp(plasma_enum_t trans) {
    return PlasmaNoTrans == trans ? blas::Op::NoTrans : (PlasmaTrans == trans ? blas::Op::Trans : blas::Op::ConjTrans);
}

static inline blas::Side
s_p2bpp(plasma_enum_t side) {
    return PlasmaLeft == side ? blas::Side::Left : blas::Side::Right;
}

static inline blas::Uplo
u_p2bpp(plasma_enum_t uplo) {
    return PlasmaUpper == uplo ? blas::Uplo::Upper : blas::Uplo::Lower;
}

static inline blas::Diag
d_p2bpp(plasma_enum_t diag) {
    return PlasmaUnit == diag ? blas::Diag::Unit : blas::Diag::NonUnit;
}

static std::atomic<int> Counter_Gemm;

extern "C" void
plasma_core_dgemm(plasma_enum_t transa, plasma_enum_t transb, int m, int n, int k, double alpha, const double *A, int lda, const double *B, int ldb, double beta, double *C, int ldc) {
    ++Counter_Gemm;
    blas::gemm(blas::Layout::ColMajor, t_p2bpp(transa), t_p2bpp(transb), m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}

static std::atomic<int> Counter_Trsm;

extern "C" void
plasma_core_dtrsm(plasma_enum_t side, plasma_enum_t uplo, plasma_enum_t trans, plasma_enum_t diag, int m, int n, double alpha, const double *A, int lda, double *B, int ldb) {
    ++Counter_Trsm;
    blas::trsm(blas::Layout::ColMajor, s_p2bpp(side), u_p2bpp(uplo), t_p2bpp(trans), d_p2bpp(diag), m, n, alpha, A, lda, B, ldb);
}

void
drndset(int m, int n, double *A, int A_ld, int seed) {
    double rcp = 1.0-1.0/RAND_MAX;

    srand(seed);

    for (int j = 0; j < n; ++j) 
        for (int i = 0; i < m; ++i) 
            A[i + j * A_ld] = rand() * rcp;
}

int
main(int argc, char *argv[]) {
    int n, nrhs, nb, ib, A_ld, B_ld, *piv, info;
    double *A, *B;

    n = 1000;
    nrhs = 1;

    for (int i = 1; i < argc && argv[i]; ++i)
        if (strncmp("--n=", argv[i], 2+1+1) != 0 || sscanf(argv[i]+2+1+1, "%d", &n) <= 0 || n < 1)
            n = 1000;

    A_ld = n;
    B_ld = n;

    A = (double *)malloc( sizeof *A * A_ld * n );
    B = (double *)malloc( sizeof *B * B_ld * nrhs );
    piv = (int *)malloc( sizeof *piv * n );

    drndset(n, n,    A, A_ld, 1313);
    drndset(n, nrhs, B, B_ld, 1313);

    plasma_init();

    plasma_get(PlasmaNb, &nb);
    plasma_get(PlasmaIb, &ib);

    Counter_Gemm = 0;
    Counter_Trsm = 0;
    double t = -wtime();
    info = plasma_dgesv(n, nrhs, A, A_ld, piv, B, B_ld);
    t += wtime();
    int cntgemm = Counter_Gemm;
    int cnttrsm = Counter_Trsm;

    plasma_finalize();

    free(piv);
    free(B);
    free(A);

    printf("n=%d nrhs=%d t=%g gflop/s=%g", n, nrhs, t, lapack::Gflop<double>::gesv(n, nrhs) / t);
    printf(" nb=%d ib=%d gemm=%d trsm=%d", nb, ib, cnttrsm, cntgemm);
    printf("\n");

    return 0;
}

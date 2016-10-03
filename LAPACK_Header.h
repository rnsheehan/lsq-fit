#ifndef LAPACK_HEADER_H
#define LAPACK_HEADER_H

// This header is needed to call functions that are defined in the LAPACK DLLs

// DPOTRF computes the Cholesky factorization of a real symmetric positive definite matrix A.
extern "C" void dpotrf_(const char *uplo, const int *N, double *A, const int *lda, int *info);

// DPOTRS solves a system of linear equations A*X = B with a symmetric positive definite matrix A using the Cholesky factorization
// A = U**T*U or A = L*L**T computed by DPOTRF.
extern "C" void dpotrs_(const char *uplo, const int *N, const int *Nrhs, double *A, const int *lda, double *B, const int *ldb, int *info);

// DPOSV computes the solution to a real system of linear equations A * X = B, where A is an N-by-N symmetric positive definite matrix and X and B
// are N-by-NRHS matrices.
extern "C" void dposv_(const char *uplo, const int *N, const int *Nrhs, double *A, const int *lda, double *B, const int *ldb, int *info);

#endif
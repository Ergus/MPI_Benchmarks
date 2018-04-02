#include "util.h"

void check_factorization(const size_t n,
                         double A1[n][n], double A2[n][n],
                         const size_t lda, const double eps)
{
	size_t i, j, k;
	const size_t len = n * n * sizeof(double);
	double *Residual = malloc( len );
	double *L1       = malloc( len );
	double *L2       = malloc( len );
	memset( L1, 0, len );
	memset( L2, 0, len );

	const double alpha = 1.0;
	LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'A', n, n, &A1[0][0], lda, Residual, n);
	LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'L', n, n, &A2[0][0], lda, L1, n);
	LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'L', n, n, &A2[0][0], lda, L2, n);

	cblas_dtrmm(
		CblasColMajor,
		CblasRight,
		CblasLower,
		CblasTrans,
		CblasNonUnit,
		n, n,
		alpha,
		L1, n,
		L2, n
		);

	// Compute the Residual A - L'*L
	// daxpy: y := a*x + y  (with alpha = -1, x:= L'*L
	cblas_daxpy( n*n, -1.0, L2, 1, Residual, 1 );

	const double Rnorm = LAPACKE_dlange(LAPACK_COL_MAJOR, 'I', n, n, Residual, n);
	const double Anorm = LAPACKE_dlange(LAPACK_COL_MAJOR, 'I', n, n, &A1[0][0], n);

	const int info_factorization = isnan(Rnorm/(Anorm*n*eps)) || (Rnorm/(Anorm*n*eps) > 60.0);

	printf("===================================\n");
	printf("Checking the Cholesky Factorization \n");
	printf("-- ||L'L-A||_oo/(||A||_oo.n.eps) = %e \n", Rnorm / (Anorm * n * eps));
	if ( info_factorization) printf("-- Factorization is suspicious... ! \n");
	else printf("-- Factorization is CORRECT ! \n");
	printf("===================================\n");

	free(Residual);
	free(L1);
	free(L2);
}

void oss_potrf(const size_t bsize, double A[bsize][bsize])
{
	int error = LAPACKE_dpotrf(
		LAPACK_COL_MAJOR,
		'L',
		bsize,
		&A[0][0],
		bsize
		);

	assert(error == 0);
}


void oss_trsm(const size_t bsize, double A[bsize][bsize], double B[bsize][bsize])
{
	cblas_dtrsm(
		CblasColMajor,
		CblasRight,
		CblasLower,
		CblasTrans,
		CblasNonUnit,
		bsize, bsize, 1.0,
		&A[0][0], bsize,
		&B[0][0], bsize
		);
}

void oss_syrk(const size_t bsize,
              double A[bsize][bsize],
              double B[bsize][bsize])
{
	cblas_dsyrk(
		CblasColMajor,
		CblasLower,
		CblasNoTrans,
		bsize, bsize, -1.0,
		&A[0][0], bsize,
		1.0,
		&B[0][0], bsize
		);
}

void oss_gemm(const size_t bsize,
              double A[bsize][bsize],
              double B[bsize][bsize],
              double C[bsize][bsize])
{
	cblas_dgemm(
		CblasColMajor,
		CblasNoTrans,
		CblasTrans,
		bsize, bsize, bsize, -1.0,
		&A[0][0], bsize,
		&B[0][0], bsize,
		1.0,
		&C[0][0], bsize
		);
}


void fill_block_ij(const size_t nblocks, const size_t bsize,
                   double matrix[nblocks][nblocks][bsize][bsize],
                   const size_t i, const size_t j)
{
	const size_t seed = i * nblocks + j;
	size_t k, l;
	struct drand48_data status;       	// using re-entrant version for rng
	srand48_r(seed, &status);
	double rnd1, rnd2;

	for (k = 0; k < bsize; ++k) {
		for (l = 0; l < bsize; ++l) {
			drand48_r(&status, &rnd1);
			drand48_r(&status, &rnd2);

			const double val = rnd1 * rnd2;
			matrix[i][j][k][l] = val;
			matrix[j][i][l][k] = val;
		}
	}
}


void fill_block_ii(const size_t nblocks, const size_t bsize,
                   double matrix[nblocks][nblocks][bsize][bsize],
                   const size_t i)
{
	const size_t seed = i * nblocks + i, ld = nblocks * bsize;
	size_t k, l;
	double rnd1, rnd2;
	struct drand48_data status;     	// using re-entrant version for rng

	srand48_r(seed, &status);

	for (k = 0; k < bsize; ++k) {
		drand48_r(&status, &rnd1);
		drand48_r(&status, &rnd2);
		matrix[i][i][k][k] = rnd1 * rnd2 + ld;

		for (l = k+1; l < bsize; ++l) {
			drand48_r(&status, &rnd1);
			drand48_r(&status, &rnd2);

			const double val = rnd1 * rnd2;
			matrix[i][i][k][l] = val;
			matrix[i][i][l][k] = val;
		}
	}
}

void initialize_matrix_blocked(const size_t nblocks, const size_t bsize,
                               double matrix[nblocks][nblocks][bsize][bsize])
{
	size_t i, j;
	for (i = 0; i < nblocks; ++i) {
		fill_block_ii(nblocks, bsize, matrix, i);

		for (j = i + 1; j < nblocks; ++j)
			fill_block_ij(nblocks, bsize, matrix, i, j);
	}
}

void blocked2flat(size_t nblocks, size_t bsize, size_t ld,
                       double blocked[nblocks][nblocks][bsize][bsize],
                       double flat[nblocks*bsize][nblocks*bsize])
{
	size_t i, j, ii, jj;

	for (i = 0; i < nblocks; ++i)
		for (j = 0; j < nblocks; ++j)
			for (ii = 0; ii < bsize; ++ii)
				for (jj = 0; jj < bsize; ++jj)
					flat[i * bsize + ii][j * bsize + jj] = blocked[i][j][ii][jj];
}

void write_matrix(const char filename[64], size_t ld, double matrix[ld][ld])
{
	size_t i, j;
	printf("Writing matrix in file %s\n", filename);
	FILE *fd = fopen(filename, "w");
	assert(fd);

	for (i = 0; i < ld; ++i) {
		for(j = 0; j < ld; ++j){
			fprintf(fd,"%lf ", matrix[i][j]);
		}
		fprintf(fd,"\n");
	}
	fclose(fd);
}

//! Write the blocked matrix reshaping, I'll keep it binary to allow parallel writing later
void write_matrix_blocked(const char filename[64],
                          const size_t nblocks, const size_t bsize,
                          double matrix[nblocks][nblocks][bsize][bsize])
{
	size_t i, j, k ,l;
	printf("Writing blocked matrix in file %s\n", filename);
	FILE *fd = fopen(filename, "w");
	assert(fd);

	for (i = 0; i < nblocks; ++i) {          // block row
		for (k = 0; k < bsize; ++k) {        // line row within block
			for (j = 0; j < nblocks; ++j) {  // block column
				for(l=0; l < bsize; ++l) {
					fprintf(fd,"%lf ", matrix[i][j][k][l]);
				}
			}
			fprintf(fd,"\n");
		}
	}
	fclose(fd);
}

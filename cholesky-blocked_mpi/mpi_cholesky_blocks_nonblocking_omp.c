#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <assert.h>
#include <stdbool.h>
#include <sys/time.h>
#include <mpi.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <limits.h>

#define HALFBYTE (CHAR_BIT >> 1)
#define GETTAG(row, col) { (row << sizeof(row) * HALFBYTE) | col }

#ifdef HAVE_MKL_H
#   include <mkl.h>
#else
#   include <cblas.h>
#   include <lapacke.h>
#endif

enum blas_cmach_type {
	blas_base      = 151,
	blas_t         = 152,
	blas_rnd       = 153,
	blas_ieee      = 154,
	blas_emin      = 155,
	blas_emax      = 156,
	blas_eps       = 157,
	blas_prec      = 158,
	blas_underflow = 159,
	blas_overflow  = 160,
	blas_sfmin     = 161
};

static inline double BLAS_dpow_di(double x, int n)
{
	double rv = 1.0;

	if (n < 0) {
		n = -n;
		x = 1.0 / x;
	}

	for (; n; n >>= 1, x *= x) {
		if (n & 1) rv *= x;
	}

	return rv;
}

static inline double BLAS_dfpinfo(enum blas_cmach_type cmach)
{
	const double b = 2.0;
	const int    t = 53, m = -1021; // l = 1024, 
	const double eps = BLAS_dpow_di(b, -t);
	const double r   = BLAS_dpow_di(b, m-1);

	//double o = ((1.0 - eps) * BLAS_dpow_di(b, l-1)) * b;

	switch (cmach) {
	case blas_eps:   return eps;
	case blas_sfmin: return r;
	default:
		fprintf(stderr, "%s %d %d %d\n", __func__, -1, cmach, 0);
		abort();
		break;
	}

	return 0.0;
}

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
	LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'A', n, n, (void *) A1, lda, Residual, n);
	LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'L', n, n, (void *) A2, lda, L1, n);
	LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'L', n, n, (void *) A2, lda, L2, n);

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
		(void *) A,
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
		(void *) A, bsize,
		(void *) B, bsize
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
		(void *) A, bsize,
		1.0,
		(void *) B, bsize
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
		(void *) A, bsize,
		(void *) B, bsize,
		1.0,
		(void *) C, bsize
		);
}


void task_fill_block_ij(const size_t nblocks, const size_t bsize,
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

//! This is to initialize diagonal blocks (assuming 1 block/task)
void task_fill_block_ii(const size_t nblocks, const size_t bsize,
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
		task_fill_block_ii(nblocks, bsize, matrix, i);

		for (j = i + 1; j < nblocks; ++j)
			task_fill_block_ij(nblocks, bsize, matrix, i, j);
	}
}

void task_blocked2flat(size_t nblocks, size_t bsize, size_t ld,
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

void printblock(size_t bsize, double A[bsize][bsize]) {
	for (size_t k = 0; k < bsize; ++k) {        // line row within block
		for(size_t l=0; l < bsize; ++l) {
			printf("%.3lf ", A[k][l]);
		}
		printf("\n");
	}
	printf("\n");
}

void printbuffer(size_t lblocks, size_t bsize, double A[lblocks][bsize][bsize]) {
	for (size_t k = 0; k < bsize; ++k) {        // line row within block
		for(size_t block = 0; block < lblocks; ++block){
			for(size_t l=0; l < bsize; ++l) {
				printf("%.3lf ", A[block][k][l]);
			}
			printf(" ");
		}
		printf("\n");
	}
	printf("\n");
}

void printfull(size_t lblocks,size_t nblocks, size_t bsize, double A[lblocks][nblocks][bsize][bsize]) {
	for (size_t i = 0; i < lblocks; ++i) {          // block row
		for (size_t k = 0; k < bsize; ++k) {        // line row within block
			for (size_t j = 0; j < nblocks; ++j) {  // block column
				for(size_t l=0; l < bsize; ++l) {
					printf("%.3lf ", A[i][j][k][l]);
				}
				printf(" ");
			}
			printf("\n");
		}
		printf("\n");
	}
}

void cholesky(size_t nblocks, size_t bsize, size_t rank, size_t wsize,
              double A[nblocks / wsize][nblocks][bsize][bsize])
{
	assert(nblocks % wsize == 0);

	size_t j, k, Ri, Rl, region, dest;

	const size_t bsize2 = bsize * bsize;

	const size_t lblocks = nblocks / wsize;     // local blocks
	const size_t first_block = lblocks * rank;       // first block

	const size_t doubles_region = lblocks * bsize2;
	const size_t nrequests = (wsize - rank -1);

	double (*buffer)[bsize][bsize] = NULL;
	double (*buffer0)[lblocks][bsize][bsize] = NULL;
	MPI_Request *requests = NULL;
	MPI_Request send_request;

	if (rank) {
		buffer0 = malloc(2 * lblocks * bsize2 * sizeof(double)); //buffer[i:nblocks]
		buffer = malloc((wsize - rank) * lblocks * bsize2 * sizeof(double)); //buffer[i:nblocks]
		requests = malloc(nrequests * sizeof(MPI_Request));

		MPI_Irecv((void *) buffer0[0], doubles_region, MPI_DOUBLE,
			          0, rank, MPI_COMM_WORLD, &requests[0]);

		for (region = rank + 1; region < wsize; ++region) {
			const size_t local_first_block = region * lblocks -  first_block;

			MPI_Irecv((void *) buffer[local_first_block], doubles_region, MPI_DOUBLE,
			          0, region, MPI_COMM_WORLD, &requests[region - rank]);
		}
	}

	// Wait data until my turn
	for(Ri = 0; Ri < first_block; ++Ri) { // Process rows from previous processes until get here
		const size_t next_Ri =  Ri + 1;
		const int source = next_Ri / lblocks;
		const int tag = next_Ri << 16;

		const size_t buff0index = Ri%2;

		MPI_Wait(&requests[0], MPI_STATUS_IGNORE);

		if (next_Ri < first_block) {
			const size_t buff0nextindex = next_Ri % 2;
			MPI_Irecv((void *) buffer0[buff0nextindex], doubles_region, MPI_DOUBLE,
			          source, tag | rank, MPI_COMM_WORLD, &requests[0]);
		}

		#pragma omp parallel
		{
			#pragma omp for
			for (j = 0; j < lblocks; ++j) { // local triangle
				for (k = 0; k < j; ++k)
					oss_gemm(bsize, buffer0[buff0index][j], buffer0[buff0index][k], A[k][first_block + j]);

				oss_syrk(bsize, buffer0[buff0index][j], A[j][first_block + j]);
			}

			#pragma omp for
			for (region = rank + 1; region < wsize; ++region) {
				const size_t localindex = region - rank;
				const int rtag = tag | region;
				const size_t local_first_block = region * lblocks -  first_block;

				MPI_Wait(&requests[localindex], MPI_STATUS_IGNORE);

				for (j = 0; j < lblocks; ++j) {
					for (k = 0; k < lblocks; ++k)
						oss_gemm(bsize, buffer[local_first_block + j],
						         buffer0[buff0index][k], A[k][region*lblocks + j]);
				}

				if (next_Ri < first_block)
					MPI_Irecv((void *) buffer[local_first_block], doubles_region, MPI_DOUBLE,
					          source, rtag, MPI_COMM_WORLD, &requests[localindex]);
			}
		}
	}

	// Process my rows
	for (Rl = 0; Rl < lblocks; ++Rl) { // Process my blocks

		oss_potrf(bsize, A[Rl][first_block + Rl]);      // Diagonal Block Factorization

		// regions first to send them
		#pragma omp parallel private(region, j, k, dest)
		{
			#pragma omp for
			for (region = rank + 1; region < wsize; ++region) {
				for (j = 0; j < lblocks; ++j) {
					oss_trsm(bsize, A[Rl][first_block + Rl], A[Rl][region*lblocks + j]);
				}
				for (dest = rank + 1; dest <= region; ++dest) {
					int rtag = (first_block + Rl) << 16 | region;
					MPI_Isend(A[Rl][region*lblocks], doubles_region, MPI_DOUBLE,
					          dest, rtag, MPI_COMM_WORLD, &send_request);
					MPI_Request_free(&send_request);
				}
			}

			// rest of first row (triangle section)
			#pragma omp for
			for (j = Rl + 1; j < lblocks; ++j)
				oss_trsm(bsize, A[Rl][first_block + Rl], A[Rl][first_block + j]);

			#pragma omp for
			for (j = Rl + 1; j < lblocks; ++j) { // local triangle
				for (k = Rl + 1; k < j; ++k)
					oss_gemm(bsize, A[Rl][first_block + j],
					         A[Rl][first_block + k], A[k][first_block + j]);
				oss_syrk(bsize, A[Rl][first_block + j], A[j][first_block + j]);
			}

			#pragma omp for  // absolute indices here in j!!!!!
			for (j = first_block + lblocks; j < nblocks; ++j) {
				for (k = Rl + 1; k < lblocks; ++k)
					oss_gemm(bsize, A[Rl][j], A[Rl][first_block + k], A[k][j]);
			}
		}
	}

	if (rank) {
		free(buffer);
		free(requests);
		}

	MPI_Barrier(MPI_COMM_WORLD);
}

void usage(int argc, char **argv) {
	printf("\nUsage: ./%s <size> <block_size> <check> [out]\n", argv[0]);
	printf("\tsize  : Matrix' order (size x size)\n");
	printf("\tblock  : Blocking factor for the matrix.\n");
	printf("\tcheck  : Whether to check the factorization's result (1/0).\n");
	printf("\tout   : [optional] The file name prefix where the matrix will be written.\n");
	printf("size %% block_size = 0\n");
	printf("Different block sizes may generate different matrixes depending on the test, use the 'out' option with caution.\n\n");
}

int main(int argc, char **argv)
{
	setbuf(stdout, NULL); // Do not buffer prints

	int rank, wsize;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &wsize);

	if (argc < 4) {
		if(!rank)
			usage(argc, argv);
		MPI_Abort(MPI_COMM_WORLD, MPI_ERR_ARG);
	}

	const size_t ld    = strtoul(argv[1], NULL, 10); // Matrix size
	const size_t bsize = strtoul(argv[2], NULL, 10); // Block size
	const bool   check = strtoul(argv[3], NULL,  2); // Factorization

	assert(ld % bsize == 0);
	assert((ld / bsize) % wsize == 0);

	//========= End Command Line ===============
	const size_t elements = ld * ld;
	const size_t lelements = ld * ld / wsize;
	const size_t nblocks = ld / bsize;           // Total number of blocks
	struct timeval start, stop;
	MPI_Request req;

	//=============== Arrays ===================
	double (*matrix)[nblocks][bsize][bsize] = NULL;
	double (*original)[ld] = NULL;
	double (*factorized)[ld] = NULL;

	if (!rank) {
		printf("%-20s -> %s\n" ,"BENCHMARK", argv[0]);
		printf("%-20s -> %lu\n","SIZE"     , ld);
		printf("%-20s -> %lu\n","BSIZE"    , bsize);
		printf("%-20s -> %d\n" ,"CHECK"    , check);

		//======= Allocate matrices ===============
		matrix = malloc(elements * sizeof(double));

		initialize_matrix_blocked(nblocks, bsize, matrix);

		MPI_Iscatter((double *)matrix, lelements, MPI_DOUBLE,
		             MPI_IN_PLACE, lelements, MPI_DOUBLE,
		             0, MPI_COMM_WORLD, &req);

		if (check) {
			original = malloc(elements * sizeof(double));
			factorized = malloc(elements * sizeof(double));

			task_blocked2flat(nblocks, bsize, ld, matrix, original);
			write_matrix_blocked("orig.txt", nblocks, bsize, matrix);
			write_matrix("orig_flat.txt", ld, original);
		}

		MPI_Wait(&req, MPI_STATUS_IGNORE);

		// ===========================================
		printf("Executing the factorization...\n");
		gettimeofday(&start, NULL);
		cholesky(nblocks, bsize, rank, wsize, matrix);
		gettimeofday(&stop, NULL);
		// ===========================================

		MPI_Gather(MPI_IN_PLACE, lelements, MPI_DOUBLE,
		           (double *) matrix, lelements, MPI_DOUBLE,
		           0, MPI_COMM_WORLD);

		double elapsed = 1000000.0 * (stop.tv_sec - start.tv_sec);
		elapsed += stop.tv_usec - start.tv_usec;

		double gflops = (ld * ld * ld) / (elapsed * 3.0e+3);
		printf(  "%-20s -> %f\n"  , "PERFORMANCE(GFlops)", gflops);
		printf(  "%-20s -> %f\n"  , "TIME(s)"            , elapsed);

		//======== Check if set =====================
		if (check) {
			task_blocked2flat(nblocks, bsize, ld, matrix, factorized);

			printf("Writting after\n");
			write_matrix_blocked("fact.txt", nblocks, bsize, matrix);
			write_matrix("fact_flat.txt", ld, factorized);

			printf("Checking the correctness of the factorization...\n");
			const double EPS = BLAS_dfpinfo(blas_eps);
			check_factorization(ld, original, factorized, ld, EPS);

			free(factorized);
			free(original);
		}

	} else {
		matrix = malloc(lelements * sizeof(double));
		MPI_Iscatter(NULL, 0, 0,
		            matrix, lelements, MPI_DOUBLE,
		            0, MPI_COMM_WORLD, &req);
		MPI_Wait(&req, MPI_STATUS_IGNORE);

		cholesky(nblocks, bsize, rank, wsize, matrix);

		MPI_Gather(matrix, lelements, MPI_DOUBLE,
		           NULL, 0, 0,
		           0, MPI_COMM_WORLD);
	}

	free(matrix);
	MPI_Finalize();
	}

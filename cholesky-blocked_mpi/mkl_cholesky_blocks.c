
#include "util.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <mkl_scalapack.h>
#include "mkl.h"


typedef MKL_INT MDESC[ 9 ];

void usage(int argc, char **argv) {
	printf("\nUsage: ./%s <size> <block_size> <check> [out]\n", argv[0]);
	printf("\tsize  : Matrix' order (size x size)\n");
	printf("\tblock  : Blocking factor for the matrix.\n");
	printf("\tcheck  : Whether to check the factorization's result (1/0).\n");
	printf("\tout   : [optional] The file name prefix where the matrix will be written.\n");
	printf("size %% block_size = 0\n");
	printf("Different block sizes may generate different matrixes depending on the test, use the 'out' option with caution.\n\n");
}

/*==== MAIN FUNCTION =================================================*/
int main( int argc, char **argv)
{
	struct timeval t[6];
	int rank, wsize, dims[2] = {0};       // MPI vars
	int context, pos[2] = {0}, info;      // blacs vars
	// Fortran interfaces
	const char order = 'R', trans = 'N';;
	const int zero = 0, one = 1, negone = -1;
	const double zerod = 0, oned = 1;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &wsize);
	
	if (argc < 4) {
		if(!rank)
			usage(argc, argv);
		MPI_Abort(MPI_COMM_WORLD, MPI_ERR_ARG);
	}

	const int ld    = strtoul(argv[1], NULL, 10); // Matrix size
	const int bsize = strtoul(argv[2], NULL, 10); // Block size
	const bool   check = strtoul(argv[3], NULL,  2); // Factorization

	assert(ld % bsize == 0);
	assert((ld / bsize) % wsize == 0);
	assert(ld >= bsize * wsize);

	fprintf(stderr,"# Wsize: %d Rank: %d \n", rank, wsize);
	//========= End Command Line ===============
	const int nblocks = ld / bsize;           // Total number of blocks
	double *A_full = NULL, *A = NULL, *A_fact = NULL;  // Arrays variables

	//========= MPI proc grid ==================
	MPI_Dims_create(wsize, 2, dims); // dims = {rows, cols}
	assert(ld % dims[0] == 0);
	assert(ld % dims[1] == 0);

	//=========== Blacs ========================
	blacs_get_(&negone, &zero, &context);
	blacs_gridinit_(&context, "R", &dims[0], &dims[1]);
	blacs_gridinfo_(&context, &dims[0], &dims[1], &pos[0], &pos[1]);

	const int lrows = numroc_(&ld, &bsize, &pos[0], &zero, &dims[0]);
	const int lcols = numroc_(&ld, &bsize, &pos[1], &zero, &dims[1]);

	// Matrix descriptors
	MDESC descA_full, descA;

	if (!rank) {
		const size_t memsize = ld * ld * sizeof(double);
		// This is for debug, simplify letter
		double (*tmp)[nblocks][bsize][bsize] = malloc(memsize);
		A_full = malloc(memsize);
		assert(tmp && A_full);

		gettimeofday(&t[0], NULL);
		initialize_matrix_blocked(nblocks, bsize, tmp);
		gettimeofday(&t[1], NULL);
		blocked2flat(nblocks, bsize, ld, tmp, (void *) A_full);
		free(tmp);

		if (check) {
			write_matrix("orig_flat.txt", ld, (void *) A_full);
			A_fact = calloc(ld * ld, sizeof(double));
			assert(A_fact);
		}
	}

	// Distributed array in all ranks
	A = calloc(lcols * lrows, sizeof(double));

	// Compute leading dimensions
	int ldd_full = numroc_(&ld, &ld, &pos[0], &zero, &dims[0]);
	ldd_full = (ldd_full > 1 ? ldd_full : 1);

	const int ldd = ( lrows > 1 ? lrows : 1 );

	// Initialize descriptors
	descinit_(descA_full, &ld, &ld, &ld, &ld, &zero, &zero, &context, &ldd_full, &info);
	descinit_(descA, &ld, &ld, &bsize, &bsize, &zero, &zero, &context, &ldd, &info);

	// Scatter
	gettimeofday(&t[2], NULL);
	pdgeadd_(&trans, &ld, &ld, &oned, A_full, &one, &one, descA_full, &zerod, A, &one, &one, descA);

	// Cholesky HEREEEE finally!!!
	gettimeofday(&t[3], NULL);
	pdpotrf( "L", &ld, A, &one, &one, descA, &info);
	gettimeofday(&t[4], NULL);

	// Gather
	if (check) {
		pdgeadd_(&trans, &ld, &ld, &oned, A, &one, &one, descA,
		         &zerod, A_fact, &one, &one, descA_full);

		gettimeofday(&t[5], NULL);
	}
	
	free(A);                   //  Destroy arrays

	if (!rank) {
		const double elapsed = getT(t[3],t[4]);
		const double gflops = (ld * ld * ld) / (elapsed * 3.0e+3);

		printf("%-20s -> %s\n" , "BENCHMARK"		   , argv[0]);
		printf("%-20s -> %lu\n", "WSIZE"     		   , wsize);
		printf("%-20s -> %lu\n", "SIZE"     		   , ld);
		printf("%-20s -> %lu\n", "BSIZE"    		   , bsize);
		printf("%-20s -> %d\n" , "CHECK"    		   , check);
		printf("%-20s -> %lf\n", "PERFORMANCE(GFlops)" , gflops);
		printf("%-20s -> %lf\n", "TIME(init)"          , getT(t[0],t[1]));
		printf("%-20s -> %lf\n", "TIME(scatter)"       , getT(t[2],t[3]));
		printf("%-20s -> %lf\n", "TIME(cholesky)"      , elapsed);

		if (check) {
			printf("%-20s -> %lf\n", "TIME(gather)"        , getT(t[4],t[5]));
			write_matrix("fact_flat.txt", ld, (void *) A_fact);

			printf("# Checking the correctness of the factorization...\n");
			const double EPS = BLAS_dfpinfo(blas_eps);
			check_factorization(ld, (void *) A_full, (void *) A_fact, ld, EPS);

			free(A_fact);
		}
		free(A_full);
	}

	blacs_gridexit_(&context);  // Destroy process grid
	blacs_exit_(&one);            // 1 to call MPI finalize ourself

	MPI_Finalize();
	fprintf(stderr, "# Process %d Finalized\n", rank);
	return 0;
}



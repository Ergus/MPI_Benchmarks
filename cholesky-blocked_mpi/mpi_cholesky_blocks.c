#include "util.h"

void cholesky(size_t nblocks, size_t bsize, size_t rank, size_t wsize,
              double A[nblocks / wsize][nblocks][bsize][bsize])
{
	assert(nblocks % wsize == 0);

	const size_t lblocks = nblocks / wsize;     // local blocks
	const size_t first_block = lblocks * rank;       // first block

	const size_t bsize2 = bsize * bsize;
	const size_t buff_blocks = (nblocks - first_block);
	const size_t buff_size = buff_blocks * bsize2;

	const size_t num_senders = wsize - rank - 1;

	double (*buffer)[bsize][bsize] = NULL;
	MPI_Request *requests;

	if (num_senders)
		 requests = malloc(num_senders * sizeof(MPI_Request));

	if (rank)
		buffer = malloc(buff_size * sizeof(double)); //buffer[i:nblocks]

	size_t i, j, k;

	// Wait data until my turn
	for(i = 0; i < first_block; ++i) { // Process rows from previous processes until get here
		MPI_Recv((double *) buffer, buff_size, MPI_DOUBLE, i/lblocks, i,
		         MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		for (j = 0; j < lblocks; ++j) { // local triangle
			for (k = 0; k < j; ++k)
				oss_gemm(bsize, buffer[j], buffer[k], A[k][first_block + j]);

			oss_syrk(bsize, buffer[j], A[j][first_block + j]);
		}
		for (j = lblocks; j < nblocks; ++j) {
			for (k = 0; k < lblocks; ++k)
				oss_gemm(bsize, buffer[j], buffer[k], A[k][j]);
		}
	}

	// Process my rows
	for (i = 0; i < lblocks; ++i) { // Process my blocks
		const size_t lfirsti = first_block + i;

		oss_potrf(bsize, A[i][lfirsti]);      // Diagonal Block Factorization

		for (j = nblocks - 1; j > lfirsti; --j) {
			oss_trsm(bsize, A[i][lfirsti], A[i][j]);
			if (j % lblocks == 0) {
				const size_t dest = j / lblocks;
				MPI_Isend((double *)A[i][j], (nblocks - j) * bsize2, MPI_DOUBLE,
				          dest, i, MPI_COMM_WORLD, &requests[dest-1]);
				MPI_Request_free(&requests[dest-1]);
			}
		}
		for (j = lfirsti + 1; j < first_block + lblocks; ++j) { // local triangle
			for (k = lfirsti + 1; k < j; ++k)
				oss_gemm(bsize, A[i][j], A[i][k], A[k - first_block][j]);

			oss_syrk(bsize, A[i][j], A[j - first_block][j]);

		}
		for (j = first_block + lblocks; j < nblocks; ++j) {
			for (k = i + 1; k < lblocks; ++k)
				oss_gemm(bsize, A[i][j], A[i][first_block + k], A[k][j]);
		}
	}

	if (num_senders)
		free(requests);

	if (rank)
		free(buffer);

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

		    blocked2flat(nblocks, bsize, ld, matrix, original);
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
			blocked2flat(nblocks, bsize, ld, matrix, factorized);

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

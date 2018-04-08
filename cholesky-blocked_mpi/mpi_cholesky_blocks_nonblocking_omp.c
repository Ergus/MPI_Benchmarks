#include "util.h"

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
	MPI_Request requests[nrequests];

	if (rank) {
		buffer0 = malloc(2 * lblocks * bsize2 * sizeof(double)); //buffer[i:nblocks]
		buffer = malloc((wsize - rank) * lblocks * bsize2 * sizeof(double)); //buffer[i:nblocks]
		assert(buffer0 && buffer);

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

		const size_t buff0index = Ri % 2;

		MPI_Wait(&requests[0], MPI_STATUS_IGNORE);

		if (next_Ri < first_block) {
			const size_t buff0nextindex = next_Ri % 2;
			MPI_Irecv((void *) buffer0[buff0nextindex], doubles_region, MPI_DOUBLE,
			          source, tag | rank, MPI_COMM_WORLD, &requests[0]);
		}

		#pragma omp parallel private (j,k,region)
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
					MPI_Request send_request;

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
		free(buffer0);
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

	int rank, wsize, provided;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
	assert(provided == MPI_THREAD_MULTIPLE);
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
	assert(ld >= bsize * wsize);

	fprintf(stderr,"# Wsize: %d Rank: %d \n", rank, wsize);
	//========= End Command Line ===============
	const size_t elements = ld * ld;
	const size_t lelements = ld * ld / wsize;
	const size_t nblocks = ld / bsize;           // Total number of blocks
	struct timeval t[6];
	MPI_Request req;

	//=============== Arrays ===================
	double (*matrix)[nblocks][bsize][bsize] = NULL;
	double (*original)[ld] = NULL;
	double (*factorized)[ld] = NULL;

	if (!rank) {
		//======= Allocate matrices ===============
		matrix = malloc(elements * sizeof(double));
		assert(matrix);

		gettimeofday(&t[0], NULL);
		initialize_matrix_blocked(nblocks, bsize, matrix);
		gettimeofday(&t[1], NULL);

		MPI_Iscatter((double *)matrix, lelements, MPI_DOUBLE,
		             MPI_IN_PLACE, lelements, MPI_DOUBLE,
		             0, MPI_COMM_WORLD, &req);

		if (check) {
			original = malloc(elements * sizeof(double));
			factorized = malloc(elements * sizeof(double));
			assert(original && factorized);

			blocked2flat(nblocks, bsize, ld, matrix, original);
			write_matrix_blocked("orig.txt", nblocks, bsize, matrix);
			write_matrix("orig_flat.txt", ld, original);
		}

		MPI_Wait(&req, MPI_STATUS_IGNORE);
		gettimeofday(&t[2], NULL);
		// ===========================================
		printf("# ecuting the factorization...\n");
		gettimeofday(&t[3], NULL);
		cholesky(nblocks, bsize, rank, wsize, matrix);
		gettimeofday(&t[4], NULL);
		// ===========================================

		MPI_Gather(MPI_IN_PLACE, lelements, MPI_DOUBLE,
		           (double *) matrix, lelements, MPI_DOUBLE,
		           0, MPI_COMM_WORLD);

		gettimeofday(&t[5], NULL);

		const double elapsed = getT(t[3],t[4]);
		const double gflops = (ld * ld * ld) / (elapsed * 3.0e+3);

		printf("%-20s -> %s\n" , "BENCHMARK"           , argv[0]);
		printf("%-20s -> %lu\n", "WSIZE"     		   , wsize);
		printf("%-20s -> %lu\n", "SIZE"    			   , ld);
		printf("%-20s -> %lu\n", "BSIZE"   			   , bsize);
		printf("%-20s -> %d\n" , "CHECK"   			   , check);
		printf("%-20s -> %lf\n", "PERFORMANCE(GFlops)" , gflops);
		printf("%-20s -> %lf\n", "TIME(init)"          , getT(t[0],t[1]));
		printf("%-20s -> %lf\n", "TIME(scatter)"       , getT(t[1],t[2]));
		printf("%-20s -> %lf\n", "TIME(cholesky)"      , elapsed);
		printf("%-20s -> %lf\n", "TIME(gather)"        , getT(t[4],t[5]));

		//======== Check if set =====================
		if (check) {
			blocked2flat(nblocks, bsize, ld, matrix, factorized);

			printf("# Writing after\n");
			write_matrix_blocked("fact.txt", nblocks, bsize, matrix);
			write_matrix("fact_flat.txt", ld, factorized);

			printf("# Checking the correctness of the factorization...\n");
			const double EPS = BLAS_dfpinfo(blas_eps);
			check_factorization(ld, original, factorized, ld, EPS);

			free(factorized);
			free(original);
		}

	} else {
		matrix = malloc(lelements * sizeof(double));
		assert(matrix);

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
	fprintf(stderr, "# Process %d Finalized\n", rank);
	}

/*
 * Copyright (C) 2021  Jimmy Aguilar Mena
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "matvec-mpi.h"

envinfo _env;

void Initialize(int *argc, char ***argv,
                const size_t dim, const size_t lthreads)
{
	MPI_Init(argc, argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &(_env.rank));
	MPI_Comm_size(MPI_COMM_WORLD, &(_env.worldsize));

	myassert(lthreads > 0);
	myassert(dim >= _env.worldsize);                 // more rows than task size
	modcheck(dim, _env.worldsize);                   // we need to split exactly

	_env.lthreads = lthreads;		                 // threads per process
	_env.tthreads = lthreads * _env.worldsize;		 // Total threads
	_env.first_local_thread = _env.rank * lthreads;  // first local thread id
	omp_set_dynamic(0);                              // Disable dynamic teams
	omp_set_num_threads(lthreads);                   // parallel regions size

	_env.dim = dim;					                 // dimension
	_env.ldim = dim / _env.worldsize;                // rows for local process

	_env.IprintA = (_env.rank == 0);						// prints A
	_env.IprintB = (imin(1,_env.worldsize-1) == _env.rank); // prints B
	_env.IprintC = (imin(2,_env.worldsize-1) == _env.rank); // prints C

	// test that the cpuset >= number of local threads
	cpu_set_t mask;
	int ret = sched_getaffinity(0, sizeof(mask), &mask);
	myassert(ret == 0);

	_env.cpu_count = CPU_COUNT(&mask);
	myassert(_env.cpu_count >= lthreads);
}

void Finalize()
{
	MPI_Finalize();
}

void init(double * const __restrict__ array,
          const size_t rows, const size_t cols)
{
	const size_t fullsize = rows * cols;

	#pragma omp parallel
	{
		size_t i;

		struct drand48_data drand_buf;
		srand48_r(_env.first_local_thread + omp_get_thread_num(), &drand_buf);
		double x;

		#pragma omp for
		for(i = 0; i < fullsize; ++i) {
			drand48_r(&drand_buf, &x);
			array[i] = x;
		}
	}
}


void __print(const double * const mat,
             const size_t rows, const size_t cols,
             const char name[64])
{
	if ((_env.rank < 3)
	    && ((_env.IprintA && (strcmp(name,"A") == 0)) ||
		    (_env.IprintB && (strcmp(name,"B") == 0)) ||
		    (_env.IprintC && (strcmp(name,"C") == 0)))) {

		printf("Printing %s in process %d\n",name,_env.rank);

		char filename[256];
		sprintf(filename,"%s.mat", name);
		FILE *fp = fopen(filename, "w+");
		myassert(fp);

		fprintf(fp, "# name: %s\n", name);
		fprintf(fp, "# type: matrix\n");
		fprintf(fp, "# rows: %lu\n", rows);
		fprintf(fp, "# columns: %lu\n", cols);

		for (size_t i = 0; i < rows; ++i) {
			for(size_t j = 0; j < cols; ++j) {
				fprintf(fp, "%3.8lf ", mat[i * cols + j]);
			}
			fprintf(fp,"\n");
		}
		fclose(fp);
	}
}

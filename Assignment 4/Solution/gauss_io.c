#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

int main(int argc, char** argv) {

	char matrix_name[200], vector_name[200], solution_name[200];
	int rows, columns,rhs_rows;
	int size, rank;
	double total_time, io_time = 0, setup_time, kernel_time, mpi_time = 0;
	double total_start, io_start, setup_start, kernel_start, mpi_start;
	MPI_File fh;
	MPI_Status status;   

	if( argc != 2 ) { 
		perror("The base name of the input matrix and vector files must be given\n"); 
		exit(-1);
	}

	int print_x = 0;

	sprintf(matrix_name,   "%s.mat.pre", argv[1]);
	sprintf(vector_name,   "%s.vec.pre", argv[1]);
	sprintf(solution_name, "%s.sol.pre", argv[1]);

	MPI_Init(&argc, &argv); 
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	if(rank == 0){
		printf("Solving the Ax=b system with Gaussian Elimination:\n");
		printf("READ:  Matrix   (A) file name: \"%s\"\n", matrix_name);
		printf("READ:  RHS      (b) file name: \"%s\"\n", vector_name);
		printf("WRITE: Solution (x) file name: \"%s\"\n", solution_name);
	}

	total_start = MPI_Wtime();

	int row, column, index, i;
	io_start = MPI_Wtime();
	MPI_File_open(MPI_COMM_WORLD, matrix_name , MPI_MODE_RDONLY,  MPI_INFO_NULL, &fh);

	MPI_File_read(fh, &rows, 1,MPI_INT, &status);
	MPI_File_read(fh, &columns, 1,MPI_INT, &status);

	if(rows != columns) {
		perror("Only square matrices are allowed\n");
		MPI_Abort(MPI_COMM_WORLD, -1);MPI_File_open(MPI_COMM_WORLD, matrix_name , MPI_MODE_RDONLY,  MPI_INFO_NULL, &fh);
	} 

	if(rows % size != 0) {
		perror("The matrix should be divisible by the number of processes\n");
		MPI_Abort(MPI_COMM_WORLD, -1);
	}    

	int local_block_size = rows / size;

	// Reading the local block, which each rank read seperately //
	double *matrix_local_block = (double *) malloc((local_block_size * rows) * sizeof(double));
	// filename, offset, buf, count, type, status
	MPI_File_read_at(fh ,rank*local_block_size*rows*sizeof(double) + 2*sizeof(int), matrix_local_block, local_block_size * rows ,MPI_DOUBLE, &status);
	//
	MPI_File_close(&fh);

	MPI_File_open(MPI_COMM_WORLD, vector_name , MPI_MODE_RDONLY,  MPI_INFO_NULL, &fh);
	MPI_File_read(fh, &rhs_rows, 1,MPI_INT, &status);

	if(rhs_rows != rows){
		perror("RHS rows must match the sizes of A");
		MPI_Abort(MPI_COMM_WORLD, -1);
	}

	double* rhs_local_block  = (double *)malloc(rows * sizeof(double));
	MPI_File_read_at(fh ,rank*local_block_size*sizeof(double) + 1*sizeof(int),rhs_local_block,local_block_size * rows ,MPI_DOUBLE, &status);
	MPI_File_close(&fh);

	io_time += MPI_Wtime()-io_start;
	setup_start = MPI_Wtime();

	int process, column_pivot;

	double tmp, pivot;
	double *pivots = (double *) malloc((local_block_size + (rows * local_block_size) + 1) * sizeof(double));
	double *local_work_buffer = (double *) malloc(local_block_size * sizeof(double));
	double *accumulation_buffer = (double *) malloc(local_block_size * 2 * sizeof(double));
	double *solution_local_block = (double *) malloc(local_block_size * sizeof(double));

	setup_time += MPI_Wtime() - setup_start;
	kernel_start = MPI_Wtime();

	for(process = 0; process < rank; process++) {
		mpi_start = MPI_Wtime();
		MPI_Recv(pivots, (local_block_size * rows + local_block_size + 1), MPI_DOUBLE, process, process, MPI_COMM_WORLD, &status);
		mpi_time += MPI_Wtime() - mpi_start;

		for(row = 0; row < local_block_size; row++){
			column_pivot = ((int) pivots[0]) * local_block_size + row;
			for (i = 0; i < local_block_size; i++){
				index = i * rows;
				tmp = matrix_local_block[index + column_pivot];
				for (column = column_pivot; column < columns; column++){
					matrix_local_block[index + column] -=  tmp * pivots[(row * rows) + (column + local_block_size + 1)];
				}
				rhs_local_block[i] -= tmp * pivots[row + 1];
				matrix_local_block[index + column_pivot] = 0.0;
			}
		}
	}

	for(row = 0; row < local_block_size; row++){
		column_pivot = (rank * local_block_size) + row;
		index = row * rows;
		pivot = matrix_local_block[index + column_pivot];
		assert(pivot!= 0);

		for (column = column_pivot; column < columns; column++){
			matrix_local_block[index + column] = matrix_local_block[index + column]/pivot; 
			pivots[index + column + local_block_size + 1] = matrix_local_block[index + column];
		}

		local_work_buffer[row] = (rhs_local_block[row])/pivot;
		pivots[row+1] =  local_work_buffer[row];

		for (i = (row + 1); i < local_block_size; i++) {
			tmp = matrix_local_block[i*rows + column_pivot];
			for (column = column_pivot+1; column < columns; column++){
				matrix_local_block[i*rows+column] -=  tmp * pivots[index + column + local_block_size + 1];
			}
			rhs_local_block[i] -= tmp * local_work_buffer[row];
			matrix_local_block[i * rows + row] = 0;
		}
	}

	for (process = (rank + 1); process < size; process++) {
		pivots[0] = (double) rank;
		mpi_start = MPI_Wtime();
		MPI_Send( pivots, (local_block_size * rows + local_block_size + 1), MPI_DOUBLE, process, rank, MPI_COMM_WORLD);
		mpi_time += MPI_Wtime() - mpi_start;
	} 

	for (process = (rank + 1); process<size; ++process) {
		mpi_start = MPI_Wtime();
		MPI_Recv( accumulation_buffer, (2 * local_block_size), MPI_DOUBLE, process, process, MPI_COMM_WORLD, &status); 
		mpi_time += MPI_Wtime() - mpi_start;

		for (row  = (local_block_size - 1); row >= 0; row--) {
			for (column = (local_block_size - 1);column >= 0; column--) {
				index = (int) accumulation_buffer[column];
				local_work_buffer[row] -= accumulation_buffer[local_block_size + column] * matrix_local_block[row * rows + index];
			}
		}
	}

	for (row = (local_block_size - 1); row >= 0; row--) {
		index = rank * local_block_size + row;
		accumulation_buffer[row] = (double) index;
		accumulation_buffer[local_block_size+row] = solution_local_block[row] = local_work_buffer[row];
		for (i = (row - 1); i >= 0; i--){
			local_work_buffer[i] -= solution_local_block[row] * matrix_local_block[ (i * rows) + index];
		}
	}

	for (process = 0; process < rank; process++){
		mpi_start = MPI_Wtime();
		MPI_Send( accumulation_buffer, (2 * local_block_size), MPI_DOUBLE, process, rank, MPI_COMM_WORLD); 
		mpi_time += MPI_Wtime() - mpi_start;
	}

	kernel_time += MPI_Wtime() - kernel_start;
	//  Output //

	MPI_File_open(MPI_COMM_WORLD, solution_name , MPI_MODE_WRONLY | MPI_MODE_CREATE,  MPI_INFO_NULL, &fh);
	double double_rows=rows;
	MPI_File_write_at(fh,0,&double_rows,1,MPI_DOUBLE,&status);
	MPI_File_write_at(fh,rank*local_block_size*sizeof(double) + 1*sizeof(double),solution_local_block,local_block_size,MPI_DOUBLE,&status);
	MPI_File_close(&fh);

	total_time += MPI_Wtime() - total_start;

	printf("[R%02d] Times: IO: %f; Setup: %f; Compute: %f; MPI: %f; Total: %f;\n", 
			rank, io_time, setup_time, kernel_time, mpi_time, total_time);

	free(matrix_local_block);
	free(rhs_local_block);
	free(pivots);
	free(local_work_buffer);
	free(accumulation_buffer);
	free(solution_local_block);

	MPI_Finalize(); 
	return 0;
}


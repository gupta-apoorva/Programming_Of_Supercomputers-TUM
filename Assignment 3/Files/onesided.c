#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

int main(int argc, char** argv) {
	char matrix_name[200], vector_name[200], solution_name[200];
	MPI_Win window, window2, window3, window4;
	int *columns;
	int *rows;
	int size, rank;
	double **matrix_2d_mapped, *matrix_1D_mapped, *rhs, *solution;
	double total_time, io_time = 0, setup_time, kernel_time, mpi_time = 0;
	double total_start, io_start, setup_start, kernel_start, mpi_start;
	FILE *matrix_file, *vector_file, *solution_file;
	MPI_Status status;     

	if( argc != 2 ) { 
		perror("The base name of the input matrix and vector files must be given\n"); 
		exit(-1);
	}

	int print_a = 0;
	int print_b = 0;
	int print_x = 0;

	sprintf(matrix_name,   "%s.mat", argv[1]);
	sprintf(vector_name,   "%s.vec", argv[1]);
	sprintf(solution_name, "%s.sol", argv[1]);
	
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
	MPI_Win_allocate(sizeof(int),sizeof(int),MPI_INFO_NULL,MPI_COMM_WORLD,&columns,&window);
	MPI_Win_allocate(sizeof(int),sizeof(int),MPI_INFO_NULL,MPI_COMM_WORLD,&rows,&window2);

	int row, column, index;
	if(rank == 0) {
		io_start = MPI_Wtime();
		if ((matrix_file = fopen (matrix_name, "r")) == NULL) {
			perror("Could not open the specified matrix file");
			MPI_Abort(MPI_COMM_WORLD, -1);
		}

		fscanf(matrix_file, "%d %d", rows, columns);     
		if(*rows != *columns) {
			perror("Only square matrices are allowed\n");
			MPI_Abort(MPI_COMM_WORLD, -1);
		}  	
		if(*rows % size != 0) {
			perror("The matrix should be divisible by the number of processes\n");
			MPI_Abort(MPI_COMM_WORLD, -1);
		}  	

		matrix_2d_mapped = (double **) malloc(*rows * sizeof(double *));
		for(row = 0; row < *rows; row++){
			matrix_2d_mapped[row] = (double *) malloc(*rows * sizeof(double));
			for(column = 0; column < *columns; column++){
				fscanf(matrix_file, "%lf", &matrix_2d_mapped[row][column]);
			}
		}
		fclose(matrix_file);

		if ((vector_file = fopen (vector_name, "r")) == NULL){
			perror("Could not open the specified vector file");
			MPI_Abort(MPI_COMM_WORLD, -1);
		}

		int rhs_rows;
		fscanf(vector_file, "%d", &rhs_rows);     
		if(rhs_rows != *rows){
			perror("RHS rows must match the sizes of A");
			MPI_Abort(MPI_COMM_WORLD, -1);
		}

		MPI_Alloc_mem(*rows*sizeof(double),MPI_INFO_NULL,&rhs);
		for (row = 0; row < *rows; row++){
			fscanf(vector_file, "%lf",&rhs[row]);
		}
		fclose(vector_file); 
		io_time += MPI_Wtime() - io_start;

		MPI_Alloc_mem(*rows*(*rows)*sizeof(double),MPI_INFO_NULL,&matrix_1D_mapped);
		index = 0;
		for(row=0; row < *rows; row++){
			for(column=0; column < *columns; column++){
				matrix_1D_mapped[index++] = matrix_2d_mapped[row][column];
			}
		}
		solution = (double *) malloc (*rows * sizeof(double));
	}

	setup_start = MPI_Wtime();

	int i;
	MPI_Win_fence(0,window);
	MPI_Win_fence(0,window2);
	MPI_Get(columns,1,MPI_INT,0,0,1,MPI_INT,window);
	MPI_Get(rows,1,MPI_INT,0,0,1,MPI_INT,window2);
	MPI_Win_fence(0,window);
	MPI_Win_fence(0,window2);
	int local_block_size = *rows / size;
	
	int process, column_pivot;
	double tmp, pivot;
	
	double *matrix_local_block;
 	MPI_Alloc_mem(local_block_size * (*rows) * sizeof(double),MPI_INFO_NULL,&matrix_local_block);
	
	double *rhs_local_block;
	MPI_Alloc_mem(local_block_size * sizeof(double),MPI_INFO_NULL,&rhs_local_block);
	
	double* pivots;	
	double *accumulation_buffer;	
	double *local_work_buffer = (double *) malloc(local_block_size * sizeof(double));
	double *solution_local_block = (double *) malloc(local_block_size * sizeof(double));

	if (rank ==0){
   	MPI_Win_create(matrix_1D_mapped,*rows*(*rows),sizeof(double),MPI_INFO_NULL,MPI_COMM_WORLD,&window3);
   	MPI_Win_create(rhs,*rows,sizeof(double),MPI_INFO_NULL,MPI_COMM_WORLD,&window4);
	}
	else{
	MPI_Win_create(matrix_1D_mapped,0,sizeof(double),MPI_INFO_NULL,MPI_COMM_WORLD,&window3);
   	MPI_Win_create(rhs,0,sizeof(double),MPI_INFO_NULL,MPI_COMM_WORLD,&window4);
	}

	MPI_Win_fence(0,window3);
	MPI_Win_fence(0,window4);
	MPI_Get(matrix_local_block,*rows*local_block_size,MPI_DOUBLE,0,rank*(*rows)*local_block_size,*rows*local_block_size,MPI_DOUBLE,window3);
	MPI_Get(rhs_local_block,local_block_size,MPI_DOUBLE,0,rank*local_block_size,local_block_size,MPI_DOUBLE,window4);
	MPI_Win_fence(0,window3);
	MPI_Win_fence(0,window4);

	setup_time = MPI_Wtime() - setup_start;
	kernel_start = MPI_Wtime();
	
	MPI_Alloc_mem((local_block_size + (*rows * local_block_size) + 1)*sizeof(double),MPI_INFO_NULL,&pivots);
	MPI_Alloc_mem(2*local_block_size*sizeof(double),MPI_INFO_NULL,&accumulation_buffer);	

	MPI_Win window5, window6, window7;
	MPI_Group comm_group, comm_group2, comm_group3, group, group2, group3;
	MPI_Comm_group(MPI_COMM_WORLD, &comm_group);
	MPI_Comm_group(MPI_COMM_WORLD, &comm_group2);
	MPI_Comm_group(MPI_COMM_WORLD, &comm_group3);
	MPI_Win_create(pivots,*rows*local_block_size + local_block_size +1,sizeof(double),MPI_INFO_NULL,MPI_COMM_WORLD,&window5);
	MPI_Win_create(accumulation_buffer,2*local_block_size,sizeof(double),MPI_INFO_NULL,MPI_COMM_WORLD,&window6);

	if (rank ==0){
		MPI_Alloc_mem(*rows*sizeof(double),MPI_INFO_NULL,&solution);
		MPI_Win_create(solution,*rows,sizeof(double),MPI_INFO_NULL,MPI_COMM_WORLD,&window7);
	}
	else{
		MPI_Win_create(solution,0,sizeof(double),MPI_INFO_NULL,MPI_COMM_WORLD,&window7);
	}

	int from_whom_1 = 0;

	for(process = 0; process < rank; process++) {
		mpi_start = MPI_Wtime();
		MPI_Group_incl(comm_group,1,&process,&group);
		MPI_Win_start(group,0,window5);
		MPI_Get(pivots,*rows*local_block_size+local_block_size +1, MPI_DOUBLE, from_whom_1 ,0,*rows*local_block_size+local_block_size+1,MPI_DOUBLE,window5);
		MPI_Win_complete(window5);
		from_whom_1 = from_whom_1 +1;
		mpi_time += MPI_Wtime() - mpi_start;

		for(row = 0; row < local_block_size; row++){
			column_pivot = ((int) pivots[0]) * local_block_size + row;
			for (i = 0; i < local_block_size; i++){
				index = i * (*rows);
				tmp = matrix_local_block[index + column_pivot];
				for (column = column_pivot; column < *columns; column++){
					matrix_local_block[index + column] -=  tmp * pivots[(row * (*rows)) + (column + local_block_size + 1)];
				}
				rhs_local_block[i] -= tmp * pivots[row + 1];
				matrix_local_block[index + column_pivot] = 0.0;
			}
		}
	}

	for(row = 0; row < local_block_size; row++){
		column_pivot = (rank * local_block_size) + row;
		index = row * (*rows);
		pivot = matrix_local_block[index + column_pivot];
		assert(pivot!= 0);

		for (column = column_pivot; column < *columns; column++){
			matrix_local_block[index + column] = matrix_local_block[index + column]/pivot; 
			pivots[index + column + local_block_size + 1] = matrix_local_block[index + column];
		}

		local_work_buffer[row] = (rhs_local_block[row])/pivot;
		pivots[row+1] =  local_work_buffer[row];

		for (i = (row + 1); i < local_block_size; i++) {
			tmp = matrix_local_block[i*(*rows) + column_pivot];
			for (column = column_pivot+1; column < *columns; column++){
				matrix_local_block[i*(*rows)+column] -=  tmp * pivots[index + column + local_block_size + 1];
			}
			rhs_local_block[i] -= tmp * local_work_buffer[row];
			matrix_local_block[i * (*rows) + row] = 0;
		}
	}

	pivots[0] = (double) rank;
	mpi_start = MPI_Wtime();
	for (process = (rank + 1); process < size; process++) {
		MPI_Group_incl(comm_group,1,&process,&group);
		MPI_Win_post(group,0,window5);
		MPI_Win_wait(window5);
	}

	mpi_time += MPI_Wtime() - mpi_start;

	int from_whom = size -1 ;

	for (process = (rank + 1); process<size; ++process) {
		mpi_start = MPI_Wtime();
		MPI_Group_incl(comm_group2,1,&from_whom,&group2);
		MPI_Win_start(group2,0,window6);
		MPI_Get(accumulation_buffer,2*local_block_size, MPI_DOUBLE, from_whom ,0,2*local_block_size,MPI_DOUBLE,window6);
		MPI_Win_complete(window6);

		from_whom = from_whom - 1;
		mpi_time += MPI_Wtime() - mpi_start;

		for (row  = (local_block_size - 1); row >= 0; row--) {
			for (column = (local_block_size - 1);column >= 0; column--) {
				index = (int) accumulation_buffer[column];
				local_work_buffer[row] -= accumulation_buffer[local_block_size + column] * matrix_local_block[row * (*rows) + index];
			}
		}
	}

	for (row = (local_block_size - 1); row >= 0; row--) {
		index = rank * local_block_size + row;
		accumulation_buffer[row] = (double) index;
		accumulation_buffer[local_block_size+row] = solution_local_block[row] = local_work_buffer[row];
		for (i = (row - 1); i >= 0; i--){
			local_work_buffer[i] -= solution_local_block[row] * matrix_local_block[ (i * (*rows)) + index];
		}
	}

	for (process = 0; process < rank; process++){
		MPI_Group_incl(comm_group2,1,&process,&group2);
		MPI_Win_post(group2,0,window6);
		MPI_Win_wait(window6);
	}

	int* rank_zero = (int*) malloc(sizeof(int));
	*rank_zero = 0;

	if(rank == 0) {
		for(i = 0; i < local_block_size; i++){
			solution[i] = solution_local_block[i];
		}
		for(i = 1; i < size; i++){
			MPI_Group_incl(comm_group3,1,&i,&group3);
			MPI_Win_post(group3,0,window7);
			MPI_Win_wait(window7);
		}
		mpi_time += MPI_Wtime() - mpi_start;
	} else {
		mpi_start = MPI_Wtime();
		MPI_Group_incl(comm_group3,1,rank_zero,&group3);
		MPI_Win_start(group3,0,window7);
		MPI_Put(solution_local_block,local_block_size, MPI_DOUBLE, 0 ,rank*local_block_size,local_block_size,MPI_DOUBLE,window7);
		MPI_Win_complete(window7);
		mpi_time += MPI_Wtime() - mpi_start;
	}
	
	kernel_time = MPI_Wtime() - kernel_start;

	if (rank == 0) {
		io_start = MPI_Wtime();
		if ((solution_file = fopen(solution_name, "w+")) == NULL) {
			perror("Could not open the solution file");
			MPI_Abort(MPI_COMM_WORLD, -1);
		}

		fprintf(solution_file, "%d\n", *rows);
		for(i = 0; i < *rows; i++) {
			fprintf(solution_file, "%f ", solution[i]);
		}
		fprintf(solution_file, "\n");
		fclose(solution_file);
		io_time += MPI_Wtime() - io_start;

		if(print_a){
			printf("\nSystem Matrix (A):\n");
			for (row = 0; row < *rows; row++) {
				for (column = 0; column < *columns; column++){
					printf("%4.1f ", matrix_2d_mapped[row][column]);
				}
				printf("\n");
			}
		}

		if(print_b){
			printf("\nRHS Vector (b):\n");
			for (row = 0; row < *rows; row++) {
				printf("%4.1f\n", rhs[row]);
			}
		}

		if(print_x){
			printf("\n\nSolution Vector (x):\n");
			for(row = 0; row < *rows; row++){
				printf("%4.4f\n",solution[row]);
			}
		}
	}

	total_time = MPI_Wtime() - total_start;

	printf("[R%02d] Times: IO: %f; Setup: %f; Compute: %f; MPI: %f; Total: %f;\n", 
			rank, io_time, setup_time, kernel_time, mpi_time, total_time);

	if(rank == 0){
		for(i = 0; i < *rows; i++){
			free(matrix_2d_mapped[i]);
		}
		free(matrix_2d_mapped);
		free(rhs);
		free(solution);
	}
	MPI_Win_free(&window);
	MPI_Win_free(&window2);
	MPI_Win_free(&window3);
	MPI_Win_free(&window4);
	MPI_Win_free(&window5);
	MPI_Win_free(&window6);
	MPI_Win_free(&window7);
	free(matrix_local_block);
	free(rhs_local_block);
	free(pivots);
	free(local_work_buffer);
	free(accumulation_buffer);
	free(solution_local_block);

	MPI_Finalize(); 
	return 0;
}

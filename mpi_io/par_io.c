#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

int main(int argc, char** argv) {

	char matrix_name[200],matrix_name_par[200];
	int rows, columns, size, rank;
	double **matrix_2d_mapped, *matrix_1D_mapped;
	FILE *matrix_file;
	MPI_File fh;
	MPI_Status status;   

	double io_time, total_time, par_io_time;  

	sprintf(matrix_name,   "%s.mat", argv[1]);
	sprintf(matrix_name_par,   "%s.mat.pre", argv[1]);

	MPI_Init(&argc, &argv); 
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	io_time = MPI_Wtime();

	int row, column, index;
		if(rank == 0) {
			matrix_file = fopen (matrix_name, "r");
			fscanf(matrix_file, "%d %d", &rows, &columns); 

			matrix_2d_mapped = (double **) malloc(rows * sizeof(double *));
			for(row = 0; row < rows; row++){
				matrix_2d_mapped[row] = (double *) malloc(rows * sizeof(double));
				for(column = 0; column < columns; column++){
					fscanf(matrix_file, "%lf", &matrix_2d_mapped[row][column]);
				}
			}
			fclose(matrix_file);

			matrix_1D_mapped = (double *) malloc(rows * rows * sizeof(double));
			index = 0;
			for(row=0; row < rows; row++){
				for(column=0; column < columns; column++){
					matrix_1D_mapped[index++] = matrix_2d_mapped[row][column];
				}
			}
		}
	io_time = MPI_Wtime() - io_time;
		total_time = MPI_Wtime();
	MPI_Bcast(&rows,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&columns,1,MPI_INT,0,MPI_COMM_WORLD);	

	double local_block_size = rows/size;



	double *matrix_local_block = (double *) malloc(local_block_size * rows * sizeof(double));

	MPI_Scatter(matrix_1D_mapped,local_block_size * rows,MPI_DOUBLE,matrix_local_block,local_block_size * rows,MPI_DOUBLE,0,MPI_COMM_WORLD);
	total_time = MPI_Wtime() - total_time + io_time;

	par_io_time = MPI_Wtime();

	MPI_File_open(MPI_COMM_WORLD, matrix_name_par , MPI_MODE_RDONLY,  MPI_INFO_NULL, &fh);

	MPI_File_read(fh, &rows, 1,MPI_INT, &status);
	MPI_File_read(fh, &columns, 1,MPI_INT, &status);
	MPI_File_read_at(fh ,rank*local_block_size*rows*sizeof(double) + 2*sizeof(int),matrix_local_block,local_block_size * rows ,MPI_DOUBLE, &status);
	MPI_File_close(&fh);

	double total_time_max, par_io_time_max;
	par_io_time = MPI_Wtime() - par_io_time;
	MPI_Reduce(&total_time,&total_time_max,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
	MPI_Reduce(&par_io_time,&par_io_time_max,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

	if (rank ==0)
	printf("%d %f %f\n",size,total_time_max,par_io_time_max );

	MPI_Finalize(); 
	return 0;
}



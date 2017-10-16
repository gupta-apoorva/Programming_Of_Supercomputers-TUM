#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char** argv) {
	
	MPI_File fh; 
	FILE *file;
	MPI_Status  status; 
	int rank,nprocs;
	char file_name[200];
	int rows;
	double* buffer;
	int row=0,index=0;

	sprintf(file_name,   "%s.post", argv[1]);

	MPI_Init(&argc, &argv); 
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs); 

	printf("Postprocessing of the solution file:\n");
	printf("write:Solution file name: \"%s\"\n", file_name);

	MPI_File_open(MPI_COMM_WORLD,argv[1] , MPI_MODE_RDONLY,  MPI_INFO_NULL, &fh); 
	MPI_File_read(fh, &rows, 1,MPI_INT, &status);
	buffer = (double*)malloc(rows*sizeof(double));
	MPI_File_read(fh, buffer, rows ,MPI_DOUBLE, &status);
	MPI_File_close(&fh);

	file = fopen(file_name,"w");
	fprintf(file, "%d  ", rows);
	for (row = 0; row<rows;row++){
		fprintf(file, "%lf  ", buffer[index]);
    		index = index + 1;
	}

		printf("Solution written in ASCII format\n");

	fclose(file);
	free(buffer);

  MPI_Finalize(); 
  return 0;
}


//The file will use only one process to convert.


#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char** argv) {

  MPI_File fh; 
  FILE *matrix_file, *vector_file;
  MPI_Status  status; 
  int rank,nprocs;
  int rows,columns,rows_vector;
  double* matrix_1d_mapped;
  double* vector_mapped;


  char matrix_name[200], vector_name[200],matrix_name_out[200], vector_name_out[200];
  sprintf(matrix_name,   "%s.mat", argv[1]);
  sprintf(vector_name,   "%s.vec", argv[1]);
  sprintf(matrix_name_out,   "%s.mat.pre", argv[1]);
  sprintf(vector_name_out,   "%s.vec.pre", argv[1]);


  MPI_Init(&argc, &argv); 
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs); 

 
  printf("Preprocesssing of the matrix and vector:\n");
  printf("READ:  Matrix   (A) file name: \"%s\"\n", matrix_name);
  printf("READ:  RHS      (b) file name: \"%s\"\n\n", vector_name);


  if ((matrix_file = fopen (matrix_name, "r")) == NULL) {
    perror("Could not open the specified matrix file");
    MPI_Abort(MPI_COMM_WORLD, -1);
  }

  fscanf(matrix_file, "%d %d", &rows,&columns);     
  if(rows != columns) {
    perror("Only square matrices are allowed\n");
    MPI_Abort(MPI_COMM_WORLD, -1);
  }

  int row,column,index;

  matrix_1d_mapped = (double *) malloc(rows *rows* sizeof(double));
  index = 0;
  for(row = 0; row < rows; row++){
    for(column = 0; column < columns; column++){
      fscanf(matrix_file, "%lf", &matrix_1d_mapped[index]);
      index +=1;
    }
  }

  int open_matrix;
  open_matrix = MPI_File_open(MPI_COMM_WORLD, matrix_name_out, MPI_MODE_WRONLY | MPI_MODE_CREATE,  MPI_INFO_NULL, &fh); 
  MPI_File_write(fh, &rows, 1,MPI_INT, &status);
  MPI_File_write(fh, &columns, 1,MPI_INT, &status);
  MPI_File_write(fh ,matrix_1d_mapped,rows*columns,MPI_DOUBLE,&status);
  
  fclose(matrix_file);
  MPI_File_close(&fh);

  printf("Matrix file done : %s\n",matrix_name_out);


  vector_file = fopen (vector_name, "r");
  fscanf(vector_file, "%d", &rows_vector);    
  vector_mapped = (double*) malloc(rows*sizeof(double));
  for (row = 0; row < rows; row++){
    fscanf(vector_file, "%lf",&vector_mapped[row]);
  }


  int open_vector;
  open_vector = MPI_File_open(MPI_COMM_WORLD, vector_name_out, MPI_MODE_WRONLY | MPI_MODE_CREATE,  MPI_INFO_NULL, &fh); 
  MPI_File_write(fh, &rows_vector, 1,MPI_INT, &status);
  MPI_File_write(fh ,vector_mapped,rows,MPI_DOUBLE,&status);
  fclose(vector_file);
  MPI_File_close(&fh);

  printf("Vector file done : %s\n",vector_name_out);


  free(matrix_1d_mapped); 
  free(vector_mapped);
  MPI_Finalize(); 
  return 0;
}

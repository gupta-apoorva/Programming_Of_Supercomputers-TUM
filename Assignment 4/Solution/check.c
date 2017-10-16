#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

int main(int argc, char** argv) {

	char solution_name_p[200];
	char solution_name_s[200];
	int rows;
	int size, rank;
	FILE * sol_s_file;
	MPI_File fh;
	MPI_Status status;   
	int i;

	if( argc != 2 ) { 
		perror("The base name of the input matrix and vector files must be given\n"); 
		exit(-1);
	}


	sprintf(solution_name_p, "../ge_data/size%sx%s.sol.pre", argv[1],argv[1]);
	sprintf(solution_name_s, "../ge_data/size%sx%s.sol", argv[1],argv[1]);
	char *p;		
	rows = strtol(argv[1], &p, 10);


	MPI_Init(&argc, &argv); 
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	if(rank == 0){
		printf("Verifying the solution with serial IO result:\n");
		printf("Parallel Solution (x) file name: \"%s\"\n", solution_name_p);
		printf("Serial Solution (x) file name: \"%s\"\n", solution_name_s);
	}



	double* sol_p = (double *) malloc( rows*sizeof(double));
	double* sol_s = (double *) malloc( rows*sizeof(double));

	// READING sol parllel
	MPI_File_open(MPI_COMM_WORLD, solution_name_p , MPI_MODE_RDONLY,  MPI_INFO_NULL, &fh);
	MPI_File_read(fh, sol_p, rows, MPI_DOUBLE, &status);
	MPI_File_close(&fh);


	// READING sol 
	if ((sol_s_file = fopen (solution_name_s, "r")) == NULL){
		perror("Could not open the specified sol_serial file");
		MPI_Abort(MPI_COMM_WORLD, -1);
	}

	for (i = 0; i < rows; i++){
		fscanf(sol_s_file, "%lf",&sol_s[i]);
	}
	fclose(sol_s_file);

	int flag=0;
	for (i=0;i<rows;i++){
		if ( fabs(sol_p[i]-sol_s[i]) >1e-5 ){
			flag=1;
			printf("ERROR ! sol_p: %f, sol_s: %f\n",sol_p[i],sol_s[i]);
		}	
	}	
	if (flag==0){
		printf("The solution is correct !! NICE JOB !!\n");
	}






	//printf("[R%02d] Times: IO: %f; Setup: %f; Compute: %f; MPI: %f; Total: %f;\n", 
	//		rank, io_time, setup_time, kernel_time, mpi_time, total_time);

	free(sol_p);
	free(sol_s);

	MPI_Finalize(); 
	return 0;
}


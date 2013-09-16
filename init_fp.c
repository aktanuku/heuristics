#include<stdlib.h>
#include"get_input.h"
#include"init_fp.h"


//function allocates memory for three different fingerprint types;
void mem_alloc_fp(float bin_size, float cutoff, int num_type_A, int num_type_B, int type, fingerprint *F){
	int i=0;
	int j=0;
	switch(type)
	{
		case 1:
			(F->rows)=3;					//for FR1/2 rows = 3 (one for AA, one for AB, one for BB)
			F->cols=(int)cutoff/bin_size;			//number of colums is always the number of cuttoff/size of bins
			F->fp = malloc((F->rows)*sizeof(float *));
			for(i=0; i<F->rows; i++)
				(F->fp)[i]=malloc(F->cols*sizeof(float));
			

			for(i=0; i<(F->rows);i++){			//initilize all elements to -1
				for(j=0; j<(F->cols); j++){
					(F->fp)[i][j]=-1;
				}
			}
			
		break;

		case 2:
			F->rows=num_type_A+num_type_B;			//for QE1 the number of rows =number of atoms
			F->cols=(int)cutoff/bin_size;
			F->fp = malloc((F->rows)*sizeof(float *));			
			for(i=0; i<F->rows; i++)
				(F->fp)[i]=malloc(F->cols*sizeof(float));
		
			for(i=0; i<(F->rows);i++){			//initilize all elements to -1
				for(j=0; j<(F->cols); j++){
					(F->fp)[i][j]=-1;
				}
			}

		break;

		case 3:
			F->rows=2*(num_type_A+num_type_B);		//for QE2 each atom in the unit cell has two diff-
			F->cols=(int)cutoff/bin_size;			//-erent fingerprints associated with it
			F->fp = malloc((F->rows)*sizeof(float *));
			for(i=0; i<F->rows; i++)
				(F->fp)[i]=malloc(F->cols*sizeof(float));
		
			
			for(i=0; i<(F->rows);i++){			//initilize all elements to -1
				for(j=0; j<(F->cols); j++){
					(F->fp)[i][j]=-1;
				}
			}


		break;
				
	}

}
//free fingerprint memory
void mem_free_fp(fingerprint *F){
	int i=0;
	for(i=0; i<F->rows; i++)
		free((F->fp)[i]);
	free((F->fp));			

						


}




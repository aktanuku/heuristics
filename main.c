#include<stdio.h>
#include"get_input.h"
#include<string.h>
#include<stdlib.h>
#include"init_fp.h"
#include"calc_bonds.h"
#include"calc_hueristics.h"
#include"SmoothPeak.h"





int main(int argc, char **argv){

if(argc != 5){
	printf("Enter four arguments in the form: POSCAR cutoff binwidth sigma (smoothing value)\n");
	return 0;
}


//get relavent parameters from user
char* filename;
filename=argv[1];

float cutoff=atof(argv[2]);
//printf("%f\n",cutoff);
float bin_width=atof(argv[3]);

float sigma=atof(argv[4]);


cell* A;	
A=(cell*)malloc(sizeof(cell));

//initialize cell
getUserInput(filename,A);	


//calculate unit cell volume
float volume = ((*A).base1[1]*(*A).base2[2]-(*A).base1[2]*(*A).base2[1])*(*A).base3[0]  		
                          +((*A).base1[2]*(*A).base2[0]-(*A).base1[0]*(*A).base2[2])*(*A).base3[1]
                          +((*A).base1[0]*(*A).base2[1]-(*A).base1[1]*(*A).base2[0])*(*A).base3[2];


//allocate memory for fingerprints FR/FR2 QE1 and QE2 respectively
fingerprint* F_FR12;
F_FR12=(fingerprint*)malloc(sizeof(fingerprint));

fingerprint* F_QE1;
F_QE1=(fingerprint*)malloc(sizeof(fingerprint));

fingerprint* F_QE2;
F_QE2=(fingerprint*)malloc(sizeof(fingerprint));


//variables to store LI and NNCE indices
float LIcurr=0;
float NNCEcurr= 0;

//calculate required periodic images of unit cell for cutoff distance
int N1, N2, N3;
calc_expansion(cutoff, A,&N1,&N2,&N3);

//allocate memory for fingerprints
mem_alloc_fp(bin_width, cutoff, A->num_type_A, A->num_type_B, 1, F_FR12);
mem_alloc_fp(bin_width, cutoff, A->num_type_A, A->num_type_B, 2, F_QE1);
mem_alloc_fp(bin_width, cutoff, A->num_type_A, A->num_type_B, 3, F_QE2);

//calculate bonds for fingerprints and NNCE/LI hueristics
bond_driver(N1,N2, N3,*A, F_FR12, F_QE1, F_QE2,&LIcurr,&NNCEcurr, cutoff, bin_width,sigma );


//calculate hueristics
float FR2=calc_FR2( F_FR12,*A,  cutoff/bin_width,  bin_width,  volume);
float FR=calc_FR( F_FR12,*A,  cutoff/bin_width,  bin_width,  volume);
float QE1=calc_QE1(F_QE1, *A, cutoff/bin_width, bin_width,  volume);
float QE2=calc_QE2(F_QE2, *A, cutoff/bin_width, bin_width,  volume);
float LI=calc_LI(*A,LIcurr);
float NNCE=calc_NNCE(*A,NNCEcurr);

//print results
printf("FR,FR2,QE1,QE2,NNCE,LI: %f %f %f %f %f %f\n",FR,FR2, QE1,QE2, NNCE, LI);

//clean up
mem_free_fp(F_FR12);
mem_free_fp(F_QE1);
mem_free_fp(F_QE2);
free(F_FR12);
free(F_QE1);
free(F_QE2);

freecell(A);

free(A);

return 0;

}

#ifndef CALC_HUERISTICS_H
#define CALC_HUERISTICS_H
#include"calc_bonds.h"
#include"get_input.h"

typedef struct{
	int type1;	//element species of 1st atom (0=A , 1=B)
	int type2;	//element species of 2nd atom
	int id_1;	//index in unit cell of 1st atom
	int id_2;	//index in unit cell of 2nd atom
	float dist;
} bond;

//updates LI index everytime we change lattice site
float update_LI(bond *arr,float LI,int n);
//update nearest neirghbor CE index every time we change lattice site
float update_NNCE(bond *arr, float NNCE, int n);


//update FR and FR2 fingerprints
void update_FR_FR2(cell C, int i, int j, int N1, int N2, int N3, float length, fingerprint *F, float volume, float bin_width,int num_bins,float sigma);

//update QE1 fingerprits
void  update_QE1(cell C, int i, int j, int N1, int N2, int N3, float length, fingerprint *F,float volume,float bin_width,int num_bins,float sigma);

//update QE2 fingerprint
void  update_QE2(cell C, int i, int j, int N1, int N2, int N3, float length, fingerprint *F,float volume,float bin_width,int num_bins, float sigma);

//claculates QE2 once fingerprint is compelte
double calc_QE2(fingerprint* F, cell C, int num_bins, float bin_width, float volume);
//calculates QE1 once fingerprint is compelte
double calc_QE1(fingerprint* F, cell C, int num_bins, float bin_width, float volume);
//calculates FR once fingerprint is compelte
float calc_FR(fingerprint* F,cell C, int num_bins, float bin_width, float volume);
	

//claculates FR2 once fingerprint is complete
float calc_FR2(fingerprint* F,cell C, int num_bins, float bin_width, float volume);
//calculates likelihood index after all bonds from all sites ahve been indexed
float calc_LI(cell C, float LIcurr);
//calculates NNCE once all bonds from all sites have been indexed
float calc_NNCE(cell C, float NNCEcurr);

#endif




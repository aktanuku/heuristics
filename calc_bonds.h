#ifndef CALC_BONDS_H
#define CALC_BONDS_H
#include<math.h>
#include"get_input.h"
#include"init_fp.h"
#include"calc_hueristics.h"
#include"SmoothPeak.h"
//calculate the number of unit cell repetitions needed 
//for the desired cutoff distance
void calc_expansion(float cutoff, cell *a, int *N1, int *N2, int *N3);


//calculate distance between two atoms
float calc_dist(int A1, int A2, int A3, int i1, int i2, cell C);


//iterates through all atoms within cutoff distance 
//calculates bonds

void bond_driver(int N1, int N2, int N3, cell C,fingerprint* F_FR_12, fingerprint* F_QE1, fingerprint* F_QE2,
	float* LIcurr,float *NNCEcurr,float cutoff,float bin_width ,float sigma);


//sorts bonds in order of distance for NNCE and LI indexes
void sort_bonds(bond* arr, int n);

#endif

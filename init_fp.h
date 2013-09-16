#ifndef INIT_FP_H
#define INIT_FP_H
#include"get_input.h"

// fingerprint objects
typedef struct{
	int rows;
	int cols;
	float **fp;
	int type; //type=1: FR/FR1 fingerprints, type=2: QE1 fingerprint, type =3: QE2 fingerprint
}fingerprint;


void mem_alloc_fp(float bin_size, float cutoff, int num_type_A, int num_type_B, int type, fingerprint *fp);
void mem_free_fp(fingerprint *fp);

#endif

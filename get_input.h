#ifndef GET_INPUT_H
#define GET_INPUT_H


typedef struct
{
	char* title;			//title of POSCAR
	float mult;			//unit cell multiplier
	float base1[3];			//cell vector a
	float base2[3];			//cell vector b
	float base3[3];			//cell vector c
		
	int num_type_A;			//number of elements of type A
	int num_type_B;			//number of elements of type B
	
	float **coordinates;		//atomic coordinates

}cell;

void getUserInput(char* a,cell* C);

#endif

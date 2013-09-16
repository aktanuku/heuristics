#include "get_input.h"
#include<string.h>
#include<stdio.h>
#include<stdlib.h>
void getUserInput(char* a, cell *C){
	
	//open file
	
	FILE *fp;
	fp = fopen(a,"r");
	char str[100];
	float x,y,z;

	
	//read and store title
	fgets(str,100,fp);
	sscanf(str,"%s",&(C->title));
	//read and store multiplier
	fgets(str,100,fp);
	sscanf(str,"%f",&(C->mult));
	
	
	
	//read and store unit cell vectors
	fgets(str,100,fp);
	sscanf(str,"%f %f %f", &x, &y, &z);
	C->base1[0]=x*(C->mult);(C->base1[1])=y*(C->mult); (C->base1[2])=z*(C->mult);

	fgets(str,100,fp);
	sscanf(str,"%f %f %f", &x, &y, &z);
	(C->base2)[0]=x*(C->mult);(C->base2)[1]=y*(C->mult); (C->base2)[2]=z*(C->mult);
	
	fgets(str,100,fp);
	sscanf(str,"%f %f %f", &x, &y, &z);
	(C->base3)[0]=x*(C->mult);(C->base3)[1]=y*(C->mult); (C->base3)[2]=z*(C->mult);
	
		

	//read and store number of elements of each type	
	fgets(str,100,fp);
	sscanf(str,"%d %d",&(C->num_type_A),&(C->num_type_B));
	
	//read and ignore direct coordinate
	char* Direct;
	fgets(str,100,fp);
	sscanf(str, "%s",&Direct);
	

	//allocate memory for coordinates
	C->coordinates = malloc((C->num_type_A+C->num_type_B)*sizeof(float *));
	int i=0;
	for(i=0; i<(C->num_type_A+C->num_type_B); i++){
		(C->coordinates)[i]=malloc(3*sizeof(float));
	}
	
		

	//read coordinates from file
	for(i=0;i<(C->num_type_A+C->num_type_B); i++){
		
		fgets(str,100,fp);
		sscanf(str,"%f %f %f",&x,&y,&z);
		(C->coordinates)[i][0]=x;	(C->coordinates)[i][1]=y;	(C->coordinates)[i][2]=z;

	}

}


//free unit cell coordinatess after use
void freecell(cell *C){
	int i=0;
	for(i=0; i<C->num_type_A+C->num_type_B; i++){
		free((C->coordinates)[i]);
	}
	free(C->coordinates);

}



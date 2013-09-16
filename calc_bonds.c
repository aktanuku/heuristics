#include"calc_bonds.h"
#include"get_input.h"
#include"calc_hueristics.h"
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"SmoothPeak.h"









//calclate how many periodic images we need to take of each atom in the unit cell
void calc_expansion(float cutoff, cell*a,int *N1,int *N2,int *N3){
		
		//N1 is the expansion along the a cell vector direction
		//N2 is the expansion along the b cell vector direction
		//N3 is the expansion along the c cell vector direction
		//all of which are measured in multiples of the unit cell vector length in each direction
		//initialize them to 1
		*N1=1; *N2=1; *N3=1;

		//calculate length of cell vectors
		float norm1, norm2, norm3;
		norm1 =sqrt(((a->base1)[0])*((a->base1)[0])+((a->base1)[1])*((a->base1)[1])+((a->base1)[2])*((a->base1)[2]));
		norm2 =sqrt(((a->base2)[0])*((a->base2)[0])+((a->base2)[1])*((a->base2)[1])+((a->base2)[2])*((a->base2)[2]));
		norm3 =sqrt(((a->base3)[0])*((a->base3)[0])+((a->base3)[1])*((a->base3)[1])+((a->base3)[2])*((a->base3)[2]));
	

	//multipl each  length by integer values until we excede 2*cuttoff
	float dist1=norm1, dist2=norm2, dist3=norm3;

	while(dist1<=2*cutoff){
		*N1=*N1+1;
		dist1=*N1*norm1;
	}

	while(dist2<=2*cutoff){
		*N2=*N2+1;
		dist2=*N2*norm2;
	}
	while(dist3<=2*cutoff){
		*N3=*N3+1;
		dist3=*N3*norm3;
	}
	
	*N1=*N1+1;
	*N2=*N2+1;
	*N3=*N3+1;

}


float calc_dist(int N1,int N2,int N3,int i1 ,int i2,cell C ){
	

	//convert first coordinate from fractional to cartesian coordinates
	float x[3];
	x[0]=C.coordinates[i1][0]*C.base1[0]+C.coordinates[i1][1]*C.base2[0]+C.coordinates[i1][2]*C.base3[0];
	x[1]=C.coordinates[i1][0]*C.base1[1]+C.coordinates[i1][1]*C.base2[1]+C.coordinates[i1][2]*C.base3[1];	
	x[2]=C.coordinates[i1][0]*C.base1[2]+C.coordinates[i1][1]*C.base2[2]+C.coordinates[i1][2]*C.base3[2];

	
	//convert shift cell to cartesian coordinates
	float A1; float A2; float A3;
	A1= N1*C.base1[0]+N2*C.base2[0]+N3*C.base3[0];
	A2= N1*C.base1[1]+N2*C.base2[1]+N3*C.base3[1];
	A3= N1*C.base1[2]+N2*C.base2[2]+N3*C.base3[2];
	
	//convert second coordinate to cartesian coordinates
 	
	float y[3];
	y[0]=C.coordinates[i2][0]*C.base1[0]+C.coordinates[i2][1]*C.base2[0]+C.coordinates[i2][2]*C.base3[0];
	y[1]=C.coordinates[i2][0]*C.base1[1]+C.coordinates[i2][1]*C.base2[1]+C.coordinates[i2][2]*C.base3[1];	
	y[2]=C.coordinates[i2][0]*C.base1[2]+C.coordinates[i2][1]*C.base2[2]+C.coordinates[i2][2]*C.base3[2];


	//add shift to second coordinate
	y[0]=y[0]+A1;
	y[1]=y[1]+A2;
	y[2]=y[2]+A3;

	//calculate distance
	float dist = sqrt((x[0]-y[0])*(x[0]-y[0])+(x[1]-y[1])*(x[1]-y[1])+(x[2]-y[2])*(x[2]-y[2]));

	return(dist);

}




void bond_driver(int N1, int N2, int N3, cell C,fingerprint* F_FR_12, fingerprint* F_QE1, fingerprint* F_QE2,
						float* LIcurr,float* NNCEcurr,float cutoff,float bin_width,float sigma ){
	
	float LI;
	float NNCE;
	//indexes
	int a=-N1;
	int b=-N2;
	int c=-N3;
	int i=0;
	int j=0;
	
	float bondlength=0;				//bond length
	
	bond* bond_arr;					//array to store bonds
	bond* temp_ptr;					//temporary pointer in case memory allocation fails	
	bond temp_bond;					//temporary bond to store bonds in array

	int n=0; 					//index to keep track of number of bonds
	int nmax=100;

	int num_bins=cutoff/bin_width;		//number of bins

	//calculate cell volume
	float volume = (C.base1[1]*C.base2[2]-C.base1[2]*C.base2[1])*C.base3[0]  		
                          +(C.base1[2]*C.base2[0]-C.base1[0]*C.base2[2])*C.base3[1]
                          +(C.base1[0]*C.base2[1]-C.base1[1]*C.base2[0])*C.base3[2];

	
	//the first three loops iterate over cell vector multipliers
	//the last two iterate over atoms in the unit cell
	while(i<(C.num_type_A+C.num_type_B)){//iterate over all sites in the lattice
		bond_arr=(bond*)malloc(100*sizeof(bond));
		nmax=100;			
		a=-N1;
		while(a<=N1){//these next three indexes iterate over periodic images of atoms
	
			b=-N2;
			while(b<=N2){ 
			
				c=-N3;
				while(c<=N3){
				
					j=0;
					while(j<(C.num_type_A+C.num_type_B)){
					
						//if the same atom do nothing
					       if(calc_dist(a,b,c,i,j,C)!=0){
							
						//calculate bond length							
						bondlength=calc_dist(a,b,c,i,j,C);

						if(bondlength<=cutoff){					
						//store bond
						//check type of first atom
						temp_bond.type1=((i+1)>C.num_type_A);
						//check type of second atom
						temp_bond.type2=((j+1)>C.num_type_A);
						//store lattice indexes of both atoms
						temp_bond.id_1=i;
						temp_bond.id_2=j;
						//sote bondlength
						temp_bond.dist=bondlength;
						//store bond
						bond_arr[n]=temp_bond;
						n=n+1;
						//update fingerprints					
						update_FR_FR2( C, i,  j,  a, b, c, bondlength,F_FR_12 , volume, bin_width, num_bins,sigma);
						update_QE1( C, i,  j,  a, b, c, bondlength,F_QE1 , volume, bin_width, num_bins,sigma);
						update_QE2( C, i,  j,  a, b, c, bondlength,F_QE2 , volume, bin_width, num_bins,sigma);
	
						//if out of memory allocate twice as much
						if(n==(nmax-2)){
							temp_ptr=realloc(bond_arr,sizeof(bond)*2*nmax);
							if(temp_ptr==NULL){
								printf("Memory allocation failed\n");
								}
							else{
								bond_arr=temp_ptr;
								}						
							nmax=2*nmax;
							}
						}
						
						}
					j=j+1;
					}
				c=c+1;	
				}
			b=b+1;
			}
		a=a+1;
		}
	
	//sort bonds and update LI/NNCE index
	sort_bonds(bond_arr,n);
	LI=update_LI(bond_arr,*LIcurr,n);
	*LIcurr=LI;
	NNCE = update_NNCE(bond_arr,*NNCEcurr,n);
	*NNCEcurr=NNCE;
	//
	free(bond_arr);
	n=0;
	

	i=i+1;
	}

}



void sort_bonds(bond* arr,int n){
//selection sort

	float min_length=arr[0].dist;
	int i=0;
	int j=i;
	int mindex=0;
	bond tempbond;
	for(i=0; i<n; i++){
		mindex=i;
		min_length=arr[i].dist;

		for(j=i;j<n;j++){
			if(arr[j].dist<min_length){
				min_length=arr[j].dist;
				mindex=j;
			}
		}

		tempbond=arr[mindex];
		arr[mindex]=arr[i];
		arr[i]=tempbond;
	}

}


				
		
















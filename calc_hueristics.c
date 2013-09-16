#include"get_input.h"
#include"init_fp.h"
#include"calc_bonds.h"
#include"calc_hueristics.h"
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include"SmoothPeak.h"


float update_LI(bond *arr, float LI, int n){
	int i=0;			//index


	float L=0;			//variable for sumation
	int shell=0;			//keeptrack of neighbor level
	float current=arr[0].dist;	//
	int m=0;			//keep track of bonds in neighbor level
	int k=0;			//keep track of total number of bonds
	while(shell<7){
	
	float like=0;
	m=0;
	//while in current neighor level
	while(((fabs(arr[k].dist-current)<.0001))){
		
		//subtract unlike bonds and add like bonds
		if(arr[k].type2!=arr[k].type1)
			like=like-1.;
		else
			like=like+1.;
		k=k+1;
		m=m+1;
	}
	current=arr[k].dist;

	shell=shell+1;
	//divide by the number of bonds and add to previous value
	L=like/m+L;


	
	}
	//average over all neighbor levels
	L=L/7;
	//add to previous result
	L=L+LI;
	
return(L);
}


float update_NNCE(bond *arr, float NNCE, int n){
	int i=0;

	
	float L=0;					//return value
	int shell=0;					//neighbor level
	float current=arr[0].dist;			//bondlength at current shell
	int m=0;					
	int k=0;
	
	//measure correlation up to nearest neighbor level
	while(shell<1){
	
	float like=0;			//variable for summation
	m=0;
	
	//while we are in the nearest neighbor shell,
	while(((fabs(arr[k].dist-current)<.0001))){
		
		//subtract unlike bonds and add like bonds
		if(arr[k].type2!=arr[k].type1)
			like=like-1.;
		else
			like=like+1.;
		k=k+1;
		m=m+1;
	}

	
	current=arr[k].dist;

	shell=shell+1;
	
		
	L=like;


	
	}

	//add to previous NNCE
	L=L+NNCE;
	
return(L);
}


void  update_FR_FR2(cell C, int i, int j, int N1, int N2, int N3, float length, fingerprint *F,float volume,float bin_width,int num_bins, float sigma)
{


	if(((N1==0)&&(N2==0)&&(N3==0))){
		if(j>i){
			//determine type of each atom 0= type A 1=type B 
			//I know the variable names are kind of comfusing here
			//typeA can be either A or B as can typeB
			int typeA=((i+1)>C.num_type_A);
			int typeB=((j+1)>C.num_type_A);
			//determin component of the fingerprint function to index in (see else staement for
			//explination
			int findex = typeA+typeB;
			int NA, NB;
			//for each type, determine the total number of atoms
			if(typeA==0)
				NA=C.num_type_A;
			else
				NA=C.num_type_B;
			
			if(typeB==0)
				NB=C.num_type_A;
			else
				NB=C.num_type_B;

			//determine PCF amplitude addition
			float fing = 1/(4*PI*length*length*NA*NB*2*bin_width*(1/volume));
						
			//for AA and BB we need to multiply by 2
			if(typeA==typeB)
				fing=2*fing;
			smoothPeak(fing,length,bin_width,num_bins,F, findex ,sigma);

			
				

			}
		

		}		
	else{
	
			//determine type of each atom in the bond
			int typeA=((i+1)>C.num_type_A);
			int typeB=((j+1)>C.num_type_A);
			//component of fingerprint function that this bond belonds to
			//typeA=0 or 1 type B =0 or 1 
			//if type A=0 type B=0 then 0+0 =0 (we end up in the first component
			//if type A=0 and type B=1 or vice versa 0+1 =1 (we end up in the second component)
			//if type A=1 and type B=1 we end up with 1+1 =2 the third component
			int findex = typeA+typeB;
			int NA, NB;
			
			//determine number of atoms of each type of bond
			if(typeA==0)
				NA=C.num_type_A;
			else
				NA=C.num_type_B;
			
			if(typeB==0)
				NB=C.num_type_A;
			else
				NB=C.num_type_B;

			//determine  amplitude of PCF addition
			float fing = 1/(4*PI*length*length*NA*NB*2*bin_width*(1/volume));
			
			//AA BB bonds must be counted twice
			if(typeA==typeB)
				fing=2*fing;
			//smooth bond and record
			smoothPeak(fing,length,bin_width,num_bins,F, findex ,sigma);
	
			

			
			

		}

}
 
void  update_QE1(cell C, int i, int j, int N1, int N2, int N3, float length, fingerprint *F,float volume,float bin_width,int num_bins, float sigma){
		
	//if the pair of atoms is inside the unit cell, make sure the bonds are only being counted once
	if(((N1==0)&&(N2==0)&&(N3==0))){
	

	if(j>i){
			
			//determine index of fingerprint matrix where to place the atom
			int findex = i;
			//calculate index amplitude of PCF addition
			float fing = 1/(4*PI*length*length*(C.num_type_A+C.num_type_B)*bin_width*(1/volume));
						
			//smooth bond and record
			smoothPeak(fing,length,bin_width,num_bins,F, findex ,sigma);

			
				

			}
		

		}		
	else{
	
			//index of fingerprint componenet (corresponds to atom sites)
			int findex = i;
			//amplitude of PCF contribution
			float fing = 1/(4*PI*length*length*(float)(C.num_type_A+C.num_type_B)*bin_width*(1/volume));
					//smooth and record
			smoothPeak(fing,length,bin_width,num_bins,F, findex ,sigma);
	
			

			
			

		}			

}




void  update_QE2(cell C, int i, int j, int N1, int N2, int N3, float length, fingerprint *F,float volume,float bin_width,int num_bins,float sigma){
	
	//record type of each atom
	int typeA=((i+1)>C.num_type_A);
	int typeB=((j+1)>C.num_type_A);
	int findex;		
	
	//if site A and B have different types, index in odd slot, 
	//if A and B have hte same type, index in an even slot
	if(typeA==typeB)		
		findex = 2*i;
	else
		findex = 2*i+1;
		
	
	//calculate fingerprint amplitude			
	float fing = 1/(4*PI*length*length*(float)(C.num_type_B+C.num_type_A)*bin_width*(1/volume));
	
	//smooth peak and accumulate
	smoothPeak(fing,length,bin_width,num_bins,F, findex ,sigma);

		
	
	

}

float calc_FR(fingerprint* F,cell C, int num_bins, float bin_width, float volume){
	

	int ii=0;
	int jj=0;
	//NA*NB for all combinations NA and NB
	float numtype[]={C.num_type_A*C.num_type_A, C.num_type_A*C.num_type_B, C.num_type_B*C.num_type_B};
	//sum of NA*NB combinations
	float Twab = C.num_type_A*C.num_type_A+ C.num_type_A*C.num_type_B+ C.num_type_B*C.num_type_B;
	
	float order=0;
	float ord=0;
	//calculate \int F^2 dR
	for ( ii=0;ii<F->rows;ii++)
	{
		ord =0;
		for ( jj=0;jj<F->cols; jj++)
		{
			ord= ord+(pow((F->fp)[ii][jj],2))*bin_width;
			} 
			ord=ord*(numtype[ii]/Twab);
		order = ord+order;
	}
		
	order = order*1/(pow(volume,1/3.));
		
return order;

}


float calc_FR2(fingerprint* F,cell C, int num_bins, float bin_width, float volume){
	

	int ii=0;
	int jj=0;
	//NA*NB for each type combo A and B
	float numtype[]={C.num_type_A*C.num_type_A, C.num_type_A*C.num_type_B, C.num_type_B*C.num_type_B};
	//sum of all NA*NB combonations
	float Twab = C.num_type_A*C.num_type_A+ C.num_type_A*C.num_type_B+ C.num_type_B*C.num_type_B;
	
	float order=0;
	float ord=0;

	//calculate \int F^2R^2 dR
	for ( ii=0;ii<F->rows;ii++)
	{
		ord =0;
		for ( jj=0;jj<F->cols; jj++)
		{	
			ord= ord+(pow((F->fp)[ii][jj],2)*pow(jj*bin_width,2))*bin_width;
			} 
			//normalize by number of atoms
			ord=ord*(numtype[ii]/Twab);
		order = ord+order;
	}
	//divide by volume and multiply by atoms
	order = order*(1/(volume/(C.num_type_A+C.num_type_B)));
		
return order;

}

double cosdist(float* a, float *b, int len){
//calculate cosine distance between two vectors
	int i=0;
	
	double a_norm=0;
	double b_norm=0;
	double D=0;
	
	
	for(i=0; i<len; i++){
		// calcualte norm first vector
		a_norm=a_norm+a[i]*a[i];
		//calculate norm of second vector
		b_norm=b_norm+b[i]*b[i];
		//calculate dot product between the two vectors
		D=D+a[i]*b[i];
	
	}
	//divide dot product by norm
	D=D/sqrt(a_norm*b_norm);
	//shift result to lie between 0 and 1
	D=(1.-D)/2.;
return D;

}



double calc_QE1(fingerprint* F, cell C, int num_bins, float bin_width, float volume){
	

//i am sacraficing memory/speed in the QE calculations in service of readability/accuracy
int i=0;
int j=0;
double QE1=0;
//Step 1: allocate separate two dimensional arrays for type A atoms and type B atoms
	float ** A_mat;
	float ** B_mat;
	//allocate base pointer for array for type A atoms
	A_mat =  malloc(C.num_type_A*sizeof(float*));
	//allocate base pointer for array for type B
	B_mat =  malloc(C.num_type_B*sizeof(float*));
	//allocate memory for each index
	for(i=0; i<C.num_type_A; i++)
		A_mat[i]=malloc((F->cols)*sizeof(float));
	//allocate memory for each index
	for(i=0; i<C.num_type_B; i++)
		B_mat[i]=malloc((F->cols)*sizeof(float));
	

	
//step 2: copy data into respective matrices
	//type A atoms
	for(i=0; i<C.num_type_A; i++){
		for(j=0; j<num_bins; j++){	
		A_mat[i][j]=(F->fp)[i][j];
		}
	}

	//type B atoms
	for(i=0; i<C.num_type_B; i++){
		for(j=0; j<num_bins; j++){
		B_mat[i][j]=(F->fp)[C.num_type_A+i][j];
		}
	}

//step 3: calculate all possible cosine distances for atoms with type A indices:

	int number_of_combinations_A=0; 	//keep track of number of pairs of type A atoms
	double dist=0;
	double S_contrib_A=0;		//quasi entropy contribution
	for(i=0; i<C.num_type_A; i++){
		for(j=i; j<C.num_type_A; j++){
		if(j>i){
		dist=cosdist(A_mat[i],A_mat[j],F->cols);	//distance between the ith and jth fingerprint component	
		S_contrib_A= S_contrib_A+(1.-dist)*log(1.-dist);//calculate contribtuion and add to total
		
		number_of_combinations_A=number_of_combinations_A+1;
		}
		}
	}
		

//step 4: calculate all possible distances between sites of type B
	int number_of_combinations_B=0; 	//keep track of number of pairs of type A atoms
	dist=0;
	double S_contrib_B=0;		//quasi entropy contribution
	for(i=0; i<C.num_type_B; i++){
		for(j=i; j<C.num_type_B-1; j++){
		if(j>i){
		
		dist=cosdist(B_mat[i],B_mat[j],(F->cols));		//distance between the ith and jth fingerprint component
		S_contrib_B= S_contrib_B+(1.-dist)*log(1.-dist);//calculate contribtuion and add to total
	
		number_of_combinations_B=number_of_combinations_B+1;
		}
		}
	}
if(number_of_combinations_A==0){
	number_of_combinations_A=1;
	}
if(number_of_combinations_B==0){
	number_of_combinations_B=1;
	}
//step 5: calculate QE
QE1=-(((float)C.num_type_A)/((float)C.num_type_A+(float)C.num_type_B))*(((float)S_contrib_A)/((float)number_of_combinations_A))
    -((float)C.num_type_B/((float)C.num_type_A+(float)C.num_type_B))*((float)S_contrib_B/(float)number_of_combinations_B);



//step 6: clean up

	//free memory

	//free memory in each index of array for elements of type A
	for(i=0; i<C.num_type_A; i++)
		free(A_mat[i]);
	//free memory for each index of array for elements of type B
	for(i=0; i<C.num_type_B; i++)
		free(B_mat[i]);
	//free base pointer for A matrix
	free(A_mat);
	//free base pointer for B matrix
	free(B_mat);


return QE1;

}




double calc_QE2(fingerprint* F, cell C, int num_bins, float bin_width, float volume){
	

//i am sacraficing memory/speed in the QE calculations in service of readability/accuracy
int i=0;
int j=0;
double QE2=0;
//Step 1: allocate separate two dimensional arrays for type A atoms(like and unlike) and type B atoms (like and unlike)
	float ** A_mat_S;
	float ** B_mat_S;
	float ** A_mat_D;
	float ** B_mat_D;
//**************for this seciton of code anything with an S prefix at the end of it 
//		is meant to denote LIKE interactions anything with a D 
//		at the end of it is meant to denote unlike interactions
//
//
	//allocate base pointer for array for type A atoms
	A_mat_S =  malloc(C.num_type_A*sizeof(float*));
	//allocate base pointer for array for type B
	B_mat_S =  malloc(C.num_type_B*sizeof(float*));
	//allocate memory for each index
	
	//allocate base pointer for array for type A atoms
	A_mat_D =  malloc(C.num_type_A*sizeof(float*));
	//allocate base pointer for array for type B
	B_mat_D =  malloc(C.num_type_B*sizeof(float*));

	for(i=0; i<C.num_type_A; i++)
		A_mat_S[i]=malloc((F->cols)*sizeof(float));
	//allocate memory for each index
	for(i=0; i<C.num_type_B; i++)
		B_mat_S[i]=malloc((F->cols)*sizeof(float));
	for(i=0; i<C.num_type_A; i++)
		A_mat_D[i]=malloc((F->cols)*sizeof(float));
	//allocate memory for each index
	for(i=0; i<C.num_type_B; i++)
		B_mat_D[i]=malloc((F->cols)*sizeof(float));


	
//step 2: copy data into respective matrices
	//type A atoms (like interactions)
	for(i=0; i<C.num_type_A; i++){
		for(j=0; j<num_bins; j++){	
		A_mat_S[i][j]=(F->fp)[2*i][j];
		}
	}

	//type like bonds for type B atoms 
	for(i=0; i<C.num_type_B; i++){
		for(j=0; j<num_bins; j++){
		B_mat_S[i][j]=(F->fp)[2*C.num_type_A+2*i][j];
		}
	}

	//type unlike bonds for type A atoms
	for(i=0; i<C.num_type_A; i++){
		for(j=0; j<num_bins; j++){	
		A_mat_D[i][j]=(F->fp)[2*i+1][j];
		}
	}

	//type unlike bonds for type B atoms
	for(i=0; i<C.num_type_B; i++){
		for(j=0; j<num_bins; j++){
		B_mat_D[i][j]=(F->fp)[2*C.num_type_A+2*i+1][j];
		}
	}

//step 3: calculate all possible cosine distances for atoms with type A indices (like bonds):

	int number_of_combinations_A_S=0; 	//keep track of number of pairs of type A atoms
	double dist=0;
	double S_contrib_A_S=0;		//quasi entropy contribution
	for(i=0; i<C.num_type_A; i++){
		for(j=i; j<C.num_type_A; j++){
		if(j>i){
		dist=cosdist(A_mat_S[i],A_mat_S[j],F->cols);	//distance between the ith and jth fingerprint component	
//		printf("%f %d %d \n",(1-dist)*log(1-dist),i,j);
		S_contrib_A_S= S_contrib_A_S+(1.-dist)*log(1.-dist);//calculate contribtuion and add to total
		number_of_combinations_A_S=number_of_combinations_A_S+1;
		}
		}
	}
		
//step 4: calculate all possible cosine distances for atoms with type A (unlike bonds)

	int number_of_combinations_A_D=0; 	//keep track of number of pairs of type A atoms
	
	double S_contrib_A_D=0;		//quasi entropy contribution
	for(i=0; i<C.num_type_A; i++){
		for(j=i; j<C.num_type_A; j++){
		if(j>i){
		dist=cosdist(A_mat_D[i],A_mat_D[j],F->cols);	//distance between the ith and jth fingerprint component	
//		printf("%f %d %d \n",(1-dist)*log(1-dist),i,j);
		S_contrib_A_D= S_contrib_A_D+(1.-dist)*log(1.-dist);//calculate contribtuion and add to total
		number_of_combinations_A_D=number_of_combinations_A_D+1;
		}
		}
	}
	
//step 5: calculate all possible distances between sites of type B (like bonds)
	int number_of_combinations_B_S=0; 	//keep track of number of pairs of type A atoms
	dist=0;
	double S_contrib_B_S=0;		//quasi entropy contribution
	for(i=0; i<C.num_type_B; i++){
		for(j=i; j<C.num_type_B-1; j++){
		if(j>i){
		dist=cosdist(B_mat_S[i],B_mat_S[j],(F->cols));		//distance between the ith and jth fingerprint component
//		printf("%f %d %d \n",(1-dist)*log(1-dist),i,j);
		S_contrib_B_S= S_contrib_B_S+(1.-dist)*log(1.-dist);//calculate contribtuion and add to total
		number_of_combinations_B_S=number_of_combinations_B_S+1;
		}
		}
	}
//step 6: calculate all possible distances between sites of type B (unlike bonds)
	int number_of_combinations_B_D=0; 	//keep track of number of pairs of type A atoms
	dist=0;
	double S_contrib_B_D=0;		//quasi entropy contribution
	for(i=0; i<C.num_type_B; i++){
		for(j=i; j<C.num_type_B-1; j++){
		if(j>i){
		dist=cosdist(B_mat_D[i],B_mat_D[j],(F->cols));		//distance between the ith and jth fingerprint component
//		printf("%f %d %d \n",(1-dist)*log(1-dist),i,j);
	
		S_contrib_B_D= S_contrib_B_D+(1.-dist)*log(1.-dist);//calculate contribtuion and add to total
		number_of_combinations_B_D=number_of_combinations_B_D+1;
		}
		}
	}


if(number_of_combinations_A_S==0){
	number_of_combinations_A_S=1;
	}
if(number_of_combinations_B_S==0){
	number_of_combinations_B_S=1;
	}
if(number_of_combinations_A_D==0){
	number_of_combinations_A_D=1;
	}
if(number_of_combinations_B_D==0){
	number_of_combinations_B_D=1;
	}



//step 7: calculate QE2
QE2=-((float)C.num_type_A/((float)C.num_type_A+(float)C.num_type_B))*((float)S_contrib_A_S/(float)number_of_combinations_A_S)
    -((float)C.num_type_A/((float)C.num_type_A+(float)C.num_type_B))*((float)S_contrib_A_D/(float)number_of_combinations_A_D)
    -((float)C.num_type_B/((float)C.num_type_A+(float)C.num_type_B))*((float)S_contrib_B_S/(float)number_of_combinations_B_S)
    -((float)C.num_type_B/((float)C.num_type_A+(float)C.num_type_B))*((float)S_contrib_B_D/(float)number_of_combinations_B_D);




//step 6: clean up

	//free memory

	//free memory in each index of array for elements of type A
	for(i=0; i<C.num_type_A; i++)
		free(A_mat_S[i]);
	//free memory for each index of array for elements of type B
	for(i=0; i<C.num_type_B; i++)
		free(B_mat_S[i]);
	
	//free memory in each index of array for elements of type A
	for(i=0; i<C.num_type_A; i++)
		free(A_mat_D[i]);
	//free memory for each index of array for elements of type B
	for(i=0; i<C.num_type_B; i++)
		free(B_mat_D[i]);




	//free base pointers
	free(A_mat_S);
	free(B_mat_S);
	free(A_mat_D);
	free(B_mat_D);


return QE2;

}

//average LI over all sites
float calc_LI(cell C, float LIcurr){
	return(LIcurr/((float)(C.num_type_A+C.num_type_B)));
}

//average NNCE over all sites
float calc_NNCE(cell C, float NNCEcurr){
	printf("NNCEcurr: %f\n", NNCEcurr);
	printf("C.num_type_A: %d\n", C.num_type_A);
	printf("C.num_type_B: %d\n", C.num_type_B);
	printf("returning: %f\n", NNCEcurr /((float)(C.num_type_A+C.num_type_B)));
	return(NNCEcurr /((float)(C.num_type_A+C.num_type_B)));
}








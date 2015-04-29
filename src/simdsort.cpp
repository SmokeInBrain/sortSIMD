#include <ctype.h>
#include <cstdio>
#include <iostream>
#include <cstdlib>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <unistd.h>
#include <getopt.h>
#include <fcntl.h>
#include <stdbool.h>
#include <time.h>
#include <smmintrin.h>
#include <cstring>
#include <climits>
#include <float.h>
#include <cmath>
#include "../include/minheap.h"

typedef struct Variables
{
   char *input;
   char *output;
   int numElements;
   int debug;
   int level;
} Var;


void print_help(){
	printf("Example: executable_name options [ arguments ...]\n");
    printf("    -h  --help             Shows this help\n");
    printf("    -i  --input        	   The name of the input binary file\n");
    printf("    -o  --output           The name of the output binary file\n");
    printf("    -N  --numElements      Number of elements in the input file\n");
    printf("    -d  --debug            If debug==0, no printing at all. Else if debug ==1, print the final sequence, one element each line.\n");
    printf("    -L  --level            Level of the binarian tree to split the list\n");
}

bool esPotenciaDeDos(long x)
{
    return ((x != 0) && !(x & (x - 1)));
}

Var Getoptions(int argc, char **argv){
	Var variables;
	variables.input = (char*)malloc(sizeof(char)*1);
	variables.output = (char*)malloc(sizeof(char)*1);
	variables.numElements = 0;
	variables.debug = 0;
	variables.level=0;

	int next_op;
	const char* const short_op = "hi:o:N:d:L:";

	const struct option long_op[]=
	{
		{"help",0,NULL,'h'},
		{"input",1,NULL,'i'},
		{"output",1,NULL,'o'},
		{"numElements",1,NULL,'N'},
		{"debug",1,NULL,'d'},
		{"level",1,NULL,'L'},
		{ NULL, 0, NULL, 0 }
	};


	if(argc==1){
		printf("ERROR. THE PROGRAM HAS BEEN EXECUTED WITHOUT PARAMETERS OR OPTIONS");
		print_help();
  		exit(EXIT_SUCCESS);
	}

	while(1){
		next_op = getopt_long(argc, argv, short_op, long_op, NULL);
		if(next_op==-1){
			break;
		}

		switch(next_op){
			case 'h':
				print_help();
				exit(EXIT_SUCCESS);
			case 'i':
				strcpy(variables.input,optarg);
				break;
			case 'o':
				strcpy(variables.output,optarg);
				break;
			case 'N':
				if(esPotenciaDeDos(atoi(optarg))){
					variables.numElements = atoi(optarg);
				}else{
					printf("The number of elements in the list is not power of 2\n");
					exit(EXIT_SUCCESS);
				}
				break;
			case 'd':
				variables.debug = atoi(optarg);
				break;
			case 'L':
				if(atoi(optarg)>=0){
					variables.level = atoi(optarg);
				}else{
					printf("The tree level has to be greater or equal to 0\n");
					exit(EXIT_SUCCESS);
				}
				break;
			case '?':
				print_help();
				exit(1);
			case -1:
				break;
			default:
				abort();
		}
	}

	return variables;
}



__m128* minMaxNetwork(__m128 A, __m128 B, __m128 C, __m128 D){
	__m128 amicmi, amacma, bmidmi, bmadma, amicmibmidmi_min, amicmibmidmi_max, amacmabmadma_min, amacmabmadma_max, acbdminmax_min, acbdminmax_max;
	__m128 *array = (__m128*)malloc(4*sizeof(__m128));
	amicmi = _mm_min_ps(A,C);
 	amacma = _mm_max_ps(A,C);
 	bmidmi = _mm_min_ps(B,D);
 	bmadma = _mm_max_ps(B,D);

 	amicmibmidmi_min = _mm_min_ps(amicmi, bmidmi);//A
 	amicmibmidmi_max = _mm_max_ps(amicmi, bmidmi);//FUERA
 	amacmabmadma_min = _mm_min_ps(amacma, bmadma);//FUERA
 	amacmabmadma_max = _mm_max_ps(amacma, bmadma);//D

 	acbdminmax_min = _mm_min_ps(amicmibmidmi_max, amacmabmadma_min);//B
 	acbdminmax_max = _mm_max_ps(amicmibmidmi_max, amacmabmadma_min);//C

 	array[0] = amicmibmidmi_min;
 	array[1] = acbdminmax_min;
 	array[2] = acbdminmax_max;
 	array[3] = amacmabmadma_max;

 	return array;

}

__m128* shuffleNetwork(__m128 A, __m128 B, __m128 C, __m128 D){
	__m128 auxA, auxB, auxC, auxD;
	__m128 auxA2, auxB2, auxC2, auxD2;
	__m128 *array = (__m128*)malloc(4*sizeof(__m128));
	auxA = _mm_shuffle_ps(A, B, _MM_SHUFFLE(1,0,1,0));
 	auxB = _mm_shuffle_ps(A, B, _MM_SHUFFLE(3,2,3,2));
 	auxC = _mm_shuffle_ps(C, D, _MM_SHUFFLE(1,0,1,0));
 	auxD = _mm_shuffle_ps(C, D, _MM_SHUFFLE(3,2,3,2));

 	auxA2 = _mm_shuffle_ps(auxA, auxC, _MM_SHUFFLE(2,0,2,0));
 	auxB2 = _mm_shuffle_ps(auxA, auxC, _MM_SHUFFLE(3,1,3,1));
 	auxC2 = _mm_shuffle_ps(auxB, auxD, _MM_SHUFFLE(2,0,2,0));
 	auxD2 = _mm_shuffle_ps(auxB, auxD, _MM_SHUFFLE(3,1,3,1));

 	array[0]=auxA2;
 	array[1]=auxB2;
 	array[2]=auxC2;
 	array[3]=auxD2;

 	return array;


}
__m128* BMN(__m128 A, __m128 B){
	__m128 *array = (__m128*)malloc(2*sizeof(__m128));
	__m128 Baux, As, Bs;
	__m128 ABmin, ABmax, r1AB, r11AB, r2AB, r22AB;

	Baux = _mm_shuffle_ps(B, B, _MM_SHUFFLE(0,1,2,3));

	As = _mm_shuffle_ps(A,A, _MM_SHUFFLE(3,1,2,0));
 	Bs = _mm_shuffle_ps(Baux,Baux, _MM_SHUFFLE(3,1,2,0));

 	ABmin = _mm_min_ps(As,Bs);
 	ABmax = _mm_max_ps(As,Bs);

 	r1AB = _mm_shuffle_ps(ABmin, ABmax, _MM_SHUFFLE(2,0,2,0));
 	r11AB = _mm_shuffle_ps(r1AB, r1AB, _MM_SHUFFLE(3,1,2,0));

 	r2AB = _mm_shuffle_ps(ABmin, ABmax, _MM_SHUFFLE(3,1,3,1));
 	r22AB = _mm_shuffle_ps(r2AB, r2AB, _MM_SHUFFLE(3,1,2,0));

 	ABmin = _mm_min_ps(r11AB,r22AB);
 	ABmax = _mm_max_ps(r11AB,r22AB);

 	r1AB = _mm_shuffle_ps(ABmin, ABmax, _MM_SHUFFLE(1,0,1,0));
 	r11AB = _mm_shuffle_ps(r1AB, r1AB, _MM_SHUFFLE(3,1,2,0));

 	r2AB = _mm_shuffle_ps(ABmin, ABmax, _MM_SHUFFLE(3,2,3,2));
 	r22AB = _mm_shuffle_ps(r2AB, r2AB, _MM_SHUFFLE(3,1,2,0));

 	ABmin = _mm_min_ps(r11AB,r22AB);
 	ABmax = _mm_max_ps(r11AB,r22AB);

 	r1AB = _mm_shuffle_ps(ABmin, ABmax, _MM_SHUFFLE(1,0,1,0));
 	r11AB = _mm_shuffle_ps(r1AB, r1AB, _MM_SHUFFLE(3,1,2,0));

 	r2AB = _mm_shuffle_ps(ABmin, ABmax, _MM_SHUFFLE(3,2,3,2));
 	r22AB = _mm_shuffle_ps(r2AB, r2AB, _MM_SHUFFLE(3,1,2,0));

 	array[0]=r11AB;
 	array[1]=r22AB;

 	return array;
}

__m128* SIMDMerge(__m128 A, __m128 B, __m128 C, __m128 D){
 
    __m128 *array = (__m128*)malloc(4*sizeof(__m128));
 	__m128 *SMerge1, *SMerge2, *SMerge3;
 	__m128 mask, result, result1, tmp_c, tmp_c2;
 	SMerge1=BMN(A, C);
 	//SMerge1[0] Primera secuencia de salida
 	tmp_c = _mm_shuffle_ps(B, B, _MM_SHUFFLE(0,0,0,0));
 	tmp_c2 = _mm_shuffle_ps(D, D, _MM_SHUFFLE(0,0,0,0));
 	mask = _mm_cmpgt_ps(tmp_c, tmp_c2);
 	result = _mm_blendv_ps(B, D, mask);
 	SMerge2 = BMN(SMerge1[1], result);
    //SMerge2[0] Segunda secuencia de salida
 	mask = _mm_cmpeq_ps(result, B);
 	result1 = _mm_blendv_ps(B, D, mask);
 	SMerge3 = BMN(SMerge2[1], result1);
 	//SMerge3[0] Tercera secuencia de salida
 	//SMerge3[1] Cuarta salida
 	array[0]=SMerge1[0];
 	array[1]=SMerge2[0];
 	array[2]=SMerge3[0];
 	array[3]=SMerge3[1];
 	return array;		
}

float *SIMDSort(float *entrada){
	__m128 a,b,c,d;__m128 Af,Bf,Cf,Df, auxA2, auxB2, auxC2, auxD2;
 	__m128 Ai,Bi,Ci,Di;
 	float p1[4]__attribute__((aligned(16)))= {entrada[0],entrada[1],entrada[2],entrada[3]};
 	float p2[4]__attribute__((aligned(16)))= {entrada[4],entrada[5],entrada[6],entrada[7]};
 	float p3[4]__attribute__((aligned(16)))= {entrada[8],entrada[9],entrada[10],entrada[11]};
 	float p4[4]__attribute__((aligned(16)))= {entrada[12],entrada[13],entrada[14],entrada[15]};
 	float resultado1[4]__attribute__((aligned(16))) = {0.0,0.0,0.0,0.0};
 	float resultado2[4]__attribute__((aligned(16))) = {0.0,0.0,0.0,0.0};
 	float resultado3[4]__attribute__((aligned(16))) = {0.0,0.0,0.0,0.0};
 	float resultado4[4]__attribute__((aligned(16))) = {0.0,0.0,0.0,0.0};
 	/////////////////////////////////////REGISTERS//////////////////////////////////7
 	a = _mm_load_ps(p1);
 	b = _mm_load_ps(p2);
 	c = _mm_load_ps(p3);d = _mm_load_ps(p4);
 	////////////////////////////////////MIN MAX NETWORK//////////////////////////////////////
 	__m128* minMaxArray;
 	minMaxArray = minMaxNetwork(a,b,c,d);
 	////////////////////////////////FLOAT CONVERSION////////////////////////////////////
 	Af = minMaxArray[0];
 	Bf = minMaxArray[1];Cf = minMaxArray[2];
 	Df = minMaxArray[3];
 	free(minMaxArray);
 	////////////////////////////////SHUFFLE NETWORK//////////////////////////////////////
 	__m128* shuffleArray;
 	shuffleArray = shuffleNetwork(Af,Bf,Cf,Df);
 	auxA2 = shuffleArray[0];
 	auxB2 = shuffleArray[1];
 	auxC2 = shuffleArray[2];
 	auxD2 = shuffleArray[3];
 	free(shuffleArray);
 	/////////////////////////////////BMN//////////////////////////////////////////
 	__m128* BMNarray1;
 	__m128* BMNarray2;
 	BMNarray1=BMN(auxA2,auxB2);
 	BMNarray2=BMN(auxC2,auxD2);

 	auxA2 = BMNarray1[0];
 	auxB2 = BMNarray1[1];
 	auxC2 = BMNarray2[0];
 	auxD2 = BMNarray2[1];
 	free(BMNarray1);
 	free(BMNarray2);
 	/////////////////////////////////////////SIMD MERGE///////////////////////////////////////////////////////////////////
 	__m128 *SIMDMergeArray;
 	SIMDMergeArray = SIMDMerge(auxA2, auxB2, auxC2, auxD2);
 	/////////////////////////////////////////SAVE RESULTS///////////////////////////////////////////////////////////////
 	Ai = SIMDMergeArray[0];
 	Bi = SIMDMergeArray[1];
 	Ci = SIMDMergeArray[2];
 	Di = SIMDMergeArray[3];
 	free(SIMDMergeArray);
 	/////////////////STORING///////////////////
 	_mm_store_ps(resultado1 , Ai );
 	_mm_store_ps(resultado2 , Bi );
 	_mm_store_ps(resultado3 , Ci );
 	_mm_store_ps(resultado4 , Di );
 	/////////////DATA TO MEMORY//////////////
 	entrada[0]=resultado1[0];
 	entrada[1]=resultado1[1];
 	entrada[2]=resultado1[2];
 	entrada[3]=resultado1[3];
 	entrada[4]=resultado2[0];
 	entrada[5]=resultado2[1];
 	entrada[6]=resultado2[2];
 	entrada[7]=resultado2[3];
 	entrada[8]=resultado3[0];
 	entrada[9]=resultado3[1];
 	entrada[10]=resultado3[2];
 	entrada[11]=resultado3[3];
 	entrada[12]=resultado4[0];
 	entrada[13]=resultado4[1];
 	entrada[14]=resultado4[2];
 	entrada[15]=resultado4[3];

 	return entrada;
}

float *kWayMergeBruteForce(float **entrada, int numList){
	int *indices = (int *)malloc(numList*sizeof(int));
	float *returnedArray = (float *)malloc(numList*16*sizeof(float));
	int j, i;
	for(i=0;i<numList;i++){
		indices[i]=0;
	}

	bool termino = false;
	j=0;
	while(!termino) {

		float valor_minimo  = FLT_MAX;
		int indice_minimo = -1; 

		for (i = 0; i < numList; ++i)
		{
			if (indices[ i ] < 16) {
				if(entrada[ i ][ indices[ i ] ] <= valor_minimo) {
					valor_minimo = entrada[ i ][ indices[ i ] ];
					indice_minimo = i;
				}
			}
		}

		if (indice_minimo == -1) {
			break;
		}
		indices[ indice_minimo ]++;
		
		returnedArray[j] = valor_minimo;
		j++;
	}
	return returnedArray;
}

void printArray(float *finalArray, int numElements, int debug, char *output){
	//printf("debug %s ", output);
	if(debug==0){
		int i;
		int fd = open(output, O_CREAT | O_RDWR);
		//printf("open %d ", fd);
		for (i = 0; i < numElements; ++i)
		{
			write(fd, &finalArray[i], sizeof(float));
		}
		close(fd);
	}else{
		int i;	
		int suma=0;
		for (i = 0; i < numElements; ++i)
		{
			suma = suma + 1;
			printf("%f\n", finalArray[i]);
		}
		printf("Total: %d\n", suma);
	}
	free(finalArray);
}

int main (int argc, char **argv){
	Var variables;
	variables = Getoptions(argc, argv);
	int level = variables.level;
	int numth = (int)pow(2, level);
	int i,n, mytid, numths;
	float **bufferGeneral = (float **)malloc((variables.numElements/16)*sizeof(float *));
	for (i = 0; i < variables.numElements/16; i++) {
		bufferGeneral[i] = (float *)malloc(16*sizeof(float));	
	}
	int fd1;

  	fd1 = open(variables.input, O_RDONLY);
  	float j;
  	int k=0;
  	int p=0;
  	while (read(fd1, &j, sizeof(float)) != 0) { 
  		if(k%16==0 && k!=0){
  			p++;
  			k=0;
  		}
    		bufferGeneral[p][k] = j;
    		k++;
  	}

 	close(fd1);
 	
 	float **finalArrayParallel = (float **)malloc((numth)*sizeof(float*));

 	int numero_de_elementos;
 	omp_set_num_threads(numth);
 	#pragma omp parallel private(mytid, numths, n)
 	{
 		int j;
 		mytid = omp_get_thread_num();
 		numths = omp_get_num_threads();
 		n = (variables.numElements/16)/numths;
 		int menor = INT_MAX;
 		int mayor = 0;
 		#pragma omp for schedule(static, n)
 		for(j=0;j<variables.numElements/16;j++){
 	   		bufferGeneral[j] = SIMDSort(bufferGeneral[j]);
 	   		if( j < menor){
 	   			menor = j;
 	   		}
 	   		if( mayor < j){
 	   			mayor = j;
 	   		}
 		}


 		numero_de_elementos = mayor - menor + 1;
		finalArrayParallel[mytid] = mergeKArrays(bufferGeneral + menor, numero_de_elementos, 16);
 	}//FIN PARALELIZACIÓN 1

 	int numthnew = numth/2;
 	level--;
 	while(level!=0 && numthnew!=1){
 		omp_set_num_threads(numthnew);
 		#pragma omp parallel private(mytid, numths, n)
 		{
 			mytid = omp_get_thread_num();
 			numths = omp_get_num_threads();
 			int elementsmenores = (variables.numElements/numths)*mytid;
 			int elementsmayores = (variables.numElements/numths)*mytid + (variables.numElements/numths)  - 1;
 			int numero_de_elementos = elementsmayores - elementsmenores +1;
 			printf("Nivel: %d, numero hebras: %d con menor=%d y mayor = %d, total = %d\n", level, numths, elementsmenores, elementsmayores, numero_de_elementos);
 		}
 		level--;
 		numthnew = numthnew/2;
 	}

 	//float *finalArray = (float*)malloc((variables.numElements)*sizeof(float));
 
 	//finalArray = mergeKArrays(finalArrayParallel, numth, variables.numElements/numth);

	//printArray(finalArray, variables.numElements, variables.debug, variables.output);

	//printArray(finalArrayParallel[3], variables.numElements/numth,variables.debug, variables.output);
	//printf("Numero hebras iniciales: %d\n", numth);

 	


	



	return 0;
}


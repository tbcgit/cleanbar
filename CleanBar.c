//----------------------------------------------------------------
//  CLEANBAR.C  : BARCODES AND SINGLE CELL SEQUENCING       
//----------------------------------------------------------------
//  1- Barcodes are searched in order
//  2- We use the direct and reverse_complementary sequences
//  3- We search for consecutive BARCODES.
//  4- We generate output FASTQ file with with 4 BARCODES.
//  5- Also with 2 or 3 Barcodes. 
//  6- Different folders for reads with 4 BARCODES and 2 or 3.
//----------------------------------------------------------------
//  7- We always show both senses of the READ
//  8- We clean the READS from the output files with 4 BARCODES
//----------------------------------------------------------------
//  9- We generate a single output file for 0 or 1 BARCODE.
//     And also for 2 or 3 barcodes.
//
//----------------------------------------------------------------  
//----------------------------------------------------------------
//  Project Manager: Maria Dzunkova
//----------------------------------------------------------------
//  Programmed by: Vicente Arnau.  20-VII-2024        
//----------------------------------------------------------------

//----------------------------------------------------------------
// Compile as :  gcc  -O2 -o CB  CleanBar.c
//----------------------------------------------------------------
// Execute as :  CB  <options>  bar_file  fastq_file
//----------------------------------------------------------------
// Examples of use:
// CB   barcodes.txt  Atrandi_1k.fq
// CB   barcodes.txt  Atrandi_100.fq
// CB   -l 55  -s 99 barcodes.txt  Atrandi_100.fq
// CB   barcodes.txt  Atrandi_2.fastq
//----------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//======>>  Constants that depend on the file with the BARCODES:
#define MAXBAR    24  //--> Number of BARCODES per position
#define BARSIZE    8  //--> BARCODE SIZE
#define LINKSIZE   4  //--> LINK SIZE
//==============================================================

				  
#define CAB_4     320  //--> large value 
                       //--> (greater than CAB+BARSIZE+LINKSIZE)
#define SS_100    100

//=====> ASCII character code: "#"  "." "-"
#define SOS       35  //--> ASCII code
#define PUNTO     46  //--> ASCII code
#define MENOS     45  //--> ASCII code


#define SS          32  //--> Size of BARCODE LABELS
//=====> longitud maxima ETI = 370
#define SSS        600
//=====> longitud maxima SEQ = 320000 nt
#define SSS4  20000000     //--> 40Mnt


//=======>>>>>>>  These folders must already be created:
#define FF_4BAR   "res_4barcodes/"
#define FF_23BAR  "res_23barcodes/"

//=======>>>>>>>  Non-output files (Ficheros de NO salida):
#define No_file   "      ------      "


//------->> Predefined functions (at the end of the file):
 int busca(char *bar, char *sss, int pos, int CAB);
 int inv_com(char *ss1, char *ss2);


int main( int argc, char *argv[])
{
 FILE   *fbar, *ffastq, *fout, *fest;
 FILE   *ff10, *ff23;
 FILE   *facu; //------>> para acumular las reads
 
 char   *ficha, *ficha_base, caracter;
 char   *ff_stats[SS_100], *ff_sum[SS_100];
 char   ff_aux_base[SS_100], ff_aux[SS_100]; //--> for output files
 
 char   *linea, *auxss, *auxll;
 int    i, j, k, contador;
 int    encontrados_A, encontrados_B, encontrados_C, encontrados_D;
 int    pos_D, pos_C, pos_B, pos_A, dd, tengo_bar;
 char   *linea1, *linea2, *linea2_inv, *linea3, *linea4;
 char   *linea2_clean, *linea4_clean;
 char   *lin_dir;
 int    len1, len2, len3, len4, maxlen1, maxlen2;
 int    cuantos=0, cc_dir=0, cc_ic=0;
 int    encontrados=0, completos=0, lim=0;

 char   *nomff1, *nomff2; //----> name of the two output files
 int    pos_max_dir, pos_max_ic;
 int    corte_dir, corte_ic;
 int    pos_D_d, pos_C_d, pos_B_d, pos_A_d;
 
 int    i_barD=-1, min_dd=0;
 int    ii, jj, kk, aux, aux2, pos_bus=0;
 int    indice_D, indice_C, indice_B, indice_A;
 int    uno_cero=0;
 
  float  faux, faux1, faux2;

 char  ss_com[BARSIZE+1]; 
//------------------->> You can change SS to BARSIZE+1:
 char  etiAA[MAXBAR][SS], SSA[MAXBAR][BARSIZE+1];
 char  etiBB[MAXBAR][SS], SSB[MAXBAR][BARSIZE+1];
 char  etiCC[MAXBAR][SS], SSC[MAXBAR][BARSIZE+1];
 char  etiDD[MAXBAR][SS], SSD[MAXBAR][BARSIZE+1];
//----------------------------

//-------------------->> default values
//---------> number of predefined reads with screen output:
 int  LLL_VEO= 2000; 
 //-------->> number of nucleotides analyzed at the start of the read
 int  HEAD   = 88; 
 int  CAB    = 80; //---> (HEAD-BARSIZE)
 


 printf (  "===================================================");
 printf ("\n=     CleanBar : Single Cell data analysis        =");
 printf ("\n===================================================");
 printf ("\n=  Vicente Arnau & Maria Dzunkova . 20-XII-2024   =");
 printf ("\n===================================================\n");


 if (argc == 1) {
   printf("\nSINTAX: ./CB  <options>  BARCODES_File  FASTQ_File ");
   printf(" >  screen_File ");
   printf("\n\nYou can use this: ./CB  --help \n");
   fflush(stdin); getchar();  exit(0);
   }


 if (argc == 2){
	// printf("\n ##%s##", argv[1]); 
	 
	aux=1; aux2=1;
	aux= strcmp ( argv[1], "-help");
	aux2= strcmp ( argv[1], "--help");
	if ((aux==0) || (aux2==0)){
		printf("\nSINTAX: ./CB  <options>  BARCODES_File  FASTQ_File");
		 printf(" >  screen_output_File \n");
		printf("\n<options>: ");
		printf("\n -l\t: Number of nt parsed at the start and the end of a read");
		printf("\n -s\t: Number of reads showed on the screen");	
		printf("\n   \t: screen_output_File is optional \n");
		fflush(stdin); getchar();  exit(0);
		}
	else {
		printf("\nYou can use this: ./CB  --help");
		fflush(stdin); getchar();  exit(0);
		}		
	}
	 
//==============================================================
//========  I reserve dynamic memory for STRINGS:

 auxll= (char *) calloc (SSS, sizeof(char));
 auxss= (char *) calloc (SSS, sizeof(char));
 
 linea = (char *) calloc (SSS, sizeof(char));
 ficha = (char *) calloc (SSS, sizeof(char));
 ficha_base = (char *) calloc (SSS, sizeof(char));
 nomff1 = (char *) calloc (SSS, sizeof(char));
 nomff2 = (char *) calloc (SSS, sizeof(char));
  
 linea1= (char *) calloc (SSS, sizeof(char));
 linea2= (char *) calloc (SSS4, sizeof(char));
 linea3= (char *) calloc (SS,  sizeof(char)); 
 linea4= (char *) calloc (SSS4, sizeof(char)); 
 
 linea2_inv  = (char *) calloc (SSS4, sizeof(char));
 linea2_clean= (char *) calloc (SSS4, sizeof(char));
 linea4_clean= (char *) calloc (SSS4, sizeof(char));
 
 lin_dir   = (char *) calloc (SSS, sizeof(char));
//==============================================================
	 
	 
//------------->> NO OPTIONS:	 
 if (argc == 3){
	 
	strcpy(ficha, argv[1]);

	if ( (fbar=fopen(ficha,"r"))==NULL) {
		printf ("\n Problems with  %s??", ficha);
		fflush(stdin); getchar();  exit(0);
		}
	printf ("\n BARCODES File:  %s", ficha);

	strcpy(ficha, argv[2]);
	if ( (ffastq=fopen(ficha,"r"))==NULL) {
		printf ("\n Problems with  %s??", ficha);
		fflush(stdin); getchar();  exit(0);
		}

	printf ("\n FASTQ File   :  %s\n", ficha);
	strcpy(ficha_base, ficha);
	printf ("\n==================================================");

	}

 else {//----->> WHIT OPTIONS:
     if (argc == 5){
		aux=1;
		aux=  strcmp ( argv[1], "-s");
		if (aux==0){ //---> we modify the lines per screen that we display
		    LLL_VEO= atoi(argv[2]);
			printf ("\nLines per screen = %d\n", LLL_VEO);
			}
		else {
			aux=1;
			aux=  strcmp ( argv[1], "-l");
			if (aux==0){ //---> we modify head (CAB+SIZEBAR)
				HEAD= atoi(argv[2]);
				printf ("\nSequence length analyzed = %d\n", HEAD);
				CAB = HEAD - BARSIZE;
				}			
		}

		
		//---> 	Now I read the data files:
		strcpy(ficha, argv[3]);

		if ( (fbar=fopen(ficha,"r"))==NULL) {
			printf ("\n Problems with  %s??", ficha);
			fflush(stdin); getchar();  exit(0);
			}
		printf ("\n BARCODES File:  %s", ficha);

		strcpy(ficha, argv[4]);
		if ( (ffastq=fopen(ficha,"r"))==NULL) {
			printf ("\n Problems with  %s??", ficha);
			fflush(stdin); getchar();  exit(0);
			}
		printf ("\n FASTQ File   :  %s\n", ficha);
		strcpy(ficha_base, ficha);
		printf ("\n==================================================");
	 
		}
		else if (argc == 7){//------->>  Ultimo IF-ELSE
		
			//--> Primer opcional:
			aux=1;
			aux=  strcmp ( argv[1], "-s");
			if (aux==0){ //---> we modify the lines per screen that we display
		    	LLL_VEO= atoi(argv[2]);
				printf ("\nLines per screen = %d\n", LLL_VEO);
				}
			else {
				aux=1;
				aux=  strcmp ( argv[1], "-l");
				if (aux==0){ //---> we modify head (CAB+SIZEBAR)
					HEAD= atoi(argv[2]);
					printf ("\nSequence length analyzed = %d\n", HEAD);
					CAB = HEAD - BARSIZE;
					}			
				}
		
			//--> Segundo opcional:
			aux=1;
			aux=  strcmp ( argv[3], "-s");
			if (aux==0){ //---> we modify the lines per screen that we display
		    	LLL_VEO= atoi(argv[4]);
				printf ("\nLines per screen = %d\n", LLL_VEO);
				}
			else {
				aux=1;
				aux=  strcmp ( argv[3], "-l");
				if (aux==0){ //---> we modify head (CAB+SIZEBAR)
					HEAD= atoi(argv[4]);
					printf ("\nSequence length analyzed = %d\n", HEAD);
					CAB = HEAD - BARSIZE;
					}			
				}
		
			//---> 	Now I read the data files:
			strcpy(ficha, argv[5]);

			if ( (fbar=fopen(ficha,"r"))==NULL) {
				printf ("\n Problems with  %s??", ficha);
				fflush(stdin); getchar();  exit(0);
				}
			printf ("\n BARCODES File:  %s", ficha);

			strcpy(ficha, argv[6]);
			if ( (ffastq=fopen(ficha,"r"))==NULL) {
				printf ("\n Problems with  %s??", ficha);
				fflush(stdin); getchar();  exit(0);
				}

			printf ("\n FASTQ File   :  %s\n", ficha);
			strcpy(ficha_base, ficha);
			printf ("\n==================================================");
		
			}//---->> del último IF-ELSE
	}
 
// printf("\n ... press a key <CR> ");
// fflush(stdin); getchar();


//============>> We generate the name of the output files:

//---> Generating "_stats.txt" file:

 ff_aux_base[0]=0; //--> inicialized

 aux= strlen(ficha_base);  // printf("\n length of %s = %d", ficha_base, aux);
 i=0;
 do {
    ff_aux_base[i]= ficha_base[i]; i++;
 }
 while ((ficha_base[i-1] != PUNTO) && (i<aux));
 ff_aux_base[i-1]=0; //--> end of string, no dot	
	
 strcpy(ff_aux, ff_aux_base);	
 strcat(ff_aux, "_stats.txt");

 fest=fopen(ff_aux, "w");
 printf("\n Stats file   : %s", ff_aux);
 
//--------------------------------------
//---> Generating "_summary.txt" file:

 strcpy(ff_aux, ff_aux_base);	
 strcat(ff_aux, "_summary.txt");

 fout=fopen(ff_aux, "w");
 printf("\n Summary file : %s", ff_aux);
 
//--------------------------------------------------
//------>> Ficheros de salida con 0+1 o 2+3 BARCODES:
 
 strcpy(ff_aux, ff_aux_base);	
 strcat(ff_aux, "_1_0_bar.fq");
 
 printf("\n File with 1 or 0 barcode : %s", ff_aux);
  
 ff10=fopen(ff_aux,"a+");
//------- 
 strcpy(ff_aux, ff_aux_base);	
 strcat(ff_aux, "_3_2_bar.fq");
 
  printf("\n File with 3 or 2 barcode : %s", ff_aux);
 ff23=fopen(ff_aux,"a+");
 
 printf ("\n==================================================\n");
  
  
// fflush(stdin); getchar();
//********  Leo datos de BARCODES.TXT y guardo en arrays  *********

//====>> BARCODES A  ======================
 caracter= getc(fbar);
 if (caracter!=SOS)
   {printf("\n We started wrong: %c ", caracter); exit(0);}

 auxss=fgets(linea, SSS, fbar);
 printf("#%s", linea);

//---------------------------
 for (i=0; i<MAXBAR; i++){
	fscanf(fbar,"%s %s", etiAA[i], SSA[i]);
    printf("\n%2d %s %s", i, etiAA[i], SSA[i]);	 
    }
 printf("\n------------------------------------------- \n");	 
 
//====>> BARCODES B  ====================== 
 caracter= getc(fbar);
 if (caracter!=SOS){
 	 auxss=fgets(linea, SSS, fbar);
     caracter= getc(fbar);
	 }

 if (caracter!=SOS)
   {printf("\n Bad second set of barcodes: %c ", caracter); exit(0);}

 auxss=fgets(linea, SSS, fbar);
 printf("#%s", linea);
//---------------------------
 for (i=0; i<MAXBAR; i++){
	fscanf(fbar,"%s %s", etiBB[i], SSB[i]);
    printf("\n%2d %s %s", i, etiBB[i], SSB[i]);	 
    }
 printf("\n------------------------------------------- \n");

//====>> BARCODES C  ======================
 caracter= getc(fbar);
 if (caracter!=SOS){
 	 auxss=fgets(linea, SSS, fbar);
     caracter= getc(fbar);
	 }

 auxss=fgets(linea, SSS, fbar);
 printf("#%s", linea);
//---------------------------
 for (i=0; i<MAXBAR; i++){
	fscanf(fbar,"%s %s", etiCC[i], SSC[i]);
    printf("\n%2d %s %s", i, etiCC[i], SSC[i]);	 
    }
 printf("\n------------------------------------------- \n");

//====>> BARCODES D  ====================== 
  caracter= getc(fbar);
 if (caracter!=SOS){
 	 auxss=fgets(linea, SSS, fbar);
     caracter= getc(fbar);
	 }

 auxss=fgets(linea, SSS, fbar);
 printf("#%s", linea);
//---------------------------
 for (i=0; i<MAXBAR; i++){
	fscanf(fbar,"%s %s", etiDD[i], SSD[i]);
    printf("\n%2d %s %s", i, etiDD[i], SSD[i]);	 
    }
    
 printf("\n------------------------------------------- \n");
//  fflush(stdin); getchar();

//==============================================================
//=========    READ   THE  FASTQ FILE     
//==============================================================

 maxlen1=0; maxlen2=0;
 contador=0; uno_cero=0;
 lim=0;
 encontrados_A=0; encontrados_B=0; encontrados_C=0; encontrados_D=0;
 	
 while (!feof(ffastq)) {
	
	linea1[0]=0; linea2[0]=0, linea3[0]=0, linea4[0]=0;
	
	//--->> PRIMERA linea de fichero FASTQ
	
	auxss= fgets(linea1, SSS, ffastq);
	len1=strlen(linea1); linea1[len1-1]=0;
	if(len1>2){
		if (len1>maxlen1) maxlen1=len1;
		}
		
	//--->> SEGUNDA linea de fichero FASTQ
	auxss= fgets(linea2, SSS4-1, ffastq);
	len2=strlen(linea2);  linea2[len2-1]=0;
	if(len1>2){
		if (len2>maxlen2) maxlen2=len2;
		}

	//--->> TERCERA linea. La leo y paso de ella:
	auxss= fgets(linea3, SS, ffastq);
	len3= strlen(linea3);  linea3[len3-1]=0;
   
	//--->> CUARTA  linea: La leo y no la uso:
	auxss= fgets(linea4, SSS4-1, ffastq);
	len4 = strlen(linea4);  linea4[len4-1]=0;

	// printf("\n--------------------------------");
    //fflush(stdin); getchar();

    
    cuantos=0; cc_dir=0;
	if (len2>CAB){//====================>>  TENGO LINEA2  <<===========
	
		strcpy(nomff1, "ff_"); //--->  inicio de nombre de fichero
		
		cc_dir=0;
	    contador++; lim++;
		fprintf(fout, "%s\tDirect\t%5d", linea1, contador); 
	  
	  
		//-------------->>  BARCODE_D :  <<--------------------------
		
		if (lim <= LLL_VEO) //-->> por pantalla las primeras LLL_VEO reads
			printf("\n---------%6d - dir ----------------", contador);
			
	    tengo_bar=0;
	    i_barD=0; dd= -1;  min_dd= CAB_4; //--> UN NUMERO GRANDE.
		
		
	    for (i=0; i<MAXBAR; i++){
			
		   //--> Funcion de busqueda de BARCODES:
	       dd= busca(SSD[i], linea2, 0, CAB);  
			
		   if (dd>=0){
				if(dd<min_dd){
					min_dd=dd;
					i_barD= i;
					//printf("[%3d]-", min_dd);
					}
				}			   
		   }//---->> del for(i)	
		//printf("\n");	
					
			
		   if ((min_dd>=0) && (min_dd!=CAB_4)){
				// printf("\n-------> linea : %6d", contador);
				if (lim <= LLL_VEO){
				   printf("\n%d\t -->[BARCODE_D]\n", min_dd);
				   for (k=0; k<(CAB+BARSIZE); k++) printf("%c", linea2[k]);	
				   printf ("\n");
				
				   for (k=0; k<min_dd; k++) printf(" ");
				   printf ("%s\t", SSD[i_barD]);
				   printf (" --> %s\n", etiDD[i_barD]);
				   }
				fprintf(fout, "\t%s", etiDD[i_barD]);
				
				cc_dir++;
				encontrados_D++;
				// tengo_bar=1; //---> SALGO: he encontrado BARCODE-D
				// fflush(stdin); getchar();
				}


		if ((min_dd>=0) && (min_dd!=CAB_4))  
			  dd= min_dd;
		else  dd= -1;

		pos_D= dd;
		
		if (dd<0) fprintf(fout, "\t...\t%3d", dd);

		else fprintf(fout, "\t%3d", pos_D);
	
		if (dd<0) strcat(nomff1, "____");
		else{
			strcat(nomff1, etiDD[i_barD]); 
			strcat(nomff1, "_");
			}	
		
		//-------------->>  BARCODE_C :  
	    tengo_bar=0; i=0; dd= -1;  
		
		if (pos_D>=0) pos_bus= pos_D + 8;
		else          pos_bus= 0;
		
	    do{

			//-------->>  Trozo de linea2 con lo que compara el BARCODE:
			 kk=0;
			for(k=pos_bus; k < (pos_bus+BARSIZE); k++, kk++) {
				ss_com[kk]= linea2[k];
				}
			ss_com[BARSIZE]	= 0;  //--> final de string !!!
		
			
			jj=0;
			do{
				aux= strcmp(ss_com, SSC[jj]);
				if (aux==0) tengo_bar=1;
				jj++;
				
			  }
			while ((jj<MAXBAR) && (tengo_bar==0) );
			
		   //dd = pos_bus;
			indice_C = jj-1;
			
		   if (tengo_bar == 1){
			   
			   dd= pos_bus;
			   
				//printf("\n-------> linea : %6d", contador);
				if (lim <= LLL_VEO){
				   printf("\n%d\t --> [BARCODE_C]\n", dd);
				   for (k=0; k < (CAB+BARSIZE); k++) printf("%c", linea2[k]);
				   printf ("\n");
				
				   for (k=0; k<dd; k++) printf(" ");
				   printf ("%s\t", SSC[indice_C]);
				   printf (" --> %s\n", etiCC[indice_C]);
				   }
				fprintf(fout, "\t%s", etiCC[indice_C]);
				
				cc_dir++;
				encontrados_C++;
				// fflush(stdin); getchar();
				}
				
		   pos_bus++;
	       }
	    while ((tengo_bar==0) && (pos_bus<CAB));
		
		
		pos_C= dd;
		
		if (dd<0) fprintf(fout, "\t...\t%3d", dd);
		else fprintf(fout, "\t%3d", pos_C);
		// fprintf(fout, "\n");
				
		if (dd<0) strcat(nomff1, "___");
		else{
			strcat(nomff1, etiCC[indice_C]); 
			strcat(nomff1, "_");
			}		
		
	//-------------->>  BARCODE_B :  
	    tengo_bar=0; i=0; dd= -1;  
		
		if (pos_C>=0) pos_bus= pos_C + 8;
		else          pos_bus= 0;
		
	    do{
			//-------->>  Trozo de linea2 con lo que compara el BARCODE:
			 kk=0;
			for(k=pos_bus; k < (pos_bus+BARSIZE); k++, kk++) {
				ss_com[kk]= linea2[k];
				}
			
			ss_com[BARSIZE]	= 0;  //--> final de string !!!
			
			jj=0;
			do{
				aux= strcmp(ss_com, SSB[jj]);
				if (aux==0) tengo_bar=1;
				jj++;
			  }
			while ((jj<MAXBAR) && (tengo_bar==0) );
			
			indice_B = jj-1;
			
		   if (tengo_bar == 1){
			   
			   dd= pos_bus;
			   
				if (lim <= LLL_VEO){
				   printf("\n%d\t --> [BARCODE_B]\n", dd);
				   for (k=0; k< (CAB+BARSIZE); k++) printf("%c", linea2[k]);
				   printf ("\n");
				
				   for (k=0; k<dd; k++) printf(" ");
				   printf ("%s\t", SSB[indice_B]);
				   printf (" --> %s\n", etiBB[indice_B]);
				   }
				fprintf(fout, "\t%s", etiBB[indice_B]);
				
				cc_dir++;
				encontrados_B++;
				// fflush(stdin); getchar();
				}
				
		   pos_bus++;
	       }
	    while ((tengo_bar==0) && (pos_bus<CAB));
		
		
		pos_B= dd;
		
		if (dd<0) fprintf(fout, "\t...\t%3d", dd);
		else fprintf(fout, "\t%3d", pos_B);
		// fprintf(fout, "\n");
				
		if (dd<0) strcat(nomff1, "___");
		else{
			strcat(nomff1, etiBB[indice_B]); 
			strcat(nomff1, "_");
			}		
				
	
		//-------------->>  BARCODE_A :  
	    tengo_bar=0; i=0; dd= -1;  
		
		if (pos_B>=0) pos_bus= pos_B + 8;
		else          pos_bus= 0;
		
	    do{
			//-------->>  Trozo de linea2 con lo que compara el BARCODE:
			kk=0;
			for(k=pos_bus; k < (pos_bus+BARSIZE); k++, kk++) {
				ss_com[kk]= linea2[k];
				}
			
			ss_com[BARSIZE]	= 0;  //--> final de string !!!
			
			jj=0;
			do{
				aux= strcmp(ss_com, SSA[jj]);
				if (aux==0) tengo_bar=1;
				jj++;
			  }
			while ((jj<MAXBAR) && (tengo_bar==0) );
			
			indice_A = jj-1;
			
		   if (tengo_bar == 1){
			   
			   dd= pos_bus;
			   
				if (lim <= LLL_VEO){
				   printf("\n%d\t --> [BARCODE_A]\n", dd);
				   for (k=0; k< (CAB+BARSIZE); k++) printf("%c", linea2[k]);
				   printf ("\n");
				
				   for (k=0; k<dd; k++) printf(" ");
				   printf ("%s\t", SSA[indice_A]);
				   printf (" --> %s\n", etiAA[indice_A]);
				   }
				fprintf(fout, "\t%s", etiAA[indice_A]);
				
				cc_dir++;
				encontrados_A++;
				}
				
		   pos_bus++;
	       }
	    while ((tengo_bar==0) && (pos_bus<CAB));
		
		
		pos_A= dd;
		
		if (dd<0) fprintf(fout, "\t...\t%3d", dd);
		else fprintf(fout, "\t%3d", pos_A);
		// fprintf(fout, "\n");
				
		// if (dd<0) strcat(nomff1, "___");
		if (dd<0) strcat(nomff1, "__");
		else{
			strcat(nomff1, etiAA[indice_A]); 
			}		
				

		
					
		//----->> Posicion en DIR para limpiar las reads:			
		pos_max_dir=0;			
		if (cc_dir==4) pos_max_dir= pos_A;
		else {
			if (pos_D > pos_max_dir) pos_max_dir=pos_D;
			if (pos_C > pos_max_dir) pos_max_dir=pos_C;
			if (pos_B > pos_max_dir) pos_max_dir=pos_B;
			if (pos_A > pos_max_dir) pos_max_dir=pos_A;
			}
					
		//----> Copio valores:
		pos_D_d= pos_D; pos_C_d= pos_C; pos_B_d= pos_B; pos_A_d= pos_A;  	
		
		
		strcat(nomff1, ".fq");  //--->> siempre con *.fq
		
		// if (cc_dir==4) 	fflush(stdin); getchar();
				
	//===================================================================		
	//========>> Search in reverse complementary (INVERSA_COMPLEMENTARIA)
	
		
		//--------->> REVERSE COMPLEMENTARY:
		
  		k = inv_com(linea2_inv, linea2);
		strcpy(nomff2, "ff_"); //--->  inicio de nombre de fichero
			
		cc_ic= 0;
		fprintf(fout, "\tRev_com\t%5d", contador); 
	  
	  
	  
	  //-------------->>  BARCODE_D :  NUEVA FUNCION!!  
  if (lim <= LLL_VEO) //-->> por pantalla las primeras LLL_VEO reads
	printf("\n---------%6d - rc  ----------------", contador);
			
	tengo_bar=0;
	i_barD=0; dd= 0;  min_dd= CAB_4; //--> UN NUMERO GRANDE.
		
		
	    for (i=0; i<MAXBAR; i++){
			
		   //--> Funcion de busqueda de BARCODES:
	       dd = busca(SSD[i], linea2_inv, 0, CAB);  
			
		   if (dd>=0){
				if(dd<min_dd){
					min_dd=dd;
					i_barD= i;
					}
				}			   
		   }//---->> del for(i)	
			
			
		   if ((min_dd>=0) && (min_dd!=CAB_4)){
				// printf("\n-------> linea : %6d", contador);
				if (lim <= LLL_VEO){
				   printf("\n%d\t -->[BARCODE_D]\n", min_dd);
				   for (k=0; k<(CAB+BARSIZE); k++) printf("%c", linea2_inv[k]);	
				   printf ("\n");
				
				   for (k=0; k<min_dd; k++) printf(" ");
				   printf ("%s\t", SSD[i_barD]);
				   printf (" --> %s\n", etiDD[i_barD]);
				   }
				fprintf(fout, "\t%s", etiDD[i_barD]);
				
				cc_ic++;
				encontrados_D++;
				// tengo_bar=1; //---> SALGO: he encontrado BARCODE-D				
				// fflush(stdin); getchar();
				}


		if ((min_dd>=0) && (min_dd!=CAB_4))  
			  dd= min_dd;
		else  dd= -1;


		pos_D= dd;
		
		if (dd<0) fprintf(fout, "\t...\t%3d", dd);
		else fprintf(fout, "\t%3d", pos_D);
		
		if (dd<0) strcat(nomff2, "____");
				
		else{
			strcat(nomff2, etiDD[i_barD]); 
			strcat(nomff2, "_");
			}	

		
		//-------------->>  BARCODE_C :  
	    tengo_bar=0; i=0; dd= -1;  
		
		if (pos_D>=0) pos_bus= pos_D + 8;
		else          pos_bus= 0;
		
	    do{
			//-------->>  Trozo de linea2_inv con lo que compara el BARCODE:
			kk=0;
			for(k=pos_bus; k < (pos_bus+BARSIZE); k++, kk++) {
				ss_com[kk]= linea2_inv[k];
				}
			
			ss_com[BARSIZE]	= 0;  //--> final de string !!!
			
			jj=0;
			do{
				aux= strcmp(ss_com, SSC[jj]);
				if (aux==0) tengo_bar=1;
				jj++;
			  }
			while ((jj<MAXBAR) && (tengo_bar==0) );
			
			indice_C = jj-1;
			
		   if (tengo_bar == 1){
			   
			   dd= pos_bus;
			   
				if (lim <= LLL_VEO){
				   printf("\n%d\t --> [BARCODE_C]\n", dd);
				   for (k=0; k<(CAB+BARSIZE); k++) printf("%c", linea2_inv[k]);
				   printf ("\n");
				
				   for (k=0; k<dd; k++) printf(" ");
				   printf ("%s\t", SSC[indice_C]);
				   printf (" --> %s\n", etiCC[indice_C]);
				   }
				fprintf(fout, "\t%s", etiCC[indice_C]);
				
				cc_ic++;
				encontrados_C++;
				// fflush(stdin); getchar();
				}
				
		   pos_bus++;
	       }
	    while ((tengo_bar==0) && (pos_bus<CAB));
		
		
		pos_C= dd;
		
		if (dd<0) fprintf(fout, "\t...\t%3d", dd);
		else fprintf(fout, "\t%3d", pos_C);
		// fprintf(fout, "\n");
				
		if (dd<0) strcat(nomff2, "___");
		else{
			strcat(nomff2, etiCC[indice_C]); 
			strcat(nomff2, "_");
			}
			

 
	//-------------->>  BARCODE_B :  
	    tengo_bar=0; i=0; dd= -1;  
		
		if (pos_C>=0) pos_bus= pos_C + 8;
		else          pos_bus= 0;
		
	    do{
			//-------->>  Trozo de linea2_inv con lo que compara el BARCODE:
			kk=0;
			for(k=pos_bus; k < (pos_bus+BARSIZE); k++, kk++) {
				ss_com[kk]= linea2_inv[k];
				}
			
			ss_com[BARSIZE]	= 0;  //--> final de string !!!
			
			jj=0;
			do{
				aux= strcmp(ss_com, SSB[jj]);
				if (aux==0) tengo_bar=1;
				jj++;
			  }
			while ((jj<MAXBAR) && (tengo_bar==0) );
			
			indice_B = jj-1;
			
		   if (tengo_bar == 1){
			   
			   dd= pos_bus;
			   
				if (lim <= LLL_VEO){
				   printf("\n%d\t --> [BARCODE_B]\n", dd);
				   for (k=0; k<(CAB+BARSIZE); k++) printf("%c", linea2_inv[k]);
				   printf ("\n");
				
				   for (k=0; k<dd; k++) printf(" ");
				   printf ("%s\t", SSB[indice_B]);
				   printf (" --> %s\n", etiBB[indice_B]);
				   }
				fprintf(fout, "\t%s", etiBB[indice_B]);
				
				cc_ic++;
				encontrados_B++;
				// fflush(stdin); getchar();
				}
				
		   pos_bus++;
	       }
	    while ((tengo_bar==0) && (pos_bus<CAB));
		
		
		pos_B= dd;
		
		if (dd<0) fprintf(fout, "\t...\t%3d", dd);
		else fprintf(fout, "\t%3d", pos_B);
				
		if (dd<0) strcat(nomff2, "___");
		else{
			strcat(nomff2, etiBB[indice_B]); 
			strcat(nomff2, "_");
			}		
		

		//-------------->>  BARCODE_A :  
	    tengo_bar=0; i=0; dd= -1;  
		
		if (pos_B>=0) pos_bus= pos_B + 8;
		else          pos_bus= 0;
		
	    do{
			//-------->>  Trozo de linea2_inv con lo que compara el BARCODE:
			kk=0;
			for(k=pos_bus; k < (pos_bus+BARSIZE); k++, kk++) {
				ss_com[kk]= linea2_inv[k];
				}
			
			ss_com[BARSIZE]	= 0;  //--> final de string !!!
			
			jj=0;
			do{
				aux= strcmp(ss_com, SSA[jj]);
				if (aux==0) tengo_bar=1;
				jj++;
			  }
			while ((jj<MAXBAR) && (tengo_bar==0) );
			
			indice_A = jj-1;
			
		   if (tengo_bar == 1){
			   
			   dd= pos_bus;
			   
				if (lim <= LLL_VEO){
				   printf("\n%d\t --> [BARCODE_A]\n", dd);
				   for (k=0; k<(CAB+BARSIZE); k++) printf("%c", linea2_inv[k]);
				   printf ("\n");
				
				   for (k=0; k<dd; k++) printf(" ");
				   printf ("%s\t", SSA[indice_A]);
				   printf (" --> %s\n", etiAA[indice_A]);
				   }
				fprintf(fout, "\t%s", etiAA[indice_A]);
				
				cc_ic++;
				encontrados_A++;
				// fflush(stdin); getchar();
				}
				
		   pos_bus++;
	       }
	    while ((tengo_bar==0) && (pos_bus<CAB));
		
		
		pos_A= dd;
		
		if (dd<0) fprintf(fout, "\t...\t%3d\n", dd);
		else fprintf(fout, "\t%3d\n", pos_A);
		// fprintf(fout, "\n");
				
		// if (dd<0) strcat(nomff2, "___");
		if (dd<0) strcat(nomff2, "__");
		else{
			strcat(nomff2, etiAA[indice_A]); 
			//strcat(nomff2, "_"); //--> Sobra al final
			}		
											
								
		//----->> Posicion en Reverse_Compl para limpiar las reads:			
		pos_max_ic=0;			
		if (cc_ic==4) pos_max_ic= pos_A;
		else {
			if (pos_D > pos_max_ic) pos_max_ic=pos_D;
			if (pos_C > pos_max_ic) pos_max_ic=pos_C;
			if (pos_B > pos_max_ic) pos_max_ic=pos_B;
			if (pos_A > pos_max_ic) pos_max_ic=pos_A;
			}						
								
								
		//------------------------------------------------------------------
		//----------->>  Resultados a FICHERO FASTQ:	
		//----------->>  Ahora tambien con 0+1 o 2+3 BARCODES:
		
		strcat(nomff2, ".fq");  //--->> siempre con *.fq
		
		
		//----------->>  4 BARCODES presentes. Resultados a FICHERO FASTQ. 
		//----------->>  Con LIMPIEZA de las reads, eliminamos BARCODES:
		
		if ( (cc_dir==4) && (cc_ic==4) ) {//------>> genero fichero de salida:
		
			//-------->> Calculo valores de corte de la READ:
			corte_dir= pos_max_dir + BARSIZE + LINKSIZE;
		    corte_ic = pos_max_ic  + BARSIZE + LINKSIZE;
		
			//-------->> Limpio read:
			
			j=0;
			for(i=corte_dir; i<(len2-corte_ic); i++){ 
			    linea2_clean[j]= linea2[i];
				linea4_clean[j]= linea4[i];
				j++;
				}
			//---> fin de linea:	
			linea2_clean[j]= 0; linea4_clean[j]= 0;

			//-------->>  Abro, escribo y cierro!!
			strcpy(lin_dir, FF_4BAR); strcat(lin_dir, nomff1);
			facu= fopen(lin_dir, "a+");
 			fprintf(facu,"%s\n%s\n%s\n%s\n", linea1, linea2_clean, linea3, linea4_clean);	
	    	fclose(facu);	
	    	encontrados++; completos++;
			
	  		}		
			
		else if ( (cc_dir==4) && (cc_ic<4) ) {
			
			//-------->> Calculo valores de corte de la READ:
			corte_dir= pos_max_dir + BARSIZE + LINKSIZE;
			corte_ic = 0;
			
			if ((cc_ic==3) && (pos_max_ic < (CAB- (BARSIZE + LINKSIZE))) )	
				corte_ic= pos_max_ic + BARSIZE + LINKSIZE;
			
			if ((cc_ic==2) && (pos_max_ic < (CAB- 2*(BARSIZE + LINKSIZE))) )	
				corte_ic= pos_max_ic + BARSIZE + LINKSIZE;		
			
			//-------->> Limpio read:
			j=0;
			for(i=corte_dir; i<(len2-corte_ic); i++){ 
			    linea2_clean[j]= linea2[i];
				linea4_clean[j]= linea4[i];
				j++;
				}
			//---> fin de linea:	
			linea2_clean[j]= 0; linea4_clean[j]= 0;

			//-------->>  Abro, escribo y cierro!!
			strcpy(lin_dir, FF_4BAR); strcat(lin_dir, nomff1);
			facu= fopen(lin_dir, "a+");
 			fprintf(facu,"%s\n%s\n%s\n%s\n", linea1, linea2_clean, linea3, linea4_clean);	
	    	fclose(facu);	
	    	encontrados++; completos++;
		    }
			
			else if ( (cc_dir<4) && (cc_ic==4) ) {
				
					//-------->> Calculo valores de corte de la READ:
					corte_dir= 0;
					corte_ic = pos_max_ic + BARSIZE + LINKSIZE;
			
					if ((cc_dir==3) && (pos_max_dir < (CAB- (BARSIZE + LINKSIZE))) )	
						corte_dir= pos_max_dir + BARSIZE + LINKSIZE;
					
					if ((cc_dir==2) && (pos_max_dir < (CAB- 2*(BARSIZE + LINKSIZE))) )	
						corte_dir= pos_max_dir + BARSIZE + LINKSIZE;
					
					//-------->> Limpio read:
					j=0;
					for(i=corte_dir; i<(len2-corte_ic); i++){ 
						linea2_clean[j]= linea2[i];
						linea4_clean[j]= linea4[i];
						j++;
						}				
					//---> fin de linea:	
					linea2_clean[j]= 0; linea4_clean[j]= 0;
				
					//-------->>  Abro, escribo y cierro!!
					strcpy(lin_dir, FF_4BAR); strcat(lin_dir, nomff2);
					facu= fopen(lin_dir, "a+");
					fprintf(facu,"%s\n%s\n%s\n%s\n", linea1, linea2_clean, linea3, linea4_clean);	
					fclose(facu);	
					encontrados++; completos++;				
				
					}


			
		//---->> si tengo mas de 1 BARCODE --> IMPRIMO EN FICHERO SIN LIMPIAR:	
		else if ((cc_dir >= cc_ic) && ( cc_dir>1) && (cc_dir<4) ){
				// strcat(nomff1, ".fq");
				strcpy(lin_dir, FF_23BAR );    strcat(lin_dir, nomff1);
				
				facu= fopen(lin_dir, "a+");
	
				fprintf(facu,"%s\n%s\n%s\n%s\n", linea1, linea2, linea3, linea4);	
				fprintf(ff23,"%s\n%s\n%s\n%s\n", linea1, linea2, linea3, linea4);	
				fclose(facu); encontrados++;
			
				}
				//----> Nueva condicion: 26-IX-2024
				
				else if ((cc_ic > cc_dir) && (cc_ic>1) && (cc_ic<4)   ){
					   // strcat(nomff2, ".fq");
					   strcpy(lin_dir, FF_23BAR);  strcat(lin_dir, nomff2);
					   facu= fopen(lin_dir, "a+");
	
					   fprintf(facu,"%s\n%s\n%s\n%s\n", linea1, linea2, linea3, linea4);	
					   fprintf(ff23,"%s\n%s\n%s\n%s\n", linea1, linea2, linea3, linea4);	
					   fclose(facu); encontrados++;
					}
			
		//---->> para cuando tengo 1 o 0 BARCODES: 	
		if ((cc_dir < 2 ) && ( cc_ic<2 )){

				fprintf(ff10,"%s\n%s\n%s\n%s\n", linea1, linea2, linea3, linea4);	
				uno_cero++;
			
				}	
			
		
		//----------->>  Estadisticas:				
		if (cc_dir >= cc_ic) cuantos =  cc_dir;
		else  cuantos = cc_ic;
				
		}//=========>> DEL IF de tengo LINEA2 > 2
		
		


   if ((cc_ic == 4) || (cc_dir==4)){
		if ((cc_ic > 1) || (cc_dir > 1)){
			
			if 	(cc_dir == 4)
			    fprintf(fest,"%s\t%d_dir\t%s\t%4d\t(#%d)\t", linea1, contador, nomff1, pos_max_dir, cc_dir);
			else  
				fprintf(fest,"%s\t%d_dir\t%s\t%4d\t(#%d)\t", linea1, contador, No_file, pos_max_dir, cc_dir);
			    
			
			//---> Imprimo nt del LINK para DIR:
			if 	(cc_dir == 4)
				fprintf(fest,"%d\t%d\t%d\t", pos_C_d-(pos_D_d+BARSIZE), pos_B_d-(pos_C_d+BARSIZE), pos_A_d-(pos_B_d+BARSIZE));
			else fprintf(fest,"\t\t\t");			
			
			if (cc_ic == 4)
			    fprintf(fest,"%d_rc \t%s\t%4d\t(#%d)\t", contador, nomff2, pos_max_ic,  cc_ic);
			else 
			    fprintf(fest,"%d_rc \t%s\t%4d\t(#%d)\t", contador, No_file, pos_max_ic,  cc_ic);
			
			//---> Imprimo nt del LINK para IC:
			if 	(cc_ic == 4)
				fprintf(fest,"%d\t%d\t%d\n", pos_C-(pos_D+BARSIZE), pos_B-(pos_C+BARSIZE), pos_A-(pos_B+BARSIZE));
			else fprintf(fest,"\t\t\n");
			}
		}




 }//-------------------->> extern WHILE  and end of loop
  
 
 printf("\n-------------------------------------------");
 printf("\nFound BARCODES-A: %4d", encontrados_A);
 printf("\nFound BARCODES-B: %4d", encontrados_B);
 printf("\nFound BARCODES-C: %4d", encontrados_C);
 printf("\nFound BARCODES-D: %4d", encontrados_D);
 printf("\n-------------------------------------------");
 printf("\nNumber of fastq lines analyzed = %7d", contador);
 faux1= (float) uno_cero*100.0; faux1=(faux1/contador);
 printf("\nReads with 0 or 1 BARCODES     = %7d\t%6.2f%%", uno_cero, faux1);
 // printf("\nReads con 2 3 o 4 BARCODES  = %7d", encontrados); 
  faux1= (float) (encontrados-completos)*100.0; faux1=(faux1/contador);
 printf("\nReads with 2 o 3 BARCODES      = %7d\t%6.2f%%", encontrados-completos, faux1);
 faux1= (float) (completos)*100.0; faux1=(faux1/contador);
 printf("\nReads with 4 BARCODES          = %7d\t%6.2f%%", completos, faux1);
 printf("\n-------------------------------------------");
 printf("\nNucleotides analyzed in header : %d", HEAD);
 printf("\nmaxlen1: %8d", maxlen1);
 printf("\nmaxlen2: %8d", maxlen2);
 printf("\n-------------------------------------------");


 fprintf(fout, "\n----------------------------------------------");
 fprintf(fout, "\nFound BARCODES-A: %d", encontrados_A);
 fprintf(fout, "\nFound BARCODES-B: %d", encontrados_B);
 fprintf(fout, "\nFound BARCODES-C: %d", encontrados_C);
 fprintf(fout, "\nFound BARCODES-D: %d", encontrados_D);
 fprintf(fout, "\n----------------------------------------------");
 fprintf(fout, "\nNucleotides analyzed in header = %d\n", HEAD);
 fprintf(fout, "\nFastq lines analyzed       = %7d", contador);
 faux1= (float) uno_cero*100.0; faux1=(faux1/contador);
 fprintf(fout, "\nReads with 0 or 1 BARCODES = %7d\t%6.2f%%", uno_cero, faux1);
 faux1= (float) (encontrados-completos)*100.0; faux1=(faux1/contador);
 fprintf(fout, "\nReads with 2 or 3 BARCODES = %7d\t%6.2f%%", encontrados-completos, faux1);
 faux1= (float) (completos)*100.0; faux1=(faux1/contador); 
 fprintf(fout, "\nReads with 4 BARCODES      = %7d\t%6.2f%%", completos, faux1);
 fprintf(fout, "\n----------------------------------------------");

 fclose(fbar); fclose(ffastq); fclose(fout); fclose(fest);
 fclose(ff10); fclose(ff23);
 
// printf("\n\n\t ... eso, eso es todo amigos.\n ");
 printf("\n\n\t ... That's all Folks.\n ");
// fflush(stdin); getchar();
 
 return (MAXBAR);
 }


//*************************************************************
//*************************************************************
//**************    DEFINED FUNCTIONS    **********************

int busca(char *bar, char *sss, int pos, int CAB)
{
 int  ii, jj, len1, len2, dd;
 int  salgo=0, igual=0, limit=0;
 
 
 len1=strlen(bar); len2=strlen(sss);
 
 if (len2>CAB) limit= CAB;
 else {
	limit= len2;
	printf("\nlimit=  %d", limit);
	fflush(stdin); getchar();
	}
  
 dd=-1;   igual=0; salgo=0;
 
 do{
	salgo=0; jj=pos; igual=0;
	for (ii=0; ii<BARSIZE; ii++, jj++){		
	   //	printf("\n%c  ", bar[ii] );
	   //	printf("%c", sss[jj] );
		
		if(bar[ii]==sss[jj]) igual++;
	}
	if (igual==BARSIZE ) salgo=1; //--> tengo coincidencia
	 
	pos++; 
		
 }while ((salgo==0) && (pos<limit));
 
 // fflush(stdin); getchar();
 
 //--> devuelvo posicion de inicio de coincidencia.
 if (igual==BARSIZE) return(pos-1);
 else            return(-1);
 
}

int inv_com(char *ss1, char *ss2)
{
 int i, j, ll, k;
 
 ll = strlen(ss2); 
 
 j=0; 
 for (i=(ll-1) ; i>=0; i-- ){
	if ( (ss2[i]=='A') || (ss2[i]=='a')) ss1[j]='T';
	else if ((ss2[i]=='T') || (ss2[i]=='t')) ss1[j]='A';
	     else if ((ss2[i]=='C') || (ss2[i]=='c')) ss1[j]='G';
		      else if ((ss2[i]=='G' || ss2[i]=='g')) ss1[j]='C';
	j++;
	}	 
 ss1[j]=0; //---> Fin de linea
 
 return(ll);
}
//*****************************************************************
//*****************************************************************


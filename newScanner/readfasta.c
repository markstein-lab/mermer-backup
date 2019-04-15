//Copyright (c) World Internet Productions, LLC, 1999
//Copyright (c) in silico Labs, LLC 2006, 2009, 2018

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#define wordsize sizeof(unsigned int)
#define YES  1
#define NO  0
#define number_of_mrnas(i) stack[i].numberOfMrnas 
#define number_of_cds(i)   stack[i].numberOfCds
#define max_Exons(i)       stack[i].maxExons
#define Mrna(i,j,k)        stack[i].mrnaData[1000*(j) + (k)]   

   FILE *dnadata;
   FILE *exceptionfile;
//   FILE *allnames;
   FILE *data;
   FILE *master;
   FILE *fastaname;
   FILE *stream;
//  FILE *xmrnastart;
//   FILE *xmrnadata; 

   void getLine (char* line, FILE *stream);
   void add(char);
   void Pop();
   void readExons();
   void closeMrna(int);
   char* readLimits(char* ptr, unsigned int*, unsigned int*);
   char* skipwhitesp(char * string);
   char* searchstr ( char *,  char *);
   char *ptr, *ptr2, *ptr3, *ptr4;
   char genome;
   char genomeM;
   char geneName[500];
   char first_name[500];
   char line[500];
   unsigned int start, stop; 
   unsigned int fileSize;
   unsigned int counter, oldcounter;
   unsigned int exceptions;
   unsigned int mrnaCounter;
   int mrnaStartCounter;
   
   int  topOfStack;
   int  DNACount, exceptionBits , excount, mrnaC;
   int  totalNames;
   int  number_of_genes;
   int  numberOfExons=0;
   int  tempExons[1000];
   
   int  eof;
   int  mostExonsPerMrna;
   int  mostSplicesPerGene;
   int  deepestStack;
   int  annot_id;
   char temp_direction;
   char geneWithMostExons[100];
   char geneWithMostMrnas[100];
   char geneMostDeeplyEmbedded[100];

   struct geneData{
      int   annot_id;
      unsigned int start, stop;
      char  direction;
      char  *geneName;
      char  **synonyms;
      int   totalNames; 
      int   *mrnaData;
      int   numberOfMrnas;
      int   numberOfCds;
      int   maxExons;
      int   maxNames;
      int   mrnaOpen;
   };
   
   struct geneData stack[100];

   struct exonData{
      int   geneNumber;
      unsigned int *mrnaInfo;
      int   maxExons;
      int   numberOfMrnas;      
   };
 
    struct exonData *geneInfo;
    int   sizeGeneInfo;
    int   genesInContig;

      	  
int main( void )
{
   FILE *names;
   char line2[500], head[500];
   char annot_name[20];
   char junk[100];
   char direction;
   char chrm_name[500];
   char full_name[500];
   char file_name[500];
   char input_file_name[500];
   char rightbr;
   unsigned int origin;
   unsigned int size;
   unsigned int contigCounter;
   int i, z, x, y, j, k, m;
   int flag, smallestMrna;
   int needToReadLine;
    
   strcpy(geneWithMostExons, "");
   strcpy(geneWithMostMrnas, "");
   strcpy(geneMostDeeplyEmbedded, "");
   mostExonsPerMrna = 0;
   mostSplicesPerGene = 0;
   deepestStack = 0;
   exceptions = 0;

   geneInfo = malloc(sizeof(struct exonData) * 16384);
   sizeGeneInfo = 16384;
        
   stream = NULL;
   fastaname = fopen("fastaname.txt", "r");
   dnadata = fopen("genome.txt", "wb");
   master = fopen("master.txt", "w");
   exceptionfile = fopen("exceptions.txt", "w");
//   allnames = fopen("xallnames.txt", "w");

   if (fastaname == NULL) {
      printf(
      "file fastaname.txt, containing the name of the fasta genome,is missing!\n");
      exit(101);
   }

// Initialize
   fscanf(fastaname, "%s", input_file_name);
   stream = fopen(input_file_name, "r");
   if (stream == NULL) {
      printf("file %s failed to open; quitting run\n", input_file_name);
      exit(101);
   }
   annot_id = 0;
   origin = 0;
   number_of_genes = 0;
   excount = 0;
   counter = 0;
   oldcounter = 0;
   mrnaC = 0;
   DNACount = 0;   //number of bits in    genome char. buffer
   exceptionBits = 0;    //number of bits in exception word buffer
   genome = 0; //genome buffer
   genomeM = 0; 
   fileSize = 0;  //bytes written to genome.txt
   totalNames = 0;
   topOfStack = 0;
   mrnaStartCounter = 0;
   
   names = fopen("xcontigs.txt", "w");
//   data = fopen("xchrdata.txt", "w");
//   xmrnadata = fopen("xmrnadata.txt", "w");
//   xmrnastart = fopen("xmrnastart.txt", "w"); 
   line[0] = 0;    //initialize line to NULL string
   getLine(line, stream);
   while (1) {
      sscanf(line, "%c%s", &rightbr, line);
      if (strlen(line) == 0) break;  //no more chromosomes
      if (rightbr != '>') {
         printf("chromosome line does not start with a >\n");
         exit (0);
      }
      fprintf(master, "%s ", line);
      strcpy(chrm_name, line);
      
      printf ("Chromosome %s starts at position %u\n", chrm_name, origin);
   
readChromosome:
         //process DNA string;
         //start by inserting an N into the genome
         add('N');
         contigCounter = 1; //gbk numbering starts at position 1.
         while (1) {
            getLine(line, stream);
            eof = feof(stream);
            if (line) for (i = 0; line[i]; i++) {
               if (line[i] != ' ') break;
            }
            if (line[i] == '>' || eof )
            break;  //found end of data for annot. file
            for (i = 0; line[i]; i++) {
               if (line[i] >= 'A') {
                  add(line[i]);
                  contigCounter++;  
               }
            }
         }
         
         printf(
              "Chrm %s has length %u, origin = %u\n", 
   			chrm_name, counter - oldcounter, oldcounter+1 );
         fprintf(names, "%s %s %u %u %u\n", 
          input_file_name, chrm_name, annot_id, oldcounter, counter-oldcounter);
         oldcounter = counter;
         annot_id++;
         if (eof) {
            printf("Found end of fasta file\n");
            break;  //essentially, goto finish
         }
 
      }		//closes the while loop
      
finish:
     add ('N');    //end chromosome with an 'N'


   //Fill out exception word with N's
   while(exceptionBits) (add ('N'));
   add(0); //make sure dnastring looks like a string, i.e. ends in a 0.
   fclose(dnadata);
   fclose(names);
   fclose(master);
   fclose (exceptionfile);
   
   master = fopen("datasize.txt", "w");
   fprintf(master, "%u %d %d %d %d\n", 
           counter, number_of_genes, annot_id, excount, totalNames);
   fprintf(master, "%u %u\n", mrnaCounter, mrnaStartCounter);
   printf("origin = %u, counter = %u, number of genes = %d, number of annotation files = %d\n",
      origin, counter, number_of_genes, annot_id);
   printf("number of exception words = %d, number of gene names = %d\n", 
      excount, totalNames);
   fclose(master);
   
   exit(100);
}

void getLine (char* line, FILE *stream) {
   int i;
   for (i = 0; i < 500; i++) {
      line[i] = getc(stream);
      if (line[i] == '\n' || line[i] == EOF || line[i] == 0) break;
   }
   line[i] = 0;
}

void add(char x) {
   int y;
   y = x & 0xDF;   //make upper case
   switch (y) {
   case 'A':
      genome <<=2;
      exceptions <<=1;
      break;
   case 'C':
      genome = 1 + (genome <<2);
      exceptions <<=1;
      break;
   case 'G':
      genome = 2 + (genome << 2);
      exceptions <<= 1;
      break;
   case 'T':
      genome = 3 + (genome << 2);
      exceptions <<= 1;
      break;
   default:
      genome <<= 2;   //encode invalid characters like 'A's
      exceptions = 1 + (exceptions << 1);  //insert a bit into exception vector
      break;
   }
   DNACount++;
   if (DNACount ==4) {   //when a genome nibble is full
      fputc(genome, dnadata);      //write it out
      fileSize++;   //add to size of dnadata file
      DNACount = 0;
      genome = 0;
   }
   exceptionBits++;
   if (exceptionBits == 32) {  //a potential word of exceptions
      if (exceptions) { //if there are exception bits
         fprintf(exceptionfile, "%d %8.8x\n", counter/32, exceptions);
         excount++; //count number of exception words
      }
      exceptionBits = 0;
      exceptions = 0;
   }
   counter++;
}





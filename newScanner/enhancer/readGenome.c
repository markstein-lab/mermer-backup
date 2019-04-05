//Copyright (c) in silico Labs, LLC 2008, 2009, 2015

#include <stdlib.h>
#include <stdio.h>
#include "timer.h"
#include "scanner.h"
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/fcntl.h>

int getpagesize();
int numberOfChromosomes;
char *chrName[100];
int readContigs(FILE *);
void readgenenames(FILE *, int);

void readGenome() {
   FILE* qnames;
   FILE* master;
   FILE* mrnastuff;
   int i, j, k, m;
   unsigned char c;
   unsigned char *DNAs;
   /*unsigned*/ char line[5000];
   /*unsigned*/ char temp[500];
   int pagesize;
   int filedes;
   struct stat sbuf;
   char genomeName[100];
 
#if (scanWidth < 5)
   unsigned char word;
#else
   unsigned short int word;
#endif
   char message[100];
   
// This first section reads in the genome itself  
//   qnames = fopen(fullName(pathName, "genome.txt"), "rb");
//   printf("File name is %s\n", fullName(pathName, "genome.txt"));
   if (utr3 | utr5 | orf | cds) strcpy(genomeName, "xmrnadata.txt");
   else strcpy(genomeName, "genome.txt");
   filedes = open(fullName(pathName, genomeName), O_RDONLY);
   pagesize = getpagesize();
   if (filedes == -1) {
      printf("Genome does not exist -- quitting\n");
      exit(0);
   }
   stat(fullName(pathName, genomeName), &sbuf);
   DNAstring = mmap(/*(caddr_t)*/0, sbuf.st_size + pagesize, PROT_READ, MAP_SHARED, filedes, 0);
   if (DNAstring == (unsigned char*) -1){
      both("File mapping did not work!\n");
      exit(1);
   }
   goto genomeHasBeenRead;

#if (scanWidth == 1 || scanWidth == 2 || scanWidth == 4 )
   DNAs = (unsigned char *)DNAstring;
   for (i = 0; i < (DNALen/4 + 8); i++) {
      DNAs[i] = fgetc(qnames);
      if(feof(qnames)) break;
   }
   
   counter = i;
#ifdef DEBUG
   sprintf(message, "compact genome took %u bytes (%u nucleotides)\n",
      counter, DNALen);
   both(message);
#endif
   fclose(qnames);
#else
   j = 0;
   k = 0;
   word = 0;
   for (i = 0; i < (DNALen/4 + 8); i++) {
      c = fgetc(qnames);
      if (j + 8 <= 2*scanWidth) {
         word = (word << 8) + c;
         j = j + 8; 
         if (j == 2*scanWidth) {
            DNAstring[k++] = word;
            j = 0;
            word = 0;
         }
      }
      else if (j + 8 > 2*scanWidth) {
         word <<= 2*scanWidth - j;
         word += (c>>(8 - 2*scanWidth + j));
         DNAstring[k++] = word;
         j = 8 - 2*scanWidth + j;
         word = c & ((1<<j) - 1);
         if (j == 2*scanWidth) {
            DNAstring[k++] = word;
            word = 0;
            j = 0;
         }
      }         
   }
   counter = k;
#ifdef DEBUG
   sprintf(message, "i = %d, k = %d, size of word = %d\n", i, k, sizeof(word));
   both (message);
   sprintf(message, "compact genome took %d bytes (%d nucleotide)\n",
      counter*sizeof(word), DNALen);
   both (message);
#endif
   fclose (qnames);
#endif

//The genome has now been read in.
//Next, we read in supporting data, that describes the chromosomes and
//tells where the genes are on the chromosomes

genomeHasBeenRead:
   master = fopen(fullName(pathName, "master.txt"), "r");
   numberOfChromosomes = 0;
   while(1) {
      line[0] = 0;                  //initialize line to empty string
      fscanf(master, "%s", line);   //read in a chromosome name
      if (strlen(line)) {           //test if a name was found
         chrName[numberOfChromosomes] = malloc(strlen(line) + 1);
         strcpy(chrName[numberOfChromosomes], line); //save chr. name
         numberOfChromosomes++;     //increase number of Chromosomes
      }
      else {
         fclose(master);
         break;
      }
   }


//next, we read in the info about chromosome arms and contigs
   arm = malloc((ncontigs+2) * sizeof(int*));
   if (arm == NULL) {
      both("Insufficient space for chromosome table\n");
      exit (0);
   }

   annotName = malloc ((ncontigs+2)* sizeof(char*));
   if (annotName == NULL) {
      both("Insufficient space for contigs table\n");
      exit(0);
   }

   chromoName = malloc ((ncontigs+2)* sizeof(char*));
   if (chromoName == NULL) {
      both ("Insufficient space for Chromosome Name Table\n");
      exit (0);
   }

   extable = malloc ((extableSize + 1) * 2 * sizeof (int));
   if (extable == NULL) {
      both("Not enough space for exception table\n");
      exit(0);
   } 

  
   qnames = fopen(fullName(pathName, "xcontigs.txt"), "r");
   k = 0;
   k = readContigs(qnames); 
   if (k != ncontigs) {
      printf(message,
            "wrong number of contings, wanted %d, found %d\n", ncontigs, k);
      both(message);
   }
   fclose(qnames);

   arm[k] = DNALen;                //end of genome
   arm[k+1] = 0;
   chromoName[k] = malloc(20);
   strcpy(chromoName[k], "end_of_genome");
   chromoName[k+1] = chromoName[k];
   annotName[k] = chromoName[k];
   annotName[k+1] = annotName[k];
   qnames = fopen(fullName(pathName, "exceptions.txt"), "r");
   for (i = 0; i < extableSize; i++) {
      getLine(line, qnames);
      sscanf(line, "%d %8x", &extable[2*i], &extable[2*i+1]);
      if (feof(qnames)) {
         printf("Exception table too short\n");
         exit(0);
      }
   }
   extable[2*i] = (DNALen/32) + 1;    //dummy last entry beyond the genome
   extable[2*i + 1] = 0xFFFFFFFF;
   fclose (qnames);



}

//Copyright (c) in silico Labs, LLC 2008, 2009, 2010, 2012, 2018

#include <stdlib.h>
#include <stdio.h>
#include "scanner.h"
#include "timer.h"
#include <float.h>
#include <string.h>

/***/
#define TIMING
/***/
typedef enum { false = 0, true = !false } bool;

void makestring(int, char*);

LOOKUPWORD* makeTables(motifDes *);

int printBoth;
int printNone;
unsigned long int DNALen;
unsigned long int chrStart(int);
unsigned int pseudoLen;
unsigned int mrnaCount;
unsigned int mapSize;
unsigned int ncontigs;
unsigned int extableSize;
unsigned int synCount;
unsigned int counter;
int numberOfTables;
unsigned long int *av;
unsigned long int *avbase;
motifDes *pv;
motifDes *pvbase;
int numberOfMatches;
LOOKUPWORD *matches;
int numberOfMotifs;
motifDes *motifs;
int *geneNumber;
unsigned char **geneName;
char** annotName;
char** chromoName;
char* chrLimit;
unsigned long int* arm;
unsigned int* extable;
ann_index *map;
mrnaID *mrnaInfo;
FILE *save;
FILE *tbl;
FILE *indexFile;
FILE *allnames;
int minGroupSize;
int windowSize;
char boolExp[5000];
char* booleanSyntax(char message[]);
void reportClusters();
void printGeneList();
void readGenome();
void makerevcompl(char *, char *, int);
void chrLengths();
char genomeName[100];
char pathName[100] = "../../";
char altPathName[] = "./";

int fromTerminal;
alist terminalInput;
int *summaryData;
int numberOfGenes;

int do_the_search(unsigned long int, unsigned long int);
int fixinput(char []);
int checkMinNum(char *line, int *minGroupSize); 
int checkBoolExpr(char *boolExp);
int checkWindowSize(char *line, int *windowSize);
void hsort2(unsigned long int[], motifDes[], unsigned long int);
void outputSummary();
int limited, utr3, utr5, intron, exon, cds, orf;
char fileName[50];
char fileName2[50];

int chrCont(char *);

extern char *colorName[26];

int mainAux (int argc, char *argv[]) {
 
   unsigned char **strings;
   int repetitions;
   int motifCount;
   char patternID;
   void uniquename(char *result);

   LOOKUPWORD *tables;
   TIMEVAL time1, time2;
   unsigned int i, j, k, errcount;
   int sum;
   unsigned int genomeLength;
   char str[] = "        ";
   char result[100];
   char output[500];
   char line[5001];
   char rev[501];

   int numberToGenerate;
   int ltdSearch;

   FILE *qnames;
   unsigned char c;
   char lookup[] = "ACGT";
   char nullstr[] = "   ";
   char filebase[50];
   char leftArrowName[50];
   char rightArrowName[50];
   char *ptr;

   int *geneList = NULL;
   int geneDistance = 0;
   int relocation;
  
   unsigned long int first, last;


normalStart:
#ifdef DEBUG
#ifdef WIN32
   printf("Compiled for Win32\n");
#endif
#endif
   result[0] = 0;
   printBoth = 1;
   printNone = 0;

   if (argc < 0) {
      result[0] = 0;
      fromTerminal = 0;
      
   }
   else {
      fromTerminal = 1;
      //Open the file for posting results
      for (i = 1; i < argc; i++) { 
         strcpy(result, argv[i]);   //read alternative result file name
         break;
      }
   }
   if (result[0]==0) {
      uniquename(result);
   }
   ltdSearch = 0;
   fileName[0] = 0;
   save = fopen("resultfiles", "r");
   if (save) fclose(save);
   else {
      mkdir("resultfiles", 0777);
   }
   strcat(fileName, "resultfiles/");
   strcat(fileName, result);
   strcpy(fileName2, fileName);
   strcat(fileName, ".txt");
   save = fopen(fileName, "w");
   strcat(fileName2, ".gff");
   tbl = fopen(fileName2, "w"); 
   
   sprintf(output, 
           "Welcome to scanner 0.4 (beta) last updated 07/26/2018\n");
   both (output);
   sprintf(output, "A log of this run will be found in file %s\n", result);
   both(output);

   
   utr3 = utr5 = intron = exon = orf = cds = 0;
   if (fromTerminal == 0) {
      sprintf(genomeName, "%s", sassoc("ON", terminalInput));
      sprintf(output, "Genome selected is %s\n", genomeName);
      both(output);
 printf("\n");

      strcat(pathName, genomeName);
      strcat(pathName, "/");
      if (strcmp(sassoc("3U", terminalInput), "on") == 0) utr3 = 1;
      if (strcmp(sassoc("5U", terminalInput), "on") == 0) utr5 = 1;
      if (strcmp(sassoc("INTRONS", terminalInput), "on") == 0) intron = 1;
      if (strcmp(sassoc("EXONS", terminalInput), "on") == 0) exon = 1;
      if (strcmp(sassoc("CDS", terminalInput), "on") == 0) cds = 1;
      if (strcmp(sassoc("ORF", terminalInput), "on") == 0) orf = 1;
      printf("<span class=\"resultfile\">Click here to download file: <A HREF=\"/resultfiles/%s\"> %s</A></span><br>\n", result, result);
   }
   else {
/**  for testing special limited searches
      exon = intron = 1;
**/  
      strcpy(pathName, altPathName);
   }
   limited = utr3 | utr5 | intron | exon | cds | orf;
   if ((utr3 | utr5 | cds | orf) & (intron | exon)) {
      intron = exon = 0;
      both("\nRESTRICTING to coding sequence regions is mutually exclusive \n");
      both("from restricting to exons or introns. Only the\n");
      both("search over CDS regions will be performed in this run\n\n");
   }

   //get size of various components of the genome
   qnames = fopen(fullName (pathName, "datasize.txt"), "r");
   if (qnames == NULL) {
      printf("No datasize.txt file found in this directory -- quitting\n");
      exit(0);
   }
   fscanf(qnames, "%lu %d %d %d %d", &DNALen, &mapSize, &ncontigs, 
          &extableSize, &synCount);
   // To allow program to work against old genome format
   mapSize = 0;

   if (utr3 | utr5 | orf | cds) { //using special mrna-only pseudogenome
      fscanf(qnames, "%u %u", &pseudoLen, &mrnaCount);
      DNALen = pseudoLen;
   }
   fclose(qnames);

//determine how many bytes the genome will take (function of how many 
//nucleotides are processed in parallel
   i = DNALen/4 + scanWidth;
   if (scanWidth & (scanWidth - 1)) {
#if (scanWidth == 3)
      i = (i * 4) / 3 + 1;
#else
      i = (i / scanWidth) * 8 + 1;
#endif 
   }
   genomeLength = i+16;
//   DNAstring = malloc(genomeLength);
#ifdef DEBUG
   sprintf(output, "original genome has %u nucleotides, %u bytes\n",
          DNALen, DNALen/4);
   both(output);
   sprintf(output, "allocated %u bytes for the genome\n", genomeLength);
   both(output);
#endif
   getTime(&time1);
   readGenome();
   getTime(&time2);
      sprintf(output, "mRNA space for your genome has %u nucleotides\n\n",
              pseudoLen);
   if (utr3 | utr5 | orf | cds) {
   } else {
      sprintf(output, "Your genome has %lu nucleotides\n\n", DNALen);
   }
   both(output);
   chrLengths();

   
#ifdef TIMING
   sprintf(output, "Time to read genome was %d milliseconds\n",
      getDiffMillisecs(&time1, &time2));
   both(output);
#endif

   avbase = malloc (50000000*sizeof(unsigned long int));
   pvbase = malloc (50000000*sizeof(motifDes*));
   motifs = malloc(500*sizeof(motifDes));


outerLoop:
   //here, we allow the user to input motifs 
   motifCount = 0;
   patternID = 'A';
   errcount = 0;

   k = 0;
   sprintf(output, "Restrict search to chrm (return for no restriction)");
   both(output);
   getLine(line, stdin);
   chrLimit = malloc(strlen(line)+1);
   strcpy(chrLimit, line);
   ltdSearch = strlen(line);
   fprintf (save, "%s\n", line);
   bool skipToSearch = false;
      
//readInputLoop:
    do{
        while(true){
           if (fromTerminal != 0){
               sprintf(output, "Type in motif %c:", patternID);
               both(output);
               getLine(line, stdin);   //read a motif from terminal
               fprintf(save, "%s\n", line);
               i = fixinput(line); //convert input to upper case
               if (i) { // improper input detected
                  printf("Improper input character %c in position %d\n", line[i-1], i);
                  continue;
               }
               if (line[0] == 0) { // null input -- end of motifs
                  both("\n");
                  skipToSearch = true;
                  break;
               }
               motifs[motifCount].motif = malloc(strlen(line) + 1);
               strcpy (motifs[motifCount].motif, line);
               motifs[motifCount].name = patternID;
               motifs[motifCount].direction = '+';
               motifs[motifCount+1].motif = malloc(strlen(line)+1);
               makerevcompl(motifs[motifCount].motif, motifs[motifCount+1].motif, 
                            strlen(line));
               motifs[motifCount+1].name = patternID;
               motifs[motifCount+1].direction = '-';
               motifs[motifCount+1].reverse = motifs[motifCount].motif;
               motifs[motifCount].reverse = motifs[motifCount+1].motif;
               motifCount += 2;
               patternID++;
               k++;
              }
          }
        if(!skipToSearch){
           for (i = 'A'; i < 'K'; i++) {
              line[0] = 'S'; line[1] = i; line[2] = 0; //construct name of field i
              //result = sassoc(line, terminalInput);
              sprintf(result, "%s", sassoc(line, terminalInput));
              if (strlen(result) == 0) continue;
              j = fixinput(result);
              if (j) {//improper input detected
                 sprintf (output, "In field %c improper input character %c in position %d\n", i, result[i-1], j);
                 both(output);
                 errcount++;
              }
              k++;
              motifs[motifCount].motif = malloc(strlen(result) + 1);
              strcpy (motifs[motifCount].motif, result);
              motifs[motifCount].name = i;
              motifs[motifCount].direction = '+';
              motifs[motifCount+1].motif = malloc(strlen(result)+1);
              makerevcompl(motifs[motifCount].motif, motifs[motifCount+1].motif, 
                           strlen(result));
              motifs[motifCount+1].name = i;
              motifs[motifCount+1].direction = '-';
              motifs[motifCount+1].reverse = motifs[motifCount].motif;
              motifs[motifCount].reverse = motifs[motifCount+1].motif;
              motifCount += 2;
           }
           //Get cluster size
           sprintf(line, "%s", sassoc("CS", terminalInput));  
           fprintf(save, "Minimum number of motifs required in a cluster: ");
           i = checkMinNum(line, &minGroupSize);  
           errcount += i;
 
           //Get boolean control expression
           sprintf(boolExp, "%s", sassoc("BC", terminalInput));
           fprintf(save, "Boolean Condition: ");
           i = checkBoolExpr(boolExp);
           errcount += i;

           //Get window size
           sprintf(line, "%s", sassoc("WS", terminalInput));
           fprintf(save, "Cluster Size: ");
           i = checkWindowSize(line, &windowSize);
           errcount += i;

 
           sprintf(line, "%s", sassoc("GN", terminalInput));
           if (strlen(line)) {
              printf("List of requested genes is %s<br>", line);
              geneList = malloc (sizeof(int) * strlen(line));
              ptr = line;
              while (strlen(ptr)) {
                 sscanf(ptr, "%s", result);
                 ptr += strlen(result);
                 if (*ptr == ',') ptr++;
                 while (*ptr == ' ') ptr++;
                 for (i = 0; i < synCount; i++) {
                     if (strcmp((char *)geneName[i], result) == 0) {
                        geneList[numberOfGenes] = geneNumber[i];
                        numberOfGenes++;
                        break;
                     }
                 }
                 if (i == synCount) {
                    printf("Gene %s was not found for this genome<br>", result);
                 }
              }      

           }
   
           sprintf(line, "%s", sassoc("GD", terminalInput));
           geneDistance = 0;
           sscanf(line, "%d", &geneDistance);
   
           if (errcount) return 0; //if errors in form, abandon search because
                                 //all error messages have been output to screen
        }
            //skipToSearch 
           if (k == 0) {
          both("Please enter at least one pattern to search for.\n");
     }while(fromTerminal);
      exit(0);
   }
   motifs[motifCount].motif = 0;
   tables = makeTables(motifs);
#ifdef DEBUG
   sprintf(output, "Number of Tables = %d\n", numberOfTables);
   both(output);
#endif

   summaryData = calloc((numberOfChromosomes + 2) * (k+1) * 2, sizeof(int));
   if (summaryData == NULL) {
      //no room for summary data -- abort run
      both ("No room for summary data -- run discontinued\n");
      exit(0);
   }


#ifdef DEBUG
   sprintf(output, "sizeof tables[] = %d, or maybe %d\n",
        sizeof(tables), sizeof (*tables));
   both(output);
   sprintf(output, "Will print %d table entries\n", 3*(1<<(2*scanWidth)));
   both(output);
   for (i = 0; i < 3*(1<<(2*scanWidth)); i++) {
      makestring(i, str);
      sprintf (output, "%5d, %s, %8x\n", i, str, tables[i]);
      both(output);
   }
#endif


   repetitions = 1;
#ifdef TIMINGx
   sprintf(output, "How many repetitions of this search?\n");
   both(output);
   scanf("%d", &repetitions);
   getchar();  //remove character that ended previous input
   fprintf(save, "%d\n", repetitions);
#endif

   getTime(&time1);
   for (i = 0; i < repetitions; i++) {
      av = avbase;
      pv = pvbase;
      if (((utr3 | utr5 | orf | cds) == 0) && (geneList != NULL)) {
         numberOfMatches = 0;
         for (j = 0; j < numberOfGenes; j++) {
             relocation = arm[map[geneList[j]].annot_number];
             numberOfMatches += do_the_search(
                map[geneList[j]].start - geneDistance + relocation,
                map[geneList[j]].stop + geneDistance + relocation);
         }
      } else if (((utr3 | utr5 | orf | cds) == 0) && (geneDistance)){
         numberOfMatches = 0;
         sum = 0;
         for (j = 0; j < mapSize; j++) {
             relocation = arm[map[j].annot_number];
             sum += (map[j].stop - map[j].start + 2*geneDistance);
             numberOfMatches += do_the_search(
                map[j].start - geneDistance + relocation,
                map[j].stop + geneDistance + relocation);
         }
  printf("Looked at %u nucleotides<br>", sum);
      } else {
         if (ltdSearch == 0) {
            numberOfMatches = do_the_search(0, DNALen);
         } else {
            j = chrCont(chrLimit);
            first = arm[j];
            k = nextCont(j);
            last = arm[k]-1;
            numberOfMatches = do_the_search(first, last);
         }
      }
   }
   hsort2(avbase, pvbase, av - avbase);
   getTime(&time2);
#ifdef TIMING
   sprintf(output, "The search took %d milliseconds, found %d matches\n",
       getDiffMillisecs(&time1, &time2), (int)(av - avbase));
   both(output);
#endif
 
   bool goToSearch = false;
   do{
       if (fromTerminal == 0) {
          goToSearch = true;
          break;
       }
   
       sprintf(output, "Enter minimum number of motifs required in a cluster: "); 
       both (output);
       getLine(line, stdin);
   }while (checkMinNum(line, &minGroupSize));

   if(!goToSearch){
     do{
       sprintf(output, "Enter boolean condition (hit return for none):\n");
       both(output);
       getLine(boolExp, stdin); 
      while(checkBoolExpr(boolExp));

    readWindowSize:
    do{
       sprintf(output, "Enter window size (max #bp in cluster): ");
       both(output);
       getLine(line, stdin);
       fprintf(save, " %s\n", line);
       windowSize = 0;
       sscanf(line, "%d", &windowSize);
       if (windowSize <= 0) {
          sprintf(output, 
            "response to this question must be a positive integer.\n");
          both(output);
          if (fromTerminal) continue;
          else errcount++;
       }
       both("\n");

       if (errcount) return 1;
     }while(fromTerminal);
   }
//readyToSearch:
   if (fromTerminal == 0) {
      both("Motifs searched for:\n");
      for (i = 0; i < motifCount; i++) {
         if (motifs[i].direction == '-') continue;
         printf("&nbsp&nbsp");
         sprintf(output, "    %c: ", motifs[i].name);
         both(output);
         if (fromTerminal==0) {
            printf("<font style='BACKGROUND-COLOR:%s'>",
               colorName[motifs[i].name - 'A']);
         }
         sprintf(output,"%s\n", motifs[i].motif);
         both(output);
         if (fromTerminal == 0) { 
            printf("</font>");
         }
      }
      sprintf(output, "Minimum cluster size = %d, window size = %d\n", 
           minGroupSize,
           windowSize);
      both(output);
      if (strlen(boolExp)) {
         sprintf(output, "Controlling boolean expression = %s\n", boolExp);
         both(output);
      }
      if (limited) {
         both("Search limited to ");
         if (utr3) both("3'utrs ");
         if (utr5) both("5'utrs ");
         if (orf)  both("ORFs ");
         if (cds)  both("CDSs ");
         if (intron) both("introns ");
         if (exon) both("exons");
         both("\n");
      }
      both("\n");
   }

   //Now, examine the hits together with the clustering requirements to
   //report the clusters that were found.
   reportClusters();
   both("\n\n");
   
   if (fromTerminal == 0) {
      printGeneList();
      return 1;
   }


   return 0;
}

void makestring(int x, char *s){ 
   int i, j;
   int len;
   char alphabet[] = "ACGT";
   len = strlen(s);
   for (i = 0; i < len; i++) {
     j = x & 3;
     s[len - i - 1] = alphabet[j];
     x>>=2;
   } 
}

int randomstring(char **array, int i, int len) {
   static long long int seed = 125*125*125*125;
   long long int product;
   int temp;
   int j; 
   unsigned char lookup[] = "ACGT";
   
   product = seed * 5;
   seed = product;
   temp = product;  
   for (j = 0; j < len; j++) { 
      array[i][j] = lookup[(temp>>(60-2*j)) & 3];
   }
   array[i][j] = 0;
   return temp;
}

void both(char stuff[]) {
   char buffer[2000];
   int i, j, n;

   if (printNone) return;
   if (fromTerminal) {
      printf("%s", stuff);
   }
   else { //communicating with an html page -- change \n to <.br>
      n = strlen(stuff);
      i = 0;
      j = 0;
      for (i = 0; i < n; i++) {
         if (stuff[i] != '\n') buffer[j] = stuff[i];
         else {
            buffer[j] = '<';
            buffer[j+1] = 'b';
            buffer[j+2] = 'r';
            buffer[j+3] = '>';
            j += 3;
         }
         j++;
         if (j > 995) {
            printf("%s", buffer);
            j = 0;
         }
      }
      if (j) {
         buffer[j] = 0;
         printf("%s", buffer);
      }   
   }
   if (printBoth) fprintf(save, "%s", stuff);
}

void uniquename(char *result) {
   struct tm {
      int tm_sec;
      int tm_min;
      int tm_hour;
      int tm_mday;
      int tm_mon;
      int tm_year;
      int tm_wday;
      int tm_yday;
      int tm_isdst;
      } *alpha;
   time_t now;
   now = time(NULL);
   alpha =(struct tm *) localtime(&now);
   sprintf (result, 
      "%d-%d-%.2d-%d\0", alpha->tm_mon+1, alpha->tm_mday,
       alpha->tm_year - 100, 
       alpha->tm_hour*3600 + alpha->tm_min * 60 + alpha->tm_sec);

}

void peter_test(alist data)
{
  char **blank ;
  terminalInput = data;
  mainAux ( -1, blank);
}

void printallsynonyms(char* resultstring, int syn_count){
   int i, j, k, m;
   int count;
   k = 0;
   count = 0;
   for (i = 1; i < primeNames[0]; i++) {
      for (j = k; j < syn_count; j++) {
         if (geneNumber[j] == primeNames[i]) {
            m = j;
            while (geneNumber[j] == geneNumber[m]) {
               if (count > 50) {
                  fprintf(save, "\n   ");
                  count = 3;
               }
               fprintf(save, "%s ", geneName[m]);
               count = count + strlen( (char *)geneName[m]) + 1;
               m++;
            }
            if (count > 0) fprintf(save, "\n");
            count = 0;
            k = (k > 20)? k - 20 : 0;
            break;
         }
      }
   }
   return;
}

void printGeneList() {
   int i, j;
   j = 0;
   fprintf(save, 
           "\nHere is a list of all genes and their synonyms"
           " mentioned in this run.\n\n");
   for (i = 0; i < mapSize; i++) {
      if (map[i].touched) {
         map[i].touched = 0;
         while(j < synCount && geneNumber[j] < i) j++;
         //Now we are at the right place in the synonym list
         while(j < synCount && geneNumber[j] == i) {
            fprintf(save, "%s\n", geneName[j]);
            j++;
         } 
      }
   }
}

int checkMinNum(char *line, int *minGroupSize) {
   char output[100];
   fprintf(save, "%s\n", line);
   *minGroupSize = 0;
   sscanf(line, "%d", minGroupSize);  //group size must be positive integer
   if (*minGroupSize <=0) {
      sprintf(output, 
         "response to Cluster Size must be a positive integer.\n");
      both(output);
      return 1;  //indicate error
   }
   return 0;  //indicate good minGroupSize input
}

int checkBoolExpr(char *boolExp) {
   char output[100];
   fprintf(save, "%s\n", boolExp);
   if (booleanSyntax(boolExp) == 0) {  //check for well-constructed boolean
      sprintf(output, "Ill-formed boolean -- reenter\n");
      both (output);
      return 1; //indicate error in boolean expression
   }
   return 0;  //indicate good boolean expression
}

int checkWindowSize(char *line, int *windowSize) {
   char output[100];
   fprintf(save, " %s\n", line);
   *windowSize = 0;
   sscanf(line, "%d", windowSize);
   if (*windowSize <= 0) {
      sprintf(output, 
        "response for window size must be a positive integer.\n");
      both(output);
      return 1;  //indicate bad window size
   }
   both("\n");
   return 0;  //indicate good window size
}

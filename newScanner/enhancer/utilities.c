//Copyright (c) in silico Labs, LLC 2008, 2009, 2105

#include <stdlib.h>
#include <stdio.h>
#include "scanner.h"
#include "timer.h"
#include <float.h>
#include <string.h>

void tableOutput(char *);
void endRow();

void getLine(char line[], FILE *f) {
   int i;
   for (i=0; i<5000; i++) {
      line[i] = getc(f);
      if (line[i] == '\n' || line[i] == EOF) break;
      if (line[i] == '\b') {
         i = (i > 1 ? i - 2 : -1);
      }

   }
   line[i] = 0;
}

int readContigs(FILE *qnames) {
   char line[1000], temp[100], temp2[100];
   int k = 0;
   int i, isize;
   while (!feof(qnames)) {
      getLine(line, qnames);   //read a line of contig info
      //Format of a contig line is name of contig file, chr. no, contig no, start of contig,
      //size of contig.  We only use the start of contig for now.
      if (line[0] == 0) break;
      sscanf(line, "%s%s%d%lu%d", temp, temp2, &i, &arm[k], &isize);
      annotName[k] = malloc(strlen(temp)+1);
      strcpy(annotName[k], temp);

           for (i = 0; i < numberOfChromosomes; i++) {
              if (strcmp(temp2, chrName[i]) == 0) {
                 chromoName[k] = chrName[i];
                 break;
              }
           }
           if (i == numberOfChromosomes) printf("Chromo name not found\n");

      k++;
   }
   return k;
}

void readgenenames(FILE *allnames, int syn_count){

   int i, n;
   char line[500], line2[1000];

   for (i = 0; i < syn_count; i++) {
      getLine(line2, allnames);
      sscanf(line2, "%[^,],%d", line, &geneNumber[i]);
      n = strlen(line);
      geneName[i] = malloc(n+1);
      strcpy((char *)geneName[i], line);
   }
}

int fixinput(char line[]) {
    int i, n;
    int j;

    n = strlen(line);
    for (i = 0; line[i]; i++) {
       line[i] &= 0xDF; //force upper case
       if (line[i] < 'A' || line[i] > 'Z') return i+1;
       switch(line[i]) {
          case 'E':
          case 'F':
          case 'I':
          case 'J':
          case 'L':
          case 'O':
          case 'P':
          case 'Q':
          case 'X':
          case 'Z':
             return(i+1);
      }
   }
   return 0;
}

void makerevcompl(char *from, char *to, int num) {
   int i;

   for (i = 0; i < num; i++) {
      switch (from[i]) {
         case 'A':
            to[num-i-1] = 'T';
            break;
         case 'B':
            to[num-i-1] = 'V';
            break;
         case 'C':
            to[num-i-1] = 'G';
            break;
         case 'D':
            to[num-i-1] = 'H';
            break;
         case 'G':
            to[num-i-1] = 'C';
            break;
         case 'H':
            to[num-i-1] = 'D';
            break;
         case 'K':
            to[num-i-1] = 'M';
            break;
         case 'M':
            to[num-i-1] = 'K';
            break;
         case 'N':
            to[num-i-1] = 'N';
            break;
         case 'R':
            to[num-i-1] = 'Y';
            break;
         case 'S':
            to[num-i-1] = 'S';
            break;
         case 'T':
         case 'U':
            to[num-i-1] = 'A';
            break;
         case 'V':
            to[num-i-1] = 'B';
            break;
         case 'W':
            to[num-i-1] = 'W';
            break;
         case 'Y':
            to[num-i-1] = 'R';
            break;
         default:
            to[num-i-1] = 0;
      }
   }
   to[num] = 0;
}
 
void cleanLine(char line[]) {
   int i;
   int j;
   while (line[i]) {
      if (line[i] == '\b') {
         j = (i > 0) ? i - 1 : 0;
         while (line[j+2]) {
            line[j] = line[j+2];
            j++;
         }
         line[j] = 0;
         i--;
      }
   }
}
          
int findchrm(unsigned int x) {
   int i;
   for (i = 1; arm[i-1] < arm[i]; i++) {
      if ( x < arm[i]) break;
   }
   
   return i-1;
}

int annotation_search(int chr, int index) {
   int i;
   unsigned int u;
   unsigned long int v, w;
   
   for (i = 0; i < mapSize; i++) {
      u = map[i].annot_number;
      v = map[i].start;
      w = map[i].stop;
      if (map[i].annot_number > chr) break ;
      if (map[i].annot_number < chr) continue;
      if (map[i].start < index) {
         if (map[i].stop >= index) return i;
         else continue;
      }
      if (i > 0 && map[i-1].annot_number == map[i].annot_number) return i-1;
      else return i;
   }
   if (i > 0 && map[i-1].annot_number != chr) return -1;
   return i-1;
}

char* fullName(char* path, char* name) {
   static char buffer[1000];
   strcpy(buffer, path);
   strcat(buffer, name);
   return buffer;
} 

unsigned long int chrStart(int chr) {
   //Find the beginning of the chromosome which contains the contig chr
   int i;
   i = chr;
   while (i > 0 && (strcmp(chromoName[i-1], chromoName[i]) == 0)) i--;
   return (arm[i]);
} 

int chrIndex(int chr) {
   //Find the index for the chromosome containing the contig chr
   unsigned int i;
   for (i = 0; i < numberOfChromosomes; i++) {
      if (chromoName[chr] == chrName[i]) return i;
   }
   both("Error - chromosome name not found\n");
   return 0;
}
 
int chrCont(char* name) {
   // find first contig for chromosome name
   int i;
   for (i = 0; i < ncontigs; i++) {
       if (strcmp(name, chromoName[i]) == 0) return i; 
   }
   printf("Chromosome name not found\n");
   return(0);
}
 
int nextCont(int c) {
   //find first contig on chromosome after the one containing c
   int i;
   for (i = c+1; i < ncontigs; i++) {
      if (strcmp(chromoName[c], chromoName[i]) == 0) continue;
      return(i);
   }
}

void outputSummary() {
   char output[1000];
   int i, j, k, sum;

   if (utr3 | utr5 | cds | orf) return;
   if (fromTerminal == 0 ) {
      printf("<table><tr><td width=800, align='center'>");
   }
   both("                    SEARCH SUMMARY\n");
   if (fromTerminal == 0) {
      printf("</td></tr></table>");
   }

   if (fromTerminal == 0) {  //start table format for output
      printf("<table cellpadding = '0' cellspacing='0' border='1'> \n");
      printf("<tr>");
   }
   tableOutput("Chr Name");  //leave space for motif identifier
   for (i = 0; i < numberOfChromosomes; i++) {
      sprintf(output, " %7s", chrName[i]);
      tableOutput(output);
   }
   sprintf(output, "   Total");
   tableOutput(output);
   endRow();
   for (i = 0; i < numberOfMotifs/2; i++) {
      sum = 0;
      //print a row of information:
      if (fromTerminal == 0) printf("<tr>");   //start new table row
      sprintf(output, " Motif %c", motifs[2*i].name);
      k = motifs[2*i].name - 'A';
      tableOutput(output);
      for (j = 0; j < numberOfChromosomes; j++) {
         sum = sum + summaryData[k*(numberOfChromosomes + 1) + j];
         sprintf(output, " %7d",
            summaryData[k*(numberOfChromosomes + 1) + j]);
         tableOutput(output);
      } 
      summaryData[k*(numberOfChromosomes + 1) + numberOfChromosomes] = sum;
      sprintf(output, " %7d", sum);
      tableOutput(output);
      endRow();
   }
   //Compute total hits in each chromosome
   for (j = 0; j <= numberOfChromosomes; j++) {
      sum = 0;
      for (i = 0; i < numberOfMotifs/2; i++) {
         k = motifs[2*i].name - 'A';
         sum = sum + summaryData[k*(numberOfChromosomes + 1) + j];
      }
      summaryData[(numberOfMotifs/2)*(numberOfChromosomes+1) + j] = sum;
   }
   if (fromTerminal == 0) printf("<tr>");
   tableOutput("  Totals");
   for (j = 0; j <= numberOfChromosomes; j++) {
      sprintf(output, " %7d",
         summaryData[(numberOfMotifs/2)*(numberOfChromosomes+1) + j]);
      tableOutput(output);
   }
   endRow();

   if (fromTerminal == 0) printf("<tr>");
   tableOutput("Clusters");
   sum = 0;
   for (j = 0; j < numberOfChromosomes; j++) {
   sum += summaryData[(1+numberOfMotifs/2)*(numberOfChromosomes+1) + j];
      sprintf(output, " %7d",
         summaryData[(1+numberOfMotifs/2)*(numberOfChromosomes+1) + j]);
      tableOutput(output);
   }
   sprintf(output, " %7d", sum);
   tableOutput(output);
   endRow();
   if (fromTerminal == 0) {
      printf("</table>");
   }
}

void tableOutput(char *output) {
   if (fromTerminal == 0) {
      printf("<td width='100' align='right'>");
      if (numberOfChromosomes > 10) {
         printf("<small>");
      }
   }
   both(output);
   if (fromTerminal == 0) {
      if (numberOfChromosomes > 10) {
         printf("</small>");
      }
      printf("</td>");
   }
}

void endRow() {
   if (fromTerminal) {
      both("\n");
   } else {
      fprintf(save, "\n");   //end of line to file output
      printf("</tr>");       //end of row to html output
   }
}

int findSplice(unsigned int x) {
    int top, bottom, mid;
    int count;
    top = mrnaCount;
    bottom = 0;
    count = 0;
    while (top > bottom) {
       count++;
       if (count > 16) {
         printf("Trouble\n");
       }
       mid = (top + bottom)/2;
       if (x == mrnaInfo[mid].mrnaStart) return mid;
       if (x < mrnaInfo[mid].mrnaStart) top = mid - 1;
       else bottom = mid + 1;
    }
    if (x < mrnaInfo[top].mrnaStart) return top-1;
    if (x >= mrnaInfo[top+1].mrnaStart) {
       printf("returnd top+1\n");
       return top+1;
    }
    return top;
}

long int chrSize(int chr) {
    int i, j;
    long int pos;
   
    for (i = 0; i < ncontigs; i++) {
        if (strcmp(chrName[chr], chromoName[i]) == 0) {
           pos = arm[i]; //start of chromosome chr
           break;
        }
    }
    for (j = i+1; j < ncontigs+1; j++) {
        if (strcmp(chrName[chr], chromoName[j]) != 0) {
           return arm[j] - pos;
        }
    }
    printf("Could not find end of chromosome %s.", chrName[chr]);
   
}

chrLengths() {

   int i;
   long int chrL;
   char line[100];

   if (fromTerminal == 0) return;
   for (i = 0; i < numberOfChromosomes; i++)  {
      chrL = chrSize(i);
      sprintf(line, "Chromosome %s has length %ld\n", chrName[i], chrL);
      both(line);
   }  
}

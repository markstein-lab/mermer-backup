//Copyright (c) in silico Labs, LLC 2008, 2009, 2015

#include <stdlib.h>
#include <stdio.h>
#include "timer.h"
#include "htmlForm.h"
#include <float.h>
#include <string.h>
/***
#define TIMING
***/

void makestring(int, char*);
void inputMissing();

alist terminalInput;


int doShowCluster (alist stuff) {
 
   unsigned char **strings;
   int number;
   FILE *datafile;
   int repetitions;
   int motifCount;
   int clusterLength;
   char patternID;
   char filename[100];
   char line[100];
   char junk[100];
   char htmlLine[10000];

   unsigned int i, j, k, n, m;
   int index;
   int numberOfHits, flag, motifNo, len, current;
   int realStart;
   int found_start;
   unsigned int genomeLength;

   unsigned char c, cmod;
   char* ptr;
   char lookup[] = "ACGT";
   char nullstr[] = "   ";

   char motifName[26];
   char *motifDNA[26];
   char motifSeq[26];
   int location[26];
   char *colorName[26] ={"cyan", "yellow", "lightgreen", \
                         "lightpink", "orange", "greenyellow", \
                         "dodgerblue", "mistyrose", "lightsalmon", "khaki"};
   char* searchstr(char*, char*);
   char* skipwhitesp(char*);

normalStart:
#ifdef DEBUG
#ifdef WIN32
   printf("Compiled for Win32\n");
#endif
#endif
   terminalInput = stuff;

   filename[0] = 0;
   sprintf(filename, "%s", sassoc("FN", terminalInput));
   if (strlen(filename) == 0) {
      printf ("Data not found<br>");
      return;
   }

   line[0] = 0;
   sprintf(line, "%s", sassoc("CN", terminalInput));
   if (strlen(line) == 0) {
      printf ("Cluster number not found<br>");
      return;
   }
   sscanf(line, "%d", &number);
   if (number < 1) {
      printf("Cluster number is less than 1 -- error<br>");
      return;
   }

   datafile = fopen(filename, "r");
   if (datafile == NULL) {
      printf("Could not open file '%s'<br>", filename);
      return;
   }
   flag = 1;
   printf("<h2>Nucleotides for Cluster %d</h2><br>\n", number);
   while (flag) { //loop to skip part of the input
      getLine(line, datafile);
      if (feof(datafile)) {
         return;
      }
      if (searchstr(line, "Motifs searched")) break;
   }

   printf("%s<br><br>", line);
   motifCount = 0;
   i = 0;
   while (flag) {
      getLine(line, datafile);
      ptr = searchstr(line, ":");
      if (ptr == 0) break;  //end of reading motif names?
      motifName[i] = *(ptr-2);  //letter of motif
      sscanf(skipwhitesp(ptr), "%s", junk);
      if (m = strlen(junk) ==0) {
         printf("error reading motifs<br>");
         return;
      }   
      index = motifName[i] - 'A';
      if (index < 0 || index >= 10) {
         printf("Bad motif name %c<br>", motifName[i]);
         return;
      }
      motifDNA[index] = malloc(m+1);
      strcpy(motifDNA[index], junk);
      motifCount++;
      printf(" %c: <FONT style='background-color:%s'>",
             motifName[i], colorName[index]);
      printf("%s</FONT><br>", motifDNA[index]);
      i++;
   }

   //Now look for the requested cluster
   while(flag) {
      getLine(line, datafile);
      if (feof(datafile)) {
         inputMissing();
         return;
      }
      if (line[0] == 'C' &&
          line[1] == 'l' && (ptr = searchstr(line, "Cluster"))) {
         sscanf(ptr, "%d", &n);
         if (number == n) break;  //found the cluster?
      }
   }

   printf("<br>%s<br>", line);
   ptr = searchstr(ptr, "length");
   sscanf(ptr, "%d", &clusterLength);   

   getLine(line, datafile);
   printf("%s<br>", line);
   ptr = (searchstr(line, "MOTIFS "));

   i = 0;
   //read sequence of letters identifying motifs found
   while(*ptr >= 'A' && *ptr <= 'J') {
      motifSeq[i] = *ptr;
      ptr+=2;
      i++;
   }
   numberOfHits = i;
   ptr = searchstr(ptr, "at positions");
   //read relative locations where hits were found
   for (i = 0; i < numberOfHits; i++) {
      sscanf (ptr, "%d", &location[i]);
      ptr = searchstr(ptr, ",");
      if (ptr == 0) break;
   }
   location[numberOfHits] = 1<<30;

   getLine(line, datafile);  //get line showing chromosome & position
   printf("%s<br><br>", line);

   while(flag) { //get to the cluster info
      getLine(line, datafile);
      c = line[0];
      c &= 0xDF; 
      if (c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'U') break;
      printf("%s<br>\n", line);
   }
   printf("<br>\n");

   //Now we are ready to print the dna of the motifs 
   motifNo = 0;
   current = 0;
   htmlLine[0] = 0;
   printf("<TT><font size='3'>\n");
   found_start = 0;
   while (flag) {
      for (i = 0; i < 65; i++) { 
          if ((found_start == 0 && line[i] >='A' && line[i] <='Z') ||
               (found_start && current - realStart == location[motifNo])) { 
             if (found_start == 0) {
                found_start = 1;
                realStart = current;
             }
             strcat (htmlLine, "<FONT style='BACKGROUND-COLOR:");
             strcat (htmlLine, colorName[motifSeq[motifNo] - 'A']);
             strcat (htmlLine, "'>");
             m = strlen(motifDNA[motifSeq[motifNo] - 'A']);
             for (j = 0; j < m; j++) {
                len = strlen(htmlLine);
                htmlLine[len] = line[i];
                htmlLine[len+1] = 0;
                if (line[i] == ' ') j--;
                else current++;
                i++;
                if (i == 65) {
                   i = 0;
                   getLine(line, datafile);
                   printf("%s<BR>\n", htmlLine);
                   htmlLine[0] = 0;
                }
             }
             i--;  //i will be bumped at the bottom of the loop
             strcat (htmlLine, "</FONT>"); 
             motifNo++;
          }
          else {  //normal processing of a character
             if (line[i] == '\n' || line[i] == 0) break;
             len = strlen(htmlLine);
             htmlLine[len] = line[i];
             htmlLine[len+1] = 0;
             if (line[i] != ' ') current++;
          }
      }
      printf("%s<BR>\n", htmlLine);
      getLine(line, datafile);
      htmlLine[0] = 0;
      if (strlen(line) == 0) return;
   }
   return(0);
}
      
 
void inputMissing() {
   printf("Data missing from input file<br>");
   return;
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





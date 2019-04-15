//Copyright (c) in silico Labs, LLC, 2015

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "scanner.h"
#include "timer.h"

void recursiveEnter (char*, int, LOOKUPWORD, LOOKUPWORD *);

LOOKUPWORD* makeTables(motifDes *motifs) {
   int i, j, k, m, maxlen;
   int len;
   int tableIndex;
   int statesPerWord;
   char **paddedMotifs;
   char lookupString[64];
   int extendedMotifLength;
   LOOKUPWORD mask, tempMask, unit;
   int paddedLength;
   TIMEVAL time1, time2;

   getTime (&time1);
   unit = 1;
   statesPerWord = wordWidth / scanWidth;
#ifdef DEBUG
   printf("wordWidth = %d, scanWidth = %d, statesPerWord = %d\n",
           wordWidth, scanWidth, statesPerWord);
   printf("Size of mask = %d\n", sizeof(mask));
#endif
   i = 0;
   while (motifs[i].motif) i++;
   numberOfMotifs = i;
   paddedMotifs = malloc(sizeof(char*) * (i + 1));

   maxlen = strlen(motifs[0].motif);
   for (i = 0; i < numberOfMotifs; i++) {
#ifdef DEBUG
      printf("motif[%d] = %s\n", i, motifs[i].motif);
#endif
      len = strlen(motifs[i].motif);
      maxlen = (len > maxlen) ? len : maxlen;  //finding length of longest motif
   }
   extendedMotifLength = maxlen + scanWidth - 1;
   numberOfTables = (extendedMotifLength + scanWidth - 1) / scanWidth;
   paddedLength = (numberOfTables + 1)*scanWidth;
#ifdef DEBUG
   printf("extendedMotifLength = %d, number of tables = %d, padded length = %d\n",
      extendedMotifLength, numberOfTables, paddedLength);
#endif
   for (i = 0; i < numberOfMotifs; i++) {
      /* For each motif, construct a padded motif with scanWidth-1 Ns padding
         the motif on the left, and a total padded motif length of 
         maxlen + 2*scanWidth - 2.w
      */
      len = strlen(motifs[i].motif);
      paddedMotifs[i] = malloc(paddedLength+1);   
      for (j = 0; j < scanWidth-1; j++) paddedMotifs[i][j] = 'N';
      for (k = 0; k < len; k++) paddedMotifs[i][j++] = motifs[i].motif[k];
      while(j < paddedLength) paddedMotifs[i][j++] = 'N';
      paddedMotifs[i][j] = 0;
#ifdef DEBUG
      printf("paddedmotifs[%d] = %s\n", i, paddedMotifs[i]);
#endif
   }
   paddedMotifs[i] = 0;  //end of list of padded motifs
   matches = calloc(numberOfTables, sizeof(LOOKUPWORD)*(1<<(2*scanWidth)));
#ifdef DEBUG
   printf("%d %d %d %d\n",
      numberOfTables, sizeof(LOOKUPWORD), (1<<(2*scanWidth)),
      numberOfTables*sizeof(LOOKUPWORD) * (1<<(2*scanWidth)));
   printf("size of matches = %d %x, %x\n",
     sizeof(matches), sizeof(matches), 
     numberOfTables * sizeof(LOOKUPWORD) * (1<<2*scanWidth));

   printf("extended motif length = %d, number of tables = %d\n",
     extendedMotifLength, numberOfTables);
#endif

   /*There should be some heuristic code to arrange the motifs in an 
     advantageous order, so that motifs of nearly the same length share
     a table, and that motifs whose leading characters are similar share
     a table.  We leave this out for now, in the interests of getting
     a minimal program to work
   */

   /* BUILD THE TABLE */
   
   for (i = 0; i < numberOfMotifs; i++) {
      j = i % statesPerWord;
      mask = (unit<<(j * scanWidth));
#ifdef DEBUG
#ifdef NORMALWORD
      printf( "i = %d, j = %d, mask = %x\n", i, j, mask);
#else
      printf( "i = %d, j = %d, mask = %llx\n", i, j, mask);
#endif
#endif
      for (k = 0; k < numberOfTables*scanWidth ; k++) {
         tempMask = mask << (k % scanWidth);
         tableIndex = ((int) k) / ((int)scanWidth); 
#ifdef DEBUG
#ifdef NORMALWORD
         printf("k = %d, tempMask = %x, tableIndex = %d\n",
                k, tempMask, tableIndex);
#else
         printf("k = %d, tempMask = %llx, tableIndex  %d\n",
                k, tempMask, tableIndex);
#endif
/***********
         printf("k mod(scanWidth) = %d, k/scanWidth = %d\n", 
                (k%scanWidth), (k/scanWidth));
*************/
#endif
         for (m = 0; m < scanWidth; m++) lookupString[m] = paddedMotifs[i][k+m];
         lookupString[m] = 0;
         recursiveEnter(lookupString, tableIndex, tempMask, matches);
      }
   }

   for (i = 0; i < numberOfTables; i++) {
      j = 0;
      for (k = 0; k < 1<<(2*scanWidth); k++) {
         if (matchTable(i, k)) j++;
      }
#ifdef DEBUG
      printf("Table %d has %d non zero entries (out of %d)\n",
         i, j, (1<<2*scanWidth));
#endif
   }

   //free temporary allocated storage
   for (i = 0; i < numberOfMotifs; i++) free (paddedMotifs[i]);
   free (paddedMotifs);
   
   getTime(&time2);
#ifdef TIMING
   printf ("Table generation took %d milliseconds.\n",
      getDiffMillisecs(&time1, &time2));
#endif
   return matches;   //return pointer to tables
}

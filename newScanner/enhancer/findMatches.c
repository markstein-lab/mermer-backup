//Copyright (c) in silico Labs, LLC 2008, 2009

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "scanner.h"
#include "timer.h"

typedef enum { false = 0, true = !false } bool;
/*************
#define nextDNA() \
   if (shift >= 0) { \
      index = (geneWord >> (shift)) & mask; \
      shift -= (2 * scanWidth); \
   } \
   else { \
      index = geneWord << (-shift); \
      geneWord = genome[i++]; \
      index = (index | (geneWord >> (wordWidth + shift))) & mask; \
      shift = wordWidth - 2*scanWidth + shift; \
   }
**********/

   int whereIndex;
   int numberOfTimesCalled ;
   int numberOfTimesMotifFound;
   int exindex;
   unsigned char lookup[] = "ACGT";
   void hsort2 (unsigned long int*, motifDes [], long int);
   int findSplice(unsigned int x);
   void outputSummary();
   int invalid(unsigned int, unsigned int);
   int counts[500];
   unsigned long int *oldav;   

int do_the_search(unsigned long int start, //where to start in the genome
                   unsigned long int stop  //where to stop search in the genome
                  ) {
    
   TIMEVAL time1, time2;
   int j, depth;
   unsigned long int i;
   unsigned long int startIndex, stopIndex;
   unsigned long int index;
   LOOKUPWORD mask0, mask1, mask2, mask3, mask4, mask5, mask6;
   LOOKUPWORD mask[500];

#ifdef STATS
   //int counts[500];
#endif
   for (i = 0; i < 500; i++) counts[i] = 0;
   unsigned int temp;
   void identifyMatches(LOOKUPWORD, unsigned long int, int);
   //unsigned long int *oldav;

   startIndex = start/scanWidth;
   stopIndex = stop/scanWidth;
   numberOfTimesCalled = 0;
   numberOfTimesMotifFound = 0;
   exindex = 0;
   oldav = av;

#ifdef STATS
   //for (i = 0; i < 500; i++) counts[i] = 0;
   goto  generaldeep;
#endif
    bool firstTime1 = true;
    bool firstTime2 = true;
    bool firstTime3 = true;
    bool firstTime4 = true;
    bool firstTime5 = true;

   void updateMasks(int num, int i){
        index = getDNA(i);
        switch(num){ //update from [0, num] masks
        case 6:
            mask6 = mask5 & matchTable(6, index);
        case 5:
            mask5 = mask4 & matchTable(5, index);
        case 4:
            mask4 = mask3 & matchTable(4, index);
        case 3:
            mask3 = mask2 & matchTable(3, index);
        case 2:
            mask2 = mask1 & matchTable(2, index);
        case 1:
            mask1 = mask0 & matchTable(1, index);
        case 0:
            mask0 = matchTable(0, index);
        }
   }
   switch (numberOfTables) { //notice no breaks in table, so case 4 will go to case 3, 2,5, and so on
      case 4: 
            i = startIndex;
            updateMasks(0, i);
            i++;
            if (i > stopIndex) return finished(counts, numberOfTables);

            firstTime1 = true;
            firstTime2 = true; 

            while(true){
                while(mask2 == 0 || firstTime2){
                    firstTime2 = false; 

                    while(mask1 == 0 || firstTime1){
                        firstTime1 = false;
                        updateMasks(1, i);
                        i++;
                        if (i > stopIndex) return finished(counts, numberOfTables);
                    }
                    updateMasks(2, i);
                    i++;
                    if (i > stopIndex) return finished(counts, numberOfTables);
                }

                updateMasks(3,i);
                if (mask3) { // we have found one or more matches
                  identifyMatches(mask3, i, numberOfTables);
                }
                i++;
                if (i > stopIndex) return finished(counts, numberOfTables);
            }
      case 3:
           i = startIndex;
           updateMasks(1, i);
           i++;
           if (i > stopIndex) return finished(counts, numberOfTables);

           firstTime1 = true; 
            while(true){
               while(firstTime1 || mask1 == 0){
                   firstTime1 = false; 
                   updateMasks(1, i);
                   i++;
                   if (i > stopIndex) return finished(counts, numberOfTables);
               }

               updateMasks(2, i);
               if (mask2) { // we have found one or more matches
                  identifyMatches(mask2, i, numberOfTables);
               }
               i++;
               if (i > stopIndex) return finished(counts, numberOfTables);
            }
      case 2: 
           i = startIndex;   //look up if assigning mask1 a value and idenitfying a match will do anything
           updateMasks(0, i);
           i++;
           if (i > stopIndex) return finished(counts, numberOfTables);

           while(true){
               updateMasks(1, i);
               if (mask1) { // we have found one or more matches
                  identifyMatches(mask1, i, numberOfTables);
               }
               i++;
               if (i > stopIndex) return finished(counts, numberOfTables);
            }
            //continue onto next case     
      case 5: 
           i = startIndex;
           updateMasks(0,i);
           i++;
           if (i > stopIndex) return finished(counts, numberOfTables);
           firstTime1 = true;
           firstTime2 = true;   
           firstTime3 = true;
           while(mask3 == 0 || firstTime3){
               firstTime3 = false; 

               while(mask2 == 0 || firstTime2){
                    firstTime2 = false;

                   while(mask1 == 0 || firstTime1){
                       firstTime1 = false; 
                       updateMasks(1,i);
                       i++;
                   if (i > stopIndex) return finished(counts, numberOfTables);
                   }
                   index = getDNA(i);
                   updateMasks(2,i);
                   i++;
                   if (i > stopIndex) return finished(counts, numberOfTables);
               }

               updateMasks(3,i);
               i++;
               if (i > stopIndex) return finished(counts, numberOfTables);
               if(mask3 == 0) continue;

               updateMasks(4,i);
               if (mask4) { // we have found one or more matches
                  identifyMatches(mask4, i, numberOfTables);
               }
               i++;
               if (i > stopIndex) return finished(counts, numberOfTables);
           }
      case 6: 
           i = startIndex;
           updateMasks(0,i);
           i++;
           if (i > stopIndex) return finished(counts, numberOfTables);

           firstTime1 = true;
           firstTime2 = true;
           firstTime3 = true;
           firstTime4 = true; 

           while(true){
               while(mask4 == 0 || firstTime4){
                   firstTime4 = false;
                   while(mask3 == 0 || firstTime3){
                       firstTime3 = false;
                       while(mask2 == 0 || firstTime2){
                           firstTime2 = false;
                           while(mask1 == 0 || firstTime1){
                               firstTime1 = false;

                               updateMasks(1,i);
                               i++;
                               if (i > stopIndex) return finished(counts, numberOfTables);
                           }
                           updateMasks(2, i);
                           i++;
                           if (i > stopIndex) return finished(counts, numberOfTables);
                       }
                       updateMasks(3,i);
                       i++;
                       if (i > stopIndex) return finished(counts, numberOfTables);
                   }
                   updateMasks(4,i);
                   i++;
                   if (i > stopIndex) return finished(counts, numberOfTables);
               }
              updateMasks(5,i);
              if (mask5) { // we have found one or more matches
                  identifyMatches(mask5, i, numberOfTables);
              }
              i++;
              if (i > stopIndex) return finished(counts, numberOfTables);
           }
      case 7: 
           i = startIndex;
           updateMasks(0,i);
           i++;
           if (i > stopIndex) return finished(counts, numberOfTables);

           firstTime1 = true;
           firstTime2 = true;
           firstTime3 = true;
           firstTime4 = true; 
           firstTime5 = true;

           while(true){
                while(mask5 == 0 || firstTime5){
                   while(mask4 == 0 || firstTime4){
                       firstTime4 = false;
                       while(mask3 == 0 || firstTime3){
                           firstTime3 = false;
                           while(mask2 == 0 || firstTime2){
                               firstTime2 = false;
                               while(mask1 == 0 || firstTime1){
                                   firstTime1 = false;

                                   updateMasks(1,i);
                                   i++;
                                   if (i > stopIndex) return finished(counts, numberOfTables);
                               }
                               updateMasks(2, i);
                               i++;
                               if (i > stopIndex) return finished(counts, numberOfTables);
                           }
                           updateMasks(3,i);
                           i++;
                           if (i > stopIndex) return finished(counts, numberOfTables);
                       }
                       updateMasks(4,i);
                       i++;
                       if (i > stopIndex) return finished(counts, numberOfTables);
                   }
                    updateMasks(5,i);
                    i++;
                    if (i > stopIndex) return finished(counts, numberOfTables);
                }
                updateMasks(6,i);
                if (mask6) { // we have found one or more matches
                    identifyMatches(mask6, i, numberOfTables);
                }
                i++;
                if (i > stopIndex) return finished(counts, numberOfTables);
           }
      case 1: 
          for(i = startIndex; i != stopIndex; i++){
               updateMasks(0, i);
               if (mask0) { // we have found one or more matches
                  identifyMatches(mask0, i, numberOfTables);
               }
          }  
          return finished(counts, numberOfTables);
      default: 
           i = startIndex;
           index = getDNA(i);
           mask[0] = matchTable(0, index);
           i++;
           if (i > stopIndex) return finished(counts, numberOfTables);
           firstTime1 = true;
           firstTime2 = true;
           firstTime3 = true;
           bool skipSection = false;
           while(true){
            if(skipSection){
               while((mask[3] == 0 || numberOfTables == 4) || firstTime3){
                    firstTime3 = false;

                   while((mask[2] == 0 || numberOfTables == 3) || firstTime2){
                    firstTime2 = false;
                    while((mask[1] == 0 | numberOfTables == 2) || firstTime1){
                        firstTime1 = false;
                        index = getDNA(i);
                        //updateMaskArrayAndIdentifyMatches(1, i, stopIndex, mask);
                   mask[1] = mask[0] & matchTable(1, index);
                   mask[0] = matchTable(0, index);
                   if (numberOfTables == 2 && mask[1]) {
                      identifyMatches(mask[1], i, numberOfTables);
                   }
                               if (i > stopIndex) return finished(counts, numberOfTables);
                        #ifdef STATS
                           counts[1]++;
                        #endif
                     }
                    index = getDNA(i);
                    //updateMaskArrayAndIdentifyMatches(2, i, stopIndex, mask);
                 mask[2] = mask[1] & matchTable(2, index);
                   mask[1] = mask[0] & matchTable(1, index);
                   mask[0] = matchTable(0, index);
                   if (numberOfTables == 3 && mask[2]) {
                      identifyMatches(mask[2], i, numberOfTables);
                   }
                   i++;
                               if (i > stopIndex) return finished(counts, numberOfTables);
 
                    #ifdef STATS
                       counts[2]++;
                    #endif
                    }
                    index = getDNA(i);
                    //updateMaskArrayAndIdentifyMatches(3, i, stopIndex, mask);
               mask[3] = mask[2] & matchTable(3, index);
               mask[2] = mask[1] & matchTable(2, index);
               mask[1] = mask[0] & matchTable(1, index);
               mask[0] = matchTable(0, index);
               if (numberOfTables == 4 && mask[3]) {
                  identifyMatches(mask[3], i, numberOfTables);
               }
               i++;
                           if (i > stopIndex) return finished(counts, numberOfTables);
                #ifdef STATS
                   counts[3]++;
                #endif
               }
               depth = 4;
            }
                   index = getDNA(i);
                   for (j = depth; j > 0; j--) {
                      mask[j] = mask[j-1] & matchTable(j, index);
                   }
                   mask[0] = matchTable(0, index);
                #ifdef STATS
                   counts[depth]++;
                #endif
                   if (mask[depth]) { //either we found a match or stack is growing
                      if (depth == (numberOfTables-1)) { //we found a match
                         identifyMatches(mask[depth], i, numberOfTables);
                         i++;
                         if (i > stopIndex) return finished(counts, numberOfTables);
                         if (mask[depth-1]) continue;
                         depth--;
                         while ((depth ) && (mask[depth-1] == 0)) depth--; 
                         if (depth < 4){ 
                            skipSection = false;
                            firstTime3 = true;
                         }
                         continue;
                      }
                      i++;
                      if (i > stopIndex) return finished(counts, numberOfTables);
                      depth++;  //increase stack depth
                      skipSection = true; 
                      continue;
                   }else{
                       i++;
                       break;
                   }
               }
           }
}

int finished(int counts[], int numberOfTables){
#ifdef STATS
   for (i = 1; i < numberOfTables; i++) {
      printf("  Level %d depth count = %d\n", i, counts[i]);
   }
#endif
#ifdef DEBUG
   printf("Number of times motif found signalled: %d\n", numberOfTimesCalled);
   printf("Number of times motifs found: %d\n", numberOfTimesMotifFound);
#endif
   return av-oldav;  
}
void identifyMatches(LOOKUPWORD mask, //state mask identifying matches
                     unsigned long int index, //index into the genome
                     int depth  //depth of stack when match found
                    ) {
   int statesPerWord = wordWidth / scanWidth;
   LOOKUPWORD extract = (1 <<scanWidth) - 1;
   unsigned long int position;
   int i, j, k, m, n, imax;
   LOOKUPWORD test, bit, testbit;
   unsigned char c;
   int length;
   int isnotok(unsigned char, unsigned char);
   int foundAtLeastOne;
   int whichMotifs;
#define bumpcount(i, j) (summaryData[i * (numberOfChromosomes+1) + j]++)
#ifdef DEBUG
   foundAtLeastOne = 0;
   numberOfTimesCalled++;
#endif
   //determine how many state fields must be examined
   imax = (numberOfMotifs<statesPerWord)? numberOfMotifs:statesPerWord;
   testbit = 1<<(scanWidth-1);
   whichMotifs = 0;
   for (i = 0; i < imax; i++) {
      test = mask & extract;
      if (test) {
         bit = testbit;
         for (j = 0; j < scanWidth; j++) {
             if (bit & mask) { //found bit corresponding to match
                position = (index - depth + 1) * scanWidth + j;
                //test genome starting at position against motif i,
                //i + statesPerWord, i + 2*statesPerWord, etc
                for (k = i; k < numberOfMotifs; k+=statesPerWord) {
                   length = strlen(motifs[k].motif);
                   for (m = 0; m < length; m++) {
                      c = DNAchar(position + m);
                      if (isnotok(motifs[k].motif[m], c)) {
                         goto nextmotif;
                      }
                   }  
                   //Test that there are no "invalid characters" in the genome
                   //where the match was found
                   if (invalid(position, length)) {
                      goto nextmotif;  //if found motif spans invalid characters
                                       //in the genome, ignore the hit
                   }    
                   
                   //confirmed a match
                   if (utr3 | utr5) {
                      m = chrIndex(
                          map[mrnaInfo[whereIndex].geneID].annot_number);
                   } else {
                      //index of where hit found 
                      m = chrIndex(findchrm(position)); 
                   }
                   n= motifs[k].name;
                   if ((whichMotifs & (1<<n)) == 0) {
                      bumpcount(n, m);
                   }
                   whichMotifs |= (1<<n);
#ifdef DEBUG
                   numberOfTimesMotifFound++;
                   foundAtLeastOne = 1;
#endif
                   *av = position;
                   *pv = motifs[k];
                   av++;
                   pv++;
nextmotif:
                  ;
                }  //for k ...
             }  // if (bit & mask) ...
            bit >>= 1; 
         }  //for j ..
      } //if (test) ... 
      extract <<= scanWidth;
      testbit <<= scanWidth;
   } //for i = ...  
#ifdef DEBUG
//   if (foundAtLeastOne) numberOfTimesMotifFound++;
#endif
}

int isnotok(unsigned char x, unsigned char y) {
    //returns true if nucleotide y IS NOT an instance of
    //IUPAC character x.
    switch(x) {
       case 'A': if (y == 'A') return 0; else return 1;
       case 'B': if (y != 'A') return 0; else return 1;
       case 'C': if (y == 'C') return 0; else return 1;
       case 'D': if (y != 'C') return 0; else return 1;
       case 'G': if (y == 'G') return 0; else return 1;
       case 'H': if (y != 'G') return 0; else return 1;
       case 'K': if (y == 'G' || y == 'T') return 0; else return 1;
       case 'M': if (y == 'A' || y == 'C') return 0; else return 1;
       case 'N': return 0;
       case 'R': if (y == 'A' || y == 'G') return 0; else return 1;
       case 'S': if (y == 'C' || y == 'G') return 0; else return 1;
       case 'T': 
       case 'U': if (y == 'T') return 0; else return 1;
       case 'V': if (y != 'T') return 0; else return 1;
       case 'W': if (y == 'A' || y == 'T') return 0; else return 1;
       case 'Y': if (y == 'C' || y == 'T') return 0; else return 1;
       default: return 1;
    };
}

int invalid(unsigned int where, unsigned int length) {
 
   int jj;
   int i;

   if (utr3 | utr5) {
/***
      for (i = 0; i < mrnaCount; i++) {
         if (where < mrnaInfo[i].mrnaStart) continue;
         if (where + length <= mrnaInfo[i+1].mrnaStart) {
            whereIndex = i;
            return 0 ;
         }
      }
      return 1;  
***/
      whereIndex = findSplice(where);
      if (where + length <= mrnaInfo[whereIndex+1].mrnaStart) {
         return 0;
      }
      else return 1;
   }
   if ((where+length)/32 < extable[2*exindex]) return 0;  //all is ok
   while(where/32 < extable[2*exindex]) exindex++;
   if (where/32 == extable[2*exindex])   {
      for (jj = 0; jj < (unsigned)length; jj++)   {
         if (extable[2*exindex+1] & (1<<(31-((where+jj)&31)))) return 1;
         if (((where + jj + 1)&31) == 0) { //finished a word of exceptions
            //maybe no more exception words overlap the pattern
            if (((where + length - 1)/32) < extable[2*exindex+2]) return 0;
            //skip pattern characters which do not overlap next
            //exception word
            while ((where + jj + 1)/32 < extable[2*exindex +2]) {
               jj += 32;
            }
            exindex++;  //continue looking with next word of exceptions
         }
      }
   }
   return 0;
}
 
int exampleOf(char sample[], char pattern[], int length) {
    int i;
    for (i = 0; i < length; i++) {
        if (isnotok(pattern[i], DNAchar(sample+i))) return 0;
    }
    return 1;
}

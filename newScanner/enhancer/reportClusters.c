// Copyright (c) in silico Labs, LLC 2008, 2009, 2012, 2015

#include <stdlib.h>
#include <stdio.h>
#include "scanner.h"
#include "timer.h"
#include <string.h>

#define PETER
#define TIMING

   unsigned long int chrStart(int);
   int counts[27][3];
   unsigned long int start, stop, realPos;
   int limitedCluster(int a);
   int foundUtr3, foundUtr5, foundIntron, foundExon;
   int foundCds, foundOrf;
   int whereIndex;
   int spliceIndex;
   int clusterNumber;
   int findSplice(unsigned int);
   int exampleOf(int, char [], int);
   char *colorName[26] ={"cyan", "yellow", "lightgreen", \
                         "lightpink", "orange", "greenyellow", \
                         "dodgerblue", "mistyrose", "lightsalmon", "khaki"};

void reportClusters() {

   unsigned long int *ax, *aq, *ay, *az, cp, *alast;
   motifDes *px, *pq, *py, *pz, *plast;
   unsigned long int boundary;
   int ref = 0;
   int findchrm(unsigned int);
   void outputSummary();
   int i, j, k, m, group_count, groupsize, sum;
   int isPalindrome;
   int trial, test, motifIndex;
   int ichr, annotation_index;
   char signature[250], longsig[500];
   char decnumber[20], hitlist[5000], line[500];
   int motifStart[500];
   char resultstring[100000];
   unsigned long int endcluster;
   void fullannot(int, unsigned int, unsigned int, char *);
   void displayUTRHit(int, int, int, int, int, char *);
   int newbool(int counts[][3], int, char boolexp[]);
   int annotation_search(int, int);
   char *clusterDNA;
   char lookup[] = "ACGT";
   int clusterSize, chCount, clIndex, motifEnd;
   TIMEVAL time1, time2;
   char searchStrings[1000];

#define PRLINE(s)   if(fromTerminal == 0) both(s); else fprintf(save,"%s",s)
#define bumpcount(i, j) (summaryData[i * (numberOfChromosomes+1) + j]++)

   group_count = 0;
   av = avbase + numberOfMatches;
   searchStrings[0] = 0;
   strcat(searchStrings, motifs[0].motif);
   for (i = 2; i < numberOfMotifs; i+=2) {
      strcat(searchStrings, ",");
      strcat(searchStrings, motifs[i].motif);
   } 

   for (trial = 0; trial < 2; trial++) {
      getTime(&time1);
      group_count = 0;
      for (ax = avbase, px = pvbase; ax < av; ax++, px++) {
         //Compute boundary limit of this potential cluster
         ichr = findchrm(*ax);
         boundary = arm[ichr+1];
         aq = ax;
         pq = px;
         groupsize = 1;
         if (ax + 1 < av || minGroupSize == 1) {
            if (*(ax+1) < *ax + windowSize || minGroupSize == 1) {
               //We may have a window!

               // if hit at *ax crosses boundary, abandon possible cluster
               if (*ax + strlen(px->motif) >= boundary) continue;

               for (ay = ax+1, py = px+1; ay < av; ay++, py++) {
                  if (*ay >= *ax + windowSize) {
                     break;    //stop when ay is outside the window
                  }
                  // abandon adding more hits if this hit crosses boundary
                  if (*ay + strlen(py->motif) >= boundary ) break;
                  groupsize++; //a hit that counts!   
                  aq = ay;
                  pq = py;
               }  //close for ay loop
       
               //Now the total window runs is [ax, ay).
               if (groupsize < minGroupSize) continue;
/*******   for hits near a gene
               if (margin) { //is hit near enough to a gene?
                  dist = locategene(*ax, margin);
                  if (dist == 0) continue;
               }
*******/
               
            //find where the annotations are
            //we already know that potential cluster does
            //not cross contig (or splice for utr search)
   
            annotation_index = annotation_search(ichr, 
                                                    *ax - arm[ichr]/* + 1*/);
   
            for (j = 0; j < 26; j++) {
               counts[j][0] = counts[j][1] = counts[j][2]= 0;
            }
            hitlist[0] = 0;  //start with empty hit list;
            alast = 0;
            for (az = ax, pz = px, i = 0; az < ay; az++, pz++){      
               isPalindrome =  exampleOf(*az, pz->reverse, strlen(pz->reverse));
               if (az != ax && *az < *(alast) + strlen( (plast)->motif)) {
                  continue;
               }
               alast = az; plast = pz;
               endcluster = *az+strlen(pz->motif)-1;
               
               signature[i] = pz->name;
               longsig[2*i] = signature[i];
               i++;
               counts[signature[i-1] - 'A'][0]++;
               if (isPalindrome) {
#ifdef PETER
                  longsig[2*i-1] = '*';
#else
                  longsig[2*i-1] = ' ';
#endif
   
                  counts[signature[i-1] - 'A'][1]++;
                  counts[signature[i-1] - 'A'][2]++;
               }
               else if (pz->direction == '+') {
#ifdef PETER
                  longsig[2*i - 1] = '+';
#else
                  longsig[2*i - 1] = ' ';
#endif
                  counts[signature[i-1] - 'A'][1]++;
               }
               else {
#ifdef PETER
                  longsig[2*i - 1] = '-';
#else
                  longsig[2*i - 1] = ' ';
#endif
                  counts[signature[i-1] - 'A'][2]++;
               }
               motifStart[i-1] = *az - *ax;
               sprintf(decnumber, "%c%ld,", longsig[2*i - 1], *az - *ax);
               strcat(hitlist, decnumber);
              
            } 
            if (i == 0) continue;  //not a cluster
            clusterSize = i;
            signature[i] = 0;
            longsig[2*i] = 0;
            hitlist[strlen(hitlist)-1] = 0;  //remove last character
            //sscanf(hitlist, "%d", &test);
            //if (test != 0) goto continueForAx;
   
            //Now we see whether cluster satisfies cluster-constraints
            if (0 == newbool(counts, 26, boolExp)) {
               continue;  //No. continue outer loop (for ax = avbase, ...)
            }
            if (utr3 | utr5 | orf | cds) {
               start = *ax;
               stop = endcluster;
               //make sure that the item is really a 
               //3'utr or 5'utr or orf or cds
               if (limitedCluster(whereIndex) == 0) continue;
            } else {
               //usual case
               start = *ax - arm[ichr] /* + 1*/;  //start of cluster
               stop = endcluster - arm[ichr] /* + 1 */; //end of cluster
   
               //printf("annotation file %s (%d), starting index is %d\n",
               //        annotName[ichr], ichr, annotation_index);
   
               if (limited) {
                  //test to see if the cluster satisfies the limiting constrants
                  if (limitedCluster(annotation_index) == 0) continue;
               }
            }   
            //Now, that we've eliminated all overlapping hits and hits
            //in the wrong direction for UTR searches, if the remaining
            //number of hits is too small, abandon the cluster.
            if (i < minGroupSize) continue;
            bumpcount((numberOfMotifs/2+1), chrIndex(ichr));
            // bump up cluster count for chromosome
            group_count++;  //we found a cluster! count it.
            clusterNumber = group_count;  //for external use
            if (trial == 0) {
               goto bookkeepping;
            }
 
            //if (count_only != 'C') 
            //The next two commented-out statements are old style output
            //The following 6 executable statements give equivalent output
            //sprintf(resultstring, "%s %s pos %d (length %d) in %s:%s\n", 
            //        longsig, hitlist, (*ax) - arm[ichr], 
            //        endcluster - (*ax) + 1, 
            //        chromoName[ichr], annotName[ichr]);
            //both(resultstring);
            //if (count_only != 'C' && gene_names == 0) 
            //    fprintf(save, "%s  %s pos %d (length %d) in %s:%s\n", 
            //            longsig, hitlist, (*ax) - arm[ichr], 
            //            endcluster - (*ax) + 1, chromo_name[ichr], 
            //            annot_name[ichr]);
            //fprintf(archive, "%x %x %s\n", 
            //                 *ax, *(ay-1) + strlen(*(py-1)) - 1, longsig);
   
            sprintf(resultstring, "\nCluster %d length %ld bp\n",
                    group_count, endcluster - (*ax) + 1);
            both(resultstring);
            sprintf(resultstring, "%s	MERMER	Cluster	%lu	%lu	.	.	.	Cluster%d search=%s\n", 
                chromoName[ichr], start, stop, clusterNumber, searchStrings);
            fprintf(tbl, "%s", resultstring);
            both(resultstring);
            sprintf(resultstring, 
               "MOTIFS %s at positions %s\n", longsig, hitlist);
            both(resultstring);
            if (utr3 | utr5 | orf | cds) { //special treatment for searches 
                               //restricted to mRNAs.
               displayUTRHit(annotation_index, whereIndex, 
                             spliceIndex, start, stop, resultstring);
            } else {
               //Normal case
               sprintf(resultstring, 
                  "On chromosome %s at position %lu (%s:%s)\n",
                  chromoName[ichr], (*ax) - chrStart(ichr), chromoName[ichr],
                  annotName[ichr]); 
            }
            both(resultstring);
               
            if (annotation_index>= 0 && (utr3 | utr5 | orf | cds) == 0) {
               //printf("annotation index = %d\n", annotation_index);
               //fprintf( save, "annotation index = %d\n", annotation_index);
              
               //get the full annotation info into resultstring
               primeNames[0] = 1;  //list of gene names involved in hit
               
               fullannot(annotation_index, 
                 *ax - arm[ichr] /*+ 1 *//*- map[annotation_index].start*/,
                 endcluster - arm[ichr]/* + 1*/ /*-map[annotation_index].start*/,
                 resultstring); 
              
               //if (count_only != 'C') {
   
               both(resultstring);
               both("\n");
               //if (gene_names) { //option to print gene name synonyms
   
               //this if-statement was used to print all synonyms after each
               //found cluster. No longer used
               //if (strlen(resultstring)) {
               //   fprintf(save, 
               //      "Here are the genes and their synonyms associated" 
               //      " with this cluster\n");
               //   printallsynonyms(resultstring, synCount);
               //}
   
               //}
               //else {
               //   fprintf(save, "%s\n\n", resultstring);
               //}
               //}
            }
            else {
               // full_results = 0;
            }
            //if (full_results && gene_names == 0) {
            //   fprintf(save, "          Window start: Annot %s, pos %d.", 
            //           annot_name[ichr], *ax - arm[ichr] - ref + 1);
            //   fprintf(save, "\n");
            //}
               
            //if (full_results) {
            //   //fetch entire cluster
if (fromTerminal) goto bookkeepping;
                 if (fromTerminal == 0) {
                    printf("<div id='content%d' style='display: none'>", 
                       clusterNumber);
                    printf("<pre><font size='3'>\n");
                    printf("<font color='grey'>");
                 }
                 start = *ax - 200;
                 stop = endcluster + 200;
                 if (!(utr3 || utr5 || orf || cds)) {
                    if (start < arm[ichr]) start = arm[ichr];
                    if (stop >= arm[ichr+1]) stop = arm[ichr+1] - 1;
                 }
                 clusterDNA = malloc(stop - start + 2);
                 for (cp = start, i = 0; cp <= stop; cp++, i++) {
                    clusterDNA[i] = DNAchar((cp)) | 0x20;  //force lower case
 //                   if (utr3 | utr5 && clusterDNA[i] == 't') {
 //                       clusterDNA[i] = 'u';
 //                   }
                 }
                 j = i;
                 clusterDNA[i] = 0;
                 for (az = ax, pz = px; az < ay; az++, pz++) {
                    if (az == ax || *az >= az[-1] + strlen(pz[-1].motif)) {
                       for (cp = *az, i = 0; 
                            i < (int)strlen(pz->motif); cp++, i++) {
                          clusterDNA[cp - start ] &= 0xDF;  //make upper case
                       }
                    }
                 }
                  
                 chCount = start - *ax;
                 clIndex = 0;
                 line[0] = 0;
                 motifEnd = j + 1;
                 motifStart[clusterSize] = j+1;
                 for (i = 0; i < j; i+=10) {
                    for (k = 0; k < 10; k++) {
                       if (   fromTerminal == 0 
                           && chCount == motifStart[clIndex]) {
                          //When encountering first motif, end grey text
                          //and make subsequent characters bold 
                          if (strlen(line)) {
                             PRLINE(line);
                             line[0] = 0;
                          }
                          if (clIndex == 0) {
                             strcat(line, "</font><b>");
                          }
                          //determine where motif info is
                          for (motifIndex = 0; motifIndex < numberOfMotifs;
                               motifIndex+=2) {
                             if (motifs[motifIndex].name == signature[clIndex]) {
                                break;
                             }
                          }
                          if (motifIndex >= numberOfMotifs) motifIndex = 0;
                          //highlight motif in appropriate color
                          strcat(line, "<FONT style='BACKGROUND-COLOR:");
                          strcat(line, colorName[motifs[motifIndex].name-'A']);
                          strcat(line, "'>");
                          motifEnd = motifStart[clIndex] + 
                             strlen(motifs[motifIndex].motif)
                             - 1;
                          clIndex++;
                          printf("%s", line);
                          line[0] = 0;
                       }
                       m = strlen(line);
                       line[m] = clusterDNA[i+k];
                       line[m+1] = 0;
                       if (    fromTerminal == 0
                            && chCount == motifEnd) {
                          //turn off the highlighting
                          if (strlen(line)) {
                             PRLINE(line);
                             line[0] = 0;
                          }
                          strcat(line, "</FONT>");
                          if (clIndex == clusterSize) {
                             //end of cluster output -- end bold 
                             //and return to grey font
                             strcat(line, "</b><font color='grey'>");
                          }
                          printf("%s", line);
                          line[0] = 0;
                       }
                       chCount++;

                       if (clusterDNA[i+k] == 0) break;
                    } //for k
                    
   
                    PRLINE(line);
                    PRLINE(" ");
                    line[0] = 0;
   
                    if (i%60 == 50) {
                       PRLINE("\n");
                    }
                 } //for i
                 fprintf(save, "\n");
                 free(clusterDNA);
                 if (fromTerminal == 0) {
                    printf("</font></pre>"); 
                    printf("</div>");
                 }
                           
                 both("\n");
   
            //}
bookkeepping:
               ax = aq;
               px = pq;
            }
         }
continueForAx: ;
      }
      if (trial == 0) {
         outputSummary();
         sprintf(line, "Number of clusters found: %d\n\n", group_count);
         both(line);
      }
#ifdef TIMING
      getTime(&time2);
      if (trial == 0) {
         sprintf(line, "Time to count hits is %d ms\n",
            getDiffMillisecs(&time1, &time2));
         printBoth = 0;  //suppress further printing to .txt file
         printNone = 1;  //suppress printing further results
      } else {
         sprintf(line, "Time to generate output is %d ms\n",
            getDiffMillisecs(&time1, &time2));
      }
      both(line);
#endif
   }
   
}

int limitedCluster(int a) {
   int altArg;
   int sum, i, j, k;
   unsigned int left, right;
   int lftt, rhtt;

   //Note: a is an index into the gene annotation array, called map.
   //First, see if there are enough hits in the correct direction
   foundUtr3 = foundUtr5 = foundIntron = foundExon = 0; //nothing found yet
   foundOrf = foundCds = 0;
   if (utr3 | utr5 | orf | cds) {
      altArg = a;
      a = mrnaInfo[a].geneID;
      k = altArg;
      while (k >= 0 && mrnaInfo[k].geneID == mrnaInfo[altArg].geneID) k--;
      k = altArg - k - 1;  //index into standard mrna tables (0-based);
      spliceIndex = k;
   }
   if (map[a].direction == '+') j = 1; else j = 2;
   for (sum = 0, i = 0; i <  26; i++) sum += counts[i][j];
   if (sum < minGroupSize) return 0; //for (ax = avbase ...

   if (cds) {
      foundCds = 1;
      goto summary;
   }
   if (exon || intron) goto testIntrons;
   left = right = 0;
   // Compute size of right UTR
   j = 2*exonCount(a,k) + 2;
   while(element(a, k, j) > element(a, k, 3)) {
      right = right + element(a, k, j+1) - element(a, k, j) + 1; 
      j -= 2;
   }
   right = right + element(a, k, j+1) - element(a, k, 3);
   // Compute size of left UTR
   j = 5;
   while (element(a, k, j) < element(a, k, 2)) {
      left = left + element(a, k, j) - element(a, k, j-1) + 1;
      j += 2;
   }   
   left = left + element(a, k, 2) - element(a, k, j-1) ;

   if (orf) {
      if (map[a].direction == '+') {
         lftt = utr5;
         rhtt = utr3;
      } else {
         lftt = utr3;
         rhtt = utr5;
      }
      if ((lftt || start - mrnaInfo[altArg].mrnaStart > left) &&
          (rhtt || mrnaInfo[altArg+1].mrnaStart - right > stop)) foundOrf = 1; 
   }

   if (utr3) {
      if (map[a].direction == '+' &&
          mrnaInfo[altArg+1].mrnaStart - start <= right) foundUtr3 = 1;
      else if (map[a].direction == '-' &&
          stop - mrnaInfo[altArg].mrnaStart + 1 <= left) foundUtr3 = 1;
   }

testUTR5:
   if (utr5) {
      if (map[a].direction == '+' &&
          stop - mrnaInfo[altArg].mrnaStart + 1 <= left) foundUtr5 = 1;
      else if (map[a].direction == '-' &&
          mrnaInfo[altArg+1].mrnaStart - start <= right) foundUtr5 = 1;
   }
   goto summary;

testIntrons:
   if (intron) {
      for (i = 0; i < map[a].numberOfmRNAs; i++) {
         for (j = 0; j < exonCount(a,i) - 1; j++) {
            if (map[a].mrnaInfo[i*rowSize(a) + 2*j + 5] < start &&
                map[a].mrnaInfo[i*rowSize(a) + 2*j + 6] > stop) {
               foundIntron = (map[a].direction == '+') ?
                              j+1 : (exonCount(a,i)) - j - 1;
               goto summary;
            }
         }
      }
   }
testExons:
   if (exon) {
      for (i = 0; i < map[a].numberOfmRNAs; i++) {
         for (j = 0; j < exonCount(a,i); j++) {
            if (map[a].mrnaInfo[i*rowSize(a) + 2*j + 4] <= start &&
                map[a].mrnaInfo[i*rowSize(a) + 2*j + 5] >= stop) {
               foundExon = (map[a].direction == '+') ? 
                           j+1 : exonCount(a,i) - j;
               goto summary;
            }
         }
      }
   }

summary:
   return foundUtr3 | foundUtr5 | foundIntron | foundExon | foundOrf | foundCds;
}

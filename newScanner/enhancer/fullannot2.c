//Copyright (c) World Internet Productions, LLC, 1999
//Copyright (c) in silico Labs, LLC, 2001
//modified by in silico Labs,LLC 10/4/01 
//Copyright (c) in silico Labs, LLC 2009
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "scanner.h"

void plusmsg(char message[], int i, int pos, int posend);
void minusmsg(char message[], int i, int pos, int posend);
char * formatNumber(int);

void makerevcompl(char *from, char *to, int num);
/***********
typedef   struct ann_index {
   int annot_number;
   unsigned int start;
   unsigned int stop;
   char direction;
   char *gene_name;
   char touched;
   } ann_index;
ann_index *map;

int mapSize;
int geneName;
int margin;
***********/

int margin;
int special;
int primeNames[100];
 
//#define recordHit(i) primeNames[primeNames[0]] = i; primeNames[0]++;
#define recordHit(i) map[i].touched = 1;
 
void outputGraphicBetween(int, int, int);
void outputGraphicInternal(int, int, int);

void fullannot (int geneindex, unsigned int pos, unsigned int posend, 
                char *results) {
   
   //map is a condensation of the annotations -- one entry per gene
   //geneindex is an index into map, presumably downstream of the found cluster
   //pos is the first position in the cluster
   //posend is the last position in the cluster
   //results is the character string that describes where the cluster was found
   
   int i, j;  
   char message[300], mess1[300], extramessage[1000];

   special = 0;
   margin = 0;   
   
   //initialize output strings to empty
   //message[0] = mess1[0] = extramessage[0] = 0;
   strcpy(message, ""); strcpy(mess1, ""); strcpy(extramessage, "");
    
   //Search for the area of the cluster only within the same chromosome
   //i.e. map[geneindex].annot_number is the chromosome number.
   for (i = geneindex; 
   i < mapSize && map[i].annot_number == map[geneindex].annot_number;
   i++) {
      if (map[i].start > posend) break; //position i is beyond cluster
      if (map[i].start <= pos && map[i].stop >= posend) {
         //it is in gene at map[i]
         i++;
         //now it is in gene at map[i-1]
         break;
      }
      if (map[i].stop > posend) {
         //in this case cluster overlaps gene at position map[i]
         break;
      }
   }
   if (i == mapSize ||
      map[i].annot_number != map[geneindex].annot_number) {
      //Cluster is after all genes or in last gene of cluster.
      if (fromTerminal == 0) {
         if (map[i-1].start <= pos && posend <= map[i-1].stop) {
            outputGraphicInternal(i-1, pos, posend);
         } else {
            outputGraphicBetween(i, pos, -posend);
         }
      }
      if (map[i-1].direction == '+') {
         plusmsg(message, i-1, pos, posend);
      }
      else {
         minusmsg(message, i-1, pos, posend);
      }
      sprintf(results, "\t%s", message);
      return;
   }

   if (i == geneindex || posend > map[i-1].stop) {
      // Normal case: cluster is between two genes
      if (fromTerminal == 0) outputGraphicBetween(i, pos, posend);
      if (i > geneindex) {
         if (map[i-1].direction == '+') plusmsg(mess1, i-1, pos, posend);
         else minusmsg (mess1, i-1, pos, posend);
      } else mess1[0] = 0;  //null string
      if (map[i].direction == '+') plusmsg(message, i, pos, posend);
      else minusmsg(message, i, pos, posend);
      if (strlen(mess1)) {
         if (special) {
            if (strlen(message)) {
               sprintf(results, "%s\n%s", mess1, message);
            }
            else
               sprintf(results, "%s", mess1);
         }
         else {
            sprintf(results, "\t%s,\n\t%s", mess1, message);
         }
      }
      else {
         if (special) {
            sprintf(results, "%s", message);
         }
         else {
             sprintf(results, "\t%s", message);
         }
      }
   }

   else {
      // previous gene overlaps cluster
      if (fromTerminal == 0) outputGraphicInternal(i-1, pos, posend);
      if (map[i-1].direction == '+') plusmsg(message, i-1, pos, posend);
      else minusmsg(message, i-1, pos, posend);
      if (special) {
         sprintf(results, "%s", message);
      }
      else {
         sprintf(results, "\t%s", message);
      }
   }

   for (j = i+1; map[j].annot_number == map[i].annot_number; j++) {
      if(map[j].start > map[i-1].stop) break;
      if (map[j].start > posend && map[j-1].start <= pos ||
         map[j].start <= pos && map[j].stop >= posend ||
         map[j].start <= posend && map[j].stop > posend)
      {
         fullannot(j, pos, posend, extramessage);
         strcat(results, "\n");
         strcat(results, extramessage);
      }
   }
}

void plusmsg(char message[], int i, int pos, int posend) {
   char plusMsg[100];
   char temp[20];
#ifdef PETER
   strcpy(plusMsg, "(+)");
#else
   strcpy(plusMsg, "");
#endif
   if (limited) {
      if (foundUtr3) strcat(plusMsg, " in the 3'utr");
      if(foundUtr5) strcat(plusMsg, " in the 5'utr");
      if(foundIntron) {
         sprintf(temp, " in intron %d", foundIntron);
         strcat(plusMsg, temp);
      }
      if(foundExon) {
         sprintf(temp, " in exon %d", foundExon);
         strcat(plusMsg, temp);
      }
      foundUtr3 = foundUtr5 = foundIntron = foundExon = 0;
   }
   if (posend < map[i].start) {
      if (special) {
         if (margin == 0 || map[i].start - pos <= margin) {
            sprintf(message, "%s", map[i].geneName);
         }
      }
      else {
         recordHit(i);
         sprintf(message, "%d bp upstr of %s%s",
            map[i].start - posend, map[i].geneName, plusMsg);
      }
   }
   else if (pos < map[i].start) {
      if (special) {
         sprintf(message, "%s", map[i].geneName);
      }
      else {
         recordHit(i);
         sprintf(message, "%d bp upstr of %s%s and overlaps",
            map[i].start - pos, map[i].geneName, plusMsg);
      }
   }
   else if (posend <= map[i].stop) {
      if (special) {
         sprintf(message, "%s", map[i].geneName);
      }
      else {
         recordHit(i);
         sprintf(message, "%d bp inside and totally within %s%s",
            pos - map[i].start, map[i].geneName, plusMsg);
      }
   } 
   else if (pos <= map[i].stop) {
      if (special) {
         sprintf(message, "%s", map[i].geneName);
      }
      else {
         recordHit(i);
         sprintf(message, "%d bp from end of %s%s and overlaps", 
            map[i].stop - pos, map[i].geneName, plusMsg); 
      }
   }
   else {
      if (special) {
         if(margin == 0 || posend - map[i].stop <= margin) {
            sprintf(message, "%s", map[i].geneName);
         }
      }
      else {
         recordHit(i);
         sprintf(message, "%d bp dnstr of %s%s",
            pos - map[i].stop, map[i].geneName, plusMsg);
      }
   }
}

void minusmsg(char message[], int i, int pos, int posend) {
   char minusMsg[100];
   char temp[20];
#ifdef PETER
   strcpy(minusMsg, "(-)");
#else
   strcpy(minusMsg, "");
#endif
   if (limited) {
      if (foundUtr3) strcat(minusMsg, " in the 3'utr");
      if(foundUtr5) strcat(minusMsg, " in the 5'utr");
      if(foundIntron) {
         sprintf(temp, " in intron %d", foundIntron);
         strcat(minusMsg, temp);
      }
      if(foundExon) {
         sprintf(temp, " in exon %d", foundExon);
         strcat(minusMsg, temp);
      }
      foundUtr3 = foundUtr5 = foundIntron = foundExon = 0;
   }
   if (posend < map[i].start) {
      if (special) {
         if (margin == 0 || map[i].start - pos <= margin) {
            sprintf(message, "%s", map[i].geneName);
         }
      }
      else {
         recordHit(i);
         sprintf(message, "%d bp dnstr of %s%s",
            map[i].start - posend, map[i].geneName, minusMsg);
      }
    }
   else if (pos < map[i].start) {
      if (special) {
         sprintf(message, "%s", map[i].geneName);
      }
      else {
         recordHit(i);
         sprintf(message, "%d bp from end of %s%s and overlaps",
            posend - map[i].start, map[i].geneName, minusMsg);
      }
   }
   else if (posend <= map[i].stop) {
      if (special) {
         sprintf(message, "%s", map[i].geneName);
      }
      else {
         recordHit(i);
         sprintf(message, "%d bp inside and totally within %s%s",
            map[i].stop - posend, map[i].geneName, minusMsg);
      }
   } 
   else if (pos <= map[i].stop) {
      if (special) {
         sprintf(message, "%s", map[i].geneName);
      }
      else {
         recordHit(i);
         sprintf(message, "%d bp upstr of %s%s and overlaps", 
            posend - map[i].stop, map[i].geneName, minusMsg); 
      }
   }
   else {
      if (special) {
         if (margin == 0 || posend - map[i].stop <= margin) {
            sprintf(message, "%s", map[i].geneName);
         }
      }
      else {
         recordHit(i);
         sprintf(message, "%d bp upstr of %s%s",
            pos - map[i].stop, map[i].geneName,  minusMsg);
      }
   }
}

void outputGraphicBetween(int i, int pos, int posend) {
    //Graphic has cluster between two genes
    int distance, leftPart, centerPart, rightPart;
    int leftStart, rightStop;
    int leftOK, rightOK;
 
    rightOK = 1;
    if (posend < 0) {
       rightOK = 0;
       leftOK = 1;
       posend = -posend;
    }
    else {
       leftOK = map[i-1].annot_number == map[i].annot_number;
    }
    distance = ((map[i].start>posend)?map[i].start:posend) - 
               ((map[i-1].stop<pos)?map[i-1].stop:pos) + 1;
    leftPart = (pos - map[i-1].stop)*300.0/distance;
    if (pos <= map[i-1].stop) {
       leftPart = -leftPart;
       leftStart = map[i-1].stop + 1;
    } else {
       leftStart = pos;
    }

    if (!leftOK) {
       //No left neighbor in the same contig
       distance = map[i].start + 1;
       leftPart = pos*300.0/distance;
       leftStart = pos;
    }
    rightPart = (map[i].start - posend)*300.0/distance;
    if (map[i].start <= posend) {
       rightPart = -rightPart;
       rightStop = map[i].start - 1;
    } else {
       rightStop = posend;
    }
   
    if (!rightOK) { //special case -- cluster at end of a contig
       distance = arm[map[i].annot_number] - arm[map[i-1].annot_number] - map[i-1].stop;
       leftPart = (pos - map[i-1].stop)*300.0/distance;
       rightPart = (distance - posend + map[i-1].stop)*300.0/distance;
       leftStart = pos;
       rightStop = posend;
    }
    
    centerPart = (300 - rightPart - leftPart);
    printf("<table cellpadding='0' cellspacing='0' border='0'>\n");
//Row 1
    printf("<tr>\n");
    printf("<td width='100' height='10' align='right'> \n"); 
    if (pos<= map[i-1].stop && leftOK) {
       printf("</td> <td width = '%d' align='right'> \n", leftPart);
    }
    if (map[i-1].direction == '+' && leftOK) {
       printf("<font size='1' align 'right'><br>3' end</font>");
    }
    printf("</td>\n");
    if (pos > map[i-1].stop || !leftOK) {
       printf("<td width='%d'></td>\n", leftPart);
    } 
    printf("<td width='%d'></td>\n", centerPart);
    printf("<td width='%d'></td>\n", rightPart);
    printf("<td width='100'>\n");
    if (map[i].direction == '-' && rightOK) {
       printf("<font size='1' align='left'><br>3' end</font>\n");
    }
    printf("</td></tr>\n");

//Row2
    printf("<tr>\n");
    printf("<td align='right'>");
    if (pos <= map[i-1].stop && leftOK) {
       printf("</td><td>\n");
    }
    if (leftOK) {
       printf("%s\n", map[i-1].geneName);
       if (map[i-1].direction == '-') {
          printf(" <img src='../graphics/leftArrow.gif' width='10' height='15' /> </td>\n");
       } else {
          printf("<big>&#8594&#124</big></td>\n");
       }
    } else printf("</td>\n");
    if(pos > map[i-1].stop || !leftOK) {
       printf("<td></td>\n");
    }
    printf("<td align='center'> %s</td>\n", 
       formatNumber(rightStop - leftStart + 1));
    if (posend < map[i].start || !rightOK) {
       printf("<td> </td>\n");
    }
    if (rightOK) {
       if (map[i].direction == '+') {
          printf("<td> <img src='../graphics/rightArrow.gif' width='10' height='15' align='left'/> \n");
       } else {
          printf("<td><big>&#124&#8592</big>\n");
       }
       printf("%s</td>\n", map[i].geneName);
    }
    else printf("<td></td>\n");
    if (posend >= map[i].start && rightOK) {
       printf("<td></td>\n");
    }
    printf("</tr>\n");

//Row 3
    printf("<tr>\n");
    if (leftOK) printf("<td height='5' bgcolor='red'/>\n");
    else printf("<td height='5' bgcolor='white'/>\n");
    if (!leftOK || pos > map[i-1].stop) { //cluster starts after left gene
       printf("<td bgcolor='aqua'/>\n");
    } else {
       printf("<td bgcolor='navy'.>\n");
    }
    printf("<td bgcolor='navy'/>\n");
    if (!rightOK || map[i].start>posend) {  //cluster ends before right gene
       printf("<td bgcolor='aqua'/>\n");
    } else {
       printf("<td bgcolor='navy'.>\n");
    }
    if (rightOK) printf("<td bgcolor='red'/>\n");
    else printf("<td bgcolor='white'/>\n");
    printf("</tr>\n");

//Row 4:
    printf("<tr>\n");
    printf("<td></td>\n");
    if (leftOK) {
       if (pos > map[i-1].stop) {
          printf("<td align='center'> %s</td>\n", formatNumber(pos - map[i-1].stop));
       } else {
          printf("<td align='center'> %s</td>\n", formatNumber(map[i-1].stop - pos+1));
       }  
    } else {
       printf("<td align='center'> %s</td>\n", formatNumber(pos));
    }
    printf("<td> </td>\n");
    if (!rightOK) {
       printf("<td align='center'> %s</td>\n", formatNumber(distance - posend + map[i-1].stop));
    }
    else {
       if (posend < map[i].start) {
          printf("<td align='center'> %s</td>\n", formatNumber(map[i].start - posend));
       } else {
          printf("<td align='center'> %s</td>\n", formatNumber(posend - map[i].start + 1));
       }
    }
    printf("<td></td>\n");
    printf("</tr>\n");

    printf("</table>\n");
    return;
}
 
#define MRNA(i, j, k)  map[i].mrnaInfo[j*(2*map[i].numberOfExons+4)+k]

void outputGraphicInternal(int i, int pos, int posend) {
    //output graphic when cluster is contained within a gene
    int distance, leftPart, centerPart, rightPart;
    int leftStart, rightStop;
    int n, k, this, sofar, used, spaces;
    int totalLength;

    distance = map[i].stop - map[i].start + 1;
    leftPart = (pos - map[i].start)*500.0 / distance;
    rightPart = (map[i].stop - posend)*500.0 / distance;
    centerPart = 500 - leftPart - rightPart;
    printf("<table cellpadding='0' cellspacing='0' border='0'>\n");

//Row 1
    printf("<tr>\n");
    printf("<td width = '250' align='left'>\n");
    if (map[i].direction == '-') {
       printf("<br><font size='1' align 'left'>3' end</font>");
    }
    printf("</td>\n");
    printf("<td width='250' align = 'right'>\n");
    if (map[i].direction == '+') {
       printf("<br><font size='1' align='right'>3' end</font>\n");
    }
    printf("</td>\n");
    printf("</tr>\n");
    
//Row 2
    printf("<tr>\n");
    if (map[i].direction == '+') {
       printf("<td> <img src='../graphics/rightArrow.gif' width='10' height='15' align='left'/>%s </td>\n",map[i].geneName);
    } else {
       printf("<td align='left'><big>&#124&#8592</big>%s </td>\n", map[i].geneName);
    }
    if (map[i].direction == '-') {
       printf("<td align='right'>%s <img src='../graphics/leftArrow.gif' width='10' height='15' /> </td>\n", map[i].geneName);
    } else {
       printf("<td align='right'>%s<big>&#8594&#124</td>\n", map[i].geneName);
    }
    printf("</tr></table>\n");  //end table for top of cartoon


   

//Row3
    for (n = 0; n < map[i].numberOfmRNAs; n++) {
       //start table for actual gene cartoons
       printf("<table cellpadding='0' cellspacing='0' border='1'><tr>\n");
       spaces = 2*MRNA(i, n, 1) - 1;
       totalLength = MRNA(i, n, spaces + 4) - MRNA(i, n, 4);
       used = 0; 
       for (k = 1; k <= spaces; k++) {
          sofar = (MRNA(i, n, 4+k) - MRNA(i, n, 4)) * 500 / totalLength;
          this = sofar - used;
          if (this == 0) this++;
          if (k & 1) { //exon
             printf("<td width='%d' height='5' bgcolor='red'></td>\n", this);
          } else { //intron
             printf("<td width='%d' height='5' bgcolor='pink'></td>\n", this);
          }
          used = used + this; 
       } 
       printf("</tr>\n");
    }
    printf("</table>\n");

//Row4
    printf("<table cellpadding='0' cellspacing='0' border='0'>");
    printf("<tr>\n");
    printf("<td width='%d' align='left'> %s</td>\n",
        leftPart, formatNumber(pos - map[i].start));
    printf("<td width='%d' height='5' bgcolor='blue'></td>\n", centerPart);
    printf("<td width='%d' align='right'> %s<\td></tr>", 
        rightPart, formatNumber(map[i].stop - posend));

//Row 5
    printf("<td rowspan='2' align='right'>%s<\td><td></td></tr>\n",
        formatNumber(posend - pos + 1));
    printf("</table>\n");
    return;
}

char * formatNumber(int x) {
   static char result[50];
   if (x < 1000) {
      sprintf(result, "%dbp", x);
   }
   else if (x < 1000000) {
      sprintf(result, "%.1fkb", x/1000.0);
   }
   else sprintf(result, "%.1fMb", x/1000000.0);
   return result;
}


//Copyright (c) in silico Labs, LLC 2009

#include <stdlib.h>
#include <stdio.h>
#include "scanner.h"
#include "timer.h"
#include <float.h>
#include <string.h>

#define UTR 0
#define CDS 1

void display(int, int);
extern char fileName[];

void displayUTRHit(int gene, int where, int splice, int start, int stop, 
                   char *output) {
    int m, realPos, length, i, k, total, untrans;
    int untransRight;
    int spaces, used, sofar, this, piece, piece2, usedDNA;


    m = chrIndex(map[gene].annot_number);
    if (map[gene].direction == '+') {
       realPos = start  - mrnaInfo[where].mrnaStart
                 + element(gene, splice, 4) -
                 map[gene].start;
    } else {
       realPos = mrnaInfo[where+1].mrnaStart - stop 
                 + map[gene].stop   
                 - element(gene, splice, 
                      3+2*exonCount(gene, splice));
    }
    sprintf(output,
       "On chromosome %s in gene %s(%c), %u bp upstream ", 
        chrName[m], map[gene].geneName, map[gene].direction, realPos); 
    both(output);
    if (foundUtr3) { 
       both ("in the 3'utr\n");
    }
    if (foundUtr5) {
       both ("in the 5'utr\n");
    }
    if (fromTerminal != 1){  //no graphic output needed

        //first, calculate the length of the gene (only exons and any 
        //untranscribed regions at either end are counted)

        spaces = 2*exonCount(gene, splice);
        length = element(gene, splice, 4) - map[gene].start ;
        untrans = length;
        for (i = 0; i < spaces; i+=2) {
           length += (element(gene, splice, i + 5) 
                      - element(gene, splice, i + 4) + 1);
        }
        untransRight = (map[gene].stop - element(gene, splice, spaces + 3));
        length += untransRight;
    

    //Row 1
        printf("<table cellpadding='0' cellspacing='0' border='0'>\n");
        printf("<td width='250' align='left'></td>\n");
        printf("<td width='250' align='right'>");
        printf("<font size='1' align='right'>3' end</font></td>\n"); 
        printf("</tr>\n");

    //Row 2
        printf("<tr>\n");
        printf("<td> <img src='../graphics/rightArrow.gif' width='10' height='15' align='left'/>%s </td>\n", map[gene].geneName);
        printf("<td align='right'>%s<big>&#8594&#124</td>\n", map[gene].geneName);
        printf("</tr></table>\n");   //end table for top of cartoon

    //Row 3
        //start actual cartoon
        printf("<table cellpadding='0' cellspacing='0' border='0'> <tr>\n");
        used = 0;

        if (map[gene].direction == '+') {
           //render a gene on the + strand
           //if mrna starts downstream of the gene start, show a blank area
           if (untrans) {
              used = (untrans) * 500 / length;
              if (used == 0) used++;
              printf ("<td width='%d' height='5'></td>\n", used);
           }
           usedDNA = untrans;
           for (k = 0; k < spaces; k+=2) {
              usedDNA = usedDNA + element(gene, splice, k+5)  
                                - element(gene, splice, k+4) + 1;
              sofar = usedDNA * 500 / length;
              this = sofar - used;
              if (this == 0) this++; 
              //test for left UTR   
              if (element(gene, splice, 2) > element(gene, splice, k + 5)) {
                 display(this, UTR);
              }
              else if (element(gene, splice, 3) < element(gene, splice, k + 4)) {
                 display(this, UTR);
              }
              else if (element(gene, splice, 2) < element(gene, splice, k + 4) &&
                       element(gene, splice, 3) > element(gene, splice, k + 5)) {
                 //exon is a CDS exon
                 display(this, CDS);
              }
              else {
                 //exon contains CDS and and UTR info
                 total = 0;
                 //compute size of left UTR piece, if any
                 piece = (element(gene, splice, 2) - element(gene, splice, k+4))
                         * 500 / length;
                 if (piece == 0) piece++;
                 if (piece < 0) piece = 0;
                 //compute size of right UTR piece, if any
                 piece2 = (element(gene, splice, k+5) - element(gene, splice, 3))
                          * 500 / length;
                 if (piece2 == 0) piece2++;
                 if (piece2 < 0) piece2 = 0;
                 if (piece) display(piece, UTR);
                 //compute size of CDS portion (make sure it is non-zero
                 while (this - piece - piece2 <= 0) this++;
                 display(this - piece - piece2, CDS);
                 if (piece2) display(piece2, UTR);
              }
              used = used + this;
           }
        }  else {
           //render a gene on the - strand
           //if mrna starts downstream of the gene start, show a blank area

           if (untransRight > 0 ) {
              used = untransRight * 500 / length;
              if (used == 0) used++;
              printf ("<td width='%d' height='5'></td>\n", used);
           }
           usedDNA = untransRight;
           for (k = spaces - 2; k >= 0; k-=2) {
              usedDNA = usedDNA + element(gene, splice, k+5)
                                - element(gene, splice, k+4) + 1;
              sofar = usedDNA*500/length;
              this = sofar - used;
              if (this == 0) this++; 
              //test for left UTR   
              if (element(gene, splice, 2) > element(gene, splice, k + 5)) {
                 display(this, UTR);
              }
              else if (element(gene, splice, 3) < element(gene, splice, k + 4)) {
                 display(this, UTR);
              }
              else if (element(gene, splice, 2) < element(gene, splice, k + 4) &&
                       element(gene, splice, 3) > element(gene, splice, k + 5)) {
                 //exon is a CDS exon
                 display(this, CDS);
              }
              else {
                 //exon contains CDS and and UTR info
                 total = 0;
                 //compute size of right UTR piece, if any
                 piece = (element(gene, splice, 2) - element(gene, splice, k+4))
                         * 500 / length;
                 if (piece == 0) piece++;
                 if (piece < 0) piece = 0;
                 //compute size of left UTR piece, if any
                 piece2 = (element(gene, splice, k+5) - element(gene, splice, 3))
                          * 500 / length;
                 if (piece2 == 0) piece2++;
                 if (piece2 < 0) piece2 = 0;
                 if (piece2) display(piece2, UTR);
                 //compute size of CDS portion (make sure it is non-zero
                 while (this - piece - piece2 <= 0) this++;
                 display(this - piece - piece2, CDS);
                 if (piece) display(piece, UTR);
              }
              used = used + this;
           }
        }
        printf("</tr></table>\n");
    //Row 4   Show the cluster:
        printf("<table cellpadding='0' cellspacing='0' border='0'><tr>");
        printf("<td width='%d'></td>\n", realPos*500/length);
        piece = (stop - start + 1)*500/length;
        if (piece < 2) piece = 2;
        printf("<td height='3' bgcolor='navy'>\n");
        printf("<A id=control%d' HREF=\"javascript:toggle('content%d');\">",
            clusterNumber, clusterNumber);
    //    printf("<A HREF='/cgi-bin/showcluster.cgi?FN=%s&CN=%d'>",
    //        fileName, clusterNumber);
        printf("<img height='5' width='%d'></A>", piece);
        printf("</td>\n");
        printf("</tr></table>\n");
    }
    map[gene].touched = 1;
    sprintf(output, "\n");

}

void display(int x, int type) {
    
    if (type == UTR) {
       printf("<td width='%d' height='5' bgcolor='lightsalmon'></td>\n", x);
    } else {
       printf("<td width='%d' height='5' bgcolor='Crimson'></td>\n", x);
    }
}

// Copyright (c) in silico Labs, LLC 2008, 2015

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "scanner.h"
void recursiveEnter(char* input, int tablenumber,
                    LOOKUPWORD mask, LOOKUPWORD *matches) {
   int i, n;
   int index;
   char u[64];
   char *s, *t;
   int patLen;

   n = strlen(input);

   for (i = 0, t = input; i < n; i++) u[i] = *t++; //copy line into u
   u[i] = 0;
   for (i = 0; i < n; i++) if (u[i] != 'N') goto work;
   //special case: mask every location
   for (i = 0; i < 1<<(2*scanWidth); i++) {
      matchTable(tablenumber, i) |= mask;
   }
   return;
work:
   for (i = 0; i < n; i++) {
      switch (u[i]) {
         case 'B':         //C or G or T
            u[i] = 'C';
            recursiveEnter(u, tablenumber, mask, matches);
            u[i] = 'G';
            recursiveEnter(u, tablenumber, mask, matches);
            u[i] = 'T';
            recursiveEnter(u, tablenumber, mask, matches);
            return;
         case 'D':         //A or G or T
            u[i] = 'A';
            recursiveEnter(u, tablenumber, mask, matches);
            u[i] = 'G';
            recursiveEnter(u, tablenumber, mask, matches);
            u[i] = 'T';
            recursiveEnter(u, tablenumber, mask, matches);
            return;
         case 'H':         //A or C or T
            u[i] = 'A';
            recursiveEnter(u, tablenumber, mask, matches);
            u[i] = 'C';
            recursiveEnter(u, tablenumber, mask, matches);
            u[i] = 'T';
            recursiveEnter(u, tablenumber, mask, matches);
            return;
         case 'K':         //G or T
            u[i] = 'G';
            recursiveEnter(u, tablenumber, mask, matches);
            u[i] = 'T';
            recursiveEnter(u, tablenumber, mask, matches);
            return;
         case 'M':         //A or C
            u[i] = 'A';
            recursiveEnter(u, tablenumber, mask, matches);
            u[i] = 'C';
            recursiveEnter(u, tablenumber, mask, matches);
            return;
         case 'N':         //A or C or G or T
            u[i] = 'A';
            recursiveEnter(u, tablenumber, mask, matches);
            u[i] = 'C';
            recursiveEnter(u, tablenumber, mask, matches);
            u[i] = 'G';
            recursiveEnter(u, tablenumber, mask, matches);
            u[i] = 'T';
            recursiveEnter(u, tablenumber, mask, matches);
            return;
         case 'R':         //A or G
            u[i] = 'A';
            recursiveEnter(u, tablenumber, mask, matches);
            u[i] = 'G';
            recursiveEnter(u, tablenumber, mask, matches);
            return;
         case 'S':         //C or G
            u[i] = 'C';
            recursiveEnter(u, tablenumber, mask, matches);
            u[i] = 'G';
            recursiveEnter(u, tablenumber, mask, matches);
            return;
         case 'U':         //allows U instead of T
            u[i] = 'T';
            recursiveEnter(u, tablenumber, mask, matches);
            return;
         case 'V':         //A or C or G
            u[i] = 'A';
            recursiveEnter(u, tablenumber, mask, matches);
            u[i] = 'C';
            recursiveEnter(u, tablenumber, mask, matches);
            u[i] = 'G';
            recursiveEnter(u, tablenumber, mask, matches);
            return;
         case 'W':         //A or T 
            u[i] = 'A';
            recursiveEnter(u, tablenumber, mask, matches);
            u[i] = 'T';
            recursiveEnter(u, tablenumber, mask, matches);
            return;
         case 'Y':         //C or T
            u[i] = 'C';
            recursiveEnter(u, tablenumber, mask, matches);
            u[i] = 'T';
            recursiveEnter(u, tablenumber, mask, matches);
            return;
         default:
            break;
      }   //closes switch
   }    //closes for (i = ...
   
   //At this point, u contains only A, C, G, or T
   //Convert to a binary number
   index = 0;
   for (i = 0; i < n; i++) {
      index <<= 2;
      switch (u[i]) {
         case 'C': index += 1; break;
         case 'G': index += 2; break;
         case 'T': index += 3; break;
         default: ; 
      } //closes switch
   }  // closes for (i = ... 
   //now we can OR in the mask into the table
   matchTable(tablenumber, index) |= mask;
   if (index >= (1<<(2*scanWidth)))
        printf("!!!Trouble in makeTable or recursiveEnter\n");
}

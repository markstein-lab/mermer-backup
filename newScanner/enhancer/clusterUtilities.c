//Copyright (c) in silico Labs, LLC  2009

#include <stdlib.h>
#include <stdio.h>
#include "scanner.h"
#include "timer.h"
#include <float.h>
#include <string.h>


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
          


char* fullName(char* path, char* name) {
   static char buffer[1000];
   strcpy(buffer, path);
   strcat(buffer, name);
   return buffer;
} 



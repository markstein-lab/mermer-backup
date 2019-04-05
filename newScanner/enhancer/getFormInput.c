//Copyright (c) 2009  insilico Labs, LLC.

#include <stdlib.h>
#include <stdio.h>
//#include "scanner.h"
#include "htmlForm.h"
#include <float.h>
#include <string.h>

int readField(char *name, int i, int inputLen);
   
alist terminalInput;
int alistLen;
char *queryString;
char *input;

void getFormInput() {
   
   char *strInputLen;
   int inputLen;
   int i;
   char name[1000];
   char data[1000];
   char c;
   alist item;
   
   inputLen = 0;
   queryString = getenv("QUERY_STRING");
   input = queryString;
   inputLen = strlen(queryString);

   i = 0;
   alistLen = 0;
   while (i < inputLen) {
      i = readField(name, i, inputLen);
      i = readField(data, i, inputLen);
      item = malloc(sizeof(struct formData)); 
      item->next = terminalInput; 
      terminalInput = item;
      item->name = malloc(strlen(name) + 2);
      item->data = malloc(strlen(data) + 2);
      strcpy(item->name, name);
      strcpy(item->data, data);
      alistLen++;
   }
   return;
}

int readField(char *name, int i, int inputLen) {
      //read a field from the standard input
      
   int j;
   char c;

   j = 0;
   while (i < inputLen) {
      sscanf(input, "%c", &c); //read an input character
      i++;
      input++;
      if (c == '=' || c == '&'|| c == 0) { //end of name
         name[j] = 0;  //end of name
         return i;
      }
      if (j > 998) continue; //parameter too long -- truncate
      if (c == '+') name[j++] = ' ';
      else if (c == '%') { //next two input chars encode a character
         sscanf(input, "%2x", &name[j]);
         j++;
         i+=2;  //we've consumed 2 more input characters
         input+=2;
      }
      else {
         name[j++] = c;
      }
   }
   //we get here when the input stream is empty.
   //finish off the current field
   name[j] = 0;
   return i; 
}

char* sassoc(char *name, alist terminalInput){ 
   int i;
   alist ptr;
   static char error[] = "NotFound";

   ptr = terminalInput;
   
   while (ptr != NULL) {
      if (strcmp(name, ptr->name) == 0) {
         //found the match.  Return pointer to value
         return ptr->data;
      }
      ptr = ptr->next;
   }
   //name not found. Return pointer to error
   return error;
}

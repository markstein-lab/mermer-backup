#include <stdlib.h>
#include <stdio.h>
#include <string.h>


void readFormInput();

typedef struct formData *alist;
   struct formData {
      char *name;
      char *data;
      alist next;
   };

extern alist terminalInput;
extern alistLen;

int main() {
    int i, n;
    char c;
    char *str;
    alist ptr;
    
    printf("Content-type: text/html\n\n\n");
    printf("<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0 Transitional//EN\">\
<html>\
<head>\
<title>DNA Match</title>\
</head>\
<BODY TEXT=#000000 LINK=#990000 VLINK=#660033 ALINK=#FFFF00>\
");
    printf("starting to process result of search screen\n");
    terminalInput = NULL;
    readFormInput();
    printf("There were %d pairs of {name,value} pairs<br>", alistLen);

    ptr = terminalInput;
    while (ptr) {
        printf("name = |%s|  value = |%s|<br>", ptr->name, ptr->data);
        ptr = ptr->next;
    }
   printf("</html>\n");
}

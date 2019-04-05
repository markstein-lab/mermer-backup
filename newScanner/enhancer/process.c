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
char* sassoc(char*, alist);
void peter_test(alist);

int main() {
    int i, n;
    char c;
    char *str, *sptr;
    alist ptr;
    
    printf("Content-type: text/html\n\n\n");
    printf("<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0 Transitional//EN\">\
<html>\
<head>\
<title>DNA Match</title>\
</head>\
<body TEXT=#000000 LINK=#990000 VLINK=#660033 ALINK=#FFFF00>\
");
    terminalInput = NULL;
    readFormInput();

    peter_test(terminalInput);

    printf("</body>\n");
   printf("</html>\n");
}

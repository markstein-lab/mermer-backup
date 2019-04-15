#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "htmlForm.h"


void readFormInput();

char* sassoc(char*, alist);
void peter_test(alist);
void getLine(char line[], FILE *f);
void copyInclude(char*);

int main() {
    int i, n;
    char c;
    char *str, *sptr;
    alist ptr;
    
    printf("Content-type: text/html\n\n\n");
    printf("<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 3.2//EN\">\
<html>\
<head>\
<title>Open Genomics : GenomeEnhancer</title>\
<link rel=\"stylesheet\" href=\"../experiment/meta/mainstyle.css\" type=\"text/css\"> \
</head>\
<body TEXT=#000000 LINK=#990000 VLINK=#660033 ALINK=#FFFF00>");

   copyInclude("../experiment/include.topleft.txt");

   printf("<div id=\"bodymain\">\n");
    terminalInput = NULL;
    readFormInput();

   printf("<script language=\"JavaScript\"> \
      function toggle(showHideDiv) { \
        var ele = document.getElementById(showHideDiv); \
        if (ele.style.display == \"block\") { \
           ele.style.display = \"none\"; \
        } \
        else { \
           ele.style.display = \"block\"; \
         } \
      } \
   </script> ");
 
    peter_test(terminalInput);
   
   printf(" </div>  <!-- close bodymain div -->\n");

//   copyInclude("../meta/include.footer.txt");

   printf("</body>\n");
   printf("</html>\n");
}
  
void copyInclude(char name[]) {
   FILE *temp;
   char line[5000];

   temp = fopen(name, "r");
   if (temp == NULL) {
      printf("File %s was not found\n<br>", name);
      return;
   }
   
   while(!feof(temp)) {
      getLine(line, temp);
      printf("%s\n", line);
   }
   fclose(temp);
   return;
}

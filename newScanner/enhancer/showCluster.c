
#include <stdlib.h>
#include <stdio.h>
#include <string.h>


void getFormInput();

typedef struct formData *alist;
   struct formData {
      char *name;
      char *data;
      alist next;
   };

extern alist terminalInput;
extern alistLen;
char* queryString;
char* sassoc(char*, alist);
void doShowCluster(alist);
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
<link rel=\"stylesheet\" href=\"/meta/mainstyle.css\" type=\"text/css\"> \
</head>\
<body TEXT=#000000 LINK=#990000 VLINK=#660033 ALINK=#FFFF00>");

   copyInclude("../meta/include.topleft.txt");

   printf("<div id=\"bodymain\">\n");
    terminalInput = NULL;
    getFormInput();
 
    doShowCluster(terminalInput);
quitting:   
   printf(" </div>  <!-- close bodymain div -->\n");

//   copyInclude("../meta/include.footer.txt");

   printf("</body>\n");
   printf("</html>\n");
}
  
void copyInclude(char name[]) {
   FILE *temp;
   char line[1000];

   temp = fopen(name, "r");
   if (temp == NULL) {
      printf("File %s was not found\n<br>", name);
      return;
   }
   
   while(!feof(temp)) {
      getLine(line, temp);
      printf("%s\n", line);
   }
   close(temp);
   return;
}

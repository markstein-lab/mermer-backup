
#include <stdlib.h>
#include <stdio.h>
#include <string.h>



int main() {
    int i, n;
    char c;
    char *str, *sptr;
    FILE *test;
    char line[1000];
    void getLine(char line[], FILE *f);

    test = fopen("../public_html/meta/include.topleft.txt", "r");
    if (test == NULL) {
       printf("File ../public_html/meta/include.toplef.txt not found\n");
       exit(0);
    }
    while (!feof(test)) {
        getLine(line, test);
        printf("%s\n", line);
    }
    close(test);
    exit(0);
    
}

void getLine(char line[], FILE *f) {
   int i;
   for (i=0; i<500; i++) {
      line[i] = getc(f);
      if (line[i] == '\n' || line[i] == EOF) break;
      if (line[i] == '\b') {
         i = (i > 1 ? i - 2 : -1);
      }

   }
   line[i] = 0;
}

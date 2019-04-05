
//Copyright (c) in silico Labs, LLC, 2009


void readFormInput();

typedef struct formData *alist;
   struct formData {
      char *name;
      char *data;
      alist next;
   };

extern alist terminalInput;
extern int alistLen;
char* sassoc(char*, alist);



//Copyright (c) in silico Labs, LLC 2008

#include <stdlib.h>
#include <stdio.h>
#include "scanner.h"
#include "timer.h"
#include <float.h>
#include <string.h>

void makestring(int, char*);

LOOKUPWORD* makeTables(motifDes *);

unsigned int DNALen;
unsigned int mapSize;
unsigned int ncontigs;
unsigned int extableSize;
unsigned int synCount;
unsigned int counter;
int numberOfTables;
unsigned int *av;
unsigned int *avbase;
motifDes *pv;
motifDes *pvbase;
int numberOfMatches;
LOOKUPWORD *matches;
int numberOfMotifs;
motifDes *motifs;
int *geneNumber;
unsigned char **geneName;
char** annotName;
char** chromoName;
unsigned int* arm;
unsigned int* extable;
ann_index *map;
FILE *save;
FILE *indexFile;
FILE *allnames;
int minGroupSize;
int windowSize;
char boolExp[200];
char* booleanSyntax(char message[]);
void reportClusters();
void printGeneList();

int do_the_search(unsigned int, unsigned int);

int main (int argc, char *argv[]) {

//Copyright (c) in silico Labs, LLC, 2008, 2009

//#define WIN32
#define NORMALWORD
#define scanWidth 4
#ifdef NORMALWORD
#define logWordWidth 5
typedef unsigned int LOOKUPWORD;
#else
#define logWordWidth 6
typedef unsigned long long int LOOKUPWORD;
#endif
#define wordWidth (1<<logWordWidth)
#include "htmlForm.h"

extern LOOKUPWORD *matches;
#define matchTable(i, j) matches[i*(1<<(2*scanWidth)) + j]

#if (scanWidth < 5)
   unsigned char *DNAstring;
#else 
   unsigned short int *DNAstring;
#endif

#if (scanWidth == 1) 
#define getDNA(x)  ((DNAstring[(unsigned)(x)/4] >> (6 - ((unsigned)(x)+(unsigned)(x))&6)) & 3)
#define DNAchar(x) lookup[(DNAstring[(unsigned)(x)/4])>>(2*(3 - ((unsigned)(x)%4))) & (3)]
#endif

#if (scanWidth == 2)
#define getDNA(x)  ((DNAstring[(unsigned)(x)/2] >> (4 - ((unsigned)(x)<<2) &4)) & 0xF)
#define DNAchar(x) lookup[(DNAstring[(unsigned)(x)/4])>>(2*(3 - ((unsigned)(x)%4))) & (3)]
#endif

#if (scanWidth > 2)
#define getDNA(x) DNAstring[x]
#define DNAchar(x) lookup[(DNAstring[(unsigned)(x)/scanWidth])>>(2*(scanWidth - 1 - ((unsigned)(x)%scanWidth))) & (3)]
#endif

typedef   struct ann_index ann_index;
struct ann_index{
   int annot_number;
   unsigned int start;
   unsigned int stop;
   int numberOfmRNAs;
   int numberOfExons;
   int *mrnaInfo;
   char direction;
   char *geneName;
   char touched;
   } ;

typedef struct mDes motifDes, *desPtr;
struct mDes {
   char *motif;
   char *reverse;
   char name;
   char direction;
   } ;
      
typedef struct mrnaID mrnaID;
struct mrnaID {
   int geneID;
   int mrnaStart;
};

void both(char []);
extern ann_index *map;
extern mrnaID *mrnaInfo;
extern unsigned long int DNALen;
extern unsigned int pseudoLen;
extern unsigned int mrnaCount;
extern unsigned int counter;
extern unsigned int mapSize;
extern unsigned int ncontigs;
extern unsigned int extableSize;
extern unsigned int synCount;
extern int numberOfTables;
extern unsigned long int *av;
extern motifDes *pv;
extern int numberOfMotifs;
extern motifDes *motifs;
extern unsigned long int *avbase;
extern motifDes *pvbase;
extern int numberOfMatches;
extern int *geneNumber;
extern unsigned char **geneName;
extern char** annotName;
extern char** chromoName;
extern unsigned long int* arm;
extern unsigned int* extable;
extern FILE *indexFile;
extern FILE *allnames;
extern FILE *save;
extern FILE *tbl;
extern int minGroupSize;
extern int windowSize;
extern char boolExp[5000];
extern int primeNames[100];
extern char pathName[];
extern int fromTerminal;
extern int numberOfGenes;
extern int printBoth;
extern int printNone;

void getLine(char line[], FILE *f);
char* fullName(char *, char *);
int findchrm(unsigned int);
int chrIndex(int);

void readFormInput();

void peter_test(alist);
extern int numberOfChromosomes;
extern char *chrName[100];
extern int *summaryData;
extern int limited;
extern int utr3;
extern int utr5;
extern int cds;
extern int orf;
extern int intron;
extern int exon;
extern int foundUtr3;
extern int foundUtr5;
extern int foundCds;
extern int foundOrf;
extern int foundIntron;
extern int foundExon;
extern int clusterNumber;
extern int printBoth;
#define rowSize(i)          (2*map[i].numberOfExons+4)
#define exonCount(a, b)     map[a].mrnaInfo[b*rowSize(a)+1]
#define element(a, b, c)    map[a].mrnaInfo[b*rowSize(a) + c]

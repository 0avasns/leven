/*****************************************************/
/*                                                   */
/* leven.c v1.0.2                                    */
/* Levenshtein distance calculator                   */
/* Athanassios Protopapas                            */
/* protopap@ilsp.gr                                  */
/* Institute for Language & Speech Processing        */
/* February 2011                                     */
/*                                                   */
/* Calculates mean Levenshtein distances for a set   */
/* of target items relative to a reference lexicon   */
/* (no cleanup or case conversion is performed;      */
/*  take care to provide pre-processed item files!)  */
/*                                                   */
/* Code based on Wikipedia algorithm                 */
/* http://en.wikipedia.org/wiki/Levenshtein_distance */
/* Thanks to Tal Yarkoni for help and examples!      */
/*                                                   */
/*****************************************************/
/* Compile with: gcc -O9 -o levencmd leven.c         */
/*****************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

#define _MAXLEN_ 50    /* maximum length of individual items (words)         */
#define _SUB 1 /* substitution cost  */
#define _DEL 1 /* deletion cost      */
#define _INS 1 /* insertion cost     */
#define _TRA 2 /* transposition cost, equal to two substitutions */
#define _TLEN 0 /* Flag to normalize distance by target length   */
#define _WLEN 0 /* Flag to normalize distance by neighbor length */
#define _WFRE 0 /* Flag to weigh distance by neighbor frequency  */

#define _NMAX_ 1000000 /* maximum number of items in lexicon or target list  */
#define _ODMAX_ 100    /* maximum number of "closest neighbors" considered   */
#define _ND_ 20        /* standard number of "closest neighbors" considered  */
#define _MINFREQ_ 0    /* minimum frequency of lexical items to consider     */ /* could be 50 */
#define _NMULT_ 20     /* maximum number of group sizes for multiple means   */

struct ldmat_struct {
  short int **d ; /* The distance calculation matric */
  int sub ; /* Substitution cost */
  int del ; /* Deletion cost     */
  int ins ; /* Insertion cost    */
  int tra ; /* Careful! this is actually the extra cost over substitution, defaults to the difference _TRA-_SUB */
} ;

struct ldmat_struct * setup_d( int sub, int del, int ins, int tra ) {
  /* Setup up and initialize the fixed components of the distance */
  /* calculation matrix, only once at the beginning, to save time */

  struct ldmat_struct *ldmat ;
  int i, j ;
  int m, n ;

  m=_MAXLEN_+1 ;
  n=_MAXLEN_+1 ;

  ldmat=(struct ldmat_struct *)malloc(sizeof(struct ldmat_struct)) ;

  ldmat->d=(short int **)malloc(m*sizeof(short int *)) ;
  for (i=0; i<m; i++) {
    ldmat->d[i]=(short int *)malloc(n*sizeof(short int)) ;
    for (j=0; j<n; j++) {
      ldmat->d[i][j]=0 ;
    }
    ldmat->d[i][0]=del*i ; /* was =i */ 
  }
  for (j=0; j<n; j++) {
    ldmat->d[0][j]=ins*j ; /* was =j */
  }
  ldmat->d[0][0]=0 ;

  ldmat->sub=sub ;
  ldmat->del=del ;
  ldmat->ins=ins ;
  ldmat->tra=tra ;

  return(ldmat) ;
}

struct neivec_struct {
  short int *OD ;    /* The vector of distances of the ND nearest neighbors  */
  float *fOD ;       /* A float version of the above, for weighted distances */
  short int ND ;     /* The number of lowest distances retained              */
  float *LDN ;       /* A vector to hold multiple averages for subsets of OD */
  char **neighbors ; /* For debugging, a list of the ND nearest neighbors    */
} ;

struct neivec_struct * setup_n( short int ND, short int ng ) {

  struct neivec_struct *n ;
  int i ;

  n = (struct neivec_struct *)malloc(sizeof(struct neivec_struct)) ;

  n->OD = (short int *)malloc(ND*sizeof(short int)) ;
  n->fOD= (float *)malloc(ND*sizeof(float)) ;
  if (ng>0) n->LDN= (float *)malloc(ng*sizeof(float)) ;
  n->ND = ND ;
  n->neighbors = (char **)malloc(ND*sizeof(char *)) ;
  for (i=0 ; i<ND; i++) n->neighbors[i]=(char *)malloc(_MAXLEN_*sizeof(char)) ;;

  return (n) ;
}

struct norm_struct {
  int tLength ;    /* Flag to normalize by target length      */
  int wLength ;    /* Flag to normalize by neighbor length    */
  int wnFreq ;     /* Flag to normalize by neighbor frequency */
  float *logfreq ; /* List of lexical item log frequencies    */
  float *nletter ; /* List of lexical item lengths            */
} ;

struct norm_struct * setup_w( int tLength, int wLength, int wnFreq, int Nlex ) {

  struct norm_struct *w ;

  w = (struct norm_struct *)malloc(sizeof(struct norm_struct)) ;

  w->tLength = tLength ;
  w->wLength = wLength ;
  w->wnFreq  = wnFreq ;
  if (wLength) w->nletter = (float *)malloc(Nlex*sizeof(float)) ;
  if (wnFreq)  w->logfreq = (float *)malloc(Nlex*sizeof(float)) ;
  
  return (w) ;

}

void show_d(struct ldmat_struct *ldmat) {
  /* For debugging purposes only */

  int i, j ;

  for (i=0; i<_MAXLEN_+1; i++) {
    for (j=0; j<_MAXLEN_+1; j++) {
      printf("%2d ",ldmat->d[i][j]);
    }
    printf("\n") ;
  }

}

void show_neighbors(struct neivec_struct *neivec) {
  /* For debugging purposes only */

  int i ;

  for (i=0; i<neivec->ND; i++) printf ("%s\t%d\n",neivec->neighbors[i],neivec->OD[i]) ;
}

int LevenshteinDistance0 ( char *s, char *t, short int **d ) {
  /* Calculate edit distance between two items       */
  /* Use fixed implicit weights sub=ins=del=1, tra=2 */
  /* This is a bit faster than the generic function  */

  int m, n, m1, n1, i, j, i1, j1 ;
  short int d10, d01, d11 ;

  m1=strlen(s) ;
  m=m1+1 ;
  n1=strlen(t) ;
  n=n1+1 ;

  for (j=1; j<n; j++) {
    j1 = j-1 ;
    for (i=1; i<m; i++) {
      i1 = i-1 ;
      if (s[i1]==t[j1]) {
	d[i][j] = d[i1][j1] ;
	}
      else {
	d10=d[i1][j] ;
	d01=d[i][j1] ;
	d11=d[i1][j1] ;
	d[i][j] = ((d10<d01) ? ((d10<d11)?d10:d11) : ((d01<d11)?d01:d11))+1 ;
      }
    }
  }

  return d[m1][n1] ;
}


int LevenshteinDistance ( char *s, char *t, short int **d,
			  int sub, int del, int ins, int tra ) {
  /* Calculate edit distance between two items  */
  /* Use variable weights for sub-del-ins-tra   */

  int m, n, m1, n1, i, j, i1, j1 ;
  int d10, d01, d11 ; /* The optimized gcc compiler apparently does not like short int additions! */
  char ss[_MAXLEN_+1], tt[_MAXLEN_+1] ; /* It is faster to keep these local and re-initialize */

  *ss=' ';
  *tt=' ';

  m1=strlen(s) ;
  m=m1+1 ;
  n1=strlen(t) ;
  n=n1+1 ;
  strcpy(ss+1,s) ;
  strcpy(tt+1,t) ;

  for (j=1; j<n; j++) {
    j1 = j-1 ;
    for (i=1; i<m; i++) {
      i1 = i-1 ;
      if (ss[i]==tt[j]) {
	d[i][j] = d[i1][j1] ;
	}
      else if ((ss[i]==tt[j1]) && (ss[i1]==tt[j])) { /* transposition */
	d10=d[i1][j]+del ;
	d01=d[i][j1]+ins ;
	d11=d[i1][j1]+tra ; /* adjusts for previous "substitution" */
	d[i][j] = (d10<d01) ? ((d10<d11)?d10:d11) : ((d01<d11)?d01:d11) ;
      }
      else {
	d10=d[i1][j]+del ;
	d01=d[i][j1]+ins ;
	d11=d[i1][j1]+sub ;
	d[i][j] = (d10<d01) ? ((d10<d11)?d10:d11) : ((d01<d11)?d01:d11) ;
      }
    }
  }

  return d[m1][n1] ;
}

int N_neigbors ( char *ts, char **lexa, int Nlex, struct ldmat_struct *ldmat ) {
  /* Return number of items with Levenshtein distance of 1 */
  /* Setting costs appropriately results in Coltheart's N  */
  /* or more inclusive measures of simple neighborhood     */ 

  int l, n, ld ;
  char *ls ;

  n=0 ;
  for (l=0; l<Nlex; l++) {
    ls=lexa[l] ;
    ld = LevenshteinDistance(ts,ls,ldmat->d,ldmat->sub,ldmat->del,ldmat->ins,ldmat->tra) ;
    if (ld==1) n++ ;
  }

  return (n) ;

}

void LDset ( char *ts, char **lexa, int Nlex, struct ldmat_struct *ldmat, struct neivec_struct *neivec ) {
  /* Fill integer array OD with the edit distance of ND nearest items to ts */

  int i, j, k, l, ld ;
  char *ls ;

  for (i=0; i<neivec->ND; i++) {
    neivec->OD[i]=_MAXLEN_ ; /* maximum difference cannot exceed maximum item length */
  }
  for (l=0; l<Nlex; l++) {
    ls=lexa[l] ;
    if (strcmp(ts,ls)!=0) {
      /* It speeds things up if parameters are passed individually to the LD loop */
      ld = LevenshteinDistance(ts,ls,ldmat->d,ldmat->sub,ldmat->del,ldmat->ins,ldmat->tra) ;
      i=0 ;
      while ((ld > neivec->OD[i]) && (i<neivec->ND)) i++ ;
      if (i<neivec->ND) {
	k=neivec->OD[i] ;
	neivec->OD[i]=ld ;
	i++ ;
      }
      for ( ; i<neivec->ND; i++) { /* continue from current i position */
	j=neivec->OD[i] ;
	neivec->OD[i]=k ;
	k=j ;
      }
    }
  }
   
}


void nLDset ( char *ts, char **lexa, int Nlex, struct ldmat_struct *ldmat, struct neivec_struct *neivec ) {
  /* A copy of LDset that saves the words along with the distances */

  int i, j, k, l, ld ;
  char *ls, ck[_MAXLEN_], cj[_MAXLEN_] ;

  for (i=0; i<neivec->ND; i++) {
    neivec->OD[i]=_MAXLEN_ ; /* maximum difference cannot exceed maximum item length */
  }
  for (l=0; l<Nlex; l++) {
    ls=lexa[l] ;
    if (strcmp(ts,ls)!=0) {
      ld = LevenshteinDistance(ts,ls,ldmat->d,ldmat->sub,ldmat->del,ldmat->ins,ldmat->tra) ;
      i=0 ;
      while ((ld > neivec->OD[i]) && (i<neivec->ND)) i++ ;
      if (i<neivec->ND) {
	k=neivec->OD[i] ;
	strcpy(ck,neivec->neighbors[i]) ;
	neivec->OD[i]=ld ;
	strcpy(neivec->neighbors[i],ls) ;
	i++ ;
      }
      for ( ; i<neivec->ND; i++) { /* continue from current i position */
	j=neivec->OD[i] ;
	strcpy(cj,neivec->neighbors[i]) ;
	neivec->OD[i]=k ;
	strcpy(neivec->neighbors[i],ck) ;
	k=j ;
	strcpy(ck,cj) ;
      }
    }
  }
   
}

void fLDset ( char *ts, char **lexa, int Nlex, struct ldmat_struct *ldmat, struct neivec_struct *neivec, struct norm_struct *wnorm ) {
  /* Fill float array fOD with the edit distance of ND nearest items to ts */
  /* Individual item distance can be weighted by item length and frequency */

  int i, l, ld ;
  float fld, fj, fk, nlettert ;
  char *ls ;

  for (i=0; i<neivec->ND; i++) {
    neivec->fOD[i]=(float)_MAXLEN_ ; /* maximum difference cannot exceed maximum item length */
  }
  if (wnorm->tLength) nlettert =  0.1*(float)strlen(ts) ; /* 0.1 is a scaling factor to retain order of magnitude. */
  for (l=0; l<Nlex; l++) {
    ls=lexa[l] ;
    if (strcmp(ts,ls)!=0) {
      ld = LevenshteinDistance(ts,ls,ldmat->d,ldmat->sub,ldmat->del,ldmat->ins,ldmat->tra) ;
      fld = (float)ld ;
      if (wnorm->tLength) fld /= nlettert ;
      if (wnorm->wLength) fld /= wnorm->nletter[l] ; 
      if (wnorm->wnFreq)  fld /= wnorm->logfreq[l] ; 
      i=0 ;
      while ((fld > neivec->fOD[i]) && (i<neivec->ND)) i++ ;
      if (i<neivec->ND) {
	fk=neivec->fOD[i] ;
	neivec->fOD[i]=fld ;
	i++ ;
      }
      for ( ; i<neivec->ND; i++) { /* continue from current i position */
	fj=neivec->fOD[i] ;
	neivec->fOD[i]=fk ;
	fk=fj ;
      }
    }
  }
   
}

float fLDcumul ( char *ts, char **lexa, int Nlex, struct ldmat_struct *ldmat, struct norm_struct *wnorm ) {
  /* Return cumulative edit distance for all items in the lexicon, added up */
  /* Individual item distance can be weighted by item length and frequency  */

  int l, ld ;
  float fld, c, nlettert ;
  char *ls ;

  c = 0.0 ;
  if (wnorm->tLength) nlettert =  0.1*(float)strlen(ts) ; /* 0.1 is a scaling factor to retain order of magnitude. */
  for (l=0; l<Nlex; l++) {
    ls=lexa[l] ;
    if (strcmp(ts,ls)!=0) {
      ld = LevenshteinDistance(ts,ls,ldmat->d,ldmat->sub,ldmat->del,ldmat->ins,ldmat->tra) ;
      fld = (float)ld ;
      if (wnorm->tLength) fld /= nlettert ;
      if (wnorm->wLength) fld /= wnorm->nletter[l] ;
      if (wnorm->wnFreq)  fld /= wnorm->logfreq[l] ;
      c += fld ;
    }
  }
  c /= (float)Nlex ; /* Normalize to something sensible */
  return(c) ;
   
}

float LDmean ( struct neivec_struct *neivec ) {
  /* Calculate mean of ND nearest-neighbor integer distances in OD */

  int i, n ;
  float m, odi ;

  m=0 ;
  n=0 ;
  for (i=0; i<neivec->ND; i++) {
    odi = (float)(neivec->OD[i]) ;
    m += odi ;
    n++ ;
  }
  m /= (float)n ;
  return (m) ;

}

float fLDmean ( struct neivec_struct *neivec ) {
  /* Calculate mean of ND nearest-neighbor float distances in fOD */

  int i, n ;
  float m, odi ;

  m=0 ;
  n=0 ;
  for (i=0; i<neivec->ND; i++) {
    odi = neivec->fOD[i] ;
    m += odi ;
    n++ ;
  }
  m /= (float)n ;
  return (m) ;
  
}

void LDmult (int ng, int dg, int Nlex, struct neivec_struct *neivec ) {
  /* Use OD to calculate mean edit distance for <ng> groups */
  /* incresasing in size in steps of <dg>, e.g., if ng=3    */
  /* and dg=5 there will be means for 5, 10, 15 neighbors   */

  int i, j, k, n ;
  float m ;

  m=0 ;
  n=0 ;
  for (j=0; j<ng; j++) {
    for (i=0; i<dg; i++) {
      k=j*dg+i ; /* increments of 5 */
      if (neivec->OD[k]<Nlex) {
	m += (float)(neivec->OD[k]) ;
	n++ ;
      }
    }
    neivec->LDN[j]=m/(float)n ;
  }

}


static const char *optString = "N:s:d:i:t:F:o:12D:flLm:h?" ;

void usage ( char *argv[] ) {
  /* Display command line parameter usage */
  printf("\nUsage: %s [-N <Nnei>] [-F <minF>:<freqFile>] [-x] [-s <sub>] [-d <del>] [-i <ins>] [-t <tra>] [-L] [-l] [-f] [-m <Ngrp>:<Nstep>] [-o <outfile>] <lexicon> [<targets>]\n\n",argv[0]) ;
  printf("\tNnei is the number of closest neighbors considered in the calculation of mean distance (default %d)\n",_ND_) ;
  printf("\t (set Nnei to 0 to add all distances, calculating an -- optionally weighted -- sum without selection)\n") ;
  printf("\t (set Nnei to -1 to receive Coltheart's N instead of mean Levenshtein distance)\n") ;
  printf("\t (set Nnei to -2 to receive number of single substitution/deletion/insertion/transposition neighbors)\n") ;
  printf("\t (if Nnei<0 any cost arguments are ignored)\n") ;
  printf("\tminF is the minimum frequency of occurrence for items in the lexicon to be considered (default %d)\n",_MINFREQ_) ;
  printf("\tfreqFile is a file with a list of frequencies of occurrence corresponding to the items in the lexicon\n") ;
  printf("\tsub is the character substitution cost (default %d)\n",_SUB) ;
  printf("\tdel is the character deletion cost (default %d)\n",_DEL) ;
  printf("\tins is the character insertion cost (default %d)\n",_INS) ;
  printf("\ttra is the character transposition cost (default %d)\n",_TRA) ;
  printf("\tThe -L switch normalizes neighbor distance by dividing over target length\n") ;
  printf("\tThe -l switch normalizes neighbor distance by dividing over neighbor length\n") ;
  printf("\tThe -f switch weights neighbor distance by dividing over neighbor log frequency (requires -F)\n") ;
  printf("\tThe -m switch calculates multiple LD means, for Ngrp sets increasing in size by Nstep\n") ;
  printf("\t (e.g., to get LD means for 10, 20, and 30 nearest neighbors, specify -m3:10)\n") ;
  printf("\toutfile is the filename to create and store the output (displayed on the terminal if not specified)\n");
  printf("\tlexicon is a file with a list of items to be considered as potential neighbors\n");
  printf("\ttargets is a file of the items for which the mean Levenshtein distance will be calculated;\n\t (the whole lexicon will be used when unspecified)\n\n");
  printf("Debugging calls:\n") ;
  printf("\n       %s [-s <sub>] [-d <del>] [-i <ins>] [-t <tra>] {-1|-2} <word1> <word2>\n\n",argv[0]) ;
  printf("\tword1 and word2 are two items for which Levenshtein distance is to be computed\n");
  printf("\t1 and 2 use different internal functions for the computation (1 fixes weights and ignores transposition)\n\n");
  printf("\n       %s  [-N <Nnei>] [-F <minF>:<freqFile>] [-s <sub>] [-d <del>] [-i <ins>] [-t <tra>] -D <word> <lexicon>\n\n",argv[0]) ;
  printf("\tword is a single target word, for which Nnei closest neighbors and corresponding distances are individually displayed\n\n") ;
  exit(-1);
}

int main( int argc, char *argv[] ) {

  int sub=_SUB, del=_DEL, ins=_INS, tra=(_TRA-_SUB), ND=_ND_ ;
  int tLength=_TLEN, wLength=_WLEN, wnFreq=_WFRE ;
  int i, t, MINFREQ=_MINFREQ_, j, n ;
  int Nlex, Ntar, freq, opt=0, Ngrp=0, Dgrp=0, totfreq=0, maxfreq=0 ; 
  int doFreq=0, allTarg=0, saveOut=0, Debug=0, ColtN=0, RDITN=0, Mult=0, Cumul=0 ;
  char word[_MAXLEN_], **lexa, **tara, *ts, *colon ;
  char *ofname=NULL, *lfname=NULL, *ffname=NULL, *tfname=NULL ;
  FILE *lexf, *tarf, *freqf, *outf ;
  float m, fnormfreq ;
  struct ldmat_struct *ldmat ;
  struct neivec_struct *neivec ;
  struct norm_struct *wnorm ;

  while ((opt = getopt(argc, argv, optString)) != -1) {
    switch(opt) {
    case 'N':
      ND=atoi(optarg);
      if (ND==0) Cumul=1 ;
      if (ND==-1) ColtN=1 ;
      if (ND==-2) RDITN=1 ;
      break ;
    case 'F':
      if ((colon=strchr(optarg,':'))==NULL) usage(argv);
      ffname=colon+1 ;
      doFreq=1 ;
      *colon='\0'; /* terminate the string to leave the number */
      MINFREQ=atoi(optarg);
      break ;
    case 's':
      sub=atoi(optarg);
      break ;
    case 'd':
      del=atoi(optarg);
      break ;
    case 'i':
      ins=atoi(optarg);
      break ;
    case 't':
      tra=atoi(optarg)-sub; 
      /* This will fail if -s is specified after -t ! */
      break ;
    case 'o':
      ofname=optarg ;
      saveOut=1 ;
      break ;
    case 'L':
      tLength=1 ;
      break ;
    case 'l':
      wLength=1 ;
      break ;
    case 'f':
      wnFreq=1 ;
      break ;
    case '1':
      Debug=1 ;
      break ;
    case '2':
      Debug=2 ;
      break ;
    case 'D':
      Debug=-1 ;
      ts=optarg ;
      break ;
    case 'm':
      if ((colon=strchr(optarg,':'))==NULL) usage(argv);
      Mult=1 ;
      *colon='\0'; /* terminate the string to leave two numbers */
      Ngrp=atoi(optarg);
      Dgrp=atoi(colon+1);
      break ;
    case 'h':
    case '?':
      usage(argv) ;
      break ;
    }
  }
  if (Debug>0) {
    ldmat=setup_d(sub,del,ins,tra);
    /* very unpolished shortcut to checking specific pairs */
    if (Debug==1) printf("%d\n",LevenshteinDistance0(argv[argc-2],argv[argc-1],ldmat->d));
    else if (Debug==2) printf("%d\n",LevenshteinDistance(argv[argc-2],argv[argc-1],ldmat->d,ldmat->sub,ldmat->del,ldmat->ins,ldmat->tra)) ;
    exit(0) ;
  }
  /* Process remaining arguments as filenames */
  if (optind>=argc) usage(argv); 
  else lfname=argv[optind++] ;
  if (optind<argc) tfname=argv[optind++] ;
  else allTarg=1 ;
  if (optind<argc) usage(argv);

  /* Read in lexicon and set up main arrays */
  lexa=(char **)malloc(_NMAX_*sizeof(char *));
  tara=(char **)malloc(_NMAX_*sizeof(char *));
  lexf=fopen(lfname,"r");
  if (lexf==NULL) { 
    printf("Could not open file %s to read lexicon\n", lfname);
    exit(-1);
  }
  if (doFreq) {
    freqf=fopen(ffname,"r");
    if (freqf==NULL) { 
      printf("Could not open file %s to read frequencies\n", ffname);
      exit(-1);
    }
  }
  for (i=0, j=0; !feof(lexf); i++) {
    if (fscanf(lexf,"%s",word)<1) break ;
    if (allTarg) {
      tara[i]=(char *)malloc(strlen(word)+1);
      strcpy(tara[i],word);
    }
    if (doFreq) {
      fscanf(freqf,"%d",&freq);
      if (freq>=MINFREQ) {
	lexa[j]=(char *)malloc(strlen(word)+1);
	strcpy(lexa[j++],word);
      }
      totfreq += freq ;
      if (freq>maxfreq) maxfreq=freq ;
    }
    else {
      lexa[j]=(char *)malloc(strlen(word)+1);
      strcpy(lexa[j++],word);
    }
  }
  fclose(lexf);
  if (doFreq) fclose(freqf);
  if (allTarg) Ntar=i ;
  Nlex=j ;

  /* A shortcut to viewing the closest ND neighbors, for debugging */
  if (Debug==-1) {
    ldmat=setup_d(sub,del,ins,tra) ;
    neivec=setup_n(ND,Ngrp) ;
    nLDset(ts,lexa,Nlex,ldmat,neivec) ; 
    show_neighbors(neivec) ;
    m=LDmean(neivec) ;
    printf("** mean of %d neighbors for <%s>: %.3f\n",ND,ts,m) ;      
    exit(0) ;
  }

  /* A bit of checking for obviously incompatible parameters */
  if (Mult) {
    if (ND <= 0) {
      printf("Cannot apply -m and -N 0/-1/-2 at the same time\n") ;
      exit(-1);
    }
    if ( Ngrp*Dgrp > Nlex ) {
      printf("Maximum number of nearest neighbors cannot exceed lexicon size\n") ;
      exit(-1);
    }
  }
  else {
    if (ND > Nlex) {
      printf("Number of nearest neighbors cannot exceed lexicon size\n") ;
      exit(-1);
    }
  }
  if ((ND > _ODMAX_) || (ND < -2)) {
    printf("Number of nearest neighbors cannot exceed %d\n",_ODMAX_) ;
    exit(-1);
  }
  /* There is no guarantee that accepted parameters will make sense! */

  wnorm=setup_w(tLength,wLength,wnFreq,Nlex) ;
  if (wnFreq) { /* Set up frequency array for weighting distances */
      if (doFreq==0) {
	printf("-f requires -F to specify frequencies file\n") ;
	exit(-1);
      }
      fnormfreq=(float)(maxfreq+1) ; /* Could be totfreq but here we control range to (0..1) */
      freqf=fopen(ffname,"r");
      if (freqf==NULL) { 
	printf("Could not open file %s to read frequencies\n", ffname);
	exit(-1);
      }
      for (i=0; i<Nlex; i++) {
	fscanf(freqf,"%d",&freq);
	if (freq>=MINFREQ) wnorm->logfreq[i] = 0.1*(1.0-log((float)freq/fnormfreq)) ;
	/* Controls frequency normalization to very large numbers for very rare items, */
	/* approaching 1.0 for very frequent items, as log(1)=0 */
	/* 0.1 is a scaling factor to retain order of magnitude. */
      }
      fclose(freqf) ;
  }
  if (wLength) { /* Set up length array for weighting distances  */
      for (i=0; i<Nlex; i++) wnorm->nletter[i] = 0.1*(float)strlen(lexa[i]) ;
      /* 0.1 is a scaling factor to retain order of magnitude. */
  }

  if (!allTarg) { /* Read separate targets file */
    tarf=fopen(tfname,"r");
    if (tarf==NULL) { 
      printf("Could not open file %s to read target words\n", tfname);
      exit(-1);
    }
    for (i=0; !feof(tarf); i++) {
      if (fscanf(tarf,"%s",word)<1) break ;
      tara[i]=(char *)malloc(strlen(word)+1);
      strcpy(tara[i],word);
    }
    Ntar=i ;
    fclose(tarf) ;
  }

  if (saveOut) {
    outf=fopen(ofname,"w");
    if (outf==NULL) { 
      printf("Could not open file %s to write output\n", ofname);
      exit(-1);
    }
  }

  if (ColtN || RDITN) { /* Options to return number of fixed-distance neighbors as opposed to distances */
    if (ColtN) ldmat=setup_d(1,2,2,2); /* set everything but substitutions greater than 1, to be rejected */
    else if (RDITN) ldmat=setup_d(1,1,1,1); /* set everything to 1, to be accepted */   
    for (t=0; t<Ntar; t++) {
      ts=tara[t] ;
      n=N_neigbors(ts,lexa,Nlex,ldmat) ;
      if (saveOut) {
	fprintf(outf,"%d\n",n) ;
	fflush(outf) ;
      }
      else printf("%s\t%d\n",ts,n) ; 
    }
  }
  else { /* Options to return mean or cumulative distances */
    ldmat=setup_d(sub,del,ins,tra) ;
    neivec=setup_n(ND,Ngrp) ;
    if (Mult) { /* Calculate distances for several different-size neighborhoods */
      for (t=0; t<Ntar; t++) {
	ts=tara[t] ;
	LDset(ts,lexa,Nlex,ldmat,neivec) ;
	LDmult(Ngrp,Dgrp,Nlex,neivec) ;
	if (saveOut) {
	  fprintf(outf,"%s",ts) ;      
	  for (j=0; j<Ngrp; j++) {
	    fprintf(outf,"\t%.3f",neivec->LDN[j]) ;
	  }
	  fprintf(outf,"\n") ;
	  fflush(outf) ;
	}
	else {
	  printf("%s",ts) ;      
	  for (j=0; j<Ngrp; j++) {
	    printf("\t%.3f",neivec->LDN[j]) ;
	  }
	  printf("\n") ;
	}
      }
    }
    else if (Cumul) { /* Add up distances for all eligible lexical items */
      for (t=0; t<Ntar; t++) {
	ts=tara[t] ;
	m=fLDcumul(ts,lexa,Nlex,ldmat,wnorm) ;
	if (saveOut) {
	  fprintf(outf,"%.3f\n",m) ;
	  fflush(outf) ;
	}
	else printf("%s\t%.3f\n",ts,m) ;      
      }
    }
    else if (wnFreq || wLength || tLength) { /* Average weighted distances, requires floats */
      for (t=0; t<Ntar; t++) {
	ts=tara[t] ;
	fLDset(ts,lexa,Nlex,ldmat,neivec,wnorm) ;
	m=fLDmean(neivec) ;
	if (saveOut) {
	  fprintf(outf,"%.3f\n",m) ;
	  fflush(outf) ;
	}
	else printf("%s\t%.3f\n",ts,m) ;      
      }
    }
    else { /* Default, the unweighted mean edit distances for a single neighborhood size */
      for (t=0; t<Ntar; t++) {
	ts=tara[t] ;
	LDset(ts,lexa,Nlex,ldmat,neivec) ;
	m=LDmean(neivec) ;
	if (saveOut) {
	  fprintf(outf,"%.3f\n",m) ;
	  fflush(outf) ;
	}
	else printf("%s\t%.3f\n",ts,m) ;      
      }
    }
  }

  if (saveOut) fclose(outf) ;
  return (0) ;

}

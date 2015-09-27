#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stddef.h>
#include <time.h>
#define MANCE 2000   /* max no. of ancestors for an individual */
#define MAXSTR 5000  /* max characters for a string */
#define MC 200       /* max no. of chiasma   */
#define MYEPS 0.5e-6
#define NR 10         /* no. of relationships considered as the null */
#define NRMLRT 11     /* no. of relationships considered in set A for
			 MLRT excluding the null relationship*/
#define PTHRESH 0.2  /* threshold for EIBD, AIBS and IBS */
#define MITER 5000   /* max no. of iteration in the EM algorithm */
#define NREP 100000   /* replicates number */
#define errorrate 0.01 /* approximate genotyping error rate */

/* for find all the relationships considered */
int findrelation(int** ped, int pedsize, int** rel);
int check_halfsib(int i, int j, int com, int ucom, int** ped, 
		  int pedsize);
int check_grandpc(int gp, int gc, int gcp, int gcup, int** ped,
		  int pedsize);
int check_avuncular(int un, int ne, int np, int nup, int** ped,
		    int pedsize);
int check_halfavuncular(int un, int ne, int np, int nup, int** ped,
			int pedsize);
int check_cousin(int i, int j, int c1p, int c2p, int c1up, int c2up,
		 int** ped, int pedsize);
int check_halfcousin(int i, int j, int c1p, int c2p, int c1up, int c2up,
		     int** ped, int pedsize);
int check_halfpluscousin(int i, int j, int com, int ucom, 
			 int** ped, int pedsize);
void build_anclist(int** ped, int pedsize, int id1, 
		   int anclist[MANCE][3], int* totalnum, int flag);
int check(int** ped, int pedsize, int id1, int id2);
int check_anclist(int** ped, int pedsize, int id2, 
		  int anclist[MANCE][3], int totalnum);
int check1_anclist(int** ped, int pedsize, int id1);
int inbred_check(int** ped, int pedsize, int id);
int recoding(int** ped, int* nomem);

/* for find genodata, and Mendelian error checking */
int findpedgeno(int** ped, int* nomem, FILE *fp_in2, 
		int nochrom, int* nomark, int** noalle,
		int** geno, char** chromnameped);
int checkMend(int** ped, int** geno, int* nomem,
	      int nochrom, int* nomark, char** chromnameped);

/* for trans. prob., cond. prob. and other prob.*/
double cp2(int f1, int m1, int f2, int m2, double* q);
double cp1(int f1, int m1, int f2, int m2, double* q);
double cp0(int f1, int m1, int f2, int m2, double* q);
void transprob(double dij, double pibdiibdj[NRMLRT+1][3+1][3+1], 
	       int reltype, int maptype);
void transproball(double dij, double pibdiibdj[NRMLRT+1][3+1][3+1], 
		  int maptype);
double pnorm(double z);

/* map functions */
double mapthetatocm(double theta, int maptype);
double mapcmtotheta(double cm, int maptype);

/* EM algorithm for estimation of p0, p1, p2 (0-1-2 IBD sharing) */
void   EMforp012(int* geno1, int* geno2, int* valid, 
		 double*** q, int nochrom, int* nomark,
		 int reltype, double phat[3+1]);

/* EIBD, AIBS, IBS statistic, mu and std calculation */
void get3statandmustd(int* geno1, int* geno2, int* valid, 
		      int commark, int reltype, int nochrom, 
		      int* nomark, double p012[3+1], 
		      double*** q, double** cenm, int maptype, 
		      double*** eeibdibd, double*** eeibdeibdibd,
		      double*** eaibsibd, double*** eaibsaibsibd,
		      double** eibsibd, double** eibsibsibd,
		      double statresult[3+1][5+1]);
void forEIBDvar(double*** eeibdeibd, double*** eeibdeibdibd, 
		int nochrom, int* nomark, int* ngeno, int*** genotype, 
		double*** q, double p012[NRMLRT+1][3+1]);
void forAIBSvar(double*** eaibsibd, double*** eaibsaibsibd, 
		int nochrom, int* nomark, int* ngeno, int*** genotype, 
		double*** q, double p012[NRMLRT+1][3+1]);
void forIBSvar(double** eibsibd, double** eibsibsibd, 
	       int nochrom, int* nomark, double*** q, int** noalle);
void get3stat(int* geno1, int* geno2, int* valid, int commark,
	      int reltype, int nochrom, int* nomark, 
	      double p012[4], double*** q,
	      double statresult[3+1]);

/* for likelihood calculation, use Markov approximation if IBD
   process is MC */
void getapprloglikelihood(int* geno1, int* geno2, int* valid,
			  double*** q, double** cenm, 
			  int nochrom, int* nomark, int maptype,
			  double loglhresult[NRMLRT+1]);

/* for simulate genodata of the relaitonships considered */
void simugeno1geno2(double*** q, double*** qcum, int nochrom, 
		    int* nomark, int **noalle, int reltype, 
		    double* chromlength, double** cenm,
		    int* simugeno1, int* simugeno2);
int crossover_new(double chromlength, double* pos, char c1, char c2,
                  char* ind);
int combine(double* p1, double* p2, double* ps, int nc1, int nc2, 
	    int ncs, char* i1, char* i2, char* is, double* pnew, 
	    char* inew);
void getibd2(double* p1,double* p2, double* pnew, int nc1,int nc2, 
	     int ncnew, char* i1,char* i2, char* inew, 
	     double chromlength, double* cenm, int nomark, 
	     int* ibd);
void getibd1(double* pnew1, double* pnew2, int ncnew1, int ncnew2,
	     char* inew1, char* inew2, double chromlength, 
	     double* cenm, int nomark, int* ibd);
void getibd0(double* pnew, int ncnew, char* inew, char c1, char c2,
	     double chromlength, double* cenm, int nomark, 
	     int* ibd);

/* for checking the errors of the data */
FILE *fp_error;

/* delete later */
/*  double checkmu[NR+1]; */

int main(int argc, char *argv[])
{

  /* -------------------------------------------------------------- */

  int    i, j, k, l, n, kk, kkk, itemp[10];
  int    flag, bugexit, numread, mpedid, mid, mmark, mpedsize;
  int mederror;
  int lineid;
  
  char   buf[MAXSTR+1], ctemp1[MAXSTR+1], ctemp2[MAXSTR+1];

  /* for pedigree info. */

  /* noped: no. of pedigrees in the pedgree file 
     nomem: pedigree id, size of the ped 
     ped: id, father id, mother id, 
     recoded fatherid, recoded motherid, gender 
     norel: no of each relationship considered 
     (1=full-sib, 2=half-sib, 3=grandparent-grandchild, 
     4=avuncular, 5=f-cousin, 6=unrelated,
     7=half-avuncular, 8=half-first-cousin, 9=half-sib+first-cousin, 
     10 = parent-offspring) 
     rel: id1 id2 reltype */
  int    pedid_now, noped, numpedid;
  int    **nomem, ***ped, **norel, ***rel;

  /* for marker info */

  /* nochrom: no. of chromosomes 
     nomark: no of markers
     noalle: no of allels 
     q: allele frequencies 
     qcum: cumulated allele frequencies 
     theta: spacing between markers in recombination fraction 
     cenm: spacing in centi Morgan 
     maptype: 1 = no-interference model, 2 = Kosambi */
  int    nochrom, *nomark, **noalle, maptype;
  double ***q, ***qcum, **theta, **cenm, *chromlength;
  char   **chromnameidx, **chromnameped;

  /* for EIBD, AIBS, IBS */

  /* p012: kinship coef.,   0, 1, 2 IBD sharing */
  static double p012[NRMLRT+1][3+1]={{0      , 0   , 0   , 0   }, 
				     {.25    , .25 , .5  , .25 }, 
				     {.125   , .5  , .5  , 0   }, 
				     {.125   , .5  , .5  , 0   },
				     {.125   , .5  , .5  , 0   }, 
				     {.0625  , .75 , .25 , 0   },
				     {0      , 1   , 0   , 0   },
				     {.0625  , .75 , .25 , 0   },
				     {.03125, .875, .125, 0   },
				     {.1875  , .375, .5  , .125},
				     {.25    , 0   , 1   , 0   },
				     {.5     , 0   , 0   , 1   }};

  /* genodata of all the members for a given pedigree */
  int    **geno;
  /* for a pair */
  int    reltype, findit, index1, index2, commark;
  int    *valid, **norelwdata;
  
  /* for variance calulation */
  int    *ngeno, ***genotype;
  double ***eeibdibd, ***eeibdeibdibd;
  double ***eaibsibd, ***eaibsaibsibd;
  double **eibsibd, **eibsibsibd;
  /* statistic, mu, std, z_obs, p-value(2-sided) */
  double statresult[3+1][5+1];

  /* for MLRT */

  int    domlrt, countbig[4+1];
  int    *simugeno1, *simugeno2, **nopairsvalid;
  double dtemp[10], phat[3+1];
  double stat_obs[3+1], stat_simu[3+1];
  double loglh_obs[NRMLRT+1], loglh_simu[NRMLRT+1];
  double mlrt_obs, mlrt_simu;
  static int countmlrt;
  
  /* for check the allele freq. */
  static double checkmarker[NR+1][2+1], checkallefreq[3+1][NR+1][5+1];
  
  FILE *fp_in1, *fp_in2;
  FILE *fp_idx, *tmpfp;
  FILE *fp_out1, *fp_out2;
  FILE *fp_out3;
  
  /* -------------------------------------------------------------- */

  if ( argc != 4 ) {
    printf("\nUsage: prest pedfile chromfiles 1 (or 2)\n\n");
    printf("       use 1 if apply only EIBD, AIBS, IBS test statistics\n");
    printf("       use 2 if further apply MLRT test\n\n");
    printf("example: prest pedfile chromfiles 1  \n\n");
    exit(0);
  }
  if ( (fp_in1 = fopen(argv[1],"r")) == 0 ) {
    printf("\nUnable to open pedfile: %s\n\n", argv[1]);
    exit(0);
  }
  if ( (fp_in2 = fopen(argv[2],"r")) == 0 ) {
    printf("\nUnable to open chromfiles: %s\n\n", argv[2]);
    exit(0);
  }
  if ( (fp_out1 = fopen("prest_out1","w")) == 0 ) {
    printf("\nUnable to open output file: prest_out1 \n\n");
    exit(0);
  }
  if ( (fp_out2 = fopen("prest_out2","w")) == 0 ) {
    printf("\nUnable to open output file: prest_out2 \n\n");
    exit(0);
  }
  if ( (fp_out3 = fopen("prest_out3","w")) == 0 ) {
    printf("\nUnable to open output file: prest_out3 \n\n");
    exit(0);
  }
  if ( (fp_error = fopen("prest_errors","w")) == 0 ) {
    printf("\nUnable to open output file: prest_errors\n\n");
    exit(0);
  }
  flag = atoi( argv[3] );
  if ( flag != 1 && flag != 2 ) {
    printf("\nIncorrect input: %s\n\n", argv[3]);
    printf("  Usage: prest pedfile chromfiles 1 (or 2)\n\n");
    printf("  use 1, if you only want to apply EIBD, AIBS, IBS\n");
    printf("  use 2, if you also want to apply MLRT\n\n");
    exit(0);
  }
  

  /* -------------------------------------------------------------- */
 
  bugexit = 0;
  maptype = 2;
  mederror = 0;
  
  srand((long)time(NULL));

  /* -------------------------------------------------------------- */

  /* read fp_in1 = pedigrees */

  /* noped: no. of the pedigrees in the file */
  noped = 0;
  pedid_now = -1;
  mpedid = mid = 1;

  lineid = 1;
  
  while ( fgets(buf, MAXSTR, fp_in1) ) {
    if (strlen(buf) >= MAXSTR - 1) {
      printf("\nWarning:  The length of the lines in file %s\n", argv[1]);
      printf("is too long (>5000).  Adjust MAXSTR and recompile.\n\n");
      exit(0);
    }
    
    numread = sscanf(buf, "%d %d %d %d %d %d", &itemp[1], &itemp[2],
		     &itemp[3], &itemp[4], &itemp[5], &itemp[6]);
    if ( (int)(log10(itemp[1])+1) > mpedid )
      mpedid = (int)(log10(itemp[1])+1);
    if ( (int)(log10(itemp[2])+1) > mid )
      mid = (int)(log10(itemp[2])+1);
    if ( numread != 6 ) {
      printf("\nInvalid %d_th line in pedfile:\n%s\n", lineid, buf);
      printf("Correct format is: pedid id f_id m_id sex affect_status\n\n");
      exit(0);
    }
    if ( itemp[1] != pedid_now ) { 
      ++noped;
      pedid_now = itemp[1];
    }
    ++lineid;
  }

  ped = (int ***)malloc((size_t) ((noped+1)*sizeof(int**)));
  nomem = (int **)malloc((size_t) ((noped+1+1)*sizeof(int*)));
  for (i = 1; i <= (noped+1); i++)
    nomem[i] = (int *)malloc((size_t) ((2+1)*sizeof(int)));
  nomem[noped+1][2] = 0;
 
  /* nomem: pedid, size of the ped*/
  rewind(fp_in1);
  pedid_now = -1;
  numpedid = 0;
   
  while ( fgets(buf, MAXSTR, fp_in1) ) {
    sscanf(buf, "%d %d %d %d %d %d", &itemp[1], &itemp[2],
	   &itemp[3], &itemp[4], &itemp[5], &itemp[6]);
    if (itemp[1] != pedid_now) { 
      ++numpedid;
      nomem[numpedid][1] = pedid_now = itemp[1];
      nomem[numpedid][2] = 1;
    }
    else {
      ++nomem[numpedid][2];
    }
  }
  
  for (i = 1; i <= noped; i++) {
    nomem[noped+1][2] += nomem[i][2];
    ped[i] = (int **)malloc((size_t) ((nomem[i][2]+1)*sizeof(int*)));
    for (j = 1; j <= nomem[i][2]; j++) {
      ped[i][j] = (int *)malloc((size_t) ((6+1)*sizeof(int)));
    }
  }

  /* ped: id, fatherid, motherid, recoded fid, recoded mid, gender */
  rewind(fp_in1);

  for (i = 1; i <= noped; i++) {
    for (j = 1; j <= nomem[i][2]; j++) {
      fscanf(fp_in1,"%d", &itemp[1]);
      fscanf(fp_in1,"%d %d %d %d", &ped[i][j][1], 
	     &ped[i][j][2], &ped[i][j][3], &ped[i][j][6]);
      fscanf(fp_in1,"%d", &itemp[6]);
    }
  }

  fclose(fp_in1);
 

  /* -------------------------------------------------------------- */

  /* recoded parents's ids and errors check */
  for (i = 1; i <= noped; i++) {
    for (j = i+1; j <= noped; j++) {
      if ( nomem[i][1] == nomem[j][1] ) {
	fprintf(fp_error, "pedigree id = %d being used twice\n\n",
		nomem[i][1]);
	++bugexit;
      }
    }
    bugexit += recoding(ped[i], nomem[i]);
  }
  if ( bugexit != 0 ) {
    printf("\nError in the pedigrees\n");
    printf("\nCheck file prest_errors file, and re-run the program\n\n");
    exit(0);
  }
  

  /* -------------------------------------------------------------- */

  /* find all the relationship types considered 
     1=full-sib, 2=half-sib, 3=grandp-c, 4=avuncular, 5=first-cousin,
     6=unrelated, 7=half-avuncular, 8=half-first-cousin, 
     9=half-sib+first-cousin, 
     10 = parent-offspring */

  rel = (int ***)malloc((size_t) ((noped+1)*sizeof(int**)));
  for (i = 1; i <= noped; i++) {
    kk = nomem[i][2]*(nomem[i][2]-1)/2;
    rel[i] = (int **)malloc((size_t) ((kk+1)*sizeof(int*)));
    for (j = 1 ; j <= kk; j++)
      rel[i][j] = (int *)malloc((size_t) ((3+1)*sizeof(int)));
  }
  norel = (int **)malloc((size_t) ((noped+1+1)*sizeof(int*)));
  for (i = 1; i <= (noped+1); i++) {
    norel[i] = (int *)malloc((size_t) ((NR+1+1)*sizeof(int)));
    for (j = 1; j <= (NR+1); j++)
      norel[i][j] = 0;
  }
  
  for (i = 1; i <= noped; i++) {
    /* input: i_th ped, size 
       return: no. of the pairs found, and
       rel: id1, id2, reltype */
    norel[i][NR+1] = findrelation(ped[i], nomem[i][2], rel[i]);
    norel[noped+1][NR+1] += norel[i][NR+1];
    for (j = 1; j <= norel[i][NR+1]; j++) {
      for (k = 1; k <= NR; k++) {
	if ( rel[i][j][3] == k )
	  ++norel[i][k];
      }
    }
  } 
  

  /* -------------------------------------------------------------- */


  /* read fp_in2 = chromfiles */

  /* nochrom: no. of chromosomes */
  nochrom = 0;

  while ( fgets(buf, MAXSTR, fp_in2) ) {
    if (strlen(buf) >= MAXSTR - 1) {
      printf("\nWarning:  The length of the lines in file %s\n", argv[2]);
      printf("is too long (>5000).  Adjust MAXSTR and recompile.\n\n");
      exit(0);
    }

    ++nochrom;
    numread = sscanf(buf, "%s %s", ctemp1, ctemp2);
    if ( numread != 2 ) {
      printf("\nInvalid line in chromfiles: %s\n%s\n", argv[2], buf);
      printf("\nCorrect format is: marker_file genotype_file\n");
      exit(0);
    }
    if ( (tmpfp = fopen(ctemp1, "r")) == 0 ) {
      printf("\nUnable to open marker_file %s\n", ctemp1);
      exit(0);
    }
    fclose(tmpfp);
    if ( (tmpfp = fopen(ctemp2, "r")) == 0 ) {
      printf("\nUnable to open genotype_file %s\n", ctemp2);
      exit(0);
    }
    fclose(tmpfp);
  }
  
  q = (double ***)malloc((size_t) ((nochrom+1)*sizeof(double**)));
  qcum = (double ***)malloc((size_t) ((nochrom+1)*sizeof(double**)));
  theta = (double **)malloc((size_t) ((nochrom+1)*sizeof(double*)));
  cenm = (double **)malloc((size_t) ((nochrom+1)*sizeof(double*)));
  noalle = (int **)malloc((size_t) ((nochrom+1)*sizeof(int*)));
  nomark = (int *)malloc((size_t) ((nochrom+1)*sizeof(int)));
  chromlength = (double *)malloc((size_t) ((nochrom+1)*sizeof(double)));
  chromnameped = (char **)malloc((size_t) ((nochrom+1)*sizeof(char*)));
  chromnameidx = (char **)malloc((size_t) ((nochrom+1)*sizeof(char*)));
  for (i = 1; i <= nochrom; i++) {
    chromnameped[i] = (char *)malloc((size_t) ((MAXSTR+1)*sizeof(char)));
    chromnameidx[i] = (char *)malloc((size_t) ((MAXSTR+1)*sizeof(char)));
    chromlength[i] = 0;
  }
  chromlength[0] = 0;
  nomark[0] = 0;
  
  /* open and read .idx files : markers info. */

  rewind(fp_in2);
    
  for (i = 1; i <= nochrom; i++) {

    /* nomark: no. of markers in each chromosome */
    fgets(buf, MAXSTR, fp_in2);
    sscanf(buf,"%s %s", chromnameidx[i], chromnameped[i]);
    fp_idx = fopen(chromnameidx[i], "r");
    fgets(buf, MAXSTR, fp_idx);
    if (strlen(buf) >= MAXSTR - 1) {
      printf("\nWarning:  The length of the lines in file %s\n", chromnameidx[i]);
      printf("is too long (>5000).  Adjust MAXSTR and recompile.\n\n");
      exit(0);
    }
    numread = sscanf(buf,"%d", &nomark[i]);
    if (numread != 1) {
      printf("\nError reading num_mark from file %s\n %s\n",
	     chromnameidx[i], buf);
      exit(0);
    }
    nomark[0] += nomark[i];

    q[i] = (double **)malloc((size_t) ((nomark[i]+1)*sizeof(double*)));
    qcum[i] = (double **)malloc((size_t) ((nomark[i]+1)*sizeof(double*)));
    theta[i] = (double *)malloc((size_t) ((nomark[i]+1)*sizeof(double)));
    cenm[i] = (double *)malloc((size_t) ((nomark[i]+1)*sizeof(double)));
    noalle[i] = (int *)malloc((size_t) ((nomark[i]+1)*sizeof(int)));

    /* noalle: no. of alleles of each marker */    
    for (j = 1; j<= nomark[i]; j++) {

      fgets(buf, MAXSTR, fp_idx);
      if (strlen(buf) >= MAXSTR - 1) {
	printf("\nWarning:  The length of the lines in file %s\n", chromnameidx[i]);
	printf("is too long (>5000).  Adjust MAXSTR and recompile.\n\n");
	exit(0);
      }
      numread = sscanf(buf,"%d %d", &itemp[1], &noalle[i][j]);
      if (numread != 2) {
	printf("\nError in file %s:\n %s\n", chromnameidx[i], buf);
	exit(0);
      }

      q[i][j] = (double *)malloc((size_t) 
				 ((noalle[i][j]+1)*sizeof(double)));
      qcum[i][j] = (double *)malloc((size_t) 
				    ((noalle[i][j]+1)*sizeof(double)));
      q[i][j][0] = qcum[i][j][0] = 0;

      /* q: allele freq. of each allele */
      for (k = 1; k <= noalle[i][j]; k++) {
	numread = fscanf(fp_idx, "%lf", &q[i][j][k]);
	if (numread != 1) {
	  printf("\nError of allele freq. in %d_th marker in file %s:\n",
		 j, chromnameidx[i]);
	  exit(0);
	}
	if (q[i][j][k] <= 0 || q[i][j][k] > 1) {
	  printf("\nError of allele freq. in %d_th marker in file %s (<=0 or > 1):\n",
		 j, chromnameidx[i]);
	  exit(0);
	}
	qcum[i][j][k] = qcum[i][j][k-1] + q[i][j][k];
      }
      fgets(buf, MAXSTR, fp_idx);
      if (strlen(buf) >= MAXSTR - 1) {
	printf("\nWarning:  The length of the lines in file %s\n", chromnameidx[i]);
	printf("is too long (>5000).  Adjust MAXSTR and recompile.\n\n");
	exit(0);
      }
      if ( fabs(qcum[i][j][noalle[i][j]]-1) > 0.01 ) {
	fprintf(fp_error," Error in %s's %d_th marker\n",
		chromnameidx[i], j);
	fprintf(fp_error, "sum of allele freq. = %f\n\n", 
		qcum[i][j][noalle[i][j]]);
	++bugexit;
      }
      qcum[i][j][noalle[i][j]] = 1;
    }

    /* theta/cenm: marker spacing in r.f. or centi-Morgan,
       important: n markers, n-1 spacings, 
       and assign 0 to the first marker */
    theta[i][1] = cenm[i][1] = 0;
    for (j = 1; j<= (nomark[i]-1); j++) {
      numread = fscanf(fp_idx,"%lf", &theta[i][j+1]);
      if ( numread != 1 ) {
	fprintf(fp_error, "Error of marker spacing in file %s:\n",
		chromnameidx[i]);
	fprintf(fp_error, "%d markers, should be (%d - 1) spacings\n\n", 
		nomark[i], nomark[i]);
	++bugexit;
      }
      if ( theta[i][j+1] > 1 ) {
	fprintf(fp_error, "Error of marker spacing in file %s:\n",
		chromnameidx[i]);
	fprintf(fp_error, "Spacing between adjacent markers should be in recombination fraction (theta), not in centi-Morgan\n\n");
	++bugexit;
      }
      if ( theta[i][j+1] < 0 ) {
	fprintf(fp_error, "Error of marker spacing in file %s:\n",
		chromnameidx[i]);
	fprintf(fp_error, " Spacing cannot be negative value\n\n");
	++bugexit;
      }
      if ( theta[i][j+1] == 0 ) 
	theta[i][j+1] = 0.0001;
      
      /* different map function */
      cenm[i][j+1] = mapthetatocm(theta[i][j+1], maptype);
      chromlength[i] += cenm[i][j+1];
    }

    /* n markers n-1 spacings, not n */
    numread = fscanf(fp_idx,"%lf", &dtemp[1]);
    if ( numread == 1 ) {
      fprintf(fp_error, "Error of marker spacing in file %s:\n",
	      chromnameidx[i]);
      fprintf(fp_error, "%d markers, should be (%d - 1) spacings\n\n", 
	      nomark[i], nomark[i]);
      ++bugexit;
    }
    /* add 10 CM for later simulation */
    chromlength[0] += chromlength[i];
    chromlength[i] += 10;

    fclose(fp_idx);
  }
      
  if ( bugexit != 0 ) {
    printf("\n Error in marker data\n");
    printf("\n Check prest_errors file, and re-run the program\n\n");
    exit(0);
  }
  

  /* -------------------------------------------------------------- */

  /* Prepare Variance of EIBD, AIBS, and IBS */

  /* ngeno: no. of possible genotypes for each marker = n(n+1)/2 
     genotype: detailed possible genotype (g1, g2) for each
     possible genotype for each marker */ 
  ngeno=(int *)malloc((size_t) ((nomark[0]+1)*sizeof(int)));
  genotype=(int ***)malloc((size_t) ((nomark[0]+1)*sizeof(int**)));

  kkk = 0;
  for (i = 1; i <= nochrom; i++) {
    for (j = 1; j <= nomark[i]; j++) {
      ngeno[kkk+j] = noalle[i][j]*(noalle[i][j]+1)/2;
      genotype[kkk+j] = (int **)malloc((size_t) 
				       ((ngeno[kkk+j]+1)*sizeof(int*)));
      for (k = 1; k <= ngeno[kkk+j]; k++) 
	genotype[kkk+j][k] = (int *)malloc((size_t) ((2+1)*sizeof(int)));
    }
    kkk += nomark[i];
  }
  
  kkk = 0;
  for (i = 1; i <= nochrom; i++) {
    for (j = 1; j <= nomark[i]; j++) {
      kk = 1;
      for (k = 1; k <= noalle[i][j]; k++) {
	for (l = k; l <= noalle[i][j]; l++) {
	  genotype[kkk+j][kk][1]=k;
	  genotype[kkk+j][kk][2]=l;
	  ++kk;
	}
      }
    }
    kkk += nomark[i];
  }

  /* for EIBD variance */
  eeibdibd = (double ***)malloc((size_t) ((NR+1)*sizeof(double**)));
  eeibdeibdibd = (double ***)malloc((size_t)
				    ((NR+1)*sizeof(double**)));
  for (i = 1; i <= NR; i++) {
    eeibdibd[i] = (double **)malloc((size_t) 
				    ((nomark[0]+1)*sizeof(double*)));
    eeibdeibdibd[i] = (double **)malloc((size_t) 
					((nomark[0]+1)*sizeof(double*)));
    for (j = 1; j <= nomark[0]; j++) {
      eeibdibd[i][j] = (double *)malloc((size_t) 
					((3+1)*sizeof(double)));
      eeibdeibdibd[i][j] = (double *)malloc((size_t) 
					    ((3+1)*sizeof(double)));
    }
  }
  forEIBDvar(eeibdibd, eeibdeibdibd, nochrom, nomark, 
	     ngeno, genotype, q, p012);
  
  /* for AIBS variance */
  eaibsibd = (double ***)malloc((size_t) ((NR+1)*sizeof(double**)));
  eaibsaibsibd = (double ***)malloc((size_t) 
				    ((NR+1)*sizeof(double**)));
  for (i = 1; i <= NR; i++) {
    eaibsibd[i] = (double **)malloc((size_t) 
				    ((nomark[0]+1)*sizeof(double*)));
    eaibsaibsibd[i] = (double **)malloc((size_t) 
					((nomark[0]+1)*sizeof(double*)));
    for (j = 1; j <= nomark[0]; j++) {
      eaibsibd[i][j] = (double *)malloc((size_t) 
					((3+1)*sizeof(double)));
      eaibsaibsibd[i][j] = (double *)malloc
	((size_t) ((3+1)*sizeof(double)));
    }
  }
  forAIBSvar(eaibsibd, eaibsaibsibd, nochrom, nomark, 
	     ngeno, genotype, q, p012);

  /* for IBS variance 
     notice that it does not depend on reltype 
     eibsibd[i][1,2,3]: E(IBS at locus i|ibd=0) E(IBS_i|ibd=1)
     E(IBS_i|ibd=2) */
  eibsibd = (double **)malloc((size_t) 
			      ((nomark[0]+1)*sizeof(double*)));
  eibsibsibd = (double **)malloc((size_t) 
				 ((nomark[0]+1)*sizeof(double*)));
  for (i = 1; i <= nomark[0]; i++) {
    eibsibd[i] = (double *)malloc((size_t) ((3+1)*sizeof(double)));
    eibsibsibd[i] = (double *)malloc((size_t) ((3+1)*sizeof(double)));
  }
  forIBSvar(eibsibd, eibsibsibd, nochrom, nomark, q, noalle);


  /* -------------------------------------------------------------- */

  /* for each pedigree, calculat EIBD, AIBS, IBS statistics, 
     mu, and std. Use Normal approximation to obtain 2-sided p-value.
     If users want to perform MLRT test, further simulate genotype
     data for each pair, and calculate empirical p-value for EIBD,
     AIBS, IBS and MLRT */

  valid = (int *)malloc((size_t) ((nomark[0]+1)*sizeof(int)));

  norelwdata = (int **)malloc((size_t) ((noped+1+1)*sizeof(int*)));
  nopairsvalid = (int **)malloc((size_t) ((noped+1+1)*sizeof(int*)));
  for (i = 1; i <= (noped+1); i++) {
    norelwdata[i] = (int *)malloc((size_t) ((NR+1+1)*sizeof(int)));
    nopairsvalid[i] = (int *)malloc((size_t) ((NR+1+1)*sizeof(int)));
    for (j = 1; j <= (NR+1); j++) {
      norelwdata[i][j] = 0;
      nopairsvalid[i][j] = 0;
    }
  }

  mpedsize = nomem[1][2];
  for (i = 2; i <= noped; i++) {
    if( nomem[i][2] > mpedsize)
      mpedsize = nomem[i][2];
  }
  geno = (int **)malloc((size_t) ((mpedsize+1)*sizeof(int*)));
  for (i = 1; i<= mpedsize; i++)
    geno[i] = (int *)malloc((size_t) ((2*nomark[0]+1)*sizeof(int)));
 
  simugeno1 = (int *)malloc((size_t) ((2*nomark[0]+1)*sizeof(int)));
  simugeno2 = (int *)malloc((size_t) ((2*nomark[0]+1)*sizeof(int)));

  /* delete later */
/*    for (k = 1; k <= NR; k++) */
/*      checkmu[k] = 0; */

  /* for each pedigree */
  for (i = 1; i <= noped; i++) {

    /* find geno data for a given pedigree */
    bugexit += findpedgeno(ped[i], nomem[i], fp_in2, nochrom,
			   nomark, noalle, geno, chromnameped);

    /* check Mendelian errors between parent-child */
    mederror += 
      checkMend(ped[i], geno, nomem[i], nochrom, nomark, chromnameped);

/*      bugexit +=  */
/*        checkMend(ped[i], geno, nomem[i], nochrom, nomark, chromnameped); */

/*      if (bugexit != 0) { */
/*        printf("\n Error in genodata\n"); */
/*        printf("\nCheck prest_errors file, and re-run the program\n\n"); */
/*        exit(0); */
/*      } */

    /* for each pair */
    for (j = 1; j <= norel[i][NR+1]; j++) {

      reltype = rel[i][j][3];
      findit = 0;
      k = 1;
      while(findit != 2) {
	if ( geno[k][0] == rel[i][j][1] ) {
	  index1 = k;
	  ++findit;
	}
	else if ( geno[k][0] == rel[i][j][2] ) {
	  index2 = k;
	  ++findit;
	}
	++k;
      }

      /* commark: no. of markers typed in both individuals */
      commark = 0;
      for (k = 1; k <= nomark[0]; k++) {
	if ( geno[index1][2*k-1] != 0 && geno[index1][2*k] != 0 &&
	     geno[index2][2*k-1] != 0 && geno[index2][2*k] != 0 ) {
	  ++commark;
	  valid[k] = 1;
	}
	else
	  valid[k] = 0;
      }
      
      /* testing only for the pair that commark != 0 */
      if (commark != 0) {
	
	domlrt = 0;

	++norelwdata[i][NR+1];
	++norelwdata[noped+1][NR+1];
	for (k = 1; k <= NR; k++) {
	  if ( reltype == k ) {
	    ++norelwdata[i][k];
	    ++norelwdata[noped+1][k];
	  }
	}
	
	/* statresult: statistic. mu, std, z_obs, p-value of 
	   EIBD, AIBS, IBS */
	get3statandmustd(geno[index1], geno[index2], valid, commark,
			 reltype, nochrom, nomark, p012[reltype], 
			 q, cenm, maptype, eeibdibd, eeibdeibdibd, 
			 eaibsibd, eaibsaibsibd, eibsibd, eibsibsibd,
			 statresult);

	mmark = (int)log10(nomark[0]) + 1;
	sprintf(ctemp1, "%%%dd %%%dd %%%dd  %%2d  %%%dd ", mpedid,
		mid, mid, mmark);
	fprintf(fp_out2, ctemp1,  
		nomem[i][1], rel[i][j][1], rel[i][j][2],
		rel[i][j][3], commark);

	/* Use EM algorithm to estimate p0, p1, and p2 
	   (prob. of 0-1-2 IBD sharing) */
	EMforp012(geno[index1], geno[index2], valid, 
		  q, nochrom, nomark, reltype, phat);
	
	/* EIBD, AIBS are not applicable for unrelated = 6 
	   use IBS for screening */
	/* EIBD is not applicable for parent-offspring = 10 and MZ
	   twins = 11, and AIBS and IBS are not powerful for
	   distinguish parent-offspring from full-sib, so use estimate
	   of p0, p1, p2 to judge*/
	if (reltype == 6) {
	  fprintf(fp_out2,"     NA  %6.4f %6.4f %6.4f       NA      NA %7.5f   ",
		  phat[1], phat[2], phat[3], statresult[3][5]);
	  if ((statresult[3][5] < PTHRESH && statresult[3][4] > 0) ||
	      phat[2] > 0.75 )
	    ++countmlrt;
	  if ( phat[2] > 0.75 ) {
	    sprintf(ctemp1, "%%%dd %%%dd %%%dd  %%2d  %%%dd ", mpedid,
		    mid, mid, mmark);
	    fprintf(fp_out3, ctemp1,  
		    nomem[i][1], rel[i][j][1], rel[i][j][2],
		    rel[i][j][3], commark);
	    fprintf(fp_out3,"  %6.4f %6.4f %6.4f  %6.4f %6.4f %6.4f\n",
		    p012[reltype][1],p012[reltype][2],p012[reltype][3],
		    phat[1], phat[2], phat[3]);
	    fflush(fp_out3);
	  }
	}
	else if (reltype == 10 || reltype == 11) {
	  fprintf(fp_out2,"     NA  %6.4f %6.4f %6.4f       NA %7.5f %7.5f   ",
		  phat[1], phat[2], phat[3],statresult[2][5], statresult[3][5]);
	  if(reltype == 10 && phat[2] < 0.9 ) {
	    ++countmlrt;
	    sprintf(ctemp1, "%%%dd %%%dd %%%dd  %%2d  %%%dd ", mpedid,
		    mid, mid, mmark);
	    fprintf(fp_out3, ctemp1,  
		    nomem[i][1], rel[i][j][1], rel[i][j][2],
		    rel[i][j][3], commark);
	    fprintf(fp_out3,"  %6.4f %6.4f %6.4f  %6.4f %6.4f %6.4f\n",
		    p012[reltype][1],p012[reltype][2],p012[reltype][3],
		    phat[1], phat[2], phat[3]);
	    fflush(fp_out3);
	  }
	  else if(reltype == 11 && phat[3] < 0.95) {
	    ++countmlrt;
	    sprintf(ctemp1, "%%%dd %%%dd %%%dd  %%2d  %%%dd ", mpedid,
		    mid, mid, mmark);
	    fprintf(fp_out3, ctemp1,  
		    nomem[i][1], rel[i][j][1], rel[i][j][2],
		    rel[i][j][3], commark);
	    fprintf(fp_out3,"  %6.4f %6.4f %6.4f  %6.4f %6.4f %6.4f\n",
		    p012[reltype][1],p012[reltype][2],p012[reltype][3],
		    phat[1], phat[2], phat[3]);
	    fflush(fp_out3);
	  }
	}
	/* for the others, use EIBD and AIBS tests, and also 
	   use p1 > 0.85 to avoid the true is parent-offspring, but
	   null is not, especially the null is full-sib */
    	else {
	  fprintf(fp_out2, " %6.4f  %6.4f %6.4f %6.4f  %7.5f %7.5f %7.5f   ", 
		  statresult[1][1], phat[1], phat[2], phat[3],
		  statresult[1][5], statresult[2][5],
		  statresult[3][5]);
	  if(statresult[1][5] < PTHRESH || statresult[2][5] < PTHRESH
	     || phat[2] > 0.75 ) {
	    ++countmlrt;
	    if( phat[2] > 0.75 ) {
	      sprintf(ctemp1, "%%%dd %%%dd %%%dd  %%2d  %%%dd ", mpedid,
		      mid, mid, mmark);
	      fprintf(fp_out3, ctemp1,  
		      nomem[i][1], rel[i][j][1], rel[i][j][2],
		      rel[i][j][3], commark);
	      fprintf(fp_out3,"  %6.4f %6.4f %6.4f  %6.4f %6.4f %6.4f\n",
		      p012[reltype][1],p012[reltype][2],p012[reltype][3],
		      phat[1], phat[2], phat[3]);
	      fflush(fp_out3);
	    }
	  }
	}
	
	/* for checking allele freq. */
	for (k = 1; k <= 3; k++) {
	  checkallefreq[k][reltype][1] += statresult[k][1];
	  checkallefreq[k][reltype][2] += statresult[k][2];
	  checkallefreq[k][reltype][3] += statresult[k][2]*statresult[k][2];
	  checkallefreq[k][reltype][4] += statresult[k][3];
	  checkallefreq[k][reltype][5] += statresult[k][3]*statresult[k][3];
	}
	checkmarker[reltype][1] += commark;
	checkmarker[reltype][2] += commark*commark;
	

	/* if user want to do MLRT test */
	if ( flag == 2 ) {
	  
	  if (reltype == 6) {
	    if ((statresult[3][5] < PTHRESH && statresult[3][4] > 0)
		|| phat[2] > 0.75)
	      domlrt = 1;
	  }
	  else if(reltype == 10 && phat[2] < 0.9 )
	    domlrt = 1;
	  else if(reltype == 11 && phat[3] < 0.95)
	    domlrt = 1;
	  else if(!(reltype == 6 || reltype == 10 || reltype == 11)) {
	    /* for the others, if use EIBD and AIBS, and if not P.O.
	       but p1 > 0.9, also do check */
	    if (statresult[1][5] < PTHRESH 
		|| statresult[2][5] < PTHRESH 
		|| phat[2] > 0.75)
	      domlrt = 1;
	  }
	  
	  /* if the pair is qualified */
	  if ( domlrt == 1 ) {

	    for (k = 1; k <= 4; k++)
	      countbig[k] = 0;

	    ++nopairsvalid[i][NR+1];
	    ++nopairsvalid[noped+1][NR+1];
	    for (k = 1; k <= NR; k++) {
	      if (reltype == k) {
		++nopairsvalid[i][k];
		++nopairsvalid[noped+1][k];
	      }
	    }

	    /* observed EIBD, AIBS, IBS statistics */ 
	    for (k = 1; k <= 3; k++)
	      stat_obs[k] = statresult[k][1];
	  
	    /* observed log( likelihood ) under each reltype */
	    getapprloglikelihood(geno[index1], geno[index2], 
				 valid, q, cenm, nochrom, nomark, maptype,
				 loglh_obs);

	    /* MLRT = max(log(likelihood_alt)) -
	       log(likelihood_null) */	      
	    if ( reltype != 1 ) 
	      dtemp[1] = loglh_obs[1];
	    else 
	      dtemp[1] = loglh_obs[2];
	    for (k = 1; k <= NRMLRT; k++) {
	      if ( k != reltype ) {
		if ( loglh_obs[k] > dtemp[1]) 
		  dtemp[1] = loglh_obs[k];
	      }
	    }
	    mlrt_obs = dtemp[1] - loglh_obs[reltype];
	    
	    /* simulate geotype data for the pair and calculate the
	       statistics */
	    for (n = 1; n <= NREP; n++) {

	      simugeno1geno2(q, qcum, nochrom, nomark, noalle, reltype, 
			     chromlength, cenm, simugeno1, simugeno2);
	    
	      get3stat(simugeno1, simugeno2, valid, 
		       commark, reltype, nochrom, nomark, 
		       p012[reltype], q, stat_simu);

	      /* delete later */
	      /*  checkmu[reltype] += stat_simu[1]; */
	   
	      getapprloglikelihood(simugeno1, simugeno2, 
				   valid, q, cenm, nochrom, nomark, maptype,
				   loglh_simu);

	      if (reltype != 1) 
		dtemp[2] = loglh_simu[1];
	      else 
		dtemp[2] = loglh_simu[2];
	      for (k = 1; k <= NRMLRT; k++) {
		if (k != reltype) {
		  if ( loglh_simu[k] > dtemp[2]) 
		    dtemp[2] = loglh_simu[k];
		}
	      }
	      mlrt_simu = dtemp[2] - loglh_simu[reltype];
	      
	    	
	      /* keep track of no. obs < simu or obs > simu  for the
		 calculation of empirical p-value */
	      for (k = 1; k <= 3; k++) {
		if (stat_simu[k] > stat_obs[k])
		  ++countbig[k];
	      }

	      if ( mlrt_simu > mlrt_obs )
		++countbig[4];
	      	   
	    }/* loop for NREP replicates simulation */
	    	
	    /* EIBD, AIBS are not applicable for unrelated = 6 */
	    if ( ! (reltype == 6 || reltype == 10 || reltype == 11)) {
	      for (k = 1; k <= 4; k++) {
		if ( countbig[k] > NREP / 2.0) {
		  fprintf(fp_out2, " %8.6f ", 2*((NREP-countbig[k])/(NREP+0.0)));
		}
		else{
		  fprintf(fp_out2, " %8.6f ", 2*(countbig[k]/(NREP+0.0)));
		}
	      }
	    }
	    else if(reltype == 6){
	      for (k = 1; k <= 2; k++) 
		fprintf(fp_out2, "       NA ");
	      for (k = 3; k <= 4; k++) {
		if ( countbig[k] > NREP / 2.0) {
		  fprintf(fp_out2, " %8.6f ", 2*((NREP-countbig[k])/(NREP+0.0)));
		}
		else{
		  fprintf(fp_out2, " %8.6f ", 2*(countbig[k]/(NREP+0.0)));
		}
	      }
	    }
	    else if(reltype == 10 || reltype == 11){
	      for (k = 1; k <= 1; k++) 
		fprintf(fp_out2, "       NA ");
	      for (k = 2; k <= 4; k++) {
		if ( countbig[k] > NREP / 2.0) {
		  fprintf(fp_out2, " %8.6f ", 2*((NREP-countbig[k])/(NREP+0.0)));
		}
		else{
		  fprintf(fp_out2, " %8.6f ", 2*(countbig[k]/(NREP+0.0)));
		}
	      }
	    }
	    
	    fprintf(fp_out2, "\n");
	    fflush(fp_out2);
	
	  } 
	  
	  /* if the pair is not qualified */
	  else {
	    for (k = 1; k <= 4; k++) 
	      fprintf(fp_out2,"       NA ");
	    fprintf(fp_out2,"\n");
	    fflush(fp_out2);
	  }
	  
	}

	/* if user does not want MLRT */	
	else {
	  fprintf(fp_out2,"\n");
	  fflush(fp_out2);
	}
	
      }
    } 
    
    printf("\npedigree %d done",nomem[i][1]);
    
  }
  printf("\n\n");

  if (mederror != 0) {
    printf("\n Prest has detected Mendelian errors\n");
    printf("   which are reported in the file prest_errors.\n");
    printf("\n In the current implementation of prest, the presence of Mendelian errors\n");
    printf("   is not taken into account in the hypothesis testing.\n");
  }
  

/*      printf("\n There are Mendelian Errors in genodata\n"); */
/*      printf("\n The testing results are not reliable\n"); */
/*      printf("\n Check prest_errors file, correct the Mendelian Errors and re-run the program\n\n"); */
/*    } */
 
  for (k = 1; k <= 3; k++) {
    for (j = 1; j <= NR; j++) {
      checkallefreq[k][j][1] = 
	checkallefreq[k][j][1]/norelwdata[noped+1][j];
      checkallefreq[k][j][2] =
	checkallefreq[k][j][2]/norelwdata[noped+1][j];
      checkallefreq[k][j][4] =
	checkallefreq[k][j][4]/norelwdata[noped+1][j];
      checkallefreq[k][j][3] = 
	sqrt((checkallefreq[k][j][3] -
	      norelwdata[noped+1][j]*checkallefreq[k][j][2]*checkallefreq[k][j][2])/(norelwdata[noped+1][j]-1));
      checkallefreq[k][j][5] = 
	sqrt((checkallefreq[k][j][5] -
	      norelwdata[noped+1][j]*checkallefreq[k][j][4]*checkallefreq[k][j][4])/(norelwdata[noped+1][j]-1));
    }
  }
  for (j = 1; j <= NR; j++) {
    checkmarker[j][1] = checkmarker[j][1]/norelwdata[noped+1][j];
    checkmarker[j][2] = sqrt((checkmarker[j][2] -
			      norelwdata[noped+1][j]*checkmarker[j][1]*checkmarker[j][1])/(norelwdata[noped+1][j]-1));
  }
  
  /* delete later */
  /*  if (flag == 2) { */
/*      for (i = 1; i <= NR; i++) { */
/*        checkmu[i] = checkmu[i]/(nopairsvalid[noped+1][i]*NREP); */
/*        printf("%f\n",checkmu[i]); */
/*      } */
/*    } */


  /* -------------------------------------------------------------- */

  /* output general information */ 

  fprintf(fp_out1, "# Number of pedigrees: \n");
  fprintf(fp_out1, "    %d\n\n", noped);

  fprintf(fp_out1, "# Total number of individuals: \n");
  fprintf(fp_out1, "    %d\n\n", nomem[noped+1][2]);

  fprintf(fp_out1, "# Number of chromosomes: \n");
  fprintf(fp_out1, "    %d\n\n", nochrom);
 
  fprintf(fp_out1, "# Total number of markers: \n");
  fprintf(fp_out1, "    %d\n\n", nomark[0]);
  
  fprintf(fp_out1, "# Number of markers on each chromosome\n\n");
  for (i = 1; i <= nochrom; i++) 
    fprintf(fp_out1, "    %s: %5d \n", chromnameidx[i], nomark[i]);
  fprintf(fp_out1, "\n");
  
  fprintf(fp_out1, "# Number of genotyped pairs found in each of the\n");
  fprintf(fp_out1, "# relationship categories analyzed by PREST: \n");
  fprintf(fp_out1, "# (Note that unrelated pairs are checked within each pedigree.)\n\n");
  
  fprintf(fp_out1, "    Total:  %d\n\n", norelwdata[noped+1][NR+1]);
  fprintf(fp_out1, "    %7d full-sib pairs\n",
	  norelwdata[noped+1][1]);
  fprintf(fp_out1, "    %7d half-sib pairs \n",
	  norelwdata[noped+1][2]);
  fprintf(fp_out1, "    %7d grandparent-grandchild pairs\n",
	  norelwdata[noped+1][3]);
  fprintf(fp_out1, "    %7d avuncular pairs\n",
	  norelwdata[noped+1][4]);
  fprintf(fp_out1, "    %7d first-cousin pairs\n",
	  norelwdata[noped+1][5]);
  fprintf(fp_out1, "    %7d unrelated pairs\n",
	  norelwdata[noped+1][6]);
  fprintf(fp_out1, "    %7d half-avuncular pairs\n",
	  norelwdata[noped+1][7]);
  fprintf(fp_out1, "    %7d half-first-cousin pairs\n",
	  norelwdata[noped+1][8]);
  fprintf(fp_out1, "    %7d half-sib-plus-first-cousin pairs\n",
	  norelwdata[noped+1][9]);
  fprintf(fp_out1, "    %7d parent-offspring pairs\n\n",
	  norelwdata[noped+1][10]);

  fprintf(fp_out1, "# Summary of genotyped pairs found in each pedigree\n\n");
  fprintf(fp_out1, "    pedid size f.sib h.sib g.p.c avun. f.cous. unrel. h.avun h.f.cous. h.sib+f.cous. p.o.   Total\n"); 
  for (i = 1; i <= noped; i++) {
    fprintf(fp_out1, " %8d %4d ", nomem[i][1], nomem[i][2]);
    for (j = 1; j <= NR; j++) 
      fprintf(fp_out1, " %4d", norelwdata[i][j]);
    fprintf(fp_out1, " %8d\n", norelwdata[i][NR+1]);
  }
  fprintf(fp_out1, "\n");
  

  if( flag == 1) {
    fprintf(fp_out1,"# Number of pairs that have p-value of either EIBD or AIBS less than 0.2\n",countmlrt);
    fprintf(fp_out1,"# (for unrelated, p-value of IBS less than 0.2).\n");
    fprintf(fp_out1,"# This is the number of pairs for which PREST would perform the MLRT.\n");
    fprintf(fp_out1,"# (The MLRT will actually be performed only if option 2 is chosen.)\n");
    fprintf(fp_out1,"# For each pair, performance of the MLRT, including 100,000\n");
    fprintf(fp_out1,"# simulated realizations takes approximately 5 minutes on \n");
    fprintf(fp_out1,"# a Sun Ultra II with 360-MHz processor.\n\n");
    fprintf(fp_out1,"    %d\n\n", countmlrt);
  }
  

  if ( flag == 2 ) {
    fprintf(fp_out1, "# Number of pairs checked by MLRT\n\n");
    fprintf(fp_out1, "    Total:  %d\n\n", nopairsvalid[noped+1][NR+1]);
    fprintf(fp_out1, "    %7d full-sib pairs\n",
	    nopairsvalid[noped+1][1]);
    fprintf(fp_out1, "    %7d half-sib pairs\n",
	    nopairsvalid[noped+1][2]);
    fprintf(fp_out1, "    %7d grandparent-grandchild pairs\n",
	    nopairsvalid[noped+1][3]);
    fprintf(fp_out1, "    %7d avuncular pairs\n",
	    nopairsvalid[noped+1][4]);
    fprintf(fp_out1, "    %7d first-cousin pairs\n",
	    nopairsvalid[noped+1][5]);
    fprintf(fp_out1, "    %7d unrelated pairs\n",
	    nopairsvalid[noped+1][6]);
    fprintf(fp_out1, "    %7d half-avuncular pairs\n",
	    nopairsvalid[noped+1][7]);
    fprintf(fp_out1, "    %7d half-first-cousin pairs\n",
	    nopairsvalid[noped+1][8]);
    fprintf(fp_out1, "    %7d half-sib-plus-first-cousin pairs\n",
	    nopairsvalid[noped+1][9]);
    fprintf(fp_out1, "    %7d parent-offspring pairs\n\n",
	    nopairsvalid[noped+1][10]);

    fprintf(fp_out1, "    In each pedigree\n\n");
    fprintf(fp_out1, "    pedid size f.sib h.sib g.p.c avun. f.cous. unrel. h.avun h.f.cous. h.sib+f.cous. p.o.   Total\n"); 
    for (i = 1; i <= noped; i++) {
      fprintf(fp_out1, " %8d %4d ", nomem[i][1], nomem[i][2]);
      for (j = 1; j <= NR; j++) 
	fprintf(fp_out1, " %4d", nopairsvalid[i][j]);
      fprintf(fp_out1, " %8d\n", nopairsvalid[i][NR+1]);
    }
    fprintf(fp_out1, "\n");
  }

  fprintf(fp_out1, "# Assuming that there are not many misclassified relative pairs in the\n");
  fprintf(fp_out1, "# data set and assuming that all the pairs are typed on roughly the\n");
  fprintf(fp_out1, "# same markers (because which markers are typed affects the null means\n");
  fprintf(fp_out1, "# of the AIBS and IBS statistics and affects the null variances of the\n");
  fprintf(fp_out1, "# EIBD, AIBS and IBS statistics), then the sample averaged observed\n");
  fprintf(fp_out1, "# statistics for each relationship type should not be significantly\n");
  fprintf(fp_out1, "# different from what is expected.  If they are significantly\n");
  fprintf(fp_out1, "# different, if may indicate that the allele frequencies are not\n");
  fprintf(fp_out1, "# accurately estimated, or that the assumptions of population\n");
  fprintf(fp_out1, "# homogeneity and Hardy-Weinberg equilibrium are not appropriate.\n");
  fprintf(fp_out1, "# It could also indicate a problem with the way the data are coded \n");
  fprintf(fp_out1, "# (see item #7 in 'Tips' file).  The following output contains \n");
  fprintf(fp_out1, "# the information of sample summary of the statistics\n\n");

  fprintf(fp_out1, "# obs      = sample averaged observed statistic\n");
  fprintf(fp_out1, "# mu(std)  = sample averaged null mean (sample std. dev. of null mean)\n");
  fprintf(fp_out1, "#            (Note that the null mean of EIBD depends only on\n");
  fprintf(fp_out1, "             the null relationship, not on the data, so for EIBD \n");
  fprintf(fp_out1, "#            we give the null mean in place of the sample averaged null mean\n");
  
  fprintf(fp_out1, "# std(std) = sample averaged standard deviation\n");
  fprintf(fp_out1, "#            (sample std. dev. of null standard deviation)\n");
  fprintf(fp_out1, "# no.mark  = sample averaged number of typed markers for a pair\n");
  fprintf(fp_out1, "#            (sample std. dev.)\n\n");

  fprintf(fp_out1, "              no.pairs no.mark(std)\n");

  if (norelwdata[noped+1][1] == 0)
    fprintf(fp_out1, "    full-sib   %7d      NA(    NA)\n",
	    norelwdata[noped+1][1]);
  else 
    fprintf(fp_out1, "    full-sib   %7d  %6.1f(%6.1f)\n",
	    norelwdata[noped+1][1], checkmarker[1][1],
	    checkmarker[1][2]);

  if (norelwdata[noped+1][2] == 0)
    fprintf(fp_out1, "    half-sib   %7d      NA(    NA)\n",
	    norelwdata[noped+1][2]);
  else
    fprintf(fp_out1, "    half-sib   %7d  %6.1f(%6.1f)\n",
	    norelwdata[noped+1][2], checkmarker[2][1],
	    checkmarker[2][2]);
  
  if (norelwdata[noped+1][3] == 0)
    fprintf(fp_out1, "    grandp-c   %7d      NA(    NA)\n",
	    norelwdata[noped+1][3]);
  else
    fprintf(fp_out1, "    grandp-c   %7d  %6.1f(%6.1f)\n",
	    norelwdata[noped+1][3], checkmarker[3][1],
	    checkmarker[3][2]);

  if (norelwdata[noped+1][4] == 0)
    fprintf(fp_out1, "   avuncular   %7d      NA(    NA)\n",
	    norelwdata[noped+1][4]);
  else
    fprintf(fp_out1, "   avuncular   %7d  %6.1f(%6.1f)\n",
	    norelwdata[noped+1][4], checkmarker[4][1],
	    checkmarker[4][2]);

  if (norelwdata[noped+1][5] == 0)
    fprintf(fp_out1, "    f-cousin   %7d      NA(    NA)\n",
	    norelwdata[noped+1][5]);
  else
    fprintf(fp_out1, "    f-cousin   %7d  %6.1f(%6.1f)\n",
	    norelwdata[noped+1][5], checkmarker[5][1],
	    checkmarker[5][2]);

  if (norelwdata[noped+1][6] == 0)
    fprintf(fp_out1, "   unrelated   %7d      NA(    NA)\n",
	    norelwdata[noped+1][6]);
  else
    fprintf(fp_out1, "   unrelated   %7d  %6.1f(%6.1f)\n",
	    norelwdata[noped+1][6], checkmarker[6][1],
	    checkmarker[6][2]);
  
  if (norelwdata[noped+1][7] == 0)
    fprintf(fp_out1, " h-avuncular   %7d      NA(    NA)\n",
	    norelwdata[noped+1][7]);
  else
    fprintf(fp_out1, " h-avuncular   %7d  %6.1f(%6.1f)\n",
	    norelwdata[noped+1][7], checkmarker[8][1],
	    checkmarker[7][2]);

  if (norelwdata[noped+1][8] == 0)
    fprintf(fp_out1, "  h-f-cousin   %7d      NA(    NA)\n",
	    norelwdata[noped+1][8]);
  else
    fprintf(fp_out1, "  h-f-cousin   %7d  %6.1f(%6.1f)\n",
	    norelwdata[noped+1][8], checkmarker[8][1],
	    checkmarker[8][2]);

  if (norelwdata[noped+1][9] == 0)
    fprintf(fp_out1, "h-sib+f-cous   %7d      NA(    NA)\n",
	    norelwdata[noped+1][9]);
  else
    fprintf(fp_out1, "h-sib+f-cous   %7d  %6.1f(%6.1f)\n",
	    norelwdata[noped+1][9], checkmarker[9][1],
	    checkmarker[9][2]);

  if (norelwdata[noped+1][10] == 0)
    fprintf(fp_out1, "pa-offspring   %7d      NA(    NA)\n",
	    norelwdata[noped+1][10]);
  else
    fprintf(fp_out1, "pa-offspring   %7d  %6.1f(%6.1f)\n\n",
	    norelwdata[noped+1][10], checkmarker[10][1],
	    checkmarker[10][2]);

  fprintf(fp_out1, "# Sample summary of EIBD \n\n");
    fprintf(fp_out1, "                   obs     mu    std(     std)\n");

  if (norelwdata[noped+1][1] == 0)
    fprintf(fp_out1, "       full-sib     NA     NA     NA(      NA)\n");
  else if (norelwdata[noped+1][1] == 1)
    fprintf(fp_out1, "       full-sib %6.4f %6.4f %6.4f(      NA)\n",
	    checkallefreq[1][1][1], checkallefreq[1][1][2], 
	    checkallefreq[1][1][4]);
  else
    fprintf(fp_out1, "       full-sib %6.4f %6.4f %6.4f(%8.6f)\n",
	    checkallefreq[1][1][1], 
	    checkallefreq[1][1][2],  
	    checkallefreq[1][1][4], checkallefreq[1][1][5]);

  if (norelwdata[noped+1][2] == 0)
    fprintf(fp_out1, "       half-sib     NA     NA     NA(      NA)\n");
  else if (norelwdata[noped+1][2] == 1)
    fprintf(fp_out1, "       half-sib %6.4f %6.4f %6.4f(      NA)\n",
	    checkallefreq[1][2][1], checkallefreq[1][2][2],
	    checkallefreq[1][2][4]);
  else
    fprintf(fp_out1, "       half-sib %6.4f %6.4f %6.4f(%8.6f)\n",
	    checkallefreq[1][2][1], 
	    checkallefreq[1][2][2], 
	    checkallefreq[1][2][4], checkallefreq[1][2][5]);

  if (norelwdata[noped+1][3] == 0)
    fprintf(fp_out1, "       grandp-c     NA     NA     NA(      NA)\n");
  else if (norelwdata[noped+1][3] == 1)
    fprintf(fp_out1, "       grandp-c %6.4f %6.4f %6.4f(      NA)\n",
	    checkallefreq[1][3][1], 
	    checkallefreq[1][3][2],checkallefreq[1][3][4]);
  else
    fprintf(fp_out1, "       grandp-c %6.4f %6.4f %6.4f(%8.6f)\n",
	    checkallefreq[1][3][1], 
	    checkallefreq[1][3][2],  
	    checkallefreq[1][3][4], checkallefreq[1][3][5]);

  if (norelwdata[noped+1][4] == 0)
    fprintf(fp_out1, "      avuncular     NA     NA     NA(      NA)\n");
  else if (norelwdata[noped+1][4] == 1)
    fprintf(fp_out1, "      avuncular %6.4f %6.4f %6.4f(      NA)\n",
	    checkallefreq[1][4][1], 
	    checkallefreq[1][4][2], checkallefreq[1][4][4]);
  else
    fprintf(fp_out1, "      avuncular %6.4f %6.4f %6.4f(%8.6f)\n",
	    checkallefreq[1][4][1], 
	    checkallefreq[1][4][2], 
	    checkallefreq[1][4][4], checkallefreq[1][4][5]);

  if (norelwdata[noped+1][5] == 0)
    fprintf(fp_out1, "       f-cousin     NA     NA     NA(      NA)\n");
  else if (norelwdata[noped+1][5] == 1)
    fprintf(fp_out1, "       f-cousin %6.4f %6.4f %6.4f(      NA)\n",
	    checkallefreq[1][5][1], 
	    checkallefreq[1][5][2],checkallefreq[1][5][4]);
  else
    fprintf(fp_out1, "       f-cousin %6.4f %6.4f %6.4f(%8.6f)\n",
	    checkallefreq[1][5][1], 
	    checkallefreq[1][5][2],
	    checkallefreq[1][5][4], checkallefreq[1][5][5]);

  fprintf(fp_out1, "      unrelated     NA     NA     NA(      NA)\n");
  
  if (norelwdata[noped+1][7] == 0)
    fprintf(fp_out1, "    h-avuncular     NA     NA     NA(      NA)\n");
  else if (norelwdata[noped+1][7] == 1)
    fprintf(fp_out1, "   h-avuncular %6.4f %6.4f %6.4f(      NA)\n",
	    checkallefreq[1][7][1], 
	    checkallefreq[1][7][2],checkallefreq[1][7][4]);
  else
    fprintf(fp_out1, "   h-avuncular  %6.4f %6.4f %6.4f(%8.6f)\n",
	    checkallefreq[1][7][1], 
	    checkallefreq[1][7][2],
	    checkallefreq[1][7][4], checkallefreq[1][7][5]);

  if (norelwdata[noped+1][8] == 0)
    fprintf(fp_out1, "     h-f-cousin     NA     NA     NA(      NA)\n");
  else if (norelwdata[noped+1][8] == 1)
    fprintf(fp_out1, "     h-f-cousin %6.4f %6.4f %6.4f(      NA)\n",
	    checkallefreq[1][8][1], 
	    checkallefreq[1][8][2],checkallefreq[1][8][4]);
  else
    fprintf(fp_out1, "     h-f-cousin %6.4f %6.4f %6.4f(%8.6f)\n",
	    checkallefreq[1][8][1], 
	    checkallefreq[1][8][2],
	    checkallefreq[1][8][4], checkallefreq[1][8][5]);

  if (norelwdata[noped+1][9] == 0)
    fprintf(fp_out1, " h-sib+f-cousin     NA     NA     NA(      NA)\n");
  else if (norelwdata[noped+1][9] == 1)
    fprintf(fp_out1, " h-sib+f-cousin %6.4f %6.4f %6.4f(      NA)\n",
	    checkallefreq[1][9][1], 
	    checkallefreq[1][9][2],checkallefreq[1][9][4]);
  else
    fprintf(fp_out1, " h-sib+f-cousin %6.4f %6.4f %6.4f(%8.6f)\n",
	    checkallefreq[1][9][1], 
	    checkallefreq[1][9][2],
	    checkallefreq[1][9][4], checkallefreq[1][9][5]);

  fprintf(fp_out1, "   pa-offspring     NA     NA     NA(      NA)\n\n");


  fprintf(fp_out1, "# Sample summary of AIBS \n\n");
  fprintf(fp_out1, "                   obs     mu(     std)    std(     std)\n");

  if (norelwdata[noped+1][1] == 0)
    fprintf(fp_out1, "       full-sib     NA     NA(      NA)     NA(      NA)\n");
  else if (norelwdata[noped+1][1] == 1)
    fprintf(fp_out1, "       full-sib %6.4f %6.4f(      NA) %6.4f(      NA)\n",
	    checkallefreq[2][1][1], checkallefreq[2][1][2],
	    checkallefreq[2][1][4]);
  else
    fprintf(fp_out1, "       full-sib %6.4f %6.4f(%8.6f) %6.4f(%8.6f)\n",
	    checkallefreq[2][1][1], 
	    checkallefreq[2][1][2], checkallefreq[2][1][3], 
	    checkallefreq[2][1][4], checkallefreq[2][1][5]);

  if (norelwdata[noped+1][2] == 0)
    fprintf(fp_out1, "       half-sib     NA     NA(      NA)     NA(      NA)\n");
  else if (norelwdata[noped+1][2] == 1)
    fprintf(fp_out1, "       half-sib %6.4f %6.4f(      NA) %6.4f(      NA)\n",
	    checkallefreq[2][2][1], 
	    checkallefreq[2][2][2],checkallefreq[2][2][4]);
  else
    fprintf(fp_out1, "       half-sib %6.4f %6.4f(%8.6f) %6.4f(%8.6f)\n",
	    checkallefreq[2][2][1], 
	    checkallefreq[2][2][2], checkallefreq[2][2][3], 
	    checkallefreq[2][2][4], checkallefreq[2][2][5]);

  if (norelwdata[noped+1][3] == 0)
    fprintf(fp_out1, "       grandp-c     NA     NA(      NA)     NA(      NA)\n");
  else if (norelwdata[noped+1][3] == 1)
    fprintf(fp_out1, "       grandp-c %6.4f %6.4f(      NA) %6.4f(      NA)\n",
	    checkallefreq[2][3][1], 
	    checkallefreq[2][3][2],checkallefreq[2][3][4]);
  else
    fprintf(fp_out1, "       grandp-c %6.4f %6.4f(%8.6f) %6.4f(%8.6f)\n",
	    checkallefreq[2][3][1], 
	    checkallefreq[2][3][2], checkallefreq[2][3][3], 
	    checkallefreq[2][3][4], checkallefreq[2][3][5]);

  if (norelwdata[noped+1][4] == 0)
    fprintf(fp_out1, "      avuncular     NA     NA(      NA)     NA(      NA)\n");
  else if (norelwdata[noped+1][4] == 1)
    fprintf(fp_out1, "      avuncular %6.4f %6.4f(      NA) %6.4f(      NA)\n",
	    checkallefreq[2][4][1], 
	    checkallefreq[2][4][2],checkallefreq[2][4][4]);
  else
    fprintf(fp_out1, "      avuncular %6.4f %6.4f(%8.6f) %6.4f(%8.6f)\n",
	    checkallefreq[2][4][1], 
	    checkallefreq[2][4][2], checkallefreq[2][4][3], 
	    checkallefreq[2][4][4], checkallefreq[2][4][5]);

  if (norelwdata[noped+1][5] == 0)
    fprintf(fp_out1, "       f-cousin     NA     NA(      NA)     NA(      NA)\n");
  else if (norelwdata[noped+1][5] == 1)
    fprintf(fp_out1, "       f-cousin %6.4f %6.4f(      NA) %6.4f(      NA)\n",
	    checkallefreq[2][5][1], 
	    checkallefreq[2][5][2],checkallefreq[2][5][4]);
  else
    fprintf(fp_out1, "       f-cousin %6.4f %6.4f(%8.6f) %6.4f(%8.6f)\n",
	    checkallefreq[2][5][1], 
	    checkallefreq[2][5][2], checkallefreq[2][5][3], 
	    checkallefreq[2][5][4], checkallefreq[2][5][5]);

  fprintf(fp_out1, "      unrelated     NA     NA(      NA)     NA(      NA)\n");

  if (norelwdata[noped+1][7] == 0)
    fprintf(fp_out1, "    h-avuncular     NA     NA(      NA)     NA(      NA)\n");
  else if (norelwdata[noped+1][7] == 1)
    fprintf(fp_out1, "    h-avuncular %6.4f %6.4f(      NA) %6.4f(      NA)\n",
	    checkallefreq[2][7][1], 
	    checkallefreq[2][7][2],checkallefreq[2][7][4]);
  else
    fprintf(fp_out1, "    h-avuncular %6.4f %6.4f(%8.6f) %6.4f(%8.6f)\n",
	    checkallefreq[2][7][1], 
	    checkallefreq[2][7][2], checkallefreq[2][7][3],
	    checkallefreq[2][7][4], checkallefreq[2][7][5]);

  if (norelwdata[noped+1][8] == 0)
    fprintf(fp_out1, "     h-f-cousin     NA     NA(      NA)     NA(      NA)\n");
  else if (norelwdata[noped+1][8] == 1)
    fprintf(fp_out1, "     h-f-cousin %6.4f %6.4f(      NA) %6.4f(      NA)\n",
	    checkallefreq[2][8][1], 
	    checkallefreq[2][8][2],checkallefreq[2][8][4]);
  else
    fprintf(fp_out1, "     h-f-cousin %6.4f %6.4f(%8.6f) %6.4f(%8.6f)\n",
	    checkallefreq[2][8][1], 
	    checkallefreq[2][8][2], checkallefreq[2][8][3],
	    checkallefreq[2][8][4], checkallefreq[2][8][5]);

  if (norelwdata[noped+1][9] == 0)
    fprintf(fp_out1, " h-sib+f-cousin     NA     NA(      NA)     NA(      NA)\n");
  else if (norelwdata[noped+1][9] == 1)
    fprintf(fp_out1, " h-sib+f-cousin %6.4f %6.4f(      NA) %6.4f(      NA)\n",
	    checkallefreq[2][9][1],  
	    checkallefreq[2][9][2],checkallefreq[2][9][4]);
  else
    fprintf(fp_out1, " h-sib+f-cousin %6.4f %6.4f(%8.6f) %6.4f(%8.6f)\n",
	    checkallefreq[2][9][1], 
	    checkallefreq[2][9][2], checkallefreq[2][9][3],
	    checkallefreq[2][9][4], checkallefreq[2][9][5]);

  if (norelwdata[noped+1][10] == 0)
    fprintf(fp_out1, "   pa-offsprin      NA     NA     NA(      NA)     NA(      NA)\n");
  else if (norelwdata[noped+1][10] == 1)
    fprintf(fp_out1, "   pa-offspring %6.4f %6.4f(      NA) %6.4f(      NA)\n",
	    checkallefreq[2][10][1], 
	    checkallefreq[2][10][2],checkallefreq[2][10][4]);
  else
    fprintf(fp_out1, "   pa-offspring %6.4f %6.4f(%8.6f) %6.4f(%8.6f)\n\n",
	    checkallefreq[2][10][1], 
	    checkallefreq[2][10][2], checkallefreq[2][10][3],
	    checkallefreq[2][10][4], checkallefreq[2][10][5]);


  fprintf(fp_out1, "# Sample summary of IBS \n\n");
  fprintf(fp_out1, "                   obs     mu(     std)    std(     std)\n");
  
  if (norelwdata[noped+1][1] == 0)
    fprintf(fp_out1, "       full-sib     NA     NA(      NA)     NA(      NA)\n");
  else if (norelwdata[noped+1][1] == 1)
    fprintf(fp_out1, "       full-sib %6.4f %6.4f(      NA) %6.4f(      NA)\n",
	    checkallefreq[3][1][1], 
	    checkallefreq[3][1][2],checkallefreq[3][1][4]);
  else
    fprintf(fp_out1, "       full-sib %6.4f %6.4f(%8.6f) %6.4f(%8.6f)\n",
	    checkallefreq[3][1][1], 
	    checkallefreq[3][1][2], checkallefreq[3][1][3], 
	    checkallefreq[3][1][4], checkallefreq[3][1][5]);

  if (norelwdata[noped+1][2] == 0)
    fprintf(fp_out1, "       half-sib     NA     NA(      NA)     NA(      NA)\n");
  else if (norelwdata[noped+1][2] == 1)
    fprintf(fp_out1, "       half-sib %6.4f %6.4f(      NA) %6.4f(      NA)\n",
	    checkallefreq[3][2][1], 
	    checkallefreq[3][2][2],checkallefreq[3][2][4]);
  else
    fprintf(fp_out1, "       half-sib %6.4f %6.4f(%8.6f) %6.4f(%8.6f)\n",
	    checkallefreq[3][2][1], 
	    checkallefreq[3][2][2], checkallefreq[3][2][3], 
	    checkallefreq[3][2][4], checkallefreq[3][2][5]);

  if (norelwdata[noped+1][3] == 0)
    fprintf(fp_out1, "       grandp-c     NA     NA(      NA)     NA(      NA)\n");
  else if (norelwdata[noped+1][3] == 1)
    fprintf(fp_out1, "       grandp-c %6.4f %6.4f(      NA) %6.4f(      NA)\n",
	    checkallefreq[3][3][1], 
	    checkallefreq[3][3][2],checkallefreq[3][3][4]);
  else
    fprintf(fp_out1, "       grandp-c %6.4f %6.4f(%8.6f) %6.4f(%8.6f)\n",
	    checkallefreq[3][3][1], 
	    checkallefreq[3][3][2], checkallefreq[3][3][3], 
	    checkallefreq[3][3][4], checkallefreq[3][3][5]);

  if (norelwdata[noped+1][4] == 0)
    fprintf(fp_out1, "      avuncular     NA     NA(      NA)     NA(      NA)\n");
  else if (norelwdata[noped+1][4] == 1)
    fprintf(fp_out1, "      avuncular %6.4f %6.4f(      NA) %6.4f(      NA)\n",
	    checkallefreq[3][4][1], 
	    checkallefreq[3][4][2], checkallefreq[3][4][4]);
  else
    fprintf(fp_out1, "      avuncular %6.4f %6.4f(%8.6f) %6.4f(%8.6f)\n",
	    checkallefreq[3][4][1], 
	    checkallefreq[3][4][2], checkallefreq[3][4][3], 
	    checkallefreq[3][4][4], checkallefreq[3][4][5]);

  if (norelwdata[noped+1][5] == 0)
    fprintf(fp_out1, "       f-cousin     NA     NA(      NA)     NA(      NA)\n");
  else if (norelwdata[noped+1][5] == 1)
    fprintf(fp_out1, "      f -cousin %6.4f %6.4f(      NA) %6.4f(      NA)\n",
	    checkallefreq[3][5][1], 
	    checkallefreq[3][5][2],checkallefreq[3][5][4]);
  else
    fprintf(fp_out1, "       f-cousin %6.4f %6.4f(%8.6f) %6.4f(%8.6f)\n",
	    checkallefreq[3][5][1], 
	    checkallefreq[3][5][2], checkallefreq[3][5][3], 
	    checkallefreq[3][5][4], checkallefreq[3][5][5]);

  if (norelwdata[noped+1][6] == 0)
    fprintf(fp_out1, "      unrelated     NA     NA(      NA)     NA(      NA)\n");
  else if (norelwdata[noped+1][6] == 1)
    fprintf(fp_out1, "      unrelated %6.4f %6.4f(      NA) %6.4f(      NA)\n",
	    checkallefreq[3][6][1], 
	    checkallefreq[3][6][2],checkallefreq[3][6][4]);
  else
    fprintf(fp_out1, "      unrelated %6.4f %6.4f(%8.6f) %6.4f(%8.6f)\n",
	    checkallefreq[3][6][1], 
	    checkallefreq[3][6][2], checkallefreq[3][6][3], 
	    checkallefreq[3][6][4], checkallefreq[3][6][5]);

  if (norelwdata[noped+1][7] == 0)
    fprintf(fp_out1, "    h-avuncular     NA     NA(      NA)     NA(      NA)\n");
  else if (norelwdata[noped+1][7] == 1)
    fprintf(fp_out1, "   h-avuncular %6.4f %6.4f(      NA) %6.4f(      NA)\n",
	    checkallefreq[3][7][1], 
	    checkallefreq[3][7][2],checkallefreq[3][7][4]);
  else
    fprintf(fp_out1, "   h-avuncular  %6.4f %6.4f(%8.6f) %6.4f(%8.6f)\n",
	    checkallefreq[3][7][1], 
	    checkallefreq[3][7][2], checkallefreq[3][7][3],
	    checkallefreq[3][7][4], checkallefreq[3][7][5]);

  if (norelwdata[noped+1][8] == 0)
    fprintf(fp_out1, "     h-f-cousin     NA     NA(      NA)     NA(      NA)\n");
  else if (norelwdata[noped+1][8] == 1)
    fprintf(fp_out1, "     h-f-cousin %6.4f %6.4f(      NA) %6.4f(      NA)\n",
	    checkallefreq[3][8][1], 
	    checkallefreq[3][8][2],checkallefreq[3][8][4]);
  else
    fprintf(fp_out1, "     h-f-cousin %6.4f %6.4f(%8.6f) %6.4f(%8.6f)\n",
	    checkallefreq[3][8][1], 
	    checkallefreq[3][8][2], checkallefreq[3][8][3],
	    checkallefreq[3][8][4], checkallefreq[3][8][5]);

  if (norelwdata[noped+1][9] == 0)
    fprintf(fp_out1, " h-sib+f-cousin     NA     NA(      NA)     NA(      NA)\n");
  else if (norelwdata[noped+1][9] == 1)
    fprintf(fp_out1, " h-sib+f-cousin %6.4f %6.4f(      NA) %6.4f(      NA)\n",
	    checkallefreq[3][9][1], 
	    checkallefreq[3][9][2],checkallefreq[3][9][4]);
  else
    fprintf(fp_out1, " h-sib+f-cousin %6.4f %6.4f(%8.6f) %6.4f(%8.6f)\n",
	    checkallefreq[3][9][1], 
	    checkallefreq[3][9][2], checkallefreq[3][9][3],
	    checkallefreq[3][9][4], checkallefreq[3][9][5]);

  if (norelwdata[noped+1][10] == 0)
    fprintf(fp_out1, "   pa-offsprin      NA     NA(      NA)     NA(      NA)\n");
  else if (norelwdata[noped+1][10] == 1)
    fprintf(fp_out1, "   pa-offspring  %6.4f %6.4f(      NA) %6.4f(      NA)\n",
	    checkallefreq[3][10][1], 
	    checkallefreq[3][10][2],checkallefreq[3][10][4]);
  else
    fprintf(fp_out1, "   pa-offspring %6.4f %6.4f(%8.6f) %6.4f(%8.6f)\n",
	    checkallefreq[3][10][1], 
	    checkallefreq[3][10][2], checkallefreq[3][10][3],
	    checkallefreq[3][10][4], checkallefreq[3][10][5]);

  /* -------------------------------------------------------------- */

  free(ped);
  free(nomem);
  free(rel);
  free(norel);
  free(norelwdata);
 
  free(q);
  free(qcum);
  free(theta);
  free(cenm);
  free(nomark);
  free(noalle);
  free(chromnameped);
  free(chromnameidx);
  free(chromlength);

  free(valid);
  free(ngeno);
  free(genotype);
  free(eeibdibd);
  free(eeibdeibdibd);
  free(eaibsibd);
  free(eaibsaibsibd);
  free(eibsibd);
  free(eibsibsibd);

  free(geno);
  free(simugeno1);
  free(simugeno2);
  free(nopairsvalid);

  fclose(fp_in2);
  fclose(fp_out1);
  fclose(fp_out2);
  fclose(fp_out3);
  fclose(fp_error);

  return 1;
  
} /* end of main */


/* for find all the relationships considered */

int findrelation(int** ped, int pedsize, int** rel)
{
  int i, j, findit;
  int nofr, inbred;
  int un, ne, np, nup;
  int gp, gc, gcp, gcup;
  int com, ucom;
  int c1p, c2p, c1up, c2up;
  int reltype;
  
  nofr = 0;
  for (i = 1; i <= (pedsize-1) ; i++) {
    for (j = (i+1); j <= pedsize; j++) {
      /* find relationship between i and j */
      findit = 0;
      
      /* parent-offspring */
      /* j is the parent of i, and j is not inbred */
      if ( ped[i][4] == j || ped[i][5] == j ) {
	if( inbred_check(ped, pedsize, j) == 0 ) {
	  reltype = 10;
	  findit = 1;
	}
      }
      if (!findit ) {
	/* i is the parent of j, and i is not inbred */
	if ( ped[j][4] == i || ped[j][5] == i ) {
	  if( inbred_check(ped, pedsize, i) == 0 ) {
	    findit = 1;
	    reltype = 10;
	  }
	}
      }
      /* not parent-offspring */
      if ( (!findit) && (! (ped[i][4] == j || ped[i][5] == j ||
			    ped[j][4] == i || ped[j][5] == i) )) {
	/* grandparent-grandchild */
	if ( findit != 1 ) {
	  /* i is the grandparent, through j'th father */
	  gp = i;
	  gc = j;
	  gcp = 4;
	  gcup = 5;
	  findit = check_grandpc(gp, gc, gcp, gcup, ped, pedsize);
	  if ( findit == 1 ) 
	    reltype = 3;
	}
	if ( findit != 1 ) {
	  /* i is the grandparent, through j'th mother */
	  gp = i;
	  gc = j;
	  gcp = 5;
	  gcup = 4;
	  findit = check_grandpc(gp, gc, gcp, gcup, ped, pedsize);
	  if ( findit == 1 ) 
	    reltype = 3;
	}
	if ( findit != 1 ) {
	  /* j is the grandparent, through i'th father */
	  gp = j;
	  gc = i;
	  gcp = 4;
	  gcup = 5;
	  findit = check_grandpc(gp, gc, gcp, gcup, ped, pedsize);
	  if ( findit == 1 ) 
	    reltype = 3;
	}  
	if ( findit != 1 ) {
	  /* j is the grandparent, through j'th mother */
	  gp = j;
	  gc = i;
	  gcp = 5;
	  gcup = 4;
	  findit = check_grandpc(gp, gc, gcp, gcup, ped, pedsize);
	  if ( findit == 1 ) 
	    reltype = 3;
	}
      	/* unrelated */
	if ( findit != 1 ) {
	  inbred = check(ped, pedsize, i, j);
	  if ( inbred == 0 ) {
	    findit = 1;
	    reltype = 6;
	  }
	}
      	/* for full-sib, half-sib, avuncular, first-cousin, 
	   half-plus-cousin relative pairs, no founder should be
	   involved */
	if ( findit != 1 && 
	     ped[i][4] != 0 && ped[i][5] != 0 && 
	     ped[j][4] != 0 && ped[j][4] != 0 ) {
	  /* full-sib 
	     parents are the same and parents are unrelated,
	     and parents are not inbred themselves*/
	  if ( (ped[i][4] == ped[j][4]) && (ped[i][5] == ped[j][5]) ) {
	    inbred = check(ped, pedsize, ped[i][4], ped[i][5]);
	    inbred = inbred + inbred_check(ped, pedsize, ped[i][4]);
	    inbred = inbred + inbred_check(ped, pedsize, ped[i][5]);
	    if ( inbred == 0 ) {
	      findit = 1;
	      reltype = 1;
	    }
	  }
	  /* avuncular e.g. uncle-nephew */
	  if ( findit != 1 ) {
	    /* if i is the uncle and sibling with j'th father */
	    un = i; 
	    ne = j;
	    np = 4;
	    nup = 5;
	    findit = check_avuncular(un, ne, np, nup, ped, pedsize);
	    if ( findit == 1 ) 
	      reltype = 4;
	  }
	  if ( findit != 1 ) {
	    /* if i is the uncle and sibling with j'th mother */
	    un = i; 
	    ne = j;
	    np = 5;
	    nup = 4;
	    findit = check_avuncular(un, ne, np, nup, ped, pedsize);
	    if ( findit == 1 ) 
	      reltype = 4;
	  }
	  if ( findit != 1 ) {
	    /* if j is the uncle and sibling with j'th father */
	    un = j; 
	    ne = i;
	    np = 4;
	    nup = 5;
	    findit = check_avuncular(un, ne, np, nup, ped, pedsize);
	    if ( findit == 1 ) 
	      reltype = 4;
	  }
	  if ( findit != 1 ) {
	    /* if j is the uncle and sibling with j'th mother */
	    un = j; 
	    ne = i;
	    np = 5;
	    nup = 4;
	    findit = check_avuncular(un, ne, np, nup, ped, pedsize);
	    if ( findit == 1 ) 
	      reltype = 4;
	  }
	  /* half-sib */
	  if ( findit != 1 ) {
	    /* same mother */
	    com = 4;
	    ucom = 5;
	    findit = check_halfsib(i, j, com, ucom, ped, pedsize);
	    if ( findit == 1 ) 
	      reltype = 2;
	  }
	  if (findit != 1) {
	    /* same father */
	    com = 5;
	    ucom = 4;
	    findit = check_halfsib(i, j, com, ucom, ped, pedsize);
	    if ( findit == 1 ) 
	      reltype = 2;
	  }

	  /* first-cousin */
	  if ( findit != 1 ) {
	    /* i's father is sibling with j's father */
	    c1p = 4;
	    c2p = 4;
	    c1up = 5;
	    c2up = 5;
	    findit = check_cousin(i, j, c1p, c2p, c1up, c2up, 
				  ped, pedsize);
	    if ( findit == 1 ) 
	      reltype = 5;
	  }
	  if ( findit != 1 ) {
	    /* i's father is sibling with j's mother */
	    c1p = 4;
	    c2p = 5;
	    c1up = 5;
	    c2up = 4;
	    findit = check_cousin(i, j, c1p, c2p, c1up, c2up, 
				  ped, pedsize);
	    if ( findit == 1 ) 
	      reltype = 5;
	  }
	  if ( findit != 1 ) {
	    /* i's mother is sibling with j's father */
	    c1p = 5;
	    c2p = 4;
	    c1up = 4;
	    c2up = 5;
	    findit = check_cousin(i, j, c1p, c2p, c1up, c2up, 
				  ped, pedsize);
	    if ( findit == 1 ) 
	      reltype = 5;
	  }
	  if ( findit != 1 ) {
	    /* i's mother is sibling with j's mother */
	    c1p = 5;
	    c2p = 5;
	    c1up = 4;
	    c2up = 4;
	    findit = check_cousin(i, j, c1p, c2p, c1up, c2up, 
				  ped, pedsize);
	    if ( findit == 1 ) 
	      reltype = 5;
	  }

	  /* half-avuncular e.g. half-uncle-nephew */
	  if ( findit != 1 ) {
	    /* if i is the uncle and sibling with j'th father */
	    un = i; 
	    ne = j;
	    np = 4;
	    nup = 5;
	    findit = check_halfavuncular(un, ne, np, nup, ped, pedsize);
	    if ( findit == 1 ) 
	      reltype = 7;
	  }
	  if ( findit != 1 ) {
	    /* if i is the uncle and sibling with j'th mother */
	    un = i; 
	    ne = j;
	    np = 5;
	    nup = 4;
	    findit = check_halfavuncular(un, ne, np, nup, ped, pedsize);
	    if ( findit == 1 ) 
	      reltype = 7;
	  }
	  if ( findit != 1 ) {
	    /* if j is the uncle and sibling with j'th father */
	    un = j; 
	    ne = i;
	    np = 4;
	    nup = 5;
	    findit = check_halfavuncular(un, ne, np, nup, ped, pedsize);
	    if ( findit == 1 ) 
	      reltype = 7;
	  }
	  if ( findit != 1 ) {
	    /* if j is the uncle and sibling with j'th mother */
	    un = j; 
	    ne = i;
	    np = 5;
	    nup = 4;
	    findit = check_halfavuncular(un, ne, np, nup, ped, pedsize);
	    if ( findit == 1 ) 
	      reltype = 7;
	  }

	  /* half-first-cousin */
	  if ( findit != 1 ) {
	    /* i's father is halfsibling with j's father */
	    c1p = 4;
	    c2p = 4;
	    c1up = 5;
	    c2up = 5;
	    findit = check_halfcousin(i, j, c1p, c2p, c1up, c2up, 
				      ped, pedsize);
	    if ( findit == 1 ) 
	      reltype = 8;
	  }
	  if ( findit != 1 ) {
	    /* i's father is halfsibling with j's mother */
	    c1p = 4;
	    c2p = 5;
	    c1up = 5;
	    c2up = 4;
	    findit = check_halfcousin(i, j, c1p, c2p, c1up, c2up, 
				      ped, pedsize);
	    if ( findit == 1 ) 
	      reltype = 8;
	  }
	  if ( findit != 1 ) {
	    /* i's mother is halfsibling with j's father */
	    c1p = 5;
	    c2p = 4;
	    c1up = 4;
	    c2up = 5;
	    findit = check_halfcousin(i, j, c1p, c2p, c1up, c2up, 
				      ped, pedsize);
	    if ( findit == 1 ) 
	      reltype = 8;
	  }
	  if ( findit != 1 ) {
	    /* i's mother is halfsibling with j's mother */
	    c1p = 5;
	    c2p = 5;
	    c1up = 4;
	    c2up = 4;
	    findit = check_halfcousin(i, j, c1p, c2p, c1up, c2up, 
				      ped, pedsize);
	    if ( findit == 1 ) 
	      reltype = 8;
	  }

	  /* half-sib plus first cousin */
	  if (findit != 1) {
	    com = 4;
	    ucom = 5;
	    findit = 
	      check_halfpluscousin(i, j, com, ucom, ped, pedsize);
	    if (findit == 1) 
	      reltype = 9;
	  }
	  if (findit != 1) {
	    com = 5;
	    ucom = 4;
	    findit = 
	      check_halfpluscousin(i, j, com, ucom, ped, pedsize);
	    if (findit == 1) 
	      reltype = 9;
	  }
	}
      }
            
      if ( findit == 1 ) {
	++nofr;
	rel[nofr][1] = ped[i][1];
	rel[nofr][2] = ped[j][1];
	rel[nofr][3] = reltype;
      }
      
    }
  }
  
  return nofr;
}

int check_halfsib(int i, int j, int com, int ucom, int** ped, 
		  int pedsize)
{
  int findit, inbred;
  
  findit = 0;
  if ( ped[i][com] == ped[j][com] ) {
    /* the other parts unrelated, and the common parent is not inbred,
       parents are unrelated */
    inbred = check(ped, pedsize, ped[i][ucom],
				   ped[j][ucom]);
    inbred += inbred_check(ped, pedsize, ped[i][com]);
    inbred += check(ped, pedsize, ped[i][4], ped[i][5]);
    inbred += check(ped, pedsize, ped[j][4], ped[j][5]);
    if ( inbred == 0 ) 
      findit = 1;
  }

  return findit;
}

int check_grandpc(int gp, int gc, int gcp, int gcup, int** ped,
		  int pedsize)
{
  int findit, inbred;
  
  findit = 0;
  if ( ped[gc][gcp] != 0 ) {
    if ( gp == ped[ped[gc][gcp]][4] || gp == ped[ped[gc][gcp]][5] ) {
      /* grandparents should not be related, and unrelated with
	 the other parent of the grandchild */
      inbred = check(ped, pedsize,
		       ped[ped[gc][gcp]][4],ped[ped[gc][gcp]][5]);
      inbred += check(ped, pedsize, 
		      ped[ped[gc][gcp]][4],ped[gc][gcup]);
      inbred += check(ped, pedsize,
		      ped[ped[gc][gcp]][5],ped[gc][gcup]);
      if (inbred == 0) 
	findit = 1;
    }
  }
  
  return findit;
}

int check_avuncular(int un, int ne, int np, int nup, int** ped,
		    int pedsize)
{
  int findit, inbred;

  findit = 0;
  /* uncles's parents are the same as the parents of nephew's related
     parent */  
  if ( (ped[ne][np] != 0 &&
	ped[un][4] == ped[ped[ne][np]][4] &&
	ped[un][5] == ped[ped[ne][np]][5]) ) {
    /* nephew's other parent should not be related with the grandparents,
       and gradparents should not be related, and should not be inbred */
    inbred = check(ped, pedsize, ped[un][4], ped[ne][nup]);
    inbred += check(ped, pedsize, ped[un][5], ped[ne][nup]);
    inbred += check(ped, pedsize, ped[un][4], ped[un][5]);
    inbred += inbred_check(ped, pedsize, ped[un][4]);
    inbred += inbred_check(ped, pedsize, ped[un][5]);
    if (inbred == 0) 
      findit = 1;
  }

  return findit;
}

int check_halfavuncular(int un, int ne, int np, int nup, int** ped,
			int pedsize)
{
  int findit, inbred;
  int unrp, unurp;
  
  findit = 0;
  /* uncles's mother are the same as the mother of nephew's related
     parent */
  unrp = 5;
  unurp = 4;
  if ( ped[ne][np] != 0 && 
       ped[un][unrp] == ped[ped[ne][np]][unrp] ){
    /* uncles's father are unrelatd with the father of of nephew's related
       parent, and they are not inbred and nephew's unrelated
       parent are not related to them */
    inbred = check(ped, pedsize, ped[un][unurp], ped[ped[ne][np]][unurp]);
    inbred += check(ped, pedsize, ped[un][unurp], ped[un][unrp]);
    inbred += check(ped, pedsize, ped[ped[ne][np]][unurp], ped[un][unrp]);
    inbred += check(ped, pedsize, ped[un][unurp], ped[ne][nup] );
    inbred += check(ped, pedsize, ped[ped[ne][np]][unurp], ped[ne][nup] );
    inbred += check(ped, pedsize, ped[un][unrp], ped[ne][nup] );
    inbred += inbred_check(ped, pedsize, ped[un][unrp]);
    inbred += inbred_check(ped, pedsize, ped[un][unurp]);
    inbred += inbred_check(ped, pedsize, ped[ped[ne][np]][unurp]);
    if (inbred == 0) 
      findit = 1;
  }
  if(findit == 0) {
    /* uncles's father are the same as the father of nephew's related
       parent */
    unrp = 4;
    unurp = 5;
    if ( ped[ne][np] != 0 && 
	 ped[un][unrp] == ped[ped[ne][np]][unrp] ){
      /* uncles's mother are unrelatd with the mother of of nephew's related
	 parent, and they are not inbred and nephew's unrelated
	 parent are not related to them */
      inbred = check(ped, pedsize, ped[un][unurp], ped[ped[ne][np]][unurp]);
      inbred += check(ped, pedsize, ped[un][unurp], ped[un][unrp]);
      inbred += check(ped, pedsize, ped[ped[ne][np]][unurp], ped[un][unrp]);
      inbred += check(ped, pedsize, ped[un][unurp], ped[ne][nup] );
      inbred += check(ped, pedsize, ped[ped[ne][np]][unurp], ped[ne][nup] );
      inbred += check(ped, pedsize, ped[un][unrp], ped[ne][nup] );
      inbred += inbred_check(ped, pedsize, ped[un][unrp]);
      inbred += inbred_check(ped, pedsize, ped[un][unurp]);
      inbred += inbred_check(ped, pedsize, ped[ped[ne][np]][unurp]);
      if (inbred == 0) 
	findit = 1;
    }
  }
  return findit;
}



int check_cousin(int i, int j, int c1p, int c2p, int c1up, int c2up,
		 int** ped, int pedsize)
{
  int findit, inbred;
  
  findit = 0;
  if ( ped[ped[i][c1p]][4] == ped[ped[j][c2p]][4] && 
       ped[ped[i][c1p]][4] != 0 &&
       ped[ped[i][c1p]][5] == ped[ped[j][c2p]][5] &&
       ped[ped[i][c1p]][5] != 0 ) {
    /* grandparents are unrelated, non-inbred, the other parents are
       not related with the grandparents */
    inbred = check(ped, pedsize, 
		   ped[ped[i][c1p]][4], ped[ped[i][c1p]][5]);
    inbred += inbred_check(ped, pedsize, ped[ped[i][c1p]][4]);
    inbred += inbred_check(ped, pedsize, ped[ped[i][c1p]][5]);
    inbred += check(ped, pedsize, 
		    ped[i][c1up], ped[ped[i][c1p]][4]);
    inbred += check(ped, pedsize, 
		    ped[i][c1up], ped[ped[i][c1p]][5]);
    inbred += check(ped, pedsize, 
		    ped[j][c2up], ped[ped[i][c1p]][4]);
    inbred += check(ped, pedsize, 
		    ped[j][c2up], ped[ped[i][c1p]][5]);
    inbred += check(ped, pedsize, 
		    ped[i][c1up], ped[j][c2up]);
    if (inbred == 0) 
      findit = 1;
  }

  return findit;
}

int check_halfcousin(int i, int j, int c1p, int c2p, int c1up, int c2up,
		     int** ped, int pedsize)
{
  int findit, inbred;
  int rp, urp;
  
  findit = 0;
  /* related grandparent is the grandmother*/
  rp = 5;
  urp = 4;
  if ( ped[ped[i][c1p]][rp] == ped[ped[j][c2p]][rp] && 
       ped[ped[i][c1p]][rp] != 0){
    /* the other grandparents are unrelated, non-inbred, and the other
       parents are not related with each other and not related with
       the three grandparents */
    inbred = check(ped, pedsize, 
		   ped[ped[i][c1p]][urp], ped[ped[j][c2p]][urp]);
    inbred = check(ped, pedsize, 
		   ped[ped[i][c1p]][urp],ped[ped[i][c1p]][rp]);
    inbred = check(ped, pedsize, 
		   ped[ped[i][c1p]][rp], ped[ped[j][c2p]][urp]);
    inbred += inbred_check(ped, pedsize, ped[ped[i][c1p]][rp]);
    inbred += inbred_check(ped, pedsize, ped[ped[i][c1p]][urp]);
    inbred += inbred_check(ped, pedsize, ped[ped[j][c2p]][urp]);
    inbred += check(ped, pedsize, 
		    ped[i][c1up], ped[j][c2up]);
    inbred += check(ped, pedsize, 
		    ped[i][c1up], ped[ped[i][c1p]][rp]);
    inbred += check(ped, pedsize, 
		    ped[i][c1up], ped[ped[i][c1p]][urp]);
    inbred += check(ped, pedsize, 
		    ped[i][c1up], ped[ped[j][c2p]][urp]);
    inbred += check(ped, pedsize, 
		    ped[j][c2up], ped[ped[i][c1p]][rp]);
    inbred += check(ped, pedsize, 
		    ped[j][c2up], ped[ped[i][c1p]][urp]);
    inbred += check(ped, pedsize, 
		    ped[j][c2up], ped[ped[j][c2p]][urp]);
    if (inbred == 0) 
      findit = 1;
  }
  /* related grandparent is the grandfather*/
  if(findit == 0) {
    rp = 4;
    urp = 5;
    if ( ped[ped[i][c1p]][rp] == ped[ped[j][c2p]][rp] && 
	 ped[ped[i][c1p]][rp] != 0){
      /* the other grandparents are unrelated, non-inbred, and the other
	 parents are not related with each other and not related with
	 the three grandparents */
      inbred = check(ped, pedsize, 
		     ped[ped[i][c1p]][urp], ped[ped[j][c2p]][urp]);
      inbred = check(ped, pedsize, 
		     ped[ped[i][c1p]][urp],ped[ped[i][c1p]][rp]);
      inbred = check(ped, pedsize, 
		     ped[ped[i][c1p]][rp], ped[ped[j][c2p]][urp]);
      inbred += inbred_check(ped, pedsize, ped[ped[i][c1p]][rp]);
      inbred += inbred_check(ped, pedsize, ped[ped[i][c1p]][urp]);
      inbred += inbred_check(ped, pedsize, ped[ped[j][c2p]][urp]);
      inbred += check(ped, pedsize, 
		      ped[i][c1up], ped[j][c2up]);
      inbred += check(ped, pedsize, 
		      ped[i][c1up], ped[ped[i][c1p]][rp]);
      inbred += check(ped, pedsize, 
		      ped[i][c1up], ped[ped[i][c1p]][urp]);
      inbred += check(ped, pedsize, 
		      ped[i][c1up], ped[ped[j][c2p]][urp]);
      inbred += check(ped, pedsize, 
		      ped[j][c2up], ped[ped[i][c1p]][rp]);
      inbred += check(ped, pedsize, 
		      ped[j][c2up], ped[ped[i][c1p]][urp]);
      inbred += check(ped, pedsize, 
		      ped[j][c2up], ped[ped[j][c2p]][urp]);
      if (inbred == 0) 
	findit = 1;
    }
  }
  
  return findit;
}

int check_halfpluscousin(int i, int j, int com, int ucom, int** ped,
			 int pedsize)
{
  int findit, inbred;
  
  findit = 0;
  if ( ped[i][com] == ped[j][com] &&
       ped[ped[i][ucom]][4] == ped[ped[j][ucom]][4] &&
       ped[ped[i][ucom]][4] != 0 &&
       ped[ped[i][ucom]][5] == ped[ped[j][ucom]][5] &&
       ped[ped[i][ucom]][5] != 0 ) {
    /* grandparents are unrelated, non-inbred, and unrelated with the
       common parent, and the common parent should not be inbred */ 
    inbred = check(ped, pedsize, 
		   ped[ped[i][ucom]][4], ped[ped[i][ucom]][5]);
    inbred += inbred_check(ped, pedsize, ped[ped[i][ucom]][4]);
    inbred += inbred_check(ped, pedsize, ped[ped[i][ucom]][5]);
    inbred += check(ped, pedsize, 
		  ped[i][com], ped[ped[i][ucom]][4]);
    inbred += check(ped, pedsize, 
		  ped[i][com], ped[ped[i][ucom]][5]);
    inbred += inbred_check(ped, pedsize, ped[i][com]);
    if (inbred == 0) 
      findit = 1;
  }

  return findit;
}


void build_anclist(int** ped, int pedsize, int id1, int
		   anclist[MANCE][3], int* totalnum, int flag)
{
  static int numinlist;

  if ( flag == 0 ) {
    numinlist = 1;
    anclist[numinlist][1] = id1;
    anclist[numinlist][2] = 0;
  }
  if ( ped[id1][4] == 0 && ped[id1][5] == 0 ) {
    *totalnum = numinlist;
    return;
  }

  ++numinlist;
  anclist[numinlist][1] = ped[id1][4];
  anclist[numinlist][2] = ped[id1][5];

  build_anclist(ped, pedsize, ped[id1][4], anclist, totalnum, 1);
  build_anclist(ped, pedsize, ped[id1][5], anclist, totalnum, 1);

  *totalnum = numinlist;

  return;

}

int check(int** ped, int pedsize, int id1, int id2)
{
  static int anclist[MANCE][3];
  int totalnum;

  if ( id1 == id2 && id1 != 0 )
    return 1;
  if ( check1_anclist(ped, pedsize, id1) ||
       check1_anclist(ped, pedsize, id2) )
    return 1;

  build_anclist(ped, pedsize, id1, anclist, &totalnum, 0);

  return check_anclist(ped, pedsize, id2, anclist, totalnum);

}

int check_anclist(int** ped, int pedsize, int id2, 
		  int anclist[MANCE][3], int totalnum)
{
  int i, c1, c2;

  if ( ped[id2][4] == 0 || ped[id2][5] == 0 ) {
    for (i = 1; i <= totalnum; ++i) {
      if ( anclist[i][1] == id2 || anclist[i][2] == id2 )
	return 1;
    }
    return 0;
  }
  
  for (i = 1; i <= totalnum; ++i) {
    if ( anclist[i][1] == ped[id2][4] || 
	 anclist[i][2] == ped[id2][5] )
      return 1;
  }

  c1 = check_anclist(ped, pedsize, ped[id2][4], anclist, totalnum);
  if ( c1 == 1 )
    return 1;
  c2 = check_anclist(ped, pedsize, ped[id2][5], anclist, totalnum);
  if ( c2 == 1 )
    return 1;

  return 0;
}

int check1_anclist(int** ped, int pedsize, int id1)
{

  int c1, c2;

  if ( ped[id1][4] == 0 && ped[id1][5] == 0 )
    return 0;
  if ( ped[id1][4] == 0 && ped[id1][5] != 0 )
    return 1;
  if ( ped[id1][5] == 0 && ped[id1][4] != 0 )
    return 1;

  c1 = check1_anclist(ped, pedsize, ped[id1][4]);
  if ( c1 == 1 )
    return 1;
  c2 = check1_anclist(ped, pedsize, ped[id1][5]);
  if ( c2 == 1 )
    return 1;

  return 0;

}

int inbred_check(int** ped, int pedsize, int id)
{

  if ( ped[id][4]==0 && ped[id][5]==0 )
    return 0;
  if ( ped[id][4]==0 && ped[id][5]!=0 )
    return 1;
  if ( ped[id][4]!=0 && ped[id][5]==0 )
    return 1;
 
  return check(ped,pedsize,ped[id][4],ped[id][5]);

}

int recoding(int** ped, int* nomem)
{
  int j, k, kkk, done, error;

  error = 0;
  for (j = 1; j <= nomem[2]; j++) {

    /* bug checking */
    done = 0;
    if ( ped[j][2] == ped[j][3] && ped[j][2] != 0 ) {
      fprintf(fp_error, "errors in pedigree %d:\n", nomem[1]);
      fprintf(fp_error, "%d'th parents have the same id\n\n",ped[j][1]);
      ++error;
    }
    for (k = j+1; k <= nomem[2]; k++) {
      if ( ped[j][1] == ped[k][1] ) {
	fprintf(fp_error, "errors in pedigree %d:\n", nomem[1]);
	fprintf(fp_error, "individual id = %d being used\n\n",ped[j][1]);
	++error;
      }
      if ( ped[j][2] == ped[k][3] && ped[j][2] != 0 ) {
	fprintf(fp_error, "errors in pedigree %d:\n", nomem[1]);
	fprintf(fp_error, "%d has two genders\n\n",ped[j][2]);
      	++error;
      }
      if ( ped[j][3] == ped[k][2] && ped[j][3] != 0 ) {
	fprintf(fp_error, "errors in pedigree %d:\n", nomem[1]);
	fprintf(fp_error, "%d has two genders\n\n",ped[j][3]);
      	++error;
      }
    }
    
    /* recoding for founders */
    if ( ped[j][2] == 0 ) { 
      ++done;
      ped[j][4] = 0;
    }
    if ( ped[j][3] == 0 ) {
      ++done;
      ped[j][5] = 0;
    }
    
    /* recoding for non-founders */
    kkk = 1;
    while ( done != 2 && kkk <= nomem[2] ) {
      if ( ped[j][2] == ped[kkk][1] ) {
	ped[j][4] = kkk;
	++done;
      }
      if ( ped[j][3] == ped[kkk][1] ) {
	ped[j][5] = kkk;
	++done;
      }
      ++kkk;
    }
    
    /* bug checking */ 
    if ( done != 2 ) {
      fprintf(fp_error, "pedigree %d is not complete:\n", nomem[1]);
      fprintf(fp_error, "%d's parents: %d and/or %d do not have their parent's ids\n\n", ped[j][1], ped[j][2], ped[j][3]);
      ++error;
    }
  }

  return error;
}  


/* for find genodata, and Mendelian error checking */

int findpedgeno(int** ped, int* nomem, FILE *fp_in2, 
		int nochrom, int* nomark, int** noalle, 
		int** geno, char** chromnameped) 
{
  int lineid;
  int j, m, k, numread, kkk;
  int findit, error, itemp[10];
  char buf[MAXSTR+1], ctemp1[MAXSTR+1], ctemp2[MAXSTR+1];
  FILE *fp_ped;
  
  error = 0;
  /* for each individual */
  for (j = 1; j <= nomem[2]; j++) {
    rewind(fp_in2);
    /* id of the individual */
    geno[j][0] = ped[j][1];
      
    kkk = 0;
    /* for each chrmosome */
    for (k = 1; k <= nochrom; k++) {
      findit = 0;
      fscanf(fp_in2,"%s %s", ctemp1, ctemp2);
      fp_ped = fopen(ctemp2, "r");
      
      lineid = 1;
      
      while( !feof(fp_ped) && !findit ) {
	
	numread = fscanf(fp_ped,"%d %d",&itemp[1], &itemp[2]);
	if ( numread != 2 ) {
	  printf("\nError reading genotype data for %d in ped %d at %s at line %d\n",
		 geno[j][0], nomem[1], chromnameped[k], lineid);
	  exit(0);
	}

	if ( itemp[2] == geno[j][0] && itemp[1] == nomem[1] ) {
	  findit = 1;
	  for (m = 1; m <= 4; m++)
	    if ( fscanf(fp_ped,"%d",&itemp[1]) != 1 ) {
	      printf("\nError reading genotype data for %d in ped %d at %s\n",
		     geno[j][0], nomem[1], chromnameped[k]);
	      exit(0);
	    }
	  
	  for (m = (kkk+1); m <= (kkk+nomark[k]); m++) {
	    if ( fscanf(fp_ped,"%d %d",&geno[j][2*m-1],&geno[j][2*m])
		 != 2) {
	      printf("\nError reading genotype data for %d in ped %d at %s\n",
		     geno[j][0], nomem[1], chromnameped[k]);
	      exit(0);
	    }
	    if ( (geno[j][2*m-1] > noalle[k][m-kkk] 
		 ||  geno[j][2*m-1] < 0) ||
		(geno[j][2*m] > noalle[k][m-kkk] 
		 ||  geno[j][2*m] < 0 )||
		(geno[j][2*m-1] == 0 && geno[j][2*m] != 0) ||
		(geno[j][2*m-1] != 0 && geno[j][2*m] == 0) ) {
	      fprintf(fp_error, "Genotyping errors:\n");
	      fprintf(fp_error, "pedigree %d, at %s's  %dth marker\n", 
		      nomem[1], chromnameped[k], (m-kkk));
	      fprintf(fp_error, "%d's genotypy is %d %d\n",
		      ped[j][1], geno[j][2*m-1], geno[j][2*m]);
	      fprintf(fp_error, "no. of alleles of this marker is %d\n\n",
		      noalle[k][m-kkk]);
	      ++error;
	    }
	  }
	  
	}
	
	fgets(buf, MAXSTR, fp_ped);
	if (strlen(buf) >= MAXSTR - 1) {
	  printf("\nWarning:  The length of the %d_th line in file %s\n", lineid, ctemp2);
	  printf("is too long (>5000).  Adjust MAXSTR and recompile.\n\n");
	  exit(0);
	}

	++lineid;
      }
                  
      if ( !findit ) {
	for (m = (kkk+1); m <= (kkk+nomark[k]); m++) {
	  geno[j][2*m-1] = 0;
	  geno[j][2*m] = 0;
	}
      }
      
      fclose(fp_ped);
      kkk += nomark[k];
    }
  }
  
  return error;
}

int checkMend(int** ped, int** geno, int* nomem,
	      int nochrom, int* nomark, char** chromnameped)
{
  int i, k, m, kkk;
  int error, findit;
  int s1, s2, f1, f2, m1, m2;
  
  error = 0;
  for (i = 1; i <= nomem[2]; i++) {

    if ( ped[i][4] != 0 && ped[i][5] != 0) {
      kkk = 0;
      for (k = 1; k <= nochrom; k++) {
	for (m = (kkk+1); m <= (kkk+nomark[k]); m++) {
	  findit = 0;
	
	  s1 = geno[i][2*m-1];
	  s2 = geno[i][2*m];
	  f1 = geno[ped[i][4]][2*m-1];
	  f2 = geno[ped[i][4]][2*m];
	  m1 = geno[ped[i][5]][2*m-1];
	  m2 = geno[ped[i][5]][2*m];
	  
	  if ( s1 != 0 && s2 != 0 ) {
	    if ( f1 != 0 && f2 != 0 && m1 != 0 && m2 != 0 ) {
	      if ( (s1 == f1 || s1 == f2) && (s2 == m1 || s2 == m2) )
		findit = 1;
	      if ( !findit ) {
		if ( (s2 == f1 || s2 == f2) && (s1 == m1 || s1 == m2) )
		  findit = 1;
		else {
		  fprintf(fp_error, "Mendelian errors:\n");
		  fprintf(fp_error, "pedigree %d at %s's  %dth marker\n", 
			  nomem[1], chromnameped[k], (m-kkk));
		  fprintf(fp_error, "%d's parents are %d %d, and genotype:", 
			  ped[i][1], ped[i][2], ped[i][3]);
		  fprintf(fp_error, " %d %d    %d %d  %d %d\n\n", 
			 s1, s2, f1, f2, m1, m2);
		  ++error;
		}
	      }
	    }
	    else if ( f1 != 0 && f2 != 0 && m1 == 0 && m2 == 0 ) {
	      if ( s1 == f1 || s1 == f2 || s2 == f1 || s2 == f2 )
		findit = 1;
	      if ( !findit ) {
		fprintf(fp_error, "Mendelian errors:\n");
		fprintf(fp_error, "pedigree %d at %s's  %dth marker\n", 
			nomem[1], chromnameped[k], (m-kkk));
		fprintf(fp_error, "%d's parents are %d %d, and genotype:", 
			ped[i][1], ped[i][2], ped[i][3]);
		fprintf(fp_error, " %d %d    %d %d  %d %d\n\n", 
			s1, s2, f1, f2, m1, m2);
		++error;
	      }
	    }
	    else if ( f1 == 0 && f2 == 0 && m1 != 0 && m2 != 0 ) {
	      if ( s1 == m1 || s1 == m2 || s2 == m1 || s2 == m2 )
		findit = 1;
	      if ( !findit ) {
		fprintf(fp_error, "Mendelian errors:\n");
		fprintf(fp_error, "pedigree %d at %s's  %dth marker\n", 
			nomem[1], chromnameped[k], (m-kkk));
		fprintf(fp_error, "%d's parents are %d %d, and genotype:\n", 
			ped[i][1], ped[i][2], ped[i][3]);
		fprintf(fp_error, " %d %d    %d %d  %d %d\n\n", 
			s1, s2, f1, f2, m1, m2);
		++error;
	      }
	    }
	  }
	  
	}
      	kkk += nomark[k];
      }
    }
  }

  return error;
}


/* for trans. prob., cond. prob. and other prob.*/

/* cp2( ) is the conditional probability of P(X|2 IBD ) 
   where X is the genotype data for the pair, 
   similar for cp1( ) and cp0( ), unordered genotype for a pair */
double cp2(int f1, int m1, int f2, int m2, double* q)
{      
  double temp;
 
  if ( f1 == m1 && m1 == f2 && f2 == m2  ) 
    temp = q[f1]*q[f1];
  else if ( (f1 == f2 && m1 == m2) || (f1 == m2 && m1 == f2) ) 
    temp = 2*q[f1]*q[m1];
  else 
    temp = 0;

  return temp;
}

double cp1(int f1, int m1, int f2, int m2, double* q)
{
  double temp;

  if ( f1 == m1 && m1 == f2 && f2 == m2 ) 
    temp = q[f1]*q[f1]*q[f1];
  else if ( (f1 == f2 && m1 == m2) || (f1 == m2 && m1 == f2) ) 
    temp = q[f1]*q[m1]*(q[f1] + q[m1]);
  else if ( (f1 == m1 && m1 == f2) || (f1 == m1 && m1 == m2) ) 
    temp = 2*q[f1]*q[f2]*q[m2];
  else if ( (f2 == m2 && m2 == f1) || (f2 == m2 && m2 == m1) ) 
    temp = 2*q[f1]*q[m1]*q[f2];
  else if ( f1 == f2 || f1 == m2 ) 
    temp = 2*q[m1]*q[f2]*q[m2];
  else if ( m1 == f2 || m1 == m2 ) 
    temp = 2*q[f1]*q[f2]*q[m2];
  else 
    temp = 0;

  return temp;
}

double cp0(int f1, int m1, int f2, int m2, double* q)
{
  double temp;

  temp = q[f1]*q[m1]*q[f2]*q[m2];
  if ( f1 == m1 && m1==f2 && f2 == m2 ) 
    temp = temp;
  else if ( f1 == m1 &&  f2 == m2 ) 
    temp = 2*temp;
  else if ( f1 == m1 || f2 == m2 || (f1 == f2 && m1 == m2) 
	  || (f1 == m2 && m1 == f2) )
    temp = 4*temp;
  else 
    temp = 8*temp;

  return temp;
}

void transprob(double dij, double pibdiibdj[NRMLRT+1][3+1][3+1], 
	       int reltype, int maptype)
{
  int i, j, k;
  double ta, phi, psi;
  
  ta = mapcmtotheta(dij, maptype);  
  phi = ta*ta+(1-ta)*(1-ta);
  psi = 2*ta*(1-ta);

  for (i = 1; i <= NRMLRT; i++) 
    for (j = 1; j <= 3; j++) 
      for (k = 1; k <= 3; k++) 
	pibdiibdj[i][j][k] = 0;
      
  /* in 0-1-2 IBD to 0-1-2 IBD order */

  /* fullsib */
  if (reltype == 1) {
    pibdiibdj[1][1][1] = pow((ta*ta+(1-ta)*(1-ta)),2);
    pibdiibdj[1][1][2] = 2*(ta*ta+(1-ta)*(1-ta))*2*ta*(1-ta);
    pibdiibdj[1][1][3] = pow((2*ta*(1-ta)),2);
    pibdiibdj[1][2][1] = (ta*ta+(1-ta)*(1-ta))*2*ta*(1-ta);
    pibdiibdj[1][2][2] = pow((ta*ta+(1-ta)*(1-ta)),2) + pow((2*ta*(1-ta)),2);
    pibdiibdj[1][2][3] = (ta*ta+(1-ta)*(1-ta))*2*ta*(1-ta);
    pibdiibdj[1][3][1] = pow((2*ta*(1-ta)),2);
    pibdiibdj[1][3][2] = 2*(ta*ta+(1-ta)*(1-ta))*2*ta*(1-ta);
    pibdiibdj[1][3][3] = pow((ta*ta+(1-ta)*(1-ta)),2);
  }
  /* half-sib */
  else if (reltype == 2) {
    pibdiibdj[2][1][1] = ta*ta+(1-ta)*(1-ta);
    pibdiibdj[2][1][2] = 2*ta*(1-ta);
    pibdiibdj[2][2][1] = 2*ta*(1-ta);
    pibdiibdj[2][2][2] = ta*ta+(1-ta)*(1-ta);
  }
  /* grandparent-grandchild */
  else if (reltype == 3) {
    pibdiibdj[3][1][1] = (1-ta);
    pibdiibdj[3][1][2] = ta;
    pibdiibdj[3][2][1] = ta;
    pibdiibdj[3][2][2] = (1-ta);
  }
  /* avuncular */
  else if (reltype ==4) {
    pibdiibdj[4][2][2] = 1-2.5*ta+4*ta*ta-2*ta*ta*ta;
    pibdiibdj[4][2][1] = 1-pibdiibdj[4][2][2];
    pibdiibdj[4][1][2] = 1-pibdiibdj[4][2][2];
    pibdiibdj[4][1][1] = pibdiibdj[4][2][2];
  }
  /* first-cousin */
  else if (reltype == 5) {
    pibdiibdj[5][2][2] = 1-4*ta+7.5*ta*ta- 6*ta*ta*ta+2*ta*ta*ta*ta;
    pibdiibdj[5][2][1] = 1-pibdiibdj[5][2][2];
    pibdiibdj[5][1][2] = (1-pibdiibdj[5][2][2])/3.0;
    pibdiibdj[5][1][1] = 1-pibdiibdj[5][1][2];
  }
  /* unrelated */
  else if (reltype == 6) 
    pibdiibdj[6][1][1] = 1;
  
  /* half-avuncular */
  else if (reltype == 7) {
    pibdiibdj[7][2][2]=(1-ta)*(1-ta)*
      (1-ta)+ta*ta*(1-ta);
    pibdiibdj[7][2][1]=1-pibdiibdj[7][2][2];
    pibdiibdj[7][1][2]=(1-pibdiibdj[7][2][2])/3.0;
    pibdiibdj[7][1][1]=1-pibdiibdj[7][1][2];
  }
  /* half-cousin */
  else if (reltype == 8) {
    pibdiibdj[8][2][2]=(1-ta)*(1-ta)*(1-ta)*(1-ta)+ta*ta*(1-ta)*(1-ta);
    pibdiibdj[8][2][1]=1-pibdiibdj[8][2][2];
    pibdiibdj[8][1][2]=(1-pibdiibdj[8][2][2])/7.0;
    pibdiibdj[8][1][1]=1-pibdiibdj[8][1][2];
  }
  /* half-sib plus first-cousin */
  else if (reltype == 9) {
    pibdiibdj[9][3][3] = pibdiibdj[5][2][2]*phi;
    pibdiibdj[9][3][1] = pibdiibdj[5][2][1]*psi;
    pibdiibdj[9][3][2] = 1-pibdiibdj[9][3][3]-pibdiibdj[9][3][1];
    
    pibdiibdj[9][1][3] = pibdiibdj[5][1][2]*psi;
    pibdiibdj[9][1][1] = pibdiibdj[5][1][1]*phi;
    pibdiibdj[9][1][2] = 1-pibdiibdj[9][1][3]-pibdiibdj[9][1][1];
    
    pibdiibdj[9][2][3] = 1.0/3*pibdiibdj[5][2][2]*psi+
      2.0/3*pibdiibdj[5][1][2]*phi;
    pibdiibdj[9][2][1] = 1.0/3*pibdiibdj[5][2][1]*phi+
      2.0/3*pibdiibdj[5][1][1]*psi;
    pibdiibdj[9][2][2] = 1-pibdiibdj[9][2][3]-pibdiibdj[9][2][1];
  }
  /* parent-offspring */
  else if (reltype == 10) 
    pibdiibdj[10][2][2] = 1;
  /* MZ twins */
  else if (reltype == 11) 
    pibdiibdj[11][3][3] = 1;

}

void transproball(double dij, double pibdiibdj[NRMLRT+1][3+1][3+1], 
		  int maptype)
{
  int i, j, k;
  double ta, phi, psi;

  ta = mapcmtotheta(dij, maptype);  
  phi = ta*ta+(1-ta)*(1-ta);
  psi = 2*ta*(1-ta);

  for (i = 1; i <= NRMLRT; i++) 
    for (j = 1; j <= 3; j++) 
      for (k = 1; k <= 3; k++) 
	pibdiibdj[i][j][k] = 0;
      
  /* in 0-1-2 IBD to 0-1-2 IBD order */

  /* fullsib */
  pibdiibdj[1][1][1] = pow((ta*ta+(1-ta)*(1-ta)),2);
  pibdiibdj[1][1][2] = 2*(ta*ta+(1-ta)*(1-ta))*2*ta*(1-ta);
  pibdiibdj[1][1][3] = pow((2*ta*(1-ta)),2);
  pibdiibdj[1][2][1] = (ta*ta+(1-ta)*(1-ta))*2*ta*(1-ta);
  pibdiibdj[1][2][2] = pow((ta*ta+(1-ta)*(1-ta)),2) + pow((2*ta*(1-ta)),2);
  pibdiibdj[1][2][3] = (ta*ta+(1-ta)*(1-ta))*2*ta*(1-ta);
  pibdiibdj[1][3][1] = pow((2*ta*(1-ta)),2);
  pibdiibdj[1][3][2] = 2*(ta*ta+(1-ta)*(1-ta))*2*ta*(1-ta);
  pibdiibdj[1][3][3] = pow((ta*ta+(1-ta)*(1-ta)),2);

  /* half-sib */
  pibdiibdj[2][1][1] = ta*ta+(1-ta)*(1-ta);
  pibdiibdj[2][1][2] = 2*ta*(1-ta);
  pibdiibdj[2][2][1] = 2*ta*(1-ta);
  pibdiibdj[2][2][2] = ta*ta+(1-ta)*(1-ta);
  
  /* grandparent-grandchild */
  pibdiibdj[3][1][1] = (1-ta);
  pibdiibdj[3][1][2] = ta;
  pibdiibdj[3][2][1] = ta;
  pibdiibdj[3][2][2] = (1-ta);

  /* avuncular */
  pibdiibdj[4][2][2] = 1-2.5*ta+4*ta*ta-2*ta*ta*ta;
  pibdiibdj[4][2][1] = 1-pibdiibdj[4][2][2];
  pibdiibdj[4][1][2] = 1-pibdiibdj[4][2][2];
  pibdiibdj[4][1][1] = pibdiibdj[4][2][2];

  /* first-cousin */
  pibdiibdj[5][2][2] = 1-4*ta+7.5*ta*ta- 6*ta*ta*ta+2*ta*ta*ta*ta;
  pibdiibdj[5][2][1] = 1-pibdiibdj[5][2][2];
  pibdiibdj[5][1][2] = (1-pibdiibdj[5][2][2])/3.0;
  pibdiibdj[5][1][1] = 1-pibdiibdj[5][1][2];
  
  /* unrelated */
  pibdiibdj[6][1][1] = 1;

  /* half-avuncular */
  pibdiibdj[7][2][2]=(1-ta)*(1-ta)*
    (1-ta)+ta*ta*(1-ta);
  pibdiibdj[7][2][1]=1-pibdiibdj[7][2][2];
  pibdiibdj[7][1][2]=(1-pibdiibdj[7][2][2])/3.0;
  pibdiibdj[7][1][1]=1-pibdiibdj[7][1][2];
  
  /* half-cousin */
  pibdiibdj[8][2][2]=(1-ta)*(1-ta)*(1-ta)*(1-ta)+ta*ta*(1-ta)*(1-ta);
  pibdiibdj[8][2][1]=1-pibdiibdj[8][2][2];
  pibdiibdj[8][1][2]=(1-pibdiibdj[8][2][2])/7.0;
  pibdiibdj[8][1][1]=1-pibdiibdj[8][1][2];

  /* half-sib plus first-cousin */
  pibdiibdj[9][3][3] = pibdiibdj[5][2][2]*phi;
  pibdiibdj[9][3][1] = pibdiibdj[5][2][1]*psi;
  pibdiibdj[9][3][2] = 1-pibdiibdj[9][3][3]-pibdiibdj[9][3][1];
  
  pibdiibdj[9][1][3] = pibdiibdj[5][1][2]*psi;
  pibdiibdj[9][1][1] = pibdiibdj[5][1][1]*phi;
  pibdiibdj[9][1][2] = 1-pibdiibdj[9][1][3]-pibdiibdj[9][1][1];
    
  pibdiibdj[9][2][3] = 1.0/3*pibdiibdj[5][2][2]*psi+
    2.0/3*pibdiibdj[5][1][2]*phi;
  pibdiibdj[9][2][1] = 1.0/3*pibdiibdj[5][2][1]*phi+
    2.0/3*pibdiibdj[5][1][1]*psi;
  pibdiibdj[9][2][2] = 1-pibdiibdj[9][2][3]-pibdiibdj[9][2][1];

  /* parent-offspring */
  pibdiibdj[10][2][2] = 1;

  /* MZ twins */
  pibdiibdj[11][3][3] = 1;

}

/* return P(Z<=z) for Z a N(0,1) RV */
double pnorm(double z)
{
  double x = 0.0;

  if (z != 0.0) {
    double y = 0.5 * fabs (z);
    if (y >= 3.0)
      x = 1.0;
    else if (y < 1.0) {
      double w = y*y;
      x = ((((((((0.000124818987 * w
        -0.001075204047) * w +0.005198775019) * w
        -0.019198292004) * w +0.059054035642) * w
        -0.151968751364) * w +0.319152932694) * w
        -0.531923007300) * w +0.797884560593) * y * 2.0;
    }
    else {
      y -= 2.0;
      x = (((((((((((((-0.000045255659 * y
        +0.000152529290) * y -0.000019538132) * y
        -0.000676904986) * y +0.001390604284) * y
        -0.000794620820) * y -0.002034254874) * y
        +0.006549791214) * y -0.010557625006) * y
        +0.011630447319) * y -0.009279453341) * y
        +0.005353579108) * y -0.002141268741) * y
        +0.000535310849) * y +0.999936657524;
    }
  }

  return (z > 0.0 ? ((x + 1.0) * 0.5) : ((1.0 - x) * 0.5));

}


/* map functions */

double mapthetatocm(double theta, int maptype)
{
  double cm;
  
  /* no-interference model */
  if ( maptype == 1 ) {
    cm = -0.5*log(1-2*theta)*100;
  }
  /* kosambi model */
  else if ( maptype == 2 ) {
    cm = 0.25*log((1+2*theta)/(1-2*theta))*100;
  }
  
  return cm;
}

double mapcmtotheta(double cm, int maptype)
{
  double theta;
  
  /* no-interference model */
  if ( maptype == 1 ) {
    theta = 0.5*(1-exp(-.02*cm));  
  }
  /* kosambi model */
  else if ( maptype == 2 ) {
    theta = 0.5*(exp(.04*cm)-1)/(exp(.04*cm)+1); 
  }
  
  return theta;
}


/* EM algorithm for estimation of p0, p1, p2 (0-1-2 IBD sharing) */

/* need to input three coefficents, initial p0, p1, and p2 values, and
   the total marker numbers */
void EMforp012(int* geno1, int* geno2, int* valid, 
	       double*** q, int nochrom, int* nomark,
	       int reltype, double phat[3+1])
{
  int i, j, kkk, commark, niter;
  int f1, m1, f2, m2;
  double *a, *b, *c;
  double p0, p1, p2;
  double p0_next, p1_next, p2_next;
  
  /* initialized the value for phat */
  p0 = p1 = p2 = 0;
  p0_next = .2;
  p1_next = .6;
  p2_next = .2;
  
  a = (double *)malloc((size_t) ((nomark[0]+1)*sizeof(double)));
  b = (double *)malloc((size_t) ((nomark[0]+1)*sizeof(double)));
  c = (double *)malloc((size_t) ((nomark[0]+1)*sizeof(double)));
  
  niter = kkk = commark = 0;
  
  for (i = 1; i <= nochrom; i++) {
    for (j = 1; j <= nomark[i]; j++) {
      if ( valid[kkk+j] != 0 ) {
	++commark;
	
	f1 = geno1[2*(kkk+j)-1];
	m1 = geno1[2*(kkk+j)];
	f2 = geno2[2*(kkk+j)-1];
	m2 = geno2[2*(kkk+j)];
	
	a[commark] = cp0(f1, m1, f2, m2, q[i][j]);
	b[commark] = cp1(f1, m1, f2, m2, q[i][j]);
	c[commark] = cp2(f1, m1, f2, m2, q[i][j]);
	
	/* for parent-offspring, use modified a b c if all of them
	   are 0 in the presence of genotyping error, so that 
	   the denomiator a[i]*p0+b[i]*p1+c[i]*p2 will not be 0 in the EM */
	if (reltype == 10 && (b[commark] < MYEPS)) {
	  b[commark] = (1-errorrate)*(1-errorrate)*b[commark]
	    +(1-(1-errorrate)*(1-errorrate))*a[commark];
	  c[commark] = (1-errorrate)*(1-errorrate)*c[commark]
	    +(1-(1-errorrate)*(1-errorrate))*a[commark];
	}
	
      }
    }
    kkk += nomark[i];
  }
  
  while( (fabs(p0_next-p0) > MYEPS || 
	  fabs(p1_next-p1) > MYEPS ) &&
	 niter <= MITER ) {
    p0 = p0_next; 
    p1 = p1_next; 
    p2 = p2_next;
    p0_next = p1_next = p2_next = 0; 
    
    for (i = 1; i <= commark; i++) {
      
      p0_next += a[i]*p0/(a[i]*p0+b[i]*p1+c[i]*p2);
      p1_next += b[i]*p1/(a[i]*p0+b[i]*p1+c[i]*p2);
      
    }
    
    p0_next = p0_next/commark;
    p1_next = p1_next/commark;
    p2_next = 1.0 - p0_next - p1_next;
    ++niter;
  } 
  
  if (niter <= MITER) {
    phat[1] = p0_next;
    phat[2] = p1_next;
    phat[3] = p2_next;
    /* to avoid print -0 */
    if(fabs(phat[3]) < MYEPS)
      phat[3] = 0.0;
  }
  else {
    phat[1] = -1;
    phat[2] = -1;
    phat[3] = -1;
  }

  free(a);
  free(b);
  free(c);
  
}
/* EIBD, AIBS, IBS statistic, mu and std calculation */

void get3statandmustd(int* geno1, int* geno2, int* valid, int commark,
		      int reltype, int nochrom, int* nomark, 
		      double p012[4],
		      double*** q, double** cenm, int maptype,
		      double*** eeibdibd, double*** eeibdeibdibd,
		      double*** eaibsibd, double*** eaibsaibsibd,
		      double** eibsibd, double** eibsibsibd,
		      double statresult[3+1][5+1])
{
  int i, j, k, l, n, kkk;
  int f1, f2, m1, m2;
  double dij, p2, p1, p0;
  double eibd, veibd, aibs, vaibs, ibs, vibs;
  double dtemp[3+1], dtemp1[3+1], pibdiibdj[NRMLRT+1][3+1][3+1];
  double **exxchr, **exchr;
  
  exxchr = (double **)malloc((size_t) ((3+1)*sizeof(double*)));
  exchr = (double **)malloc((size_t) ((3+1)*sizeof(double*)));
  for (i = 1; i <= 3; i++) {
    exxchr[i] = (double *)malloc((size_t) ((nochrom+1)*sizeof(double)));
    exchr[i] = (double *)malloc((size_t) ((nochrom+1)*sizeof(double)));
  }
  
  eibd = aibs = ibs = 0;
  veibd = vaibs = vibs = 0;
  for (i = 1; i <= 3; i++) {
    dtemp1[i] = 0;
    for (n = 1; n <= nochrom; n++) 
      exxchr[i][n] = exchr[i][n] = 0;
  }
  
  kkk = 0;
  for (n = 1; n <= nochrom; n++) {
    for (i = 1; i <= nomark[n]; i++) {
      /* statistic calculation */
      if ( valid[kkk+i] != 0 ) {  
	f1 = geno1[2*(kkk+i)-1];
	m1 = geno1[2*(kkk+i)];
	f2 = geno2[2*(kkk+i)-1];
	m2 = geno2[2*(kkk+i)];
 	/* eibd */
	p2 = cp2(f1,m1,f2,m2,q[n][i]);
	p1 = cp1(f1,m1,f2,m2,q[n][i]);
	p0 = cp0(f1,m1,f2,m2,q[n][i]); 

	/* if it is not parent-offspring nor unrelated, nor MZ twins */
	if( !(reltype == 10 || reltype == 6 || reltype == 11) )
	  eibd += (2*p012[3]*p2+p012[2]*p1)/
	    (p012[3]*p2+p012[2]*p1+p012[1]*p0);
	/* aibs and ibs */
        if ( (f1 == f2 && m1 == m2) || (f1 == m2 && m1 == f2) ) {
          aibs += p012[0]/(p012[0] + (1-p012[0])*q[n][i][f1]) + 
	    p012[0]/(p012[0] + (1-p012[0])*q[n][i][m1]);
          ibs += 2;
        }
        else if ( f1 == f2 || f1 == m2 ) {
          aibs += p012[0]/(p012[0] + (1-p012[0])*q[n][i][f1]);
          ibs += 1;
        }
        else if ( m1 == f2 || m1 == m2 ) {
          aibs += p012[0]/(p012[0] + (1-p012[0])*q[n][i][m1]);
          ibs += 1;
        }
      }
      
      /* variance and expection calculation */
      for (j = i; j <= nomark[n]; j++) {

	if (i == j && valid[kkk+i] == 1) {
	  for (k = 1; k <= 3; k++)
	    dtemp[k] = 0;
	  for (k = 1; k <= 3; k++) {
	    dtemp[1] += eeibdeibdibd[reltype][kkk+i][k]*p012[k];
            dtemp[2] += eaibsaibsibd[reltype][kkk+i][k]*p012[k];
	    dtemp[3] += eibsibsibd[kkk+i][k]*p012[k];
 	  }
	  for (k = 1; k <= 3; k++)
	    exxchr[k][n] += dtemp[k];
	}

	else if ( j > i && valid[kkk+i] == 1 && valid[kkk+j] == 1) {
	  dij = 0;
	  for (k = (i+1); k <= j; k++) 
	    dij += cenm[n][k];
	  transprob(dij, pibdiibdj, reltype, maptype);
	  for (k = 1; k <= 3; k++)
	    dtemp[k] = 0;
	  for (k = 1; k <= 3; k++) 
	    for (l = 1; l <=3; l++) {
	      dtemp[1] += eeibdibd[reltype][kkk+i][k]*
		eeibdibd[reltype][kkk+j][l]*p012[k]*
		pibdiibdj[reltype][k][l];
	      dtemp[2] += eaibsibd[reltype][kkk+i][k]*
		eaibsibd[reltype][kkk+j][l]*p012[k]*
		pibdiibdj[reltype][k][l];
	      dtemp[3] += eibsibd[kkk+i][k]*
		eibsibd[kkk+j][l]*p012[k]*
		pibdiibdj[reltype][k][l];
	    }
	  for (l = 1; l <=3; l++) 
	    exxchr[l][n] += 2*dtemp[l];
	}
      }
            
      if ( valid[kkk+i] == 1 ) {
	exchr[1][n] += 4*p012[0];
	dtemp[2] = 0;
	for (k = 1; k <= 3; k++)
          dtemp[2] += eaibsibd[reltype][kkk+i][k]*p012[k];
        exchr[2][n] += dtemp[2];
	dtemp[3] = 0;
        for (k = 1; k <= 3; k++)
          dtemp[3] += eibsibd[kkk+i][k]*p012[k];
        exchr[3][n] += dtemp[3];
      }

    }
    
    for (k =  1; k <= 3; k++)    
      dtemp1[k] += exchr[k][n];

    if (!(reltype == 10 || reltype == 6 || reltype == 11) )
      veibd += exxchr[1][n] - exchr[1][n]*exchr[1][n];
    if (!(reltype == 6) )
      vaibs += exxchr[2][n] - exchr[2][n]*exchr[2][n];
    vibs += exxchr[3][n] - exchr[3][n]*exchr[3][n];
    kkk += nomark[n];
    
  }

  statresult[1][1] = eibd/commark;
  if(reltype == 6)
    statresult[1][1] = 0;
  else if(reltype == 10)
    statresult[1][1] = 1;
  else if(reltype == 11)
    statresult[1][1] = 2;

  /* E(EIBD) = kinship-coef*4 */
  statresult[1][2] =  p012[0]*4;
  /* delete later */
  if (fabs(dtemp1[1]/commark-statresult[1][2])>0.001) {
    printf("program is wrong for EIBD\n");
    exit(0);
  }

  statresult[1][3] = sqrt(veibd/(commark*commark));
  if(reltype == 6 || reltype == 10 || reltype == 11) 
    statresult[1][3] = 0;
  
  statresult[1][4] = (statresult[1][1] - statresult[1][2])/
    statresult[1][3];
  statresult[1][5] = 2*pnorm(-fabs(statresult[1][4]));
    
  statresult[2][1] = aibs/commark;
  if(reltype == 6)
    statresult[2][1] = 0;
  
  statresult[2][2] = dtemp1[2]/commark;
  statresult[2][3] = sqrt(vaibs/(commark*commark));
  if(reltype == 6)
    statresult[2][3] = 0;

  statresult[2][4] = (statresult[2][1] - statresult[2][2])/
    statresult[2][3];
  statresult[2][5] = 2*pnorm(-fabs(statresult[2][4]));

  statresult[3][1] = ibs/commark;
  statresult[3][2] = dtemp1[3]/commark;
  statresult[3][3] = sqrt(vibs/(commark*commark));
  statresult[3][4] = (statresult[3][1] - statresult[3][2])/
    statresult[3][3];
  statresult[3][5] = 2*pnorm(-fabs(statresult[3][4]));

  for (i = 1; i <= 3; i++) {
    free(exchr[i]);
    free(exxchr[i]);
  }
  free(exchr);
  free(exxchr);
  
}



void forEIBDvar(double*** eeibdibd, double*** eeibdeibdibd, 
		int nochrom, int* nomark, int* ngeno, int*** genotype, 
		double*** q, double p012[NRMLRT+1][3+1])
{
  int i, j, k, l, m, kkk;
  int f1, f2, m1, m2;
  double p2, p1, p0;
  double eibd;
  
  kkk=0;
  for (i = 1; i <= nochrom; i++) {
    for (j = 1; j <= nomark[i]; j++) {
      
      for (k = 1; k <= NR; k++) {
	for (l = 1; l <= 3; l++)
	  eeibdibd[k][kkk+j][l] = eeibdeibdibd[k][kkk+j][l] = 0;
      }
      
      for (k = 1; k <= ngeno[kkk+j]; k++) {
	for (l = k; l <= ngeno[kkk+j]; l++) {
	  f1 = genotype[kkk+j][k][1];
	  m1 = genotype[kkk+j][k][2];
	  f2 = genotype[kkk+j][l][1];
	  m2 = genotype[kkk+j][l][2];
		  
	  p2 = cp2(f1,m1,f2,m2,q[i][j]);
	  p1 = cp1(f1,m1,f2,m2,q[i][j]);
	  p0 = cp0(f1,m1,f2,m2,q[i][j]);
	  
	  for (m = 1; m <= NR; m++) {
	    eibd = (2*p012[m][3]*p2+p012[m][2]*p1)/
	      (p012[m][3]*p2+p012[m][2]*p1+p012[m][1]*p0);
	    
	    eeibdibd[m][kkk+j][3] += eibd*p2;
	    eeibdeibdibd[m][kkk+j][3] += eibd*eibd*p2;
	    
	    eeibdibd[m][kkk+j][2] += eibd*p1;
	    eeibdeibdibd[m][kkk+j][2] += eibd*eibd*p1;
	    
	    eeibdibd[m][kkk+j][1] += eibd*p0;
	    eeibdeibdibd[m][kkk+j][1] += eibd*eibd*p0;
	  }
	}
      }
    }
    kkk += nomark[i];
  }
  
}

void forAIBSvar(double*** eaibsibd, double*** eaibsaibsibd, 
		int nochrom, int* nomark, int* ngeno, int*** genotype, 
		double*** q, double p012[NRMLRT+1][3+1])
{
  int i, j, k, l, m, kkk;
  int f1, f2, m1, m2;
  double p2, p1, p0;
  double aibs;

  kkk=0;

  for (i = 1; i <= nochrom; i++) {
    for (j = 1; j <= nomark[i]; j++) {

      for (k = 1; k <= NR; k++) {
	for (l = 1; l <= 3; l++)
	eaibsibd[k][kkk+j][l] = eaibsaibsibd[k][kkk+j][l] = 0;
      }
      
      for (k = 1; k <= ngeno[kkk+j]; k++) {
	for (l = k; l <= ngeno[kkk+j]; l++) {
	  
	  f1 = genotype[kkk+j][k][1];
	  m1 = genotype[kkk+j][k][2];
	  f2 = genotype[kkk+j][l][1];
	  m2 = genotype[kkk+j][l][2];
		  
	  p2 = cp2(f1,m1,f2,m2,q[i][j]);
	  p1 = cp1(f1,m1,f2,m2,q[i][j]);
	  p0 = cp0(f1,m1,f2,m2,q[i][j]);
	  
	  for (m = 1; m <= NR; m++) {
	    if ( (f1 == f2 && m1 == m2) || (f1 == m2 && m1 == f2) ) 
	      aibs =
		p012[m][0]/(p012[m][0]+(1-p012[m][0])*q[i][j][f1]) + 
		p012[m][0]/(p012[m][0]+(1-p012[m][0])*q[i][j][m1]);
	    else if ( f1 == f2 || f1 == m2 ) 
	      aibs = p012[m][0]/(p012[m][0]+(1-p012[m][0])*q[i][j][f1]);
	    else if ( m1 == f2 || m1 == m2 ) 
	      aibs = p012[m][0]/(p012[m][0]+(1-p012[m][0])*q[i][j][m1]);
	    else 
	      aibs = 0;
	  
	    eaibsibd[m][kkk+j][3] += aibs*p2;
	    eaibsaibsibd[m][kkk+j][3] += aibs*aibs*p2;

            eaibsibd[m][kkk+j][2] += aibs*p1;
            eaibsaibsibd[m][kkk+j][2] += aibs*aibs*p1;

            eaibsibd[m][kkk+j][1] += aibs*p0;
            eaibsaibsibd[m][kkk+j][1] += aibs*aibs*p0;
	  }
	}
      }
    }
    kkk += nomark[i];
  }
  
}

void forIBSvar(double** eibsibd, double** eibsibsibd, 
	       int nochrom, int* nomark, double*** q, int** noalle)
{
  int i, j, k, l, m, kkk;
  double dtemp, dtemp1, dtemp2, dtemp3, dtemp4;

  kkk=0;
  for (i = 1; i <= nochrom; i++) {
    for (j = 1; j <= nomark[i]; j++) {
      
      eibsibd[kkk+j][3] = 2;
      eibsibsibd[kkk+j][3] = 4;

      dtemp = 0;
      for (k = 1; k <= noalle[i][j]; k++)
        dtemp += q[i][j][k]*q[i][j][k];
      eibsibd[kkk+j][2] = dtemp*2 + (1-dtemp);
      eibsibsibd[kkk+j][2] = dtemp*4 + (1-dtemp);

      dtemp1 = dtemp2 = dtemp3 = dtemp4 = 0;
      for (k = 1; k <= noalle[i][j]; k++) {
        dtemp1 += q[i][j][k]*q[i][j][k]*q[i][j][k]*q[i][j][k];
        for (l = 1; l <= noalle[i][j];l++) {
          if (l != k) {
            dtemp2 += 2*q[i][j][k]*q[i][j][k]*q[i][j][l]*q[i][j][l];
            dtemp3 += 4*q[i][j][k]*q[i][j][k]*q[i][j][k]*q[i][j][l];
            for (m = 1; m <= noalle[i][j]; m++) {
              if (m != k && m != l)
                dtemp4 +=
		  4*q[i][j][k]*q[i][j][k]*q[i][j][l]*q[i][j][m];
            }
          }
        }
      }
      eibsibd[kkk+j][1] = 2*(dtemp1+dtemp2)+(dtemp3+dtemp4);
      eibsibsibd[kkk+j][1] = 4*(dtemp1+dtemp2)+(dtemp3+dtemp4);
    }
    kkk += nomark[i];
  }
  
}

void get3stat(int* geno1, int* geno2, int* valid, int commark,
	      int reltype, int nochrom, int* nomark, 
	      double p012[4], double*** q,
	      double statresult[3+1])
{
  int i, n, kkk;
  int f1, f2, m1, m2;
  double p0, p1, p2;
  double eibd, aibs, ibs;
  

  kkk = 0;
  eibd = aibs = ibs = 0;
  for (n = 1; n <= nochrom; n++) {
    for (i = 1; i <= nomark[n]; i++) {
      if ( valid[kkk+i] != 0 ) {  
	f1 = geno1[2*(kkk+i)-1];
	m1 = geno1[2*(kkk+i)];
	f2 = geno2[2*(kkk+i)-1];
	m2 = geno2[2*(kkk+i)];
 	/* eibd */
	p2 = cp2(f1, m1, f2, m2, q[n][i]);
	p1 = cp1(f1, m1, f2, m2, q[n][i]);
	p0 = cp0(f1, m1, f2, m2, q[n][i]); 

        /* if it is not parent-offspring nor unrelated, nor MZ twins */
        if( !(reltype == 10 || reltype == 6 || reltype == 11) )
          eibd += (2*p012[3]*p2+p012[2]*p1)/
            (p012[3]*p2+p012[2]*p1+p012[1]*p0);

	/* aibs and ibs */
        if ( (f1 == f2 && m1 == m2) || (f1 == m2 && m1 == f2) ) {
          aibs += p012[0]/(p012[0] + (1-p012[0])*q[n][i][f1]) + 
	    p012[0]/(p012[0] + (1-p012[0])*q[n][i][m1]);
          ibs += 2;
        }
        else if ( f1 == f2 || f1 == m2 ) {
          aibs += p012[0]/(p012[0] + (1-p012[0])*q[n][i][f1]);
          ibs += 1;
        }
        else if ( m1 == f2 || m1 == m2 ) {
          aibs += p012[0]/(p012[0] + (1-p012[0])*q[n][i][m1]);
          ibs += 1;
        }
      }
    }
    kkk += nomark[n];
  }
  
  statresult[1] = eibd/commark;
  if(reltype == 6)
    statresult[1] = 0;
  else if(reltype == 10)
    statresult[1] = 1;
  else if(reltype == 11)
    statresult[1] = 2;
  statresult[2] = aibs/commark;
  if(reltype == 6)
    statresult[2] = 0;
  statresult[3] = ibs/commark;
  
}

/* for likelihood calculation */
void getapprloglikelihood(int* geno1, int* geno2, int* valid,
			  double*** q, double** cenm, 
			  int nochrom, int* nomark, int maptype,
			  double loglhresult[NRMLRT+1])
{
  int i, k, n, kkk;
  int validmin, validmax;
  int f1, m1, f2, m2;
  double s0, s1, s2, h0, h1, g0, g1;
  double au0, au1, ac0, ac1, unr0;
  double s0_next, s1_next, s2_next, h0_next, h1_next, g0_next, g1_next;
  double au0_next, au1_next, ac0_next, ac1_next, unr0_next;
  double po1, po1_next, mz2, mz2_next;
  
  double ahu0, ahu1, ahc0, ahc1, ahplusc0, ahplusc1, ahplusc2;
  double ahu0_next, ahu1_next, ahc0_next, ahc1_next;
  double ahplusc0_next, ahplusc1_next, ahplusc2_next;

  double p0, p1, p2;
  double dij, pij[NRMLRT+1][3+1][3+1];

  
  kkk = 0;
  for (i = 1; i <= NRMLRT; i++)
    loglhresult[i] = 0;
  
  for (n = 1; n <= nochrom; n++) {
    /* initialization for HMM for LRT */

    /* full-sib */
    s0 = s2 = 0.25;
    s1 = 0.5;
    /* half-sib */
    h0 = h1 = 0.5;
    /* grandparent-granchild */
    g0 = g1 = 0.5;
    /* avuncular, use approximation */
    au0 = au1 = 0.5;
    /* first-cousin, use approximation */
    ac0 = 0.75;
    ac1 = 0.25;
    /* unrelated */
    unr0 = 1;
    /* half-avuncular */
    ahu0 = 0.75;
    ahu1 = 0.25;
    /* half-cousin */
    ahc0 = 0.875;
    ahc1 = 0.125;
    /* halfsib plus firstcousin */
    ahplusc0 = 0.375;
    ahplusc1 = 0.5;
    ahplusc2 = 0.125;
    /* parent-offspring */
    po1 = 1;
    /* MZ twins */
    mz2 = 1;
    
    /* get the location for the first typed marker and the last typed
       marker on each chromosome */
    k = 1;
    while( valid[kkk+k] == 0 )
      ++k;
    validmin = k;
    k = nomark[n];
    while( valid[kkk+k] == 0 )
      k = k - 1;
    validmax = k;

    /* for each chromosome */
    for (i = validmin; i <= validmax; i++) {
      if ( valid[kkk+i] == 1 ) {
        
	f1 = geno1[2*(kkk+i)-1];
        m1 = geno1[2*(kkk+i)];
        f2 = geno2[2*(kkk+i)-1];
        m2 = geno2[2*(kkk+i)];
  
	p2 = cp2(f1, m1, f2, m2, q[n][i]);
        p1 = cp1(f1, m1, f2, m2, q[n][i]);
        p0 = cp0(f1, m1, f2, m2, q[n][i]);

	if ( i != validmax) {
	  /* get the distrance between i_th marker and the next typed
	     marker */
	  k = i + 1;
	  dij = cenm[n][k];
	  while (valid[kkk+k] == 0) {
	    ++k;
	    dij += cenm[n][k];
	  }

	  /* get the transition probabilty for all the reltype */
	  transproball(dij, pij, maptype);

	  s0_next = s0*p0*pij[1][1][1] + s1*p1*pij[1][2][1] +
	    s2*p2*pij[1][3][1];
	  s1_next = s0*p0*pij[1][1][2] + s1*p1*pij[1][2][2] +
	    s2*p2*pij[1][3][2];
	  s2_next = s0*p0*pij[1][1][3] + s1*p1*pij[1][2][3] + 
	    s2*p2*pij[1][3][3];
	  h0_next = h0*p0*pij[2][1][1] + h1*p1*pij[2][2][1];
	  h1_next = h0*p0*pij[2][1][2] + h1*p1*pij[2][2][2];
	  g0_next = g0*p0*pij[3][1][1] + g1*p1*pij[3][2][1];
	  g1_next = g0*p0*pij[3][1][2] + g1*p1*pij[3][2][2];
	  au0_next = au0*p0*pij[4][1][1] + au1*p1*pij[4][2][1];
	  au1_next = au0*p0*pij[4][1][2] + au1*p1*pij[4][2][2];
	  ac0_next = ac0*p0*pij[5][1][1] + ac1*p1*pij[5][2][1];
	  ac1_next = ac0*p0*pij[5][1][2] + ac1*p1*pij[5][2][2];
	
	  unr0_next = unr0*p0;
	  	  
	  ahu0_next = ahu0*p0*pij[7][1][1] + ahu1*p1*pij[7][2][1];
	  ahu1_next = ahu0*p0*pij[7][1][2] + ahu1*p1*pij[7][2][2];
	  ahc0_next = ahc0*p0*pij[8][1][1] + ahc1*p1*pij[8][2][1];
	  ahc1_next = ahc0*p0*pij[8][1][2] + ahc1*p1*pij[8][2][2];
	  ahplusc0_next = ahplusc0*p0*pij[9][1][1] + 
	    ahplusc1*p1*pij[9][2][1] + ahplusc2*p2*pij[9][3][1];
	  ahplusc1_next = ahplusc0*p0*pij[9][1][2] + 
	    ahplusc1*p1*pij[9][2][2] + ahplusc2*p2*pij[9][3][2];
	  ahplusc2_next = ahplusc0*p0*pij[9][1][3] + 
	    ahplusc1*p1*pij[9][2][3] + ahplusc2*p2*pij[9][3][3];

	  /* for parent-offspring, and MZ twins, take into account of
	     genotyping error */
	  po1_next =
	    po1*((1-errorrate)*(1-errorrate)*p1+(1-(1-errorrate)*(1-errorrate))*p0);

	  mz2_next =
	    mz2*((1-errorrate)*(1-errorrate)*p2+(1-(1-errorrate)*(1-errorrate))*p0);

	  
	  s0 = s0_next;
	  s1 = s1_next; 
	  s2 = s2_next;
	  h0 = h0_next;
	  h1 = h1_next;
	  g0 = g0_next;
	  g1 = g1_next;
	  au0 = au0_next;
	  au1 = au1_next;
	  ac0 = ac0_next;
	  ac1 = ac1_next;
	  
	  unr0 = unr0_next;
	  
	  ahu0 = ahu0_next;
	  ahu1 = ahu1_next;
	  ahc0 = ahc0_next;
	  ahc1 = ahc1_next;
	  ahplusc0 = ahplusc0_next;
	  ahplusc1 = ahplusc1_next; 
	  ahplusc2 = ahplusc2_next;

	  po1 = po1_next;
	  mz2 = mz2_next;

	}
      }
    }
    loglhresult[1] += log(s0*p0+s1*p1+s2*p2);
    loglhresult[2] += log(h0*p0+h1*p1);
    loglhresult[3] += log(g0*p0+g1*p1);
    loglhresult[4] += log(au0*p0+au1*p1);
    loglhresult[5] += log(ac0*p0+ac1*p1);
    loglhresult[6] += log(unr0*p0);
        
    loglhresult[7] += log(ahu0*p0+ahu1*p1);
    loglhresult[8] += log(ahc0*p0+ahc1*p1);
    loglhresult[9] += log(ahplusc0*p0+ahplusc1*p1+ahplusc2*p2);
 
    loglhresult[10] +=
      log(po1*((1-errorrate)*(1-errorrate)*p1+(1-(1-errorrate)*(1-errorrate))*p0));
    
    loglhresult[11] +=
      log(mz2*((1-errorrate)*(1-errorrate)*p2+(1-(1-errorrate)*(1-errorrate))*p0));
       
    kkk += nomark[n];
  }
}


/* for simulate genodata of the relaitonships considered */

void simugeno1geno2(double*** q, double*** qcum, int nochrom, 
		    int* nomark, int **noalle, int reltype, 
		    double* chromlength, double** cenm,
		    int* simugeno1, int* simugeno2)
{
  int i, j, k, kkk;
  int *ibd1, *ibd2, ibd;
  int nc1, nc2, ncs, ncnew1, ncnew2;
  double pos1[MC+1], pos2[MC+1], poss[MC+1];
  double posnew1[MC+1], posnew2[MC+1];
  char ind1[MC+1], ind2[MC+1], inds[MC+1];
  char indnew1[MC+1], indnew2[MC+1];
  double sign;
  

  kkk = 0;
  for (i = 1; i <= nochrom; i++) {
   
    ibd1 = (int *)malloc((size_t) ((nomark[i]+1)*sizeof(int)));
    ibd2 = (int *)malloc((size_t) ((nomark[i]+1)*sizeof(int)));
    
    for (j = 1; j <= nomark[i]; j++) 
      ibd1[j] = ibd2[j] = 0;

    /* simulate full-sib */
    if ( reltype == 1 ) {
      nc1 = crossover_new((chromlength[i]), pos1, '1', '2', ind1);
      nc2 = crossover_new((chromlength[i]), pos2, '1', '2', ind2);
      getibd1(pos1, pos2, nc1, nc2, ind1, ind2, 
	      chromlength[i], cenm[i], nomark[i], ibd1);
      
      nc1 = crossover_new((chromlength[i]), pos1,'3','4', ind1);
      nc2 = crossover_new((chromlength[i]), pos2,'3','4', ind2);
      getibd1(pos1, pos2, nc1, nc2, ind1, ind2,
	      chromlength[i], cenm[i], nomark[i], ibd2);
    }

    /*simulate half-sib */
    else if ( reltype == 2 ) {
      nc1 = crossover_new((chromlength[i]), pos1, '1', '2', ind1);
      nc2 = crossover_new((chromlength[i]), pos2, '1', '2', ind2);
      getibd1(pos1, pos2, nc1, nc2, ind1, ind2,
	      chromlength[i], cenm[i], nomark[i], ibd1);
    }
  
    /* simulate grandparent-grandchild */
    else if ( reltype == 3 ) {
      nc1 = crossover_new((chromlength[i]), pos1, '1', '2', ind1);
      nc2 = crossover_new((chromlength[i]), pos2, '3', '4', ind2);
      ncs = crossover_new((chromlength[i]), poss, 'p', 'm', inds);
      ncnew1 = combine(pos1, pos2, poss, nc1, nc2, ncs,
		       ind1, ind2, inds, posnew1, indnew1);
      getibd0(posnew1, ncnew1, indnew1, '1', '2',
	      chromlength[i],cenm[i], nomark[i], ibd1);
    }

    /* simulate avuncular */
    else if ( reltype == 4 ) {
      nc1 = crossover_new((chromlength[i]), pos1, '1', '2', ind1);
      nc2 = crossover_new((chromlength[i]), pos2, '3', '4', ind2);
      ncs = crossover_new((chromlength[i]), poss, 'p', 'm', inds);
      ncnew1 = combine(pos1, pos2, poss, nc1, nc2, ncs, ind1, ind2,
		       inds, posnew1, indnew1);

      nc1 = crossover_new((chromlength[i]), pos1, '1', '2', ind1);
      nc2 = crossover_new((chromlength[i]), pos2, '3', '4', ind2);

      getibd2(pos1, pos2, posnew1, nc1, nc2, ncnew1, 
	      ind1, ind2, indnew1, chromlength[i], cenm[i], nomark[i],
	      ibd1);
    }

    /* simulate first-cousin */
    else if ( reltype == 5 ) {
      nc1 = crossover_new((chromlength[i]), pos1, '1', '2', ind1);
      nc2 = crossover_new((chromlength[i]), pos2, '3', '4', ind2);
      ncs = crossover_new((chromlength[i]), poss, 'p', 'm', inds);
      ncnew1 = combine(pos1, pos2, poss, nc1, nc2, ncs, ind1, ind2,
		       inds, posnew1, indnew1);

      nc1 = crossover_new((chromlength[i]), pos1, '1', '2', ind1);
      nc2 = crossover_new((chromlength[i]), pos2, '3', '4', ind2);
      ncs = crossover_new((chromlength[i]), poss, 'p', 'm', inds);
      ncnew2 = combine(pos1, pos2, poss, nc1, nc2, ncs, ind1, ind2,
		       inds, posnew2, indnew2);

      getibd1(posnew1, posnew2, ncnew1, ncnew2, indnew1, indnew2,
	      chromlength[i], cenm[i], nomark[i], ibd1);
    }

    /* for unrelated reltype == 6 */
    /* no crossover simulation needed */

    /* for half-avuncular */
    else if (reltype == 7) {
      nc1 = crossover_new((chromlength[i]), pos1, '1', '2', ind1);
      nc2 = crossover_new((chromlength[i]), pos2, '3', '4', ind2);
      ncs = crossover_new((chromlength[i]), poss, 'p', 'm', inds);
      ncnew1 = combine(pos1, pos2, poss, nc1, nc2, ncs, ind1, ind2,
		       inds, posnew1, indnew1);
      
      nc1 = crossover_new((chromlength[i]), pos1, '1', '2', ind1);
      
      getibd1(posnew1, pos1, ncnew1, nc1, indnew1, ind1,
	      chromlength[i], cenm[i], nomark[i], ibd1);
    }

    /* for half-cousin */
    else if (reltype == 8) {
      nc1 = crossover_new((chromlength[i]), pos1, '1', '2', ind1);
      nc2 = crossover_new((chromlength[i]), pos2, '3', '4', ind2);
      ncs = crossover_new((chromlength[i]), poss, 'p', 'm', inds);
      ncnew1 = combine(pos1, pos2, poss, nc1, nc2, ncs, ind1, ind2,
		       inds, posnew1, indnew1);

      nc1 = crossover_new((chromlength[i]), pos1, '1', '2', ind1);
      nc2 = crossover_new((chromlength[i]), pos2, '5', '6', ind2);
      ncs = crossover_new((chromlength[i]), poss, 'p', 'm', inds);
      ncnew2 = combine(pos1, pos2, poss, nc1, nc2, ncs, ind1, ind2,
		       inds, posnew2, indnew2);
      
      getibd1(posnew1, posnew2, ncnew1, ncnew2, indnew1, indnew2,
	      chromlength[i], cenm[i], nomark[i], ibd1);
    }
    
    /* for half-sib plus first-cousin */
    else if (reltype == 9) {
      
      /* half-sib part */
      nc1 = crossover_new((chromlength[i]), pos1, '1', '2', ind1);
      nc2 = crossover_new((chromlength[i]), pos2, '1', '2', ind2);
      getibd1(pos1, pos2, nc1, nc2, ind1, ind2,
	      chromlength[i], cenm[i], nomark[i], ibd1);
      /* first-cousin part */
      nc1 = crossover_new((chromlength[i]), pos1, '1', '2', ind1);
      nc2 = crossover_new((chromlength[i]), pos2, '3', '4', ind2);
      ncs = crossover_new((chromlength[i]), poss, 'p', 'm', inds);
      ncnew1 = combine(pos1, pos2, poss, nc1, nc2, ncs, ind1, ind2,
		       inds, posnew1, indnew1);

      nc1 = crossover_new((chromlength[i]), pos1, '1', '2', ind1);
      nc2 = crossover_new((chromlength[i]), pos2, '3', '4', ind2);
      ncs = crossover_new((chromlength[i]), poss, 'p', 'm', inds);
      ncnew2 = combine(pos1, pos2, poss, nc1, nc2, ncs, ind1, ind2,
		       inds, posnew2, indnew2);
 
      getibd1(posnew1, posnew2, ncnew1, ncnew2, indnew1, indnew2,
	      chromlength[i], cenm[i], nomark[i], ibd2);
    }

    /* now based on the ibd value for each marker, 
       simulate allele data*/
    for (j = 1; j <= nomark[i]; j++) {

      /* if not parent-offspring or MZ twins*/
      if( !(reltype == 10 || reltype == 11) ) {
	ibd = ibd1[j] + ibd2[j];

	if ( ibd == 0 ) {
	  
	  sign = 1.0 * rand() / RAND_MAX;
	  for (k = 1; k <= noalle[i][j]; k++) {
	    if ( sign >= qcum[i][j][k-1] && sign <= qcum[i][j][k] )
	      simugeno1[2*(kkk+j)-1] = k;
	  }
	  sign = 1.0 * rand() / RAND_MAX;
	  for (k = 1; k <= noalle[i][j]; k++) {
	    if ( sign >= qcum[i][j][k-1] && sign <= qcum[i][j][k] )
	      simugeno1[2*(kkk+j)] = k;
	  }
	  sign = 1.0 * rand() / RAND_MAX;
	  for (k = 1; k <= noalle[i][j]; k++) {
	    if ( sign >= qcum[i][j][k-1] && sign <= qcum[i][j][k] )
	      simugeno2[2*(kkk+j)-1] = k;
	  }
	  sign=1.0 * rand() / RAND_MAX;
	  for (k = 1; k <= noalle[i][j]; k++) {
	    if ( sign >= qcum[i][j][k-1] && sign <= qcum[i][j][k] )
	      simugeno2[2*(kkk+j)] = k;
	  }
	}
	
	else if ( ibd == 1 ) {
	  sign = 1.0 * rand() / RAND_MAX;
	  for (k = 1; k <= noalle[i][j]; k++) {
	    if ( sign >= qcum[i][j][k-1] && sign <= qcum[i][j][k] )
	      simugeno1[2*(kkk+j)-1] = k;
	  }
	  sign = 1.0 * rand() / RAND_MAX;
	  for (k = 1; k <= noalle[i][j]; k++) {
	    if ( sign >= qcum[i][j][k-1] && sign <= qcum[i][j][k] )
	      simugeno1[2*(kkk+j)] = k;
	  }
	  sign = 1.0 * rand() / RAND_MAX;
	  for (k = 1; k <= noalle[i][j]; k++) {
	    if ( sign >= qcum[i][j][k-1] && sign <= qcum[i][j][k] )
	      simugeno2[2*(kkk+j)-1] = k;
	  }
	  simugeno2[2*(kkk+j)] = simugeno1[2*(kkk+j)];
	}

	else {
	  sign = 1.0 * rand() / RAND_MAX;
	  for (k = 1; k <= noalle[i][j]; k++) {
	    if ( sign >= qcum[i][j][k-1] && sign <= qcum[i][j][k] )
	      simugeno1[2*(kkk+j)-1] = k;
	  }
	  sign = 1.0 * rand() / RAND_MAX;
	  for (k = 1; k <= noalle[i][j]; k++) {
	    if ( sign >= qcum[i][j][k-1] && sign <= qcum[i][j][k] )
	      simugeno1[2*(kkk+j)] = k;
	  }
	  simugeno2[2*(kkk+j)-1] = simugeno1[2*(kkk+j)-1];
	  simugeno2[2*(kkk+j)] = simugeno1[2*(kkk+j)];
	}
      }

      else if (reltype == 10) {
	/* simulate data for parent-offspring with genotypying rate for
	   each marker = errorrate */
	sign = 1.0 * rand() / RAND_MAX;
	/* if error occurr, sample as if unrelated */
	if( sign <= (1-(1-errorrate)*(1-errorrate))) {
	  sign = 1.0 * rand() / RAND_MAX;
	  for (k = 1; k <= noalle[i][j]; k++) {
	    if ( sign >= qcum[i][j][k-1] && sign <= qcum[i][j][k] )
	      simugeno1[2*(kkk+j)-1] = k;
	  }
	  sign = 1.0 * rand() / RAND_MAX;
	  for (k = 1; k <= noalle[i][j]; k++) {
	    if ( sign >= qcum[i][j][k-1] && sign <= qcum[i][j][k] )
	      simugeno1[2*(kkk+j)] = k;
	  }
	  sign = 1.0 * rand() / RAND_MAX;
	  for (k = 1; k <= noalle[i][j]; k++) {
	    if ( sign >= qcum[i][j][k-1] && sign <= qcum[i][j][k] )
	      simugeno2[2*(kkk+j)-1] = k;
	  }
	  sign=1.0 * rand() / RAND_MAX;
	  for (k = 1; k <= noalle[i][j]; k++) {
	    if ( sign >= qcum[i][j][k-1] && sign <= qcum[i][j][k] )
	      simugeno2[2*(kkk+j)] = k;
	  }
	}
	/* if no error, sample as IBD = 1 */
	else {
	  sign = 1.0 * rand() / RAND_MAX;
	  for (k = 1; k <= noalle[i][j]; k++) {
	    if ( sign >= qcum[i][j][k-1] && sign <= qcum[i][j][k] )
	      simugeno1[2*(kkk+j)-1] = k;
	  }
	  sign = 1.0 * rand() / RAND_MAX;
	  for (k = 1; k <= noalle[i][j]; k++) {
	    if ( sign >= qcum[i][j][k-1] && sign <= qcum[i][j][k] )
	      simugeno1[2*(kkk+j)] = k;
	  }
	  sign = 1.0 * rand() / RAND_MAX;
	  for (k = 1; k <= noalle[i][j]; k++) {
	    if ( sign >= qcum[i][j][k-1] && sign <= qcum[i][j][k] )
	      simugeno2[2*(kkk+j)-1] = k;
	  }
	  simugeno2[2*(kkk+j)] = simugeno1[2*(kkk+j)];
	}
      }

      else if (reltype == 11) {
	/* simulate data for MZ twins with genotypying rate for
	   each marker = errorrate */
	sign = 1.0 * rand() / RAND_MAX;
	/* if error occurr, sample as if unrelated */
	if( sign <= (1-(1-errorrate)*(1-errorrate))) {
	  sign = 1.0 * rand() / RAND_MAX;
	  for (k = 1; k <= noalle[i][j]; k++) {
	    if ( sign >= qcum[i][j][k-1] && sign <= qcum[i][j][k] )
	      simugeno1[2*(kkk+j)-1] = k;
	  }
	  sign = 1.0 * rand() / RAND_MAX;
	  for (k = 1; k <= noalle[i][j]; k++) {
	    if ( sign >= qcum[i][j][k-1] && sign <= qcum[i][j][k] )
	      simugeno1[2*(kkk+j)] = k;
	  }
	  sign = 1.0 * rand() / RAND_MAX;
	  for (k = 1; k <= noalle[i][j]; k++) {
	    if ( sign >= qcum[i][j][k-1] && sign <= qcum[i][j][k] )
	      simugeno2[2*(kkk+j)-1] = k;
	  }
	  sign=1.0 * rand() / RAND_MAX;
	  for (k = 1; k <= noalle[i][j]; k++) {
	    if ( sign >= qcum[i][j][k-1] && sign <= qcum[i][j][k] )
	      simugeno2[2*(kkk+j)] = k;
	  }
	}
	/* if no error, sample as IBD = 2 */
	else {
	  sign = 1.0 * rand() / RAND_MAX;
	  for (k = 1; k <= noalle[i][j]; k++) {
	    if ( sign >= qcum[i][j][k-1] && sign <= qcum[i][j][k] )
	      simugeno1[2*(kkk+j)-1] = k;
	  }
	  sign = 1.0 * rand() / RAND_MAX;
	  for (k = 1; k <= noalle[i][j]; k++) {
	    if ( sign >= qcum[i][j][k-1] && sign <= qcum[i][j][k] )
	      simugeno1[2*(kkk+j)] = k;
	  }
	  simugeno2[2*(kkk+j)-1] = simugeno1[2*(kkk+j)-1];
	  simugeno2[2*(kkk+j)] = simugeno1[2*(kkk+j)];
	}
      }
    
    }
    kkk += nomark[i];
  
    free(ibd1);
    free(ibd2);
  }
}

/* input the length of chromosome and position vector 
   output the number of crossover and update the position vector */
int crossover_new(double chromlength, double* pos, char c1, char c2,
                  char* ind)
{
  int  i, K, cn, cn_temp;
  double sign, d0, d_next, d_cum;
  double pos_temp[MC];

  /* first crossover K iid exp(1/10), K is uniform{1,2,3,4,5} */
  sign = 1.0 * rand() / RAND_MAX;
  if ( sign == 1.0 ) 
    K = 5;
  else 
    K = (int)floor(sign/0.2)+1;
  d0 = 0;
  for (i = 1; i <= K; i++) {
    sign = 1.0 * rand() / RAND_MAX;
    d0 += ((sign == 1.0) ? 100000 : (-0.1*log(1-sign)));
  }
  d0 = d0*100;
  /* following crossovers Gamma(5,1/10)= 5 iid exp(1/10) */
  cn_temp = 0;
  d_cum = 0; 
  d_cum += d0;
  if ( d_cum <= chromlength ) {
    ++cn_temp;
    pos_temp[cn_temp] = d_cum;
  }
  else
    pos_temp[1] = d_cum;

  while( d_cum < chromlength && cn_temp > 0 ) {
    d_next = 0;
    K = 5;
    for (i = 1; i <= K; i++) {
      sign = 1.0 * rand() / RAND_MAX;
      d_next += ((sign == 1.0) ? 100000 : (-0.1*log(1-sign)));
    }
    d_next = d_next*100;
    d_cum += d_next;
    if ( d_cum <= chromlength ) {
      ++cn_temp;
      pos_temp[cn_temp] = d_cum;
    }
  }
  pos_temp[cn_temp+1] = d_cum;

  /* with 1/2 prob. delete each crossover point independently and
     reallocate the leftover points in the position vector */
  cn = 0;
  pos[0] = 0;
  for (i = 1; i <= cn_temp; i++) {
    sign = 1.0 * rand() / RAND_MAX;
    if ( sign >= 0.5 ) {
      ++cn;
      pos[cn] = pos_temp[i];
    }
  }
  pos[cn+1] = pos_temp[cn_temp+1];
  
  /* with 1/2 prob. decide the first index, then alternate the index
     for following markers */
  sign = 1.0 * rand() / RAND_MAX;
  if ( sign >= 0.5 ) {
    for (i = 0; i <= cn; i = i+2) {
      ind[i] = c1;
      ind[i+1] = c2;
    }
  }
  else{
    for (i = 0; i <= cn; i = i+2) {
      ind[i] = c2;
      ind[i+1] = c1;
    }
  }
  return cn;
}

/* combine two crossovered strings to a new one
   based on the crossover information */
int combine(double* p1, double* p2, double* ps, int nc1, int nc2, int ncs,
            char* i1, char* i2, char* is, double* pnew, char* inew)
{
  int currs, curr1, curr2, newpos;
  char currstr;

  curr1 = curr2 = newpos = 0;
  for (currs = 1; currs <= ncs; ++currs) {
    currstr = is[currs - 1];
    if ( currstr == 'p' ) {
      while ( curr1 != nc1 && p1[curr1 + 1] <= ps[currs - 1] )
        ++curr1;
      pnew[newpos] = ps[currs - 1];
      inew[newpos] = i1[curr1];
      newpos++;
      while ( curr1 != nc1 && p1[curr1 + 1] <= ps[currs] ) {
        ++curr1;
        pnew[newpos] = p1[curr1];
        inew[newpos] = i1[curr1];
        newpos++;
      }
    }
    else {
      while ( curr2 != nc2 && p2[curr2 + 1] <= ps[currs - 1] )
        ++curr2;
      pnew[newpos] = ps[currs - 1];
      inew[newpos] = i2[curr2];
      newpos++;
      while ( curr2 != nc2 && p2[curr2 + 1] <= ps[currs] ) {
        ++curr2;
        pnew[newpos] = p2[curr2];
        inew[newpos] = i2[curr2];
        newpos++;
      }
    }
  }
  if ( is[ncs] == 'p' ) {
    while ( curr1 <= nc1 && p1[curr1] <= ps[ncs] )
      ++curr1;
    curr1 = curr1 - 1;
    pnew[newpos] = ps[ncs];
    inew[newpos] = i1[curr1];
    newpos++;
    while ( curr1 < nc1 ) {
      ++curr1;
      pnew[newpos] = p1[curr1];
      inew[newpos] = i1[curr1];
      newpos++;
    }
    pnew[newpos] = p1[nc1+1];
  }
  else {
    while ( curr2 <= nc2 && p2[curr2] <= ps[ncs] )
      ++curr2;
    curr2 = curr2 - 1;
    pnew[newpos] = ps[ncs];
    inew[newpos] = i2[curr2];
    newpos++;
    while ( curr2 < nc2 ) {
      ++curr2;
      pnew[newpos] = p2[curr2];
      inew[newpos] = i2[curr2];
      newpos++;
    }
    pnew[newpos] = p2[nc2+1];
  }
  return newpos - 1;
}

/* getibd2(): compare one crossovered strings with the other two
   crossoved strings */
void getibd2(double* p1,double* p2, double* pnew, int nc1,int nc2, 
	     int ncnew, char* i1,char* i2, char* inew, 
	     double chromlength, double* cenm, int nomark, 
	     int* ibd)
{
  int i, currind, curr1, curr2, currnew;
  double  ipos;

  currind = 1;
  curr1 = curr2 = currnew = 0;
  ipos = 0.0;
  
  for (i = 1; i <= nomark; i++) {
    while ( curr1 != nc1 && p1[curr1 + 1] <= ipos )
      ++curr1;
    while ( curr2 != nc2 && p2[curr2 + 1] <= ipos )
      ++curr2;
    while ( currnew != ncnew && pnew[currnew + 1] <= ipos )
      ++currnew;
    if ( i1[curr1] == inew[currnew] || i2[curr2] == inew[currnew] )
      ibd[currind] = 1;
    else
      ibd[currind] = 0;
    
    ++currind;

    if (i < nomark)
      ipos += cenm[i+1];
  }
}

/* getibd1(): compare two crossovered strings 
   important: ibd starts from ibd[0] */
void getibd1(double* pnew1, double* pnew2, int ncnew1, int ncnew2,
	     char* inew1, char* inew2, double chromlength, 
	     double* cenm, int nomark, int* ibd)
{
  int i, currind, curr1, curr2;
  double  ipos;

  currind = 1;
  curr1 = curr2 = 0;
  ipos = 0.0;

  for (i = 1; i <= nomark; i++) {

    while ( curr1 != ncnew1 && pnew1[curr1 + 1] <= ipos )
      ++curr1;
    while ( curr2 != ncnew2 && pnew2[curr2 + 1] <= ipos )
      ++curr2;
    if ( inew1[curr1] == inew2[curr2] )
      ibd[currind] = 1;
    else
      ibd[currind] = 0;

    ++currind;

    if (i < nomark)
      ipos += cenm[i+1];
  }
}

/* getibd0(): compare one crossovered string with one set of orignial
   two indexes */
void getibd0(double* pnew, int ncnew, char* inew, char c1, char c2,
	     double chromlength, double* cenm, int nomark, 
	     int* ibd)
{
  int i, currind, currnew;
  double  ipos;
  
  currind = 1;
  currnew = 0;
  ipos = 0.0;

  for (i = 1; i <= nomark; i++) {
    while ( currnew != ncnew && pnew[currnew + 1] <= ipos )
      ++currnew;
    if ( inew[currnew]==c1 || inew[currnew]==c2 )
      ibd[currind] = 1;
    else
      ibd[currind] = 0;
    ++currind;
      
    if (i < nomark)
      ipos += cenm[i+1];
  }
}

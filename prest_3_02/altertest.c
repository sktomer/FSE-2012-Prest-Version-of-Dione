#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stddef.h>
#include <time.h>
#define MANCE 2000   /* max no. of ancestors for an individual */
#define MAXSTR 5000   /* max characters for a string */
#define MC 200       /* max no. of chiasma   */
#define MYEPS 0.5e-6
#define NR 10         /* no. of relationships considered */
#define NRMLRT 11     /* no. of relationships considered in MLRT */
#define PTHRESH 0.2  /* threshold for EIBD, AIBS and IBS */
#define MITER 5000   /* max no. of iteration in the EM algorithm */
#define NREP 100000 /* replicates number */
#define errorrate 0.01 /* approximate genotyping error rate */

/* for find genodata, and Mendelian error checking */
int findpedgeno(int** ped, int* nomem, FILE *fp_in2, 
		int nochrom, int* nomark, int** noalle,
		int** geno, char** chromnameped);

/* for trans. prob., cond. prob. and other prob.*/
double cp2(int f1, int m1, int f2, int m2, double* q);
double cp1(int f1, int m1, int f2, int m2, double* q);
double cp0(int f1, int m1, int f2, int m2, double* q);
void transprob(double dij, double pibdiibdj[NRMLRT+1][3+1][3+1], 
	       int reltype, int maptype);
void transproball(double dij, double pibdiibdj[NRMLRT+1][3+1][3+1], 
		  int maptype);

/* map functions */
double mapthetatocm(double theta, int maptype);
double mapcmtotheta(double cm, int maptype);

/* EIBD, IABS, IBS */
void get3stat(int* geno1, int* geno2, int* valid, int commark,
              int reltype, int nochrom, int* nomark,
              double p012[4], double*** q,
              double statresult[3+1]);

/* EM algorithm for estimation of p0, p1, p2 (0-1-2 IBD sharing) */
void   EMforp012(int* geno1, int* geno2, int* valid, 
		 double*** q, int nochrom, int* nomark,
		 int reltype, double phat[3+1]);

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
  
  char   buf[MAXSTR+1], ctemp1[MAXSTR+1], ctemp2[MAXSTR+1];

  /* for pairs info */
  int nopairs;
  int **nomem, ***ped, **reltype;
  
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

  /* p012: kinship coef., 0, 1, 2 IBD sharing */
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

  int    commark, *valid, **geno;
  int    *simugeno1, *simugeno2;
  int    countbig[4+1];
  double dtemp[10], phat[3+1];
  double stat_obs[3+1], stat_simu[3+1];
  double loglh_obs[NRMLRT+1], loglh_simu[NRMLRT+1];
  double mlrt_obs, mlrt_simu;
    
  FILE *fp_in1, *fp_in2;
  FILE *fp_idx, *tmpfp;
  FILE *fp_out;


  /* -------------------------------------------------------------- */

  if ( argc != 3 ) {
    printf("\nUsage: altertest altertest_input chromfiles \n\n");
    exit(0);
  }
  if ( (fp_in1 = fopen(argv[1],"r")) == 0 ) {
    printf("\nUnable to open altertest_input: %s\n\n", argv[1]);
    exit(0);
  }
  if ( (fp_in2 = fopen(argv[2],"r")) == 0 ) {
    printf("\nUnable to open chromfiles: %s\n\n", argv[2]);
    exit(0);
  }
  if ( (fp_out = fopen("altertest_out","w")) == 0 ) {
    printf("\nUnable to open output file: altertest_out \n\n");
    exit(0);
  }
  if ( (fp_error = fopen("altertest_errors","w")) == 0 ) {
    printf("\nUnable to open output file: altertest_errors\n\n");
    exit(0);
  }

  /* -------------------------------------------------------------- */
 
  bugexit = 0;
  maptype = 2;
  
  srand((long)time(NULL));

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
	fprintf(fp_error, "Spacing between adjacent markers should be in recombination fraction (theta), not in centi-Morgan\n\n");
	++bugexit;
      }
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
  /* read altertest_infile */

  nopairs = 0;
  mpedid = mid = 1;

  while ( fgets(buf, MAXSTR, fp_in1) ) {
    ++nopairs;
    numread = sscanf(buf, "%d %d %d %d", &itemp[1], &itemp[2],
                     &itemp[3], &itemp[4]);
    
    if ( (int)(log10(itemp[1])+1) > mpedid )
      mpedid = (int)(log10(itemp[1])+1);
    if ( (int)(log10(itemp[2])+1) > mid )
      mid = (int)(log10(itemp[2])+1);
    if ( (int)(log10(itemp[3])+1) > mid )
      mid = (int)(log10(itemp[3])+1);
    if ( numread != 4 ) {
      printf("\nInvalid line in altertest_infile:\n%s\n", buf);
      printf("Correct format is: pedid id1 id2 null_relationship \n\n");
      exit(0);
    }
    if ( !(itemp[4] == 1 || itemp[4] == 2 || itemp[4] == 3 ||
	   itemp[4] == 4 || itemp[4] == 5 || itemp[4] == 6 ||
	   itemp[4] == 7 || itemp[4] == 8 || itemp[4] == 9 || 
	   itemp[4] == 10 || itemp[4] == 11) ) {

      printf("\nInvalid line in altertest_infile:\n%s\n", buf);
      printf("ALTERTEST can only test the following relationship types\n");
      printf(" 1 = full-sib,         2 = half-sib,           3 = grandparent-child,\n");
      printf(" 4 = avuncular,        5 = first-cousin,       6 = unrelated,\n");
      printf(" 7 = half-avuncular,   8 = half-first-cousin,  9 = half-sib+first-cousin,\n");
      printf("10 = parent-offspring,11 = MZ twins \n\n");
      exit(0);
    }
  }
  
  ped = (int ***)malloc((size_t) ((nopairs+1)*sizeof(int**)));
  nomem = (int **)malloc((size_t) ((nopairs+1)*sizeof(int*)));
  reltype = (int **)malloc((size_t) ((nopairs+1)*sizeof(int*)));

  for (i = 1; i <= nopairs; i++) {
    ped[i] = (int **)malloc((size_t) ((2+1)*sizeof(int*)));
    nomem[i] = (int *)malloc((size_t) ((2+1)*sizeof(int)));
    reltype[i] = (int *)malloc((size_t) ((1+1)*sizeof(int)));
    for (j = 1; j <= 2; j++) {
      ped[i][j] = (int *)malloc((size_t) ((1+1)*sizeof(int)));
    }
  }

  rewind(fp_in1);

  for (i = 1; i <= nopairs; i++) {
    fscanf(fp_in1, "%d %d %d %d",
	   &nomem[i][1], &ped[i][1][1], &ped[i][2][1], 
	   &reltype[i][1]);
    nomem[i][2] =2;
  }
  
  fclose(fp_in1);
  
  /* for each pair perform MLRT */
  valid = (int *)malloc((size_t) ((nomark[0]+1)*sizeof(int)));
  geno = (int **)malloc((size_t) ((2+1)*sizeof(int*)));
  for (i = 1; i<= 2; i++)
    geno[i] = (int *)malloc((size_t) ((2*nomark[0]+1)*sizeof(int)));
  simugeno1 = (int *)malloc((size_t) ((2*nomark[0]+1)*sizeof(int)));
  simugeno2 = (int *)malloc((size_t) ((2*nomark[0]+1)*sizeof(int)));

  /* for each candidate pair */
  printf("\n");
  
  for (i = 1; i <= nopairs; i++) {

    bugexit += findpedgeno(ped[i], nomem[i], fp_in2, nochrom,
			   nomark, noalle, geno, chromnameped);

    /* commark: no. of markers typed in both individuals */
    commark = 0;
    for (k = 1; k <= nomark[0]; k++) {
      if ( geno[1][2*k-1] != 0 && geno[1][2*k] != 0 &&
	   geno[2][2*k-1] != 0 && geno[2][2*k] != 0 ) {
	++commark;
	valid[k] = 1;
      }
      else
	valid[k] = 0;
    }
      
    /* testing only for the pair that commark != 0 */
    if (commark != 0) {
	
      mmark = (int)log10(nomark[0]) + 1;
      sprintf(ctemp1, "%%%dd %%%dd %%%dd  %%2d  %%%dd ", mpedid,
	      mid, mid, mmark);
      fprintf(fp_out, ctemp1,  
	      nomem[i][1], ped[i][1][1], ped[i][2][1],
	      reltype[i][1], commark);
      
      countbig[1] = countbig[2] = countbig[3] = countbig[4] = 0;
      
      /* observed EIBD, AIBS, and IBS under the hypothetic null */
      get3stat(geno[1], geno[2], valid,
	       commark, reltype[i][1], nochrom, nomark,
	       p012[reltype[i][1]], q, stat_obs);
        
      EMforp012(geno[1], geno[2], valid, 
		q, nochrom, nomark, reltype[i][1], phat);

      if( !(reltype[i][1] == 6 || reltype[i][1] == 10 || reltype[i][1] == 11))      
	fprintf(fp_out, " %6.4f  %6.4f %6.4f %6.4f  ", 
		stat_obs[1], phat[1], phat[2], phat[3]);
      else 
	fprintf(fp_out,"     NA  %6.4f %6.4f %6.4f  ",
		  phat[1], phat[2], phat[3]);

      /* observed log( likelihood ) under each reltype */
      getapprloglikelihood(geno[1], geno[2], 
			   valid, q, cenm, nochrom, nomark, maptype,
			   loglh_obs);
      
      /* MLRT = max(log(likelihood_alt)) -
	 log(likelihood_null) */	      
      if ( reltype[i][1] != 1 ) 
	dtemp[1] = loglh_obs[1];
      else 
	dtemp[1] = loglh_obs[2];
      for (k = 1; k <= NRMLRT; k++) {
	if ( k != reltype[i][1] ) {
	  if ( loglh_obs[k] > dtemp[1]) 
	    dtemp[1] = loglh_obs[k];
	}
      }
      mlrt_obs = dtemp[1] - loglh_obs[reltype[i][1]];
	    
      /* simulate geotype data for the pair and calculate the
	 statistics */
      for (n = 1; n <= NREP; n++) {

	simugeno1geno2(q, qcum, nochrom, nomark, noalle, reltype[i][1], 
		       chromlength, cenm, simugeno1, simugeno2);
	
	get3stat(simugeno1, simugeno2, valid,
		 commark, reltype[i][1], nochrom, nomark,
		 p012[reltype[i][1]], q, stat_simu);

	getapprloglikelihood(simugeno1, simugeno2, 
			     valid, q, cenm, nochrom, nomark, maptype,
			     loglh_simu);

	if (reltype[i][1] != 1) 
	  dtemp[2] = loglh_simu[1];
	else 
	  dtemp[2] = loglh_simu[2];
	for (k = 1; k <= NRMLRT; k++) {
	  if (k != reltype[i][1]) {
	    if ( loglh_simu[k] > dtemp[2]) 
	      dtemp[2] = loglh_simu[k];
	  }
	}
	mlrt_simu = dtemp[2] - loglh_simu[reltype[i][1]];
	
	for (k = 1; k <= 3; k++) {
	  if (stat_simu[k] > stat_obs[k])
	    ++countbig[k];
	}

	if ( mlrt_simu > mlrt_obs )
	  ++countbig[4];
	      	   
      }/* loop for NREP replicates simulation */

      if ( ! (reltype[i][1] == 6 || reltype[i][1] == 10 || reltype[i][1] == 11)) {
 	for (k = 1; k <= 4; k++) {
	  if ( countbig[k] > NREP / 2.0) {
	    fprintf(fp_out, " %8.6f ", 2*((NREP-countbig[k])/(NREP+0.0)));
	  }
	  else{
	    fprintf(fp_out, " %8.6f ", 2*(countbig[k]/(NREP+0.0)));
	  }
	}
      }
      else if(reltype[i][1] == 6){
	for (k = 1; k <= 2; k++) 
	  fprintf(fp_out, "       NA ");
	for (k = 3; k <= 4; k++) {
	  if ( countbig[k] > NREP / 2.0) {
	    fprintf(fp_out, " %8.6f ", 2*((NREP-countbig[k])/(NREP+0.0)));
	  }
	  else{
	    fprintf(fp_out, " %8.6f ", 2*(countbig[k]/(NREP+0.0)));
	  }
	}
      }
      else if(reltype[i][1] == 10 ) {
	for (k = 1; k <= 1; k++) 
	  fprintf(fp_out, "       NA ");
	for (k = 2; k <= 4; k++) {
	  if ( countbig[k] > NREP / 2.0) {
	    fprintf(fp_out, " %8.6f ", 2*((NREP-countbig[k])/(NREP+0.0)));
		}
	  else{
	    fprintf(fp_out, " %8.6f ", 2*(countbig[k]/(NREP+0.0)));
	  }
	}
      }
      /* for MZ twins, apply only the MLLR test, because the EIBD test
	 does not apply and the ABIS and IBS tests are sensible to the
	 misclassification of the genotyping error rate */
      else if(reltype[i][1] == 11 ) {
	for (k = 1; k <= 3; k++) 
	  fprintf(fp_out, "       NA ");
	for (k = 4; k <= 4; k++) {
	  if ( countbig[k] > NREP / 2.0) {
	    fprintf(fp_out, " %8.6f ", 2*((NREP-countbig[k])/(NREP+0.0)));
	  }
	  else{
	    fprintf(fp_out, " %8.6f ", 2*(countbig[k]/(NREP+0.0)));
	  }
	}
      }
      
      fprintf(fp_out, "\n");
      fflush(fp_out);
      printf("%d_th pair done\n",i);
    }
  }
  printf("\n");
  
  /* -------------------------------------------------------------- */

  free(ped);
  free(nomem);
  free(reltype);
  
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
  free(geno);
  
  free(simugeno1);
  free(simugeno2);
  
  fclose(fp_in2);
  fclose(fp_out);
  fclose(fp_error);

  return 1;
  
} /* end of main */


/* for find genodata, and Mendelian error checking */
int findpedgeno(int** ped, int* nomem, FILE *fp_in2, 
		int nochrom, int* nomark, int** noalle, 
		int** geno, char** chromnameped) 
{
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
      
      while( !feof(fp_ped) && !findit ) {
	
	numread = fscanf(fp_ped,"%d %d",&itemp[1], &itemp[2]);
	if ( numread != 2 ) {
	  printf("\nError reading genotype data for %d in ped %d at %s\n",
		 geno[j][0], nomem[1], chromnameped[k]);
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
  if ( f1 == m1 && f2 == m2 ) 
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

/* EIBD AIBS IBS */
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

/***********************************************************************
 * MEX TEMPLATE FOR DD1M-FILTER
 * ----------------------------
 * Make a copy of this file, give it a new name and do the following:
 * o Write a function (in C) that contains the state equation.
 * o In the same file write a function containing the (multiple)
 *   observation equations.
 * o Specify the name of the file in the "define variable" KALMFILE.
 *   Remember to use " " around the name.
 * o Assign XFUNC to the name of the state equation function and assign
 *   YFUNC to the name of the observation equation function (no " "!).
 * o Compile with the appropriate Matlab command:
 *   >> mex my_dd1m.c kalmlblx.o    % PC-Linux gcc compiler
 *   >> mex my_dd1m.c kalmlblcc.obj % Matlab lcc compiler
 *
 ***********************************************************************/
#define KALMFILE "demosysm.c"
#define XFUNC demotu
#define YFUNC demoobs
/***********************************************************************/

/*
 *     INCLUDE HEADERS
 */
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "mex.h"


/*
 *     DEFINES ASSOCIATED WITH MATRIX MANIPULATION 
 */
#define ON 1
#define OFF 0
#define RUNCHK ON             /* Run-time checks switched on/off. */

/* "Inline" functions.                                                 */
/* ( No run-time checks is performed, when inline functions are used ) */
#define nof_rows(ptm)                      (ptm->row)                           /* See getrows */
#define nof_cols(ptm)                      (ptm->col)                           /* See getcols */
#define vec_len(ptv)                       (ptv->row+ptv->col-1)                /* See length  */
#define get_val(ptm,row_pos,col_pos)       (ptm->mat[row_pos][col_pos])         /* See mget    */
#define put_val(ptm,row_pos,col_pos,value) ((ptm->mat[row_pos][col_pos])=value) /* See mput    */
#define rvget(ptv,element)                 (ptv->mat[0][element])               /* See vget    */
#define cvget(ptv,element)                 (ptv->mat[element][0])               /* See vget    */
#define rvput(ptv,element,value)           ((ptv->mat[0][element])=value)       /* See vput    */
#define cvput(ptv,element,value)           ((ptv->mat[element][0])=value)       /* See vput    */


/* Declaration of the "abstract" data-type. */

typedef struct {               /* Matrix structure for C library  */
	int row;               /* These are meant to be "private" */
	int col;               /* and should only be accessed via */
	double **mat;          /* the "member functions" below.   */
} matrix;


typedef struct {               /* Matrix structure for C library  */
	int row;               /* These are meant to be "private" */
	int col;               /* and should only be accessed via */
	int **mat;             /* the "member functions" below.   */
} intmatrix;

typedef struct {          /* Optional initializations        */
	matrix *init;          /* Initialization parameters       */
	int Aflag;             /* Linear state update (deterministic term)    */
	int Fflag;             /* Linear state update (stochastic term)       */
	int Cflag;             /* Linear output equation (deterministic term) */
	int Gflag;             /* Linear output equation (stochastic term)    */
	matrix *A;             /* State transition matrix         */
	matrix *C;             /* Output sensitivity matrix       */
	matrix *F;             /* Process noise coupling matrix   */
	matrix *G;             /* Observ. noise coupling matrix   */
} optpar;


typedef struct {          /* Observation data structure       */
   int    ny;             /* Dimension of observation vector  */
	int    nobs;           /* Number of observations in stream */
	int    nw;             /* Dimension of obs. noise vector   */
	matrix *y;             /* Mean of process noise            */
   intmatrix *tidx;       /* Time stamps for obs. in .y       */
	matrix *wmean;         /* Mean of measurement noise        */
	matrix *Sw;            /* Root of noise covariance matrix  */
	matrix *y0;            /* Storage of current outp. estimate*/
   matrix *y02;           /* Storage of 2 * output estimate   */
	matrix *ybar;          /* Storage of scaled output estimate*/
	matrix *Sy;            /* Root of outp. covariance matrix  */
	matrix *Sy0;           /* Rectangular root of outp. cov mat*/
	matrix *Sytmp;         /* Storage of intermediate results  */
	matrix *Syx;           /* Root of cross-covariance matrix  */
	matrix *Syw;           /* Root of cross-covariance matrix  */
	matrix *Syx2;          /* Root of cross-covariance matrix  */
	matrix *Syw2;          /* Root of cross-covariance matrix  */
	matrix *K;             /* Kalman gain                      */
	matrix *Sxlong;        /* Long root of covariance matrix   */
	matrix *hSwt;          /* Sw multiplied by h and transposed*/
	int    yidx;           /* Index to next observation        */
	int    lasttime;       /* Sample no. for last observation  */
	matrix *wtmp;          /* temp vector                      */
	matrix *wtmp2;         /* temp vector                      */
	matrix *syp;           /* temp vector                      */
	matrix *sym;           /* temp vector                      */
	int Cflag;             /* Linear output equation (deterministic term) */
	int Gflag;             /* Linear output equation (stochastic term)    */
	matrix *C;             /* Output sensitivity matrix        */
	matrix *G;             /* Observ. noise sensitivity matrix */
	int syx2_start;        /* Start index of Syx2 in Sy        */
	int syw2_start;        /* Start index of Syw2 in Sy        */
	double scal3;          /* Constant used in dd2mfilt        */
} obsstruct;


/* Declaration of the "member functions".   */
matrix *mmake( int, int );
void mfree( matrix* );
void mprint( matrix* );
void merror( char* );
int getrows( matrix* );
int getcols( matrix* );
void minit( matrix* );
void madd( matrix*, matrix*, matrix* );
void mset( matrix*, matrix*);
intmatrix *intmmake( int, int );
void intmfree( intmatrix* );


/*
 *     PROTOTYPE DECLARATION
 */
matrix* mat2sm(const mxArray*);
void sm2mat(mxArray*, matrix*);
intmatrix* mat2intsm(const mxArray*);
void intsm2mat(mxArray*, intmatrix*);
int dd1filtm(int (*xfct)(matrix*, matrix*, matrix*, matrix*, int),
            int (*yfct)(matrix*, matrix*, matrix*, int),
	    matrix*, matrix*, matrix*, matrix*, matrix*, matrix*, 
            obsstruct*, int, optpar*);

#include KALMFILE

/*********************************************************************************
 *                                                                               *
 *    dd1mc gateway routine                                                      *
 *    ---------------------                                                      *
 *                                                                               *
 *    This is a C-version of the Matlab function 'dd1m'.                         *
 *    Type 'help dd1m' from Matlab for information on                            *
 *    how to call the function.                                                  *
 *                                                                               *
 *                                                                               *
 *    Written by Magnus Norgaard                                                 *
 *    LastEditDate: Nov. 11, 2001                                                *
 *                                                                               *
 *********************************************************************************/


/*********************************************************************************
 *                                                                               *
 *                           G A T E W A Y   R O U T I N E                       *
 *                                                                               *
 *********************************************************************************/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /*
   >>>>>>>>>>>>>>>>>>           VARIABLE DECLARATIONS          <<<<<<<<<<<<<<<<<<<
   */
   matrix *xhat_data, *Smat, *xbar, *Sxbar, *Sv, *u;
   optpar *opt;
   int a, samples=0, nx, ny=0, *ptmi, streams, i;
   mxArray  *Matmatrix, *M;
   obsstruct *obs;


  /*
   >>>>>>>>>>>>>>>>      CHECK FOR PROPER NUMBER OF ARGUMENTS      <<<<<<<<<<<<<<<
   */
   if (nrhs<7 || nrhs>8)
      mexErrMsgTxt("Wrong number of input arguments");
   else if (nlhs > 2)
       mexErrMsgTxt("Too many output arguments");
 
  /*
   >>>>>>>>>>>>>>>>>     CONVERT INPUT ARGUMENTS TO SM FORMAT     <<<<<<<<<<<<<<<<
   */
  /* Convert Sxbar */
  a = 1;

  nx = mxGetN(prhs[a]);
  if (nx==0 || nx!=mxGetM(prhs[a]))
     mexErrMsgTxt("Dimension mismatch in P0");
  else
     Sxbar = mat2sm(prhs[a]);
     
  /* Convert xbar and initialize if necessary */
  a = 0;
  if (mxGetN(prhs[a])==0 || mxGetM(prhs[a])==0){
     xbar = mmake(nx,1);
     minit(xbar);
  }
  else
     xbar = mat2sm(prhs[a]);
    
  /* Convert Sv */
  a = 2;
  if (mxGetN(prhs[a])!=mxGetM(prhs[a]))
     mexErrMsgTxt("Dimension mismatch in Q");
  else
     Sv = mat2sm(prhs[a]);


  /* Convert y */
  a  = 5;
  if(!mxIsCell(prhs[a]))
     mexErrMsgTxt("Argument 'y' must be a cell array");
  if(mxGetM(prhs[a])!=1 || mxGetM(prhs[a])>mxGetN(prhs[a]))
     mexErrMsgTxt("Argument 'y' must be a 'row' vector of cells");
  streams = mxGetN(prhs[a]);          /* Observation streams */

  /* Allocate memory for array of 'obs' structures */
  obs = (obsstruct*) malloc(streams*sizeof(obsstruct));
  
  /* Extract each 'y' matrix from cell array */
  for(i=0;i<streams;i++){
     M = mxGetCell(prhs[a],i);
     obs[i].ny = mxGetN(M);
     obs[i].nobs = mxGetM(M);
     if (obs[i].ny==0 || obs[i].nobs==0)
     mexErrMsgTxt("Observation matrix is empty");
     obs[i].y = mat2sm(M);
     ny += obs[i].ny;
  }


  /* Convert timeidx */
  a = 6;
  if(!mxIsCell(prhs[a]))
     mexErrMsgTxt("Argument 'timeidx' must be a cell array");
  if(mxGetM(prhs[a])!=1 || mxGetN(prhs[a])!=streams)
     mexErrMsgTxt("'timeidx' must have the same dimension as 'y'");

  /* Extract each 'timeidx' matrix from cell array */
  for(i=0;i<streams;i++){
     M = mxGetCell(prhs[a],i);
     if (mxGetN(M)==0 || mxGetM(M)==0)
        mexErrMsgTxt("No time stamps - 'timeidx' is empty");
     if (mxGetM(M)!=obs[i].nobs)
        mexErrMsgTxt("Dimension mismatch between dimension of cells in 'y' and 'timeidx'");
     obs[i].tidx = mat2intsm(M);
     if(samples<cvget(obs[i].tidx,obs[i].nobs-1))
         samples=cvget(obs[i].tidx,obs[i].nobs-1);
  }
  
  /* Convert Sw */
  a = 3;
  if(!mxIsCell(prhs[a]))
     mexErrMsgTxt("Argument 'Sw' must be a cell array");
  if(mxGetM(prhs[a])!=1 || mxGetN(prhs[a])!=streams)
     mexErrMsgTxt("'Sw' must have the same dimension as 'y'");

  /* Extract each 'Sw' matrix from cell array */
  for(i=0;i<streams;i++){
     M = mxGetCell(prhs[a],i);
     if (mxGetN(M)!=mxGetM(M))
        mexErrMsgTxt("Cell in 'Sw' not quadratic");
     obs[i].Sw = mat2sm(M);
     obs[i].nw = mxGetN(M);
  }
  

  /* Convert u */
  a = 4;
  if (mxGetN(prhs[a])==0 || mxGetM(prhs[a])==0){
     u = mmake(1,1);
     u->row = 0;
  }
  else{
     u = mat2sm(prhs[a]);
     samples = mxGetM(prhs[a]);
  }


  /* Convert optpar */
  opt = (optpar*) malloc(sizeof(optpar)); /* Allocate mem for par structure */
  opt->Aflag = 0; opt->Fflag=0;
  if (nrhs==8){
    a  = 7;
    Matmatrix = mxGetField(prhs[a], 0, "init");
    if(Matmatrix==NULL){
       opt->init = mmake(1,1);
       minit(opt->init);
    }
    else
       opt->init = mat2sm(Matmatrix);
       
	 
	  /* Explore if linear terms are present */
	  /* Relationship between new and past states linear */
	  Matmatrix = mxGetField(prhs[a], 0, "A");
     if(Matmatrix!=NULL){
		  if(mxGetM(Matmatrix)!=nx || mxGetN(Matmatrix)!=nx)
			  mexErrMsgTxt("optpar.A has the wrong dimension");
        opt->A     = mat2sm(Matmatrix);
        opt->Aflag = 1;
     }

     /* Relationship between states and process noise is linear */
	  Matmatrix = mxGetField(prhs[a], 0, "F");
     if(Matmatrix!=NULL){
		  if(mxGetM(Matmatrix)!=nx || mxGetN(Matmatrix)!=getrows(Sv))
			  mexErrMsgTxt("optpar.F has the wrong dimension");
        opt->F     = mat2sm(Matmatrix);
        opt->Fflag = 1;
     }


	  /* Initialize C matrices if present */
     Matmatrix = mxGetField(prhs[a], 0, "C");
	  for(i=0;i<streams;i++)          /* By default, initialize Cflags=0 */
        obs[i].Cflag = 0;
	  if(Matmatrix!=NULL){
        if(!mxIsCell(Matmatrix))
           mexErrMsgTxt("Argument 'opt.C' must be a cell array");
        for(i=0;i<streams;i++){
			  M = mxGetCell(Matmatrix,i);
           if(mxGetM(M)!=0 && mxGetN(M)!=0){
				  if((mxGetM(M)!=obs[i].ny) || (mxGetN(M)!=nx))
					  mexErrMsgTxt("Wrong dimension of a matrix in 'opt.C'");
				  else{
					  obs[i].C = mat2sm(M);
					  obs[i].Cflag = 1;
				  }
			  }
		  }
	  }

	  /* Initialize G matrices if present */
     Matmatrix = mxGetField(prhs[a], 0, "G");
	  for(i=0;i<streams;i++)          /* By default, initialize Gflags=0 */
        obs[i].Gflag = 0;
	  if(Matmatrix!=NULL){
        if(!mxIsCell(Matmatrix))
           mexErrMsgTxt("Argument 'opt.G' must be a cell array");
        for(i=0;i<streams;i++){
			  M = mxGetCell(Matmatrix,i);
           if(mxGetM(M)!=0 && mxGetN(M)!=0){
				  if((mxGetM(M)!=obs[i].ny) || (mxGetN(M)!=obs[i].nw))
					  mexErrMsgTxt("Wrong dimension of a matrix in 'opt.G'");
				  else{
					  obs[i].G = mat2sm(M);
					  obs[i].Gflag = 1;
				  }
			  }
		  }
	  }
  }
  else{
     opt->init = mmake(1,1);
     minit(opt->init);
	  for(i=0;i<streams;i++)          /* By default, initialize Cflags=0 */
        obs[i].Cflag = 0;
	  for(i=0;i<streams;i++)          /* By default, initialize Gflags=0 */
        obs[i].Gflag = 0;
  }


  /* Allocate memory for output matrices */
  xhat_data = mmake(samples+1,nx);
  Smat      = mmake(samples+1,floor(0.5*(nx*(nx+1))+0.5));



  /*
   >>>>>>>>>>>>>>>>>>>>>>         CALL THE C-ROUTINE         <<<<<<<<<<<<<<<<<<<<<
   */
  dd1filtm(XFUNC,YFUNC, xhat_data, Smat, xbar, Sxbar, Sv, u, obs, streams,opt);


  /*
   >>>>>>>>>>>>>>>>>>>         CREATE OUTPUT MATICES            <<<<<<<<<<<<<<<<<<
   */
  if(nlhs>0){
     plhs[0] = mxCreateDoubleMatrix(getrows(xhat_data),nx,mxREAL);
     sm2mat(plhs[0],xhat_data);
  }

  if(nlhs>1){
     plhs[1] = mxCreateDoubleMatrix(getrows(Smat),getcols(Smat),mxREAL);
     sm2mat(plhs[1],Smat);
  }

  /*
   >>>>>>>>>>>>>>>>>>>>        FREE ARGUMENT MATRICES        <<<<<<<<<<<<<<<<<<<<<
   */
  mfree(xbar); mfree(Sxbar); mfree(Sv); mfree(u);
  if(opt->Aflag) mfree(opt->A);
  if(opt->Fflag) mfree(opt->F);
  mfree(opt->init); free(opt); free(obs);
}

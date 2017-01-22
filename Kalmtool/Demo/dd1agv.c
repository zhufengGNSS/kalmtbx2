/***********************************************************************
 * MEX TEMPLATE FOR DD1-FILTER
 * ---------------------------
 * Make a copy of this file, give it a new name and do the following:
 * o Write a function (in C) that contains the state equation.
 * o In the same file write a function containing the observation equation.
 * o Specify the name of the file in the "define variable" KALMFILE.
 *   Remember to use " " around the name.
 * o Assign XFUNC to the name of the state equation function and assign
 *   YFUNC to the name of the observation equation function (no " "!).
 * o Compile with the appropriate Matlab command:
 *   >> mex my_dd1.c kalmlblx.o    % PC/Linux, gcc
 *   >> mex my_dd1.c kalmlblcc.obj % Matlab lcc compiler
 *
 ***********************************************************************/
#define KALMFILE "agvfct.c"
#define XFUNC agvtu
#define YFUNC agvobs
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

typedef struct {          /* Matrix structure for C library  */
	int row;               /* These are meant to be "private" */
	int col;               /* and should only be accessed via */
	double **mat;          /* the "member functions" below.   */
} matrix;


typedef struct {          /* Matrix structure for C library  */
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
int dd1filt(int (*xfct)(matrix*, matrix*, matrix*, matrix*, int),
            int (*yfct)(matrix*, matrix*, matrix*, int),
	    matrix*, matrix*, matrix*, matrix*, matrix*, matrix*, 
            matrix*, matrix*, intmatrix*, optpar*);

#include KALMFILE

/*********************************************************************************
 *                                                                               *
 *    dd1c gateway routine                                                       *
 *    ------------------------                                                   *
 *                                                                               *
 *    This is a C-version of the Matlab function 'dd1'.                          *
 *    Type 'help dd1' from Matlab for information on                             *
 *    how to call the function.                                                  *
 *                                                                               *
 *                                                                               *
 *    Written by: Magnus Norgaard                                                *
 *    LastEditDate: Nov. 9, 2001                                                 *
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
   matrix *xhat_data, *Smat, *xbar, *Sxbart, *hSvt, *hSwt, *u, *y;
   intmatrix *tidx;
   optpar *opt;
   int a, samples, nx, ny, *ptmi;
   mxArray  *Matmatrix;


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
  /* Convert Sxbart */
  a = 1;
  nx = mxGetN(prhs[a]);
  if (nx==0 || nx!=mxGetM(prhs[a]))
     mexErrMsgTxt("Dimension mismatch in P0");
  else
     Sxbart = mat2sm(prhs[a]);
     
  /* Convert xbar and initialize if necessary */
  a = 0;
  if (mxGetN(prhs[a])==0 || mxGetM(prhs[a])==0){
     xbar = mmake(nx,1);
     minit(xbar);
  }
  else
     xbar = mat2sm(prhs[a]);
    
  /* Convert hSvt */
  a = 2;
  if (mxGetN(prhs[a])!=mxGetM(prhs[a]))
     mexErrMsgTxt("Dimension mismatch in Sv");
  else
     hSvt = mat2sm(prhs[a]);

  /* Convert hSwt */
  a = 3;
  if (mxGetN(prhs[a])!=mxGetM(prhs[a]))
     mexErrMsgTxt("Dimension mismatch in Sw");
  else
     hSwt = mat2sm(prhs[a]);

  /* Convert u */
  a = 6;
  if (mxGetN(prhs[a])==0 || mxGetM(prhs[a])==0){
     mexErrMsgTxt("No time stamps. 'timeidx' is empty");
  }
  else{
     tidx = mat2intsm(prhs[a]);
  }
  
  /* Convert u */
  a = 4;
  if (mxGetN(prhs[a])==0 || mxGetM(prhs[a])==0){
     samples = cvget(tidx,nof_rows(tidx)-1);
     u = mmake(1,1);
     u->row = 0;
  }
  else{
     u = mat2sm(prhs[a]);
     samples = mxGetM(prhs[a]);
  }

  /* Convert y */
  a  = 5;
  ny = mxGetN(prhs[a]);
  if (ny==0 || mxGetM(prhs[a])==0)
     mexErrMsgTxt("Observation matrix, y, is empty");
  else if (mxGetM(prhs[a])!= mxGetM(prhs[a+1]))
     mexErrMsgTxt("'y' and 'timeidx' must have the same number of rows");
  else
     y = mat2sm(prhs[a]);

  /* Convert optpar */
  opt = (optpar*) malloc(sizeof(optpar)); /* Allocate mem for par structure */
  opt->Aflag = 0; opt->Fflag = 0; opt->Cflag = 0; opt->Gflag = 0;
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
		  if(mxGetM(Matmatrix)!=nx || mxGetN(Matmatrix)!=getrows(hSvt))
			  mexErrMsgTxt("optpar.F has the wrong dimension");
        opt->F     = mat2sm(Matmatrix);
        opt->Fflag = 1;
     }
	  
     Matmatrix = mxGetField(prhs[a], 0, "C");
     if(Matmatrix!=NULL){
		  if(mxGetM(Matmatrix)!=ny || mxGetN(Matmatrix)!=nx)
			  mexErrMsgTxt("optpar.C has the wrong dimension");
        opt->C     = mat2sm(Matmatrix);
        opt->Cflag = 1;
     }

     Matmatrix = mxGetField(prhs[a], 0, "G");
     if(Matmatrix!=NULL){
		  if(mxGetM(Matmatrix)!=ny || mxGetN(Matmatrix)!=getrows(hSwt))
			  mexErrMsgTxt("optpar.G has the wrong dimension");
        opt->G     = mat2sm(Matmatrix);
        opt->Gflag = 1;
     }
  }
  else{
     opt->init = mmake(1,1);
     minit(opt->init);
}

  

  /* Allocate memory for output matrices */
  xhat_data = mmake(samples+1,nx);
  Smat      = mmake(samples+1,floor(0.5*(nx*(nx+1))+0.5));


  /*
   >>>>>>>>>>>>>>>>>>>>>>         CALL THE C-ROUTINE         <<<<<<<<<<<<<<<<<<<<<
   */
  dd1filt(XFUNC,YFUNC, xhat_data, Smat, xbar, Sxbart, hSvt, hSwt, u, y, tidx, opt);


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
  intmfree(tidx); mfree(xbar); mfree(Sxbart); mfree(hSwt); mfree(hSvt);
  mfree(u); mfree(y); mfree(opt->init);
  if(opt->Aflag) mfree(opt->A);
  if(opt->Fflag) mfree(opt->F);
  if(opt->Cflag) mfree(opt->C);
  if(opt->Gflag) mfree(opt->G);
  free(opt);
}

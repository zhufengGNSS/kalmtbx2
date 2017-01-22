/***********************************************************************
 * XYTEST
 * ------
 * Test state equation and observation equation before integrating them
 * in the filer routines. 
 * Do the following:
 * o Write a function (in C) that contains the state equation.
 * o In the same file write a function containing the observation equation.
 * o Specify the name of the file in the "define variable" KALMFILE.
 *   Remember to use " " around the name.
 * o Assign XFUNC to the name of the state equation function and assign
 *   YFUNC to the name of the observation equation function (no " "!).
 * o Compile with the appropriate Matlab command:
 *   >> mex xytest.c kalmlblx.o    % PC-Linux gcc compiler
 *   >> mex xytest.c kalmlblcc.obj % Matlab lcc compiler
 * o Call 'xytest' from Matlab as:
 *   [y,x] = xytest(x,u,ny,v,w,init);
 *   All 6 arguments must be passed. Use [] if an argument is not available.
 *
 ***********************************************************************/
#define KALMFILE "fallfct.c"
#define XFUNC falltu
#define YFUNC fallobs
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

/* Inline functions with similar output as the library functions listed below. */
/* ( No run-time checks is performed, when inline functions are used )         */

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

typedef struct {               /* Optional initializations        */
	matrix *wmean;         /* Mean of process noise           */
	matrix *vmean;         /* Mean of measurement noise       */
	matrix *init;          /* Initialization parameters       */
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



/*
 *     PROTOTYPE DECLARATION
 */
matrix* mat2sm(const mxArray*);
void sm2mat(mxArray*, matrix*);
intmatrix* mat2intsm(const mxArray*);
void intsm2mat(mxArray*, intmatrix*);

#include KALMFILE

/*********************************************************************************
 *                                                                               *
 *    xytest gateway routine                                                     *
 *    ----------------------                                                     *
 *                                                                               *
 *    This is a small mex-gateway to test C-functions written by the user.       *
 *                                                                               *
 *                                                                               *
 *    Programmed by: Magnus Norgaard                                             *
 *    LastEditDate : Jan. 14, 2000                                               *
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
   matrix *xout, *yout, *x, *u,*init, *v, *w;
   int a, ny;


  /*
   >>>>>>>>>>>>>>>>      CHECK FOR PROPER NUMBER OF ARGUMENTS      <<<<<<<<<<<<<<<
   */
   if (nrhs!=6)
      mexErrMsgTxt("Wrong number of input arguments");
   else if (nlhs > 2)
       mexErrMsgTxt("Too many output arguments");
       
   /*
   >>>>>>>>>>>>>>>>>     CONVERT INPUT ARGUMENTS TO SM FORMAT     <<<<<<<<<<<<<<<<
   */    
       
   /* Convert "x" */
   a=0;
   if (mxGetN(prhs[a])!=0 && mxGetM(prhs[a])!=0){
      x = mat2sm(prhs[a]);
      xout = mmake(mxGetM(prhs[a]),1);
   }
   else
       mexErrMsgTxt("State vector is empty");
  
   /* Convert "u" */
   a=1;
   if (mxGetN(prhs[a])!=0 && mxGetM(prhs[a])!=0)
      u = mat2sm(prhs[a]);
   else
      u = mmake(1,1);
      
   /* Make yout */
   a=2;
   if (mxGetN(prhs[a])!=1 && mxGetM(prhs[a])!=1)
      mexErrMsgTxt("ny must be scalar");
   else{
      ny = (int)(*mxGetPr(prhs[a]));
      yout = mmake(ny,1);
   }
   
   /* Convert "v" */
   a=3;
   if (mxGetN(prhs[a])!=0 && mxGetM(prhs[a])!=0)
      v = mat2sm(prhs[a]);
   else
      v = mmake(1,1);
      
   /* Convert "w" */
   a=4;
   if (mxGetN(prhs[a])!=0 && mxGetM(prhs[a])!=0)
      w = mat2sm(prhs[a]);
    else
      w = mmake(1,1);

   /* Convert "init" */
   a=5;
   if (mxGetN(prhs[a])!=0 && mxGetM(prhs[a])!=0)
      init = mat2sm(prhs[a]);
   else
      init = mmake(1,1);

   
  /* Initialize state update and observation functions */
   XFUNC(init,x,u,v,-1);
   YFUNC(init,x,w,-1);
   printf("\nInitialization performed.\n\n");
   
   printf("Evaluate functions:\n");
   XFUNC(xout,x,u,v,0);
   YFUNC(yout,x,w,0);
   printf("\nf(x,u,v):\n"); mprint(xout);
   printf("\n\ng(x,w):\n"); mprint(yout);
   
   
   XFUNC(xout,x,u,v,-2);
   YFUNC(yout,x,w,-2);
   printf("\n\nCleaning up performed.\n\n");
   
  /*
   >>>>>>>>>>>>>>>>>>>         CREATE OUTPUT MATICES            <<<<<<<<<<<<<<<<<<
   */
   
  mfree(x); mfree(u); mfree(init); mfree(v); mfree(w);
  if(nlhs>0){
     plhs[0] = mxCreateDoubleMatrix(getrows(yout),1,mxREAL);
     sm2mat(plhs[0],yout);
  }
  if(nlhs>1){
     plhs[1] = mxCreateDoubleMatrix(getrows(xout),1,mxREAL);
     sm2mat(plhs[1],xout);
  }
}       

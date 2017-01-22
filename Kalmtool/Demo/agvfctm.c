#define PMN_PI 3.141592653589793108624468950438
int agvobs(matrix *ybar, matrix *xbar, matrix *wmean, int flag);
int agvtu(matrix *xbar, matrix *xhat, matrix *u, matrix *vmean, int flag);

int agvobs(matrix *ybar, matrix *xbar, matrix *w, int flag)
{
   int i, j;
   double dum;
   double *ptm1, *ptm2;
   static int nx, ny0, ny1;
   static matrix *C;
   
   /* Initializations */
   if (flag == -1){
      nx = 5;
      ny0 = 3;
      ny1 = 2;
      C = mmake(ny0,nx); minit(C);
      put_val(C,0,0,1.0); put_val(C,1,1,1.0); put_val(C,2,2,1.0);   
      return 0;
   }
   
   /* Clean up */
   else if (flag == -2){
      mfree(C);
      return 0;
   }
   
   /* Normal call of function */
   else if (flag == 0){
      ptm1 = C->mat[0];
      for(i=0; i<ny0; i++){
         dum  = cvget(w,i);
         for(j=0, ptm2=xbar->mat[0];j<nx;j++)
	    dum += *(ptm1++)* *(ptm2++);
	 cvput(ybar,i,dum);
      }
      return 0;
   }
   
   else if (flag == 1){
      for(i=0; i<ny1; i++)
	 cvput(ybar,i,cvget(xbar,i)+cvget(w,i));
      return 0;
   }
}



int agvtu(matrix *xbar, matrix *x, matrix *u, matrix *v, int flag)
{
   int N;
   double kenc, u0, u1, s, t;
   static int nx, ny;
   static double k1;
   
   /* Initializations */
   if (flag == -1){
      kenc = 800/2/PMN_PI;
      N  = 36;
      k1 = 1/(kenc*N);
      return 0;
   }
   
   /* Clean up */
   else if (flag == -2){
      return 0;
   }
   
   /* Normal call of function */
   else{
      u0 = cvget(u,0) + cvget(v,0);   /* Add process noise */
      u1 = cvget(u,1) + cvget(v,1);
      
      s  = 0.5*k1*(cvget(x,3)*u0+cvget(x,4)*u1);
      t  = 0.5*k1*(cvget(x,3)*u0-cvget(x,4)*u1)/cvget(x,5);

      cvput(xbar,0,cvget(x,0)+s*cos(cvget(x,2)+t));
      cvput(xbar,1,cvget(x,1)+s*sin(cvget(x,2)+t));
      cvput(xbar,2,cvget(x,2)+2*t);
      cvput(xbar,3,cvget(x,3)); /* Add process noise */
      cvput(xbar,4,cvget(x,4));
      cvput(xbar,5,cvget(x,5));
      return 0;
   }
}



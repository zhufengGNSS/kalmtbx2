int fallobs(matrix *ybar, matrix *xbar, matrix *wmean, int flag);
int falltu(matrix *xbar, matrix *xhat, matrix *u, matrix *vmean, int flag);
int bodydot(matrix *xbar, matrix *xhat, matrix *u, matrix *vmean, int flag);

int fallobs(matrix *y, matrix *x, matrix *w, int flag)
{
   int i, j;
   double tmp;
   static double M, H;
   
   /* Initializations */
   if (flag == -1){
      M = rvget(y,1);
      H = rvget(y,2);
      return 0;
   }
   
   /* Clean up */
   else if (flag == -2){
      return 0;
   }
   
   /* Normal call of function */
   else{
      tmp = cvget(x,0)-H;
      cvput(y,0,sqrt(M*M+tmp*tmp)+cvget(w,0));
      return 0;
   }
}



int falltu(matrix *xbar, matrix *x, matrix *u, matrix *v, int flag)
{
   int i, nx=3;
   double *ptm1, *ptm2, *ptm3, *ptm4;
   static matrix *k, *xtmp;
   static double Ts, mfact;
   
   /* Initializations */
   if (flag == -1){
      Ts   = rvget(xbar,0);
      k    = mmake(nx,1); minit(k);
      xtmp = mmake(nx,1);
      mfact= Ts/6;
      bodydot(xbar,x,u,v,-1);
      return 0;
   }
   
   /* Clean up */
   else if (flag == -2){
      mfree(k); mfree(xtmp);
      return 0;
   }
   
   /* Normal call of function */
   else{
   
      /* Calculate k1, add k1/6 to x and prepare call of k2 */
      bodydot(k,x,u,v,0);
      ptm2=xtmp->mat[0]; ptm3=k->mat[0];
      ptm1=xbar->mat[0]; ptm4=x->mat[0]; 
      for(i=0; i<nx;i++){
         *(ptm1++) = *ptm4 + mfact* *ptm3;          /* xbar = x + k1/6   */ 
	 *(ptm2++) = *(ptm4++) + 0.5*Ts* *(ptm3++); /* xtmp = x + 0.5*k1 */ 
      }
      
      /* Calculate k2, add 2k2/6 to x and prepare call of k3 */
      bodydot(k,xtmp,u,v,0);
      ptm2=xtmp->mat[0]; ptm3=k->mat[0];
      ptm1=xbar->mat[0]; ptm4=x->mat[0]; 
      for(i=0; i<nx;i++){
         *(ptm1++) +=  2*mfact* *ptm3;              /* xbar = xbar + 2*k2/6 */ 
	 *(ptm2++) = *(ptm4++) + 0.5*Ts* *(ptm3++); /* xtmp = x + 0.5*k2    */ 
      }
      
      /* Calculate k3, add 2k3/6 to x and prepare call of k4 */
      bodydot(k,xtmp,u,v,0);
      ptm2=xtmp->mat[0]; ptm3=k->mat[0];
      ptm1=xbar->mat[0]; ptm4=x->mat[0]; 
      for(i=0; i<nx;i++){
         *(ptm1++) +=  2*mfact* *ptm3;              /* xbar = xbar + 2*k3/6 */ 
	 *(ptm2++) = *(ptm4++) + Ts* *(ptm3++);     /* xtmp = x + k3        */ 
      }
      
      /* Calculate k4, add k4/6 to x */
      bodydot(k,xtmp,u,v,0);
      ptm3=k->mat[0]; ptm1=xbar->mat[0]; 
      for(i=0; i<nx;i++){
         *(ptm1++) +=  mfact* *(ptm3++);            /* xbar = xbar + k4/6   */ 
      }
      return 0;
   }
}



int bodydot(matrix *xbar, matrix *x, matrix *u, matrix *v, int flag)
{
   double x2;
   static double gamma;
   
   /* Initializations */
   if (flag == -1){
      gamma = rvget(xbar,3);
      return 0;
   }
   
   /* Clean up */
   else if (flag == -2){
      return 0;
   }
   
   /* Normal call of function */
   else{
      x2 = cvget(x,1);
      cvput(xbar,0,-x2+cvget(v,0));
      cvput(xbar,1,-exp(-gamma*cvget(x,0))*x2*x2*cvget(x,2)+cvget(v,1));
      cvput(xbar,2,cvget(v,2));
      return 0;
   }
}

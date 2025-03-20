/* CONTENTS
void out(void);
void out_tec(void);
void memo_on(void);
*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "func.h"

#define IN(i,j)  ( (i)*(JT+1)+(j) )
#define IC(i,j)  ( (i)*JT+(j) )
#define Ilr(i,j)  ( (i)*JT+(j) )
#define Idu(i,j)  ( (i)*(JT+1)+(j) )

#define  GAM  1.4


extern double *gdf[6],*x,*y,*fllr[5],*fldu[5],*xc,*yc,
       *nlr[3],*tlr[2],*ndu[3],*tdu[2],*vol,
	   *crlr00,*crlr01,*crlr10,*crlr11,*cllr00,*cllr01,*cllr10,*cllr11,
	   *crdu00,*crdu01,*crdu10,*crdu11,*cldu00,*cldu01,*cldu10,*cldu11;

extern double *dxflr[5],*dyflr[5],*dxfdu[5],*dyfdu[5];
extern double t,dt,t_max,coord,turb,Mach,Re_ref,u_ref,mu_ref;
extern int n_stp,max_stp,inp_par,IT,JT;

extern double *wlr[2],*wdu[2],*gdf_[6];
extern double tau;
extern int i_qgd;

extern double *mu_turb_lr, *mu_turb_du;

extern double lxmin,lymin,mu_max,H_Z,fix_dt,P_ref,D_ref;





//*****************************************************
void out(void)
{
 FILE *fpw;
 int i,j,k;
  fpw=fopen("reg_o.dat","wb");
  fwrite(&IT,sizeof(int),1,fpw);
  fwrite(&JT,sizeof(int),1,fpw);
  fwrite(&t,sizeof(double),1,fpw);
  fwrite(&dt,sizeof(double),1,fpw);
   for(i=2; i < IT-2; i++)
    for(j=2; j < JT-2; j++)
     for(k=0; k < 6; k++)
      fwrite(&gdf[k][IC(i,j)],sizeof(double),1,fpw);
   fclose(fpw);  
}

//*****************************************************
void out_tec(char *name)
{
 FILE *fpw;
 int i,j,k,iprn;
 fpw=fopen(name,"w");

 fprintf(fpw," TITLE = \" NS2D-data \"  \n");
 //fprintf(fpw," VARIABLES = \"X\", \"Y\", \"D\", \"U\", \"V\", \"P\",\"MU_TURB\", \"T\", \"Gam\" \n");
 fprintf(fpw," VARIABLES = \"X\", \"Y\", \"D\", \"U\", \"V\", \"P\", \"T\", \"Gam\" \n");
 fprintf(fpw," ZONE F=BLOCK  I=%d  J=%d \n", JT-3, IT-3);
 fprintf(fpw," VARLOCATION=([3-8]=CELLCENTERED) \n");
 fprintf(fpw," AUXDATA time=\" %12.6e\" \n",t);
 
 iprn=1;
 for(i=2; i < IT-1; i++)
  for(j=2; j < JT-1; j++)
   {
    fprintf(fpw," %12.6e ",x[IN(i,j)]);
    iprn++;
    if(iprn > 6) { iprn=1; fprintf(fpw,"\n"); }
   }
 fprintf(fpw,"\n");

 iprn=1;
 for(i=2; i < IT-1; i++)
  for(j=2; j < JT-1; j++)
   {
    fprintf(fpw," %12.6e ",y[IN(i,j)]);
    iprn++;
    if(iprn > 6) { iprn=1; fprintf(fpw,"\n"); }
   }
 fprintf(fpw,"\n");

 iprn=1;
 for(i=2; i < IT-2; i++)
  for(j=2; j < JT-2; j++)
   {
    fprintf(fpw," %12.6e ",gdf[0][IC(i,j)]/D_ref);
    iprn++;
    if(iprn > 6) { iprn=1; fprintf(fpw,"\n"); }
   }
 fprintf(fpw,"\n");

 iprn=1;
 for(i=2; i < IT-2; i++)
  for(j=2; j < JT-2; j++)
   {
    fprintf(fpw," %12.6e ",gdf[1][IC(i,j)]);
    iprn++;
    if(iprn > 6) { iprn=1; fprintf(fpw,"\n"); }
   }
 fprintf(fpw,"\n");

 iprn=1;
 for(i=2; i < IT-2; i++)
  for(j=2; j < JT-2; j++)
   {
    fprintf(fpw," %12.6e ",gdf[2][IC(i,j)]);
    iprn++;
    if(iprn > 6) { iprn=1; fprintf(fpw,"\n"); }
   }
 fprintf(fpw,"\n");

 iprn=1;
 for(i=2; i < IT-2; i++)
  for(j=2; j < JT-2; j++)
   {
    fprintf(fpw," %12.6e ",gdf[3][IC(i,j)]/P_ref);
    iprn++;
    if(iprn > 6) { iprn=1; fprintf(fpw,"\n"); }
   }
 fprintf(fpw,"\n");

 /*iprn=1;
 for(i=2; i < IT-2; i++)
  for(j=2; j < JT-2; j++)
   {
    fprintf(fpw," %12.6e ", mu_ref*pow((gdf[3][IC(i,j)]/gdf[0][IC(i,j)]),0.76)); // mu_turb_arr[IC(i,j)]
    iprn++;
    if(iprn > 6) { iprn=1; fprintf(fpw,"\n"); }
  }
 fprintf(fpw,"\n");*/

 iprn=1;
    for(i=2; i < IT-2; i++)
        for(j=2; j < JT-2; j++)
        {
            fprintf(fpw," %12.4e ",gdf[4][IC(i,j)]);
            iprn++;
            if(iprn > 6) { iprn=1; fprintf(fpw,"\n"); }
        }
    fprintf(fpw,"\n");

 iprn=1;
    for(i=2; i < IT-2; i++)
        for(j=2; j < JT-2; j++)
        {
            fprintf(fpw," %12.4e ",gdf[5][IC(i,j)]);
            iprn++;
            if(iprn > 6) { iprn=1; fprintf(fpw,"\n"); }
        }
 fprintf(fpw,"\n");

 fclose(fpw);
}

//*****************************************************
void memo_on(void)
{
 int k,len;
 len=(IT+1)*(JT+1);

 x=(double*)malloc(len*sizeof(double));
 y=(double*)malloc(len*sizeof(double));
 vol=(double*)malloc(len*sizeof(double));
 xc=(double*)malloc(len*sizeof(double));
 yc=(double*)malloc(len*sizeof(double));

 crlr00=(double*)malloc(len*sizeof(double));
 crlr01=(double*)malloc(len*sizeof(double));
 crlr10=(double*)malloc(len*sizeof(double));
 crlr11=(double*)malloc(len*sizeof(double));
 cllr00=(double*)malloc(len*sizeof(double));
 cllr01=(double*)malloc(len*sizeof(double));
 cllr10=(double*)malloc(len*sizeof(double));
 cllr11=(double*)malloc(len*sizeof(double));
 crdu00=(double*)malloc(len*sizeof(double));
 crdu01=(double*)malloc(len*sizeof(double));
 crdu10=(double*)malloc(len*sizeof(double));
 crdu11=(double*)malloc(len*sizeof(double));
 cldu00=(double*)malloc(len*sizeof(double));
 cldu01=(double*)malloc(len*sizeof(double));
 cldu10=(double*)malloc(len*sizeof(double));
 cldu11=(double*)malloc(len*sizeof(double));

 mu_turb_lr=(double*)malloc(len*sizeof(double));
 mu_turb_du=(double*)malloc(len*sizeof(double));

 for(k=0; k < 2; k++)
  {
  tlr[k]=(double*)malloc(len*sizeof(double));
  tdu[k]=(double*)malloc(len*sizeof(double));
  } 

 for(k=0; k < 3; k++)
  {
  nlr[k]=(double*)malloc(len*sizeof(double));
  ndu[k]=(double*)malloc(len*sizeof(double));
  } 

 for(k=0; k < 5; k++)
  {
  fllr[k]=(double*)malloc(len*sizeof(double));
  dxflr[k]=(double*)malloc(len*sizeof(double));
  dyflr[k]=(double*)malloc(len*sizeof(double));
  fldu[k]=(double*)malloc(len*sizeof(double));
  dxfdu[k]=(double*)malloc(len*sizeof(double));
  dyfdu[k]=(double*)malloc(len*sizeof(double));
  } 
 
 for(k=0; k < 6; k++)
  {
  gdf[k]=(double*)malloc(len*sizeof(double));
if(gdf[k]==NULL) {printf("No memo gdf =%d",k);  exit(-1);}

  gdf_[k]=(double*)malloc(len*sizeof(double));
if(gdf_[k]==NULL) {printf("No memo gdf_ =%d",k);  exit(-1);}
  }
for(k=0; k < 2; k++)
  {
  wlr[k]=(double*)malloc(len*sizeof(double));
  wdu[k]=(double*)malloc(len*sizeof(double));
  }

}

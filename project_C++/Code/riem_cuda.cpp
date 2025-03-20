/* CONTENTS
    int riemi (int lom );
*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define maxit 50
#define ee    0.000001
#define estr  (1.+ee)


//#define	g  (double) 1.4
//#define	ga (double) (0.5*(g+1.)/g)
//#define	gb (double) (0.5*(g-1.)/g)
//#define	gc (double) (2./(g-1.))
//#define	gg (double) (1./g)


#define max(a,b)  (((a) > (b)) ? (a) : (b))
#define min(a,b)  (((a) < (b)) ? (a) : (b))


       
/*-------------------------------riemi--------------------------------*/
       void riemi(double *par,int *alarm, double *g)
       {
       static double ss[2][2],rr[2],cc[2],f[2];
       static double rqp[3][2],r_o,u_o,v_o,p_o,om,s1,s2,s3;
       int k,i,it;
       double rc_1,e_pmin,u_vc,p0,p1,pp,uu,ai,tmp1,tmp2,ga[2],gb[2],gc[2],gg[2];
       static double c[2],rc[2],du,ff,df;
       int lom=0;

       *alarm=0;
       om=0.;

       for (i = 0; i < 2; i++)
       {
           ga[i] = (double) (0.5*(g[i]+1.)/g[i]);
           gb[i] = (double) (0.5*(g[i]-1.)/g[i]);
           gc[i] = (double) (2./(g[i]-1.));
           gg[i] = (double) (1./g[i]);
       }

       for(i=0; i < 2; i++)
        for(k=0; k < 3; k++)
         rqp[k][i]=par[k+3*i]; 

//printf("L: %e %e %e \n",rqp[0][0],rqp[1][0],rqp[2][0]);
//printf("R: %e %e %e \n",rqp[0][1],rqp[1][1],rqp[2][1]);

       for ( i=0; i <= 1; i++ )
       {
       c[i]=sqrt(g[i]*rqp[2][i]/rqp[0][i]);
       rc[i]=rqp[0][i]*c[i];
       }
       e_pmin=min(rqp[2][0],rqp[2][1])*ee;

/*       check existence    */

       u_vc=-c[0]*gc[0]-c[1]*gc[1];
       du=rqp[1][0]-rqp[1][1];
       if(du <= u_vc)
        {
        printf(" ATT: vacuum - no solution  \n");
        *alarm=777;  
        return;
        }

/*      inicial guess  */

       p1=(rqp[2][0]*rc[1]+rqp[2][1]*rc[0]+du*rc[0]*rc[1])/(rc[0]+rc[1]);
       p1=max(p1,e_pmin);
       it=1;
/*      Newton iter  */
       do
        {
         if( it > maxit )
                { printf("ITER > MAXIT \n"); 
                  *alarm=77;
                  return;
                }
         p0=p1;
//************************************************
       df=0.;
       for ( i=0; i <= 1; i++ )
         {
          pp=p0/rqp[2][i];
          if(pp > 1.)
            {
             rc_1=1./rc[i];
             tmp1=1./sqrt(ga[i]*pp+gb[i]);
             f[i]=(p0-rqp[2][i])*tmp1*rc_1;
             tmp2=(g[i]+1.)*pp+3.*g[i]-1.;
             df=df+0.25*tmp2*(tmp1*tmp1*tmp1)*rc_1*gg[i];
             }
          else
             {
             tmp1=pow(pp,gb[i]);
             f[i]=gc[i]*c[i]*(tmp1-1.);
             df=df+c[i]*tmp1*gg[i]/p0;
             }
         }
       ff=f[0]+f[1]-du;
       s2=rqp[1][0]-f[0];
//************************************************
         p1=p0-ff/df;
         it++;
        }
        while ( fabs(p1-p0) > e_pmin );
        
/*     wave coeff. and dens   */
       pp=p0;  uu=s2;
       for ( i=0; i <= 1; i++ )
         {
         ai=2.*( (double) i -0.5);
         if(p0 < rqp[2][i]*estr)
           {
           ss[0][i]=rqp[1][i]+ai*c[i];
           cc[i]=c[i]-0.5*ai*(g[i]-1.)*(rqp[1][i]-uu);
           ss[1][i]=uu+ai*cc[i];
           rr[i]=g[i]*pp/(cc[i]*cc[i]);
           }
         else
           {
           tmp1=(g[i]+1.)*pp+(g[i]-1.)*rqp[2][i];
           tmp2=(g[i]-1.)*pp+(g[i]+1.)*rqp[2][i];
           rr[i]=rqp[0][i]*tmp1/tmp2;
           tmp1=rqp[0][i]-rr[i];
           tmp2=rqp[0][i]*rqp[1][i]-rr[i]*uu;
           ss[0][i]=tmp2/tmp1;
           ss[1][i]=ss[0][i];
           }
         }
         
/*     data on line  */

       if(lom == 1) om=ss[0][0];
       if(lom == 2) om=uu;
       if(lom == 3) om=ss[0][1];
       s1=ss[0][0];
       s3=ss[0][1];
/*           1          */
       if( om <= ss[0][0] )
         {   r_o=rqp[0][0];  u_o=rqp[1][0]; p_o=rqp[2][0];  goto rtn_o;  }
/*           2          */
       if( ( om > ss[0][0] ) && ( om <= ss[1][0]) )
         {
         tmp1=1./(g[0]+1.);
         cc[0]=(2.*c[0]+(g[0]-1.)*(rqp[1][0]-om))*tmp1;
         u_o=om+cc[0];
         r_o=rqp[0][0]*pow( (cc[0]/c[0]),gc[0] );
         p_o=cc[0]*cc[0]*r_o/g[0];
         goto rtn_o;
         }
/*          3           */
       if( (om > ss[1][0]) && (om <= s2) )
         {  r_o=rr[0]; u_o=uu; p_o=pp; goto rtn_o;  }
/*          4           */
       if( (om > s2) && (om <= ss[1][1]) )
         {  r_o=rr[1]; u_o=uu; p_o=pp;  goto rtn_o; }
/*          5           */
       if( (om > ss[1][1]) && (om <= ss[0][1]) )
         {
          tmp2=1./(g[1]+1.);
          cc[1]=(2.*c[1]+(g[1]-1.)*(om-rqp[1][1]))*tmp2;
          u_o=om-cc[1];
          r_o=rqp[0][1]*pow( (cc[1]/c[1]),gc[1] );
          p_o=cc[1]*cc[1]*r_o/g[1];
          goto rtn_o;
         }
/*        6             */
       if( om > ss[0][1] )
         { r_o=rqp[0][1]; u_o=rqp[1][1]; p_o=rqp[2][1];  goto rtn_o; }
rtn_o:
       par[6]=r_o; par[7]=u_o; par[8]=p_o;
       par[9]=s1;  par[10]=s2; par[11]=s3; 
       }

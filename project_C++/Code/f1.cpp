/*
CONTENTS
void inp_mesh(void);
void inp_gdf(void);
void deriv_lr(void);
void deriv_du(void);
void step(void);
double flow_lr(void);
double flow_du(void);
void new_gdf(void);
void new_gdf_corr(void);
*/

//void turb_lr(void);
//void turb_du(void);

//void TGAS4 (double P, double RHO, double *H);
//double GAM_cell (double PIN, double DIN);


#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include "func.h"

#define IN(i,j)   ( (i)*(JT+1)+(j) )
#define IC(i,j)   ( (i)*JT+(j) )
#define Ilr(i,j)  ( (i)*JT+(j) )
#define Idu(i,j)  ( (i)*(JT+1)+(j) )

#define  GAM  1.4
#define PI 3.141592653 
#define  e_VA  1.e-6

#define Rair   287.

double min_mod(double a, double b);

double *gdf[6], *x, *y, *fllr[5], *fldu[5], *xc, *yc,
       *nlr[3], *tlr[2], *ndu[3], *tdu[2], *vol,
       *crlr00, *crlr01, *crlr10, *crlr11, *cllr00, *cllr01, *cllr10, *cllr11,
       *crdu00, *crdu01, *crdu10, *crdu11, *cldu00, *cldu01, *cldu10, *cldu11;

double *dxflr[5], *dyflr[5], *dxfdu[5], *dyfdu[5];
double t, dt, t_max, coord, turb, Mach, Re_ref, u_ref, mu_ref, Tem_ref;
double lxmin, lymin, mu_max, H_Z, fix_dt, P_ref, D_ref;
int n_stp, max_stp, inp_par, IT, JT;

double *wlr[2], *wdu[2], *gdf_[6];
double tau;
int t_ind, gam_ind, i_qgd, i_appr, i_rk;
double RES;

double *mu_turb_lr, *mu_turb_du;
int Flag_turb;

double eforstep = 1.e-13;

char number[11] = "0123456789";
char name[40] = "fld00000.dat";
int namec, namem;

double rade, ye, pe, xue, de;

int TI;
double TW;

int global_file_counter = 1;

//*****************************************************
void update_par_file(double P_ref, double D_ref, double xue, double ye) {
    FILE *fpr = fopen("par_i.dat", "w");
    fprintf(fpr, "132 260\n");
    fprintf(fpr, "80.E-5 0 1.E-6\n");  // t_max = 80.E-5
    fprintf(fpr, "100000 100 0\n");    // nout = 100
    fprintf(fpr, "0.001 2.6E-5\n");
    fprintf(fpr, "0. 0\n");
    fprintf(fpr, "0.1\n");
    fprintf(fpr, "0\n");
    fprintf(fpr, "1\n");
    fprintf(fpr, "0\n");
    fprintf(fpr, "0\n");
    fprintf(fpr, "%lg %lg\n", P_ref, D_ref);
    fprintf(fpr, "%lg %lg %lg %lg\n", 0.001, ye, 2.114, 1.);
    fprintf(fpr, "1 293.\n");
    fprintf(fpr, "%lg\n", xue);
    fclose(fpr);
}

int main(void) {
    double P_ref_values[] = {6132, 12265, 18400, 21450, 24530};
    double D_ref_values[] = {0.07369, 0.14738, 0.22107, 0.25741, 0.29476};
    double xue_values[] = {-0.002, 0.002};
    double ye_values[] = {0.001, 0.002};

    // Перебор всех комбинаций
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            for (int k = 0; k < 2; k++) {
                for (int l = 0; l < 2; l++) {
                    update_par_file(P_ref_values[i], D_ref_values[j], xue_values[k], ye_values[l]);
                    FILE *fpr = fopen("par_i.dat", "r");
                    int lout, nout;

                    fscanf(fpr, "%d %d", &IT, &JT);
                    fscanf(fpr, "%lg %d %lg", &t_max, &t_ind, &fix_dt);
                    fscanf(fpr, "%d %d %d", &max_stp, &nout, &inp_par);
                    fscanf(fpr, "%lg %lg", &Mach, &mu_ref);
                    fscanf(fpr, "%lg %d", &coord, &Flag_turb);
                    fscanf(fpr, "%lg", &H_Z);
                    fscanf(fpr, "%d", &i_qgd);
                    fscanf(fpr, "%d", &i_appr);
                    fscanf(fpr, "%d", &i_rk);
                    fscanf(fpr, "%d", &gam_ind);
                    fscanf(fpr, "%lg %lg", &P_ref, &D_ref);
                    fscanf(fpr, "%lg %lg %lg %lg", &rade, &ye, &pe, &de);
                    fscanf(fpr, "%d %lg", &TI, &TW);
                    fscanf(fpr, "%lg", &xue);

                    fclose(fpr);

                    // Формируем уникальное имя файла

                    t = 0.;
                    dt = 1.e-17;
                    if (t_ind == 1) dt = fix_dt;

                    inp_mesh();
                    inp_gdf();
                    n_stp = 1;
                    lout = 1;
                    //namec = 1;

                    while (t < t_max) {
                        step();
                        t += dt;
                        //printf(" N_step=%d t=%e dt=%e RES/dt=%e \n",n_stp,t,dt,RES/dt);

                        if (lout > nout) {
                            out();
                            printf(" Writing fld file by max time. File number: %05d \n", global_file_counter - 1);

                            sprintf(name, "fld%05d.dat", global_file_counter);
                            global_file_counter++;  // Увеличиваем счетчик
                            // Используем основное имя файла
                            out_tec(name);
                            process_fld(name);

                            lout = 1;
                            //namec++;
                        }

                        n_stp++;
                        lout++;

                        if (n_stp > max_stp) break;
                    }

                    out();
                    printf(" Writing fld file by max time \n");
                    out_tec(name);
                    process_fld(name);
                }
            }
        }
    }

    return 0;
}

void process_fld(const char* name) {
    std::string command = "python make_picture_res.py ";
    command += name;  // Конкатенация строки
    system(command.c_str());  // Выполняем команду
}

//*****************************************************
void step(void)
{
 int i,j,k;
 double dt_lr,dt_du,dt_lr_,dt_du_,lcell,ki,bi;


// BOUNDARY CONDITION:  ghost cells

// down wall
    for(i=2; i < IT-2; i++)
        if (TI == 0)
        {
            gdf[0][IC(i, 0)] = gdf[0][IC(i, 1)] = gdf[0][IC(i, 2)];
            gdf[1][IC(i, 0)] = gdf[1][IC(i, 1)] = -gdf[1][IC(i, 2)];
            gdf[2][IC(i, 0)] = gdf[2][IC(i, 1)] = -gdf[2][IC(i, 2)];
            gdf[3][IC(i, 0)] = gdf[3][IC(i, 1)] = gdf[3][IC(i, 2)];
            gdf[4][IC(i, 0)] = gdf[4][IC(i, 1)] = gdf[4][IC(i, 2)];
            gdf[5][IC(i, 0)] = gdf[5][IC(i,1)]  = gdf[5][IC(i,2)];
        }
        else
        {
            gdf[0][IC(i, 0)] = gdf[0][IC(i, 1)] = gdf[3][IC(i, 2)] / (TW*Rair);
            gdf[1][IC(i, 0)] = gdf[1][IC(i, 1)] = -gdf[1][IC(i, 2)];
            gdf[2][IC(i, 0)] = gdf[2][IC(i, 1)] = -gdf[2][IC(i, 2)];
            gdf[3][IC(i, 0)] = gdf[3][IC(i, 1)] = gdf[3][IC(i, 2)];
            gdf[4][IC(i, 0)] = gdf[4][IC(i, 1)] = TW;
            gdf[5][IC(i, 0)] = gdf[5][IC(i, 1)] = gdf[5][IC(i,2)];
        }                        

 // right out
    for(j=0; j < JT; j++)
     for(k=0; k < 6; k++)
       gdf[k][IC(IT-2,j)]=gdf[k][IC(IT-1,j)]=gdf[k][IC(IT-3,j)];
/* right period
       {
       gdf[k][IC(IT-2,j)]=gdf[k][IC(2,j)];
       gdf[k][IC(IT-1,j)]=gdf[k][IC(3,j)];
       }
*/

// left out
    for(j=0; j < JT; j++)
     for(k=0; k < 6; k++)
       gdf[k][IC(0,j)]=gdf[k][IC(1,j)]=gdf[k][IC(2,j)];
/*  left period
       {
       gdf[k][IC(0,j)]=gdf[k][IC(IT-4,j)];
       gdf[k][IC(1,j)]=gdf[k][IC(IT-3,j)];
       }
*/
//up wall
    for(i=2; i < IT-2; i++)
    {
    if(TI==0)
     { 
        gdf[0][IC(i,JT-1)]=gdf[0][IC(i,JT-2)]= gdf[0][IC(i,JT-3)];
        gdf[1][IC(i,JT-1)]=gdf[1][IC(i,JT-2)]=-gdf[1][IC(i,JT-3)];
        gdf[2][IC(i,JT-1)]=gdf[2][IC(i,JT-2)]=-gdf[2][IC(i,JT-3)];
        gdf[3][IC(i,JT-1)]=gdf[3][IC(i,JT-2)]= gdf[3][IC(i,JT-3)];
        gdf[4][IC(i,JT-1)]=gdf[4][IC(i,JT-2)]= gdf[4][IC(i,JT-3)];
        gdf[5][IC(i,JT-1)]=gdf[5][IC(i,JT-2)]= gdf[5][IC(i,JT-3)];
     }
    else
     {
        gdf[0][IC(i,JT-1)]=gdf[0][IC(i,JT-2)]= gdf[3][IC(i,JT-3)]/(TW*Rair);
        gdf[1][IC(i,JT-1)]=gdf[1][IC(i,JT-2)]=-gdf[1][IC(i,JT-3)];
        gdf[2][IC(i,JT-1)]=gdf[2][IC(i,JT-2)]=-gdf[2][IC(i,JT-3)];
        gdf[3][IC(i,JT-1)]=gdf[3][IC(i,JT-2)]= gdf[3][IC(i,JT-3)];
        gdf[4][IC(i,JT-1)]=gdf[4][IC(i,JT-2)]= TW;
        gdf[5][IC(i,JT-1)]=gdf[5][IC(i,JT-2)]= gdf[5][IC(i,JT-3)];
     }
    }


 deriv_lr();
 deriv_du();
 if(Flag_turb == 1)  turb_lr();
 dt_lr=flow_lr();
 if(Flag_turb == 1)  turb_du();
 dt_du=flow_du();
 new_gdf();

 lcell=min(lxmin,lymin);
 ki=min(dt_lr,dt_du);
 bi=0.25*lcell*lcell/mu_max;
 dt=min(bi,ki);
 dt=H_Z*dt;
 if(t_ind==1) dt=fix_dt;

/*
 for(i=0; i < IT; i++)
  for(j=0; j < JT; j++)
	{
	for(k=0; k < 6; k++)
           if(gdf[k][IC(i,j)]!=gdf[k][IC(i,j)]) {printf("NaN: i=%d j=%d k=%d \n",i,j,k); exit(-1); }

             if(gdf[0][IC(i,j)]<=0.) {printf("d=0: i=%d j=%d k=%d \n",i,j,k); } //exit(-1); }
	     if(gdf[3][IC(i,j)]<=0.) {printf("p=0: i=%d j=%d k=%d \n",i,j,k); } //exit(-1); }
	}
*/
}

//*****************************************************
void inp_mesh(void)
{
 int i,j;
 double dx,dy;
 double t1,t2,t3,t4,t5,t6,t7;
 FILE *fpr;

 fpr=fopen("mesh.dat","r");
// inp coord 
 fscanf(fpr,"%d %d",&IT,&JT);
 memo_on();

for(i=2; i <= IT-2; i++)
  for(j=2; j <= JT-2; j++)

//  for(j=2; j <= JT-2; j++)
//   for(i=2; i <= IT-2; i++)

     fscanf(fpr,"%lg %lg",&x[IN(i,j)],&y[IN(i,j)]);

 lxmin=9.9e+9; lymin=9.9e+9;

// norms left_right
 for(i=2; i <= IT-2; i++)
  for(j=2; j < JT-2; j++)
   {
    dx=x[IN(i,j+1)]-x[IN(i,j)];
    dy=y[IN(i,j+1)]-y[IN(i,j)];
    nlr[0][Ilr(i,j)]=dy;
    nlr[1][Ilr(i,j)]=-dx;
    nlr[2][Ilr(i,j)]=sqrt(dx*dx+dy*dy);
    lxmin=min(lxmin,nlr[2][Ilr(i,j)]);
    tlr[0][Ilr(i,j)]=dx;
    tlr[1][Ilr(i,j)]=dy;
   }
// norms down_up
 for(i=2; i < IT-2; i++)
  for(j=2; j <= JT-2; j++)
   {
    dx=x[IN(i+1,j)]-x[IN(i,j)];
    dy=y[IN(i+1,j)]-y[IN(i,j)];
    ndu[0][Idu(i,j)]=-dy;
    ndu[1][Idu(i,j)]=dx;
    ndu[2][Idu(i,j)]=sqrt(dx*dx+dy*dy);
    lymin=min(lymin,ndu[2][Idu(i,j)]);
    tdu[0][Idu(i,j)]=dx;
    tdu[1][Idu(i,j)]=dy;
   }
// volumes
 for(i=2; i < IT-2; i++)
  for(j=2; j < JT-2; j++)
   {
    dx=x[IN(i+1,j+1)]-x[IN(i,j)];
    dy=y[IN(i,j+1)]-y[IN(i+1,j)];
    vol[IC(i,j)]=0.5*dx*dy;
    dx=x[IN(i,j+1)]-x[IN(i+1,j)];
    dy=y[IN(i+1,j+1)]-y[IN(i,j)];
    vol[IC(i,j)]=vol[IC(i,j)]-0.5*dx*dy;
    xc[IC(i,j)]=0.25*(x[IN(i,j)]+x[IN(i+1,j)]+x[IN(i+1,j+1)]+x[IN(i,j+1)]);
    yc[IC(i,j)]=0.25*(y[IN(i,j)]+y[IN(i+1,j)]+y[IN(i+1,j+1)]+y[IN(i,j+1)]);
   }
for( j=0; j < JT; j++)
{
 vol[IC(0,j)]=vol[IC(2,j)];
 vol[IC(1,j)]=vol[IC(2,j)];
 vol[IC(IT-1,j)]=vol[IC(IT-3,j)];
 vol[IC(IT-2,j)]=vol[IC(IT-3,j)];
}

for( i=0; i < IT; i++)
{
 vol[IC(i,0)]=vol[IC(i,2)];
 vol[IC(i,1)]=vol[IC(i,2)];
 vol[IC(i,JT-1)]=vol[IC(i,JT-3)];
 vol[IC(i,JT-2)]=vol[IC(i,JT-3)];
}

for(j=0; j < JT; j++)
  {
  xc[IC(IT-2,j)]=2.*xc[IC(IT-3,j)]-xc[IC(IT-4,j)];
  yc[IC(IT-2,j)]=2.*yc[IC(IT-3,j)]-yc[IC(IT-4,j)];
  xc[IC(IT-1,j)]=2.*xc[IC(IT-2,j)]-xc[IC(IT-3,j)];
  yc[IC(IT-1,j)]=2.*yc[IC(IT-2,j)]-yc[IC(IT-3,j)];

  xc[IC(1,j)]=2.*xc[IC(2,j)]-xc[IC(3,j)];
  yc[IC(1,j)]=2.*yc[IC(2,j)]-yc[IC(3,j)];
  xc[IC(0,j)]=2.*xc[IC(1,j)]-xc[IC(2,j)];
  yc[IC(0,j)]=2.*yc[IC(1,j)]-yc[IC(2,j)];
  }

for(i=0; i < IT; i++)
  {
  xc[IC(i,JT-2)]=2.*xc[IC(i,JT-3)]-xc[IC(i,JT-4)];
  yc[IC(i,JT-2)]=2.*yc[IC(i,JT-3)]-yc[IC(i,JT-4)];
  xc[IC(i,JT-1)]=2.*xc[IC(i,JT-2)]-xc[IC(i,JT-3)];
  yc[IC(i,JT-1)]=2.*yc[IC(i,JT-2)]-yc[IC(i,JT-3)];
                              
  xc[IC(i,1)]=2.*xc[IC(i,2)]-xc[IC(i,3)];
  yc[IC(i,1)]=2.*yc[IC(i,2)]-yc[IC(i,3)];
  xc[IC(i,0)]=2.*xc[IC(i,1)]-xc[IC(i,2)];
  yc[IC(i,0)]=2.*yc[IC(i,1)]-yc[IC(i,2)];
  }

 fclose(fpr);
}

//*****************************************************
void inp_gdf(void)
{
 int i,j,k;
 double d,u,v,p,e;
 double L_ref;
 double r;
 FILE *fpr;

 d=1.; p=1.;  v=0.;
 u_ref=u=Mach*sqrt(GAM*p/d);

 L_ref=1.;
// mu_ref=u*L_ref*d/Re_ref;
 Tem_ref=P_ref/(Rair*D_ref); 
 for(i=0; i < IT; i++)
  for(j=0; j < JT; j++)
   {
    k=IC(i,j);
    gdf[0][k]=D_ref;
    gdf[1][k]=u;
    gdf[2][k]=v;
    gdf[3][k]=P_ref;
//down
r=sqrt(xc[IC(i,j)]*xc[IC(i,j)]+(yc[IC(i,j)]-ye)*(yc[IC(i,j)]-ye));
if(r <= rade) 
{
  gdf[3][k]=pe*P_ref;
  if(yc[k] < ye-0.5*rade) gdf[3][k]=2.*pe*P_ref; 

}                     

//up
r=sqrt((xc[IC(i,j)]-xue)*(xc[IC(i,j)]-xue)+(yc[IC(i,j)]-0.024+ye)*(yc[IC(i,j)]-0.024+ye));
if(r <= rade)
{
 gdf[3][k]=pe*P_ref;
 if(yc[k] > 0.024-ye+0.5*rade) gdf[3][k]=2.*pe*P_ref; 
}

    gdf[4][k]=gdf[3][k]/(gdf[0][k]*Rair);    
    gdf[5][k]=1.4;
   }

 if(inp_par == 1)
  {
   fpr=fopen("reg_i.dat","rb");
   fread(&IT,sizeof(int),1,fpr);
   fread(&JT,sizeof(int),1,fpr);
   fread(&t,sizeof(double),1,fpr);
   fread(&dt,sizeof(double),1,fpr);
   for(i=2; i < IT-2; i++)
    for(j=2; j < JT-2; j++)
     for(k=0; k < 6; k++)
      fread(&gdf[k][IC(i,j)],sizeof(double),1,fpr);
   fclose(fpr);  
  }
// BOUNDARY CONDITION:  ghost cells

// down wall
    for(i=2; i < IT-2; i++)
        if (TI == 0)
        {
            gdf[0][IC(i, 0)] = gdf[0][IC(i, 1)] = gdf[0][IC(i, 2)];
            gdf[1][IC(i, 0)] = gdf[1][IC(i, 1)] = -gdf[1][IC(i, 2)];
            gdf[2][IC(i, 0)] = gdf[2][IC(i, 1)] = -gdf[2][IC(i, 2)];
            gdf[3][IC(i, 0)] = gdf[3][IC(i, 1)] = gdf[3][IC(i, 2)];
            gdf[4][IC(i, 0)] = gdf[4][IC(i, 1)] = gdf[4][IC(i, 2)];
            gdf[5][IC(i, 0)] = gdf[5][IC(i, 1)] = gdf[5][IC(i,2)];
        }
        else
        {
            gdf[0][IC(i, 0)] = gdf[0][IC(i, 1)] = gdf[3][IC(i, 2)] / (TW*Rair);
            gdf[1][IC(i, 0)] = gdf[1][IC(i, 1)] = -gdf[1][IC(i, 2)];
            gdf[2][IC(i, 0)] = gdf[2][IC(i, 1)] = -gdf[2][IC(i, 2)];
            gdf[3][IC(i, 0)] = gdf[3][IC(i, 1)] = gdf[3][IC(i, 2)];
            gdf[4][IC(i, 0)] = gdf[4][IC(i, 1)] = TW;
            gdf[5][IC(i,0)]=gdf[5][IC(i,1)]= gdf[5][IC(i,2)];
        }
 // right out
    for(j=0; j < JT; j++)
     for(k=0; k < 6; k++)
       gdf[k][IC(IT-2,j)]=gdf[k][IC(IT-1,j)]=gdf[k][IC(IT-3,j)];
/* right period
       {
       gdf[k][IC(IT-2,j)]=gdf[k][IC(2,j)];
       gdf[k][IC(IT-1,j)]=gdf[k][IC(3,j)];
       }
*/
// left out
    for(j=0; j < JT; j++)
     for(k=0; k < 6; k++)
       gdf[k][IC(0,j)]=gdf[k][IC(1,j)]=gdf[k][IC(2,j)];
/* left period
       {
       gdf[k][IC(0,j)]=gdf[k][IC(IT-4,j)];
       gdf[k][IC(1,j)]=gdf[k][IC(IT-3,j)];
       }
*/

//up wall
    for(i=2; i < IT-2; i++)
    {
    if(TI==0)
     { 
        gdf[0][IC(i,JT-1)]=gdf[0][IC(i,JT-2)]= gdf[0][IC(i,JT-3)];
        gdf[1][IC(i,JT-1)]=gdf[1][IC(i,JT-2)]=-gdf[1][IC(i,JT-3)];
        gdf[2][IC(i,JT-1)]=gdf[2][IC(i,JT-2)]=-gdf[2][IC(i,JT-3)];
        gdf[3][IC(i,JT-1)]=gdf[3][IC(i,JT-2)]= gdf[3][IC(i,JT-3)];
        gdf[4][IC(i,JT-1)]=gdf[4][IC(i,JT-2)]= gdf[4][IC(i,JT-3)];
        gdf[5][IC(i,JT-1)]=gdf[5][IC(i,JT-2)]= gdf[5][IC(i,JT-3)];
     }
    else
     {
        gdf[0][IC(i,JT-1)]=gdf[0][IC(i,JT-2)]= gdf[3][IC(i,JT-3)]/(TW*Rair);
        gdf[1][IC(i,JT-1)]=gdf[1][IC(i,JT-2)]=-gdf[1][IC(i,JT-3)];
        gdf[2][IC(i,JT-1)]=gdf[2][IC(i,JT-2)]=-gdf[2][IC(i,JT-3)];
        gdf[3][IC(i,JT-1)]=gdf[3][IC(i,JT-2)]= gdf[3][IC(i,JT-3)];
        gdf[4][IC(i,JT-1)]=gdf[4][IC(i,JT-2)]= TW;
        gdf[5][IC(i,JT-1)]=gdf[5][IC(i,JT-2)]= gdf[5][IC(i,JT-3)];
     }
    }

}

//*****************************************************
void deriv_lr(void)
{
 int i,j,k;
 double dx1,dy1,dx2,dy2,df1,df2;
 double DDx,DDy,DD;
 double xc1,yc1,xc2,yc2,f3,f4;

 double mum,c,gdb[6];


   for(i=2; i <= IT-2; i++)
    for(j=2; j < JT-2; j++)
     {
     xc1=xc[IC(i-1,j)];  
     yc1=yc[IC(i-1,j)];  
     xc2=xc[IC(i,j)];    
     yc2=yc[IC(i,j)];    
     dx1=xc2-xc1;
     dy1=yc2-yc1;
     dx2=x[IN(i,j+1)]-x[IN(i,j)];
     dy2=y[IN(i,j+1)]-y[IN(i,j)];
     DD=dx1*dy2-dx2*dy1;
     for(k=0; k < 5; k++)
      {
      df1=gdf[k][IC(i,j)]-gdf[k][IC(i-1,j)];
      f3=0.25*(gdf[k][IC(i-1,j-1)]+gdf[k][IC(i,j-1)]+
               gdf[k][IC(i,j)]+gdf[k][IC(i-1,j)]); 
      f4=0.25*(gdf[k][IC(i-1,j)]+gdf[k][IC(i,j)]+
               gdf[k][IC(i,j+1)]+gdf[k][IC(i-1,j+1)]); 

      df2=f4-f3;
      DDx=df1*dy2-df2*dy1;
      DDy=df2*dx1-df1*dx2;

      dxflr[k][Ilr(i,j)]=DDx/DD;
      dyflr[k][Ilr(i,j)]=DDy/DD;  
      gdb[k]=0.5*(gdf[k][IC(i,j)]+gdf[k][IC(i-1,j)]);
      }
     }
}
//*****************************************************
void deriv_du(void)
{
 int i,j,k;
 double dx1,dy1,dx2,dy2,df1,df2;
 double DDx,DDy,DD;
 double xc3,yc3,xc4,yc4,f1,f2;

 double mum,c,gdb[6];

   for(i=2; i < IT-2; i++)
    for(j=2; j <= JT-2; j++)
     {
     xc3=xc[IC(i,j-1)];
     yc3=yc[IC(i,j-1)];
     xc4=xc[IC(i,j)]  ;
     yc4=yc[IC(i,j)]  ;
     dx2=xc4-xc3;
     dy2=yc4-yc3;
     dx1=x[IN(i+1,j)]-x[IN(i,j)];
     dy1=y[IN(i+1,j)]-y[IN(i,j)];
     DD=dx1*dy2-dx2*dy1;
     for(k=0; k < 5; k++)
      {
      df2=gdf[k][IC(i,j)]-gdf[k][IC(i,j-1)];
      f1=0.25*(gdf[k][IC(i-1,j-1)]+gdf[k][IC(i,j-1)]+
               gdf[k][IC(i,j)]+gdf[k][IC(i-1,j)]); 
      f2=0.25*(gdf[k][IC(i,j-1)]+gdf[k][IC(i+1,j-1)]+
               gdf[k][IC(i+1,j)]+gdf[k][IC(i,j)]);    
      df1=f2-f1;
      DDx=df1*dy2-df2*dy1;
      DDy=df2*dx1-df1*dx2;
      dxfdu[k][Idu(i,j)]=DDx/DD;
      dyfdu[k][Idu(i,j)]=DDy/DD;  
      gdb[k]=0.5*(gdf[k][IC(i,j)]+gdf[k][IC(i,j-1)]);
      }
      
     }
}

//*****************************************************

double flow_lr(void)
{
 int i,j,k,l,lom,alarm;
 double dx,dy,db,ub,vb,pb,e,NV1,TV1,NV2,TV2;
 double t11,t12,t21,t22,q1,q2;
 double mum,mut,muef,kapm,kapt,kapef;
 double dfx[6],dfy[6],ff[6];
 double xc1,yc1,xc2,yc2;
 double par[32],rqp[3][2],r_o,u_o,v_o,p_o,om,s1,s2,s3;
 double dksi1,dksi2;
 double delf,delf1,delf2,s;
 double w1,w2;
 double dt_lr,dt1,dt2,dt_max;
 double Sb,Sl,Sr;

int n;
double wzx,wzy,rz,dv; 
double pr1,pr2;
double om0,om1,om2,b0,b1,b2,d0,d1,d2,a0,a1,a2;
double g[2];

 lom=0;  om=0.; alarm=0;
 mu_max=0.;
 dt_lr=1.e17;
               
 for(i=2; i <= IT-2; i++)
  for(j=2; j < JT-2; j++)
   {
// left  values for riem
/*
		dx = xc[IC(i-1,j)]-xc[IC(i-2,j)];
		dy = yc[IC(i-1,j)]-yc[IC(i-2,j)];
		dksi1 = sqrt(dx*dx+dy*dy);
		dx = xc[IC(i,j)]-xc[IC(i-1,j)];
		dy = yc[IC(i,j)]-yc[IC(i-1,j)];
		dksi2 = sqrt(dx*dx+dy*dy);
*/
	for( k=0; k < 5; k++)
	{	
		d0=2./3.;  d1=1./3.;
		b0=(gdf[k][IC(i,j)]-gdf[k][IC(i-1,j)])*(gdf[k][IC(i,j)]-gdf[k][IC(i-1,j)]);
		b1=(gdf[k][IC(i-1,j)]-gdf[k][IC(i-2,j)])*(gdf[k][IC(i-1,j)]-gdf[k][IC(i-2,j)]);
		a0=d0/((e_VA+b0)*(e_VA+b0));
		a1=d1/((e_VA+b1)*(e_VA+b1));
		om0=a0/(a0+a1);  om1=a1/(a0+a1);
		ff[k]=om0*(0.5*gdf[k][IC(i-1,j)]+0.5*gdf[k][IC(i,j)])+
				om1*(-0.5*gdf[k][IC(i-2,j)]+1.5*gdf[k][IC(i-1,j)]);
        }

    NV1=(nlr[0][Ilr(i,j)]*ff[1]+nlr[1][Ilr(i,j)]*ff[2])/nlr[2][Ilr(i,j)];
    TV1=(tlr[0][Ilr(i,j)]*ff[1]+tlr[1][Ilr(i,j)]*ff[2])/nlr[2][Ilr(i,j)];
    rqp[0][0]=ff[0];
    rqp[1][0]=NV1;
    rqp[2][0]=ff[3];

// right values for riem  
/*
		dx = xc[IC(i,j)]-xc[IC(i-1,j)];
		dy = yc[IC(i,j)]-yc[IC(i-1,j)];
		dksi1 = sqrt(dx*dx+dy*dy);
		dx = xc[IC(i+1,j)] - xc[IC(i,j)];
		dy = yc[IC(i+1,j)] - yc[IC(i,j)];
		dksi2 = sqrt(dx*dx+dy*dy);
*/
	for( k=0; k < 5; k++)
	{
		d0=2./3.;  d1=1./3.;
                b0=(gdf[k][IC(i,j)]-gdf[k][IC(i-1,j)])*(gdf[k][IC(i,j)]-gdf[k][IC(i-1,j)]);
		b1=(gdf[k][IC(i+1,j)]-gdf[k][IC(i,j)])*(gdf[k][IC(i+1,j)]-gdf[k][IC(i,j)]);
		a0=d0/((e_VA+b0)*(e_VA+b0));
		a1=d1/((e_VA+b1)*(e_VA+b1));
		om0=a0/(a0+a1);  om1=a1/(a0+a1);
		ff[k]=om0*(0.5*gdf[k][IC(i-1,j)]+0.5*gdf[k][IC(i,j)])+
				om1*(1.5*gdf[k][IC(i,j)]-0.5*gdf[k][IC(i+1,j)]);
	} 

    NV2=(nlr[0][Ilr(i,j)]*ff[1]+nlr[1][Ilr(i,j)]*ff[2])/nlr[2][Ilr(i,j)];
    TV2=(tlr[0][Ilr(i,j)]*ff[1]+tlr[1][Ilr(i,j)]*ff[2])/nlr[2][Ilr(i,j)];
    rqp[0][1]=ff[0];
    rqp[1][1]=NV2;
    rqp[2][1]=ff[3];

    for(l=0; l < 2; l++)
		for(k=0; k < 3; k++)
			par[k+3*l]=rqp[k][l]; 

    g[1] = gdf[5][IC(i,j)];
    g[0] = gdf[5][IC(i - 1,j)];

    riemi(par,&alarm,g);
        if(alarm )
          {
          printf("ALARM_RIEM FROM l_r : i=%d j=%d  \n",i,j);
          printf("L: %e %e %e \n",rqp[0][0],rqp[1][0],rqp[2][0]);
          printf("R: %e %e %e \n",rqp[0][1],rqp[1][1],rqp[2][1]);
          exit(-1);
          }

     r_o=par[6];  u_o=par[7] ;   p_o=par[8] ;
     s1 =par[9];  s2 =par[10];   s3 =par[11]; 
/*
	 if(r_o!=r_o) {printf("NaN in d : i=%d j=%d  \n",i,j); exit(-1); }
	 if(u_o!=u_o) {printf("NaN in u : i=%d j=%d  \n",i,j); exit(-1); }
	 if(p_o!=p_o) {printf("NaN in p : i=%d j=%d  \n",i,j); exit(-1); }
*/

 double gb;
     gb=gdf[5][IC(i-1,j)];
     v_o=TV1;
     if(s2 < 0.)
		{
			v_o=TV2;
			gb=gdf[5][IC(i,j)];
		} 

     db=r_o; pb=p_o;
     ub=(u_o*nlr[0][Ilr(i,j)]+v_o*tlr[0][Ilr(i,j)])/nlr[2][Ilr(i,j)];
     vb=(u_o*nlr[1][Ilr(i,j)]+v_o*tlr[1][Ilr(i,j)])/nlr[2][Ilr(i,j)];
     e=pb/(gb-1.)+0.5*db*(ub*ub+vb*vb);

// Eu_flux
	w1=0.; 
	w2=0.; 
     fllr[0][Ilr(i,j)]=db*(ub-w1)*nlr[0][Ilr(i,j)]+db*(vb-w2)*nlr[1][Ilr(i,j)];
     fllr[1][Ilr(i,j)]=(db*ub*(ub-w1)+pb)*nlr[0][Ilr(i,j)]+db*(vb-w2)*ub*nlr[1][Ilr(i,j)];
     fllr[2][Ilr(i,j)]=db*vb*(ub-w1)*nlr[0][Ilr(i,j)]+(db*vb*(vb-w2)+pb)*nlr[1][Ilr(i,j)];
     fllr[3][Ilr(i,j)]=(ub-w1)*(e+pb)*nlr[0][Ilr(i,j)]+(vb-w2)*(e+pb)*nlr[1][Ilr(i,j)];


// NS_flux
// !!!!!!!!!!!!!!!!!!!!
     if(Flag_turb == 1)
		{
		muef = mu_ref * (pow((pb/(db*Rair*Tem_ref)),0.76) + mu_turb_lr[Ilr(i,j)]);
		}
     else
		{
		muef=mu_ref*pow((pb/(db*Rair*Tem_ref)),0.76);
		}
 
	 t11=2.*(2.*dxflr[1][Ilr(i,j)]-dyflr[2][Ilr(i,j)])/3.;
	 t12=(dyflr[1][Ilr(i,j)]+dxflr[2][Ilr(i,j)]);
	 t22=2.*(2.*dyflr[2][Ilr(i,j)]-dxflr[1][Ilr(i,j)])/3.;

     kapef=muef/0.72;                                                                    
     kapef*=3.5;

     if(mu_max < muef) mu_max=muef;

     t11*=muef; 
     t12*=muef; 
     t21=t12;
     t22*=muef; 

 n=Ilr(i,j);
 if(coord > 0.)
	{
	t11-=2.*muef*vb/(3.*yc[IC(i,j)]);
	t22-=2.*muef*vb/(3.*yc[IC(i,j)]);
	}

     q1=kapef*dxflr[4][Ilr(i,j)];
     q2=kapef*dyflr[4][Ilr(i,j)];
     
     fllr[1][Ilr(i,j)]-=t11*nlr[0][Ilr(i,j)]+t12*nlr[1][Ilr(i,j)];
     fllr[2][Ilr(i,j)]-=t21*nlr[0][Ilr(i,j)]+t22*nlr[1][Ilr(i,j)];
     fllr[3][Ilr(i,j)]-=(ub*t11+vb*t12+q1)*nlr[0][Ilr(i,j)]+
                        (ub*t21+vb*t22+q2)*nlr[1][Ilr(i,j)];
// calc time-step
	 if(i > 2 && i < IT-2)
		{
		dksi1=ndu[2][Idu(i-1,j)];
		dksi2=ndu[2][Idu(i,j)];
		dt1=1.e17;
		if( s1 < 0.) dt1=dksi1/(fabs(s1)+eforstep);
		dt2=1.e17;
		if( s3 > 0.) dt2=dksi2/(fabs(s3)+eforstep);
		dt1=min(dt1,dt2);
		dt2=0.25*dksi1*dksi2/muef;
		dt_max=min(dt1,dt2); 
		dt_lr=min(dt_max,dt_lr);
		}


   }

return dt_lr;
}


//*****************************************************
/*
double min_mod(double a, double b)
{
 double tmp;
 tmp=0.;

 if( a*b > 0.)
  {
   tmp=a;
   if( fabs(a) > fabs(b) )  tmp=b;
  }
  return tmp;
}
*/

double min_mod(double a, double b)
// Van Albada
{
 double tmp;
 tmp=((b*b+eforstep)*a+(a*a+eforstep)*b)/(a*a+b*b+2.*eforstep);
  return tmp;
}


//*****************************************************

double  flow_du(void)
{
 int i,j,k,l,lom,alarm;
 double dx,dy,db,ub,vb,pb,e,NV1,TV1,NV2,TV2;
 double t11,t12,t21,t22,q1,q2;
 double mum,mut,muef,kapm,kapt,kapef;
 double dfx[6],dfy[6],ff[6];
 double xc3,yc3,xc4,yc4;
 double par[32],rqp[3][2],r_o,u_o,v_o,p_o,om,s1,s2,s3;
 double deta1,deta2;
 double delf,delf1,delf2,s;
 double w1,w2;
 double dt_du,dt1,dt2,dt_max;
 double Sb,Sl,Sr;

int n;
double wzx,wzy,rz,dv; 
double pr1, pr2;
double om0,om1,om2,b0,b1,b2,d0,d1,d2,a0,a1,a2;
double g[2];

 lom=0;  om=0.;
 mu_max=0.;
 dt_du=1.e17;

 for(i=2; i < IT-2; i++)
  for(j=2; j <= JT-2; j++)
   {
// left(down)  values for riem
/*
	dx=xc[IC(i,j-1)]-xc[IC(i,j-2)];
	dy=yc[IC(i,j-1)]-yc[IC(i,j-2)];
	deta1=sqrt(dx*dx+dy*dy);
	dx=xc[IC(i,j)]-xc[IC(i,j-1)];
	dy=yc[IC(i,j)]-yc[IC(i,j-1)];
	deta2=sqrt(dx*dx+dy*dy);
*/
	for( k=0; k < 5; k++)
	{ 
		d0=2./3.;  d1=1./3.;
		b0=(gdf[k][IC(i,j)]-gdf[k][IC(i,j-1)])*(gdf[k][IC(i,j)]-gdf[k][IC(i,j-1)]);
		b1=(gdf[k][IC(i,j-1)]-gdf[k][IC(i,j-2)])*(gdf[k][IC(i,j-1)]-gdf[k][IC(i,j-2)]);
		a0=d0/((e_VA+b0)*(e_VA+b0));
		a1=d1/((e_VA+b1)*(e_VA+b1));
		om0=a0/(a0+a1);  om1=a1/(a0+a1);
		ff[k]=om0*(0.5*gdf[k][IC(i,j-1)]+0.5*gdf[k][IC(i,j)])+
				om1*(-0.5*gdf[k][IC(i,j-2)]+1.5*gdf[k][IC(i,j-1)]);
        }

    NV1=(ndu[0][Idu(i,j)]*ff[1]+ndu[1][Idu(i,j)]*ff[2])/ndu[2][Idu(i,j)];
    TV1=(tdu[0][Idu(i,j)]*ff[1]+tdu[1][Idu(i,j)]*ff[2])/ndu[2][Idu(i,j)];
    rqp[0][0]=ff[0];
    rqp[1][0]=NV1;
    rqp[2][0]=ff[3];

// right (up) values for riem	
/*
		dx=xc[IC(i,j)]-xc[IC(i,j-1)];
		dy=yc[IC(i,j)]-yc[IC(i,j-1)];
		deta1=sqrt(dx*dx+dy*dy);
		dx=xc[IC(i,j+1)]-xc[IC(i,j)];
		dy=yc[IC(i,j+1)]-yc[IC(i,j)];
		deta2=sqrt(dx*dx+dy*dy);
*/
	for( k=0; k < 5; k++)
	{	
		d0=2./3.;  d1=1./3.;
		b0=(gdf[k][IC(i,j)]-gdf[k][IC(i,j-1)])*(gdf[k][IC(i,j)]-gdf[k][IC(i,j-1)]);
		b1=(gdf[k][IC(i,j+1)]-gdf[k][IC(i,j)])*(gdf[k][IC(i,j+1)]-gdf[k][IC(i,j)]);
		a0=d0/((e_VA+b0)*(e_VA+b0));
		a1=d1/((e_VA+b1)*(e_VA+b1));
		om0=a0/(a0+a1);  om1=a1/(a0+a1);
		ff[k]=om0*(0.5*gdf[k][IC(i,j-1)]+0.5*gdf[k][IC(i,j)])+
				om1*(1.5*gdf[k][IC(i,j)]-0.5*gdf[k][IC(i,j+1)]);
         }

    NV2=(ndu[0][Idu(i,j)]*ff[1]+ndu[1][Idu(i,j)]*ff[2])/ndu[2][Idu(i,j)];
    TV2=(tdu[0][Idu(i,j)]*ff[1]+tdu[1][Idu(i,j)]*ff[2])/ndu[2][Idu(i,j)];
    rqp[0][1]=ff[0];
    rqp[1][1]=NV2;
    rqp[2][1]=ff[3];
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     if(j==2)
      {
      NV1=-NV2;
      TV1=-TV2;
      rqp[0][0]=rqp[0][1];
      rqp[1][0]=-rqp[1][1];
      rqp[2][0]=rqp[2][1];
      }
      if(j==JT-2)
      {
      NV2=-NV1;
      TV2=TV1;
      rqp[0][1]=rqp[0][0];
      rqp[1][1]=-rqp[1][0];
      rqp[2][1]=rqp[2][0];
      }

       for(l=0; l < 2; l++)
        for(k=0; k < 3; k++)
         par[k+3*l]=rqp[k][l]; 

       g[1] = gdf[5][IC(i,j)];
       g[0] = gdf[5][IC(i,j-1)];

       riemi(par,&alarm,g);
        if(alarm )
          {
          printf("ALARM_RIEM FROM l_r : i=%d j=%d  \n",i,j);
          printf("L: %e %e %e \n",rqp[0][0],rqp[1][0],rqp[2][0]);
          printf("R: %e %e %e \n",rqp[0][1],rqp[1][1],rqp[2][1]);
          exit(-1);
          }
			

     r_o=par[6];  u_o=par[7] ;   p_o=par[8] ;
     s1 =par[9];  s2 =par[10];   s3 =par[11]; 
/*
	 if(r_o!=r_o) {printf("NaN in d : i=%d j=%d  \n",i,j); exit(-1); }
	 if(u_o!=u_o) {printf("NaN in u : i=%d j=%d  \n",i,j); exit(-1); }
	 if(p_o!=p_o) {printf("NaN in p : i=%d j=%d  \n",i,j); exit(-1); }
*/

	 double gb;
     gb = gdf[5][IC(i,j-1)];
     v_o=TV1;
     if(s2  < 0.) 
		{
	 	v_o=TV2;
                gb = gdf[5][IC(i,j)];
		}
     db=r_o; pb=p_o;
     ub=(u_o*ndu[0][Idu(i,j)]+v_o*tdu[0][Idu(i,j)])/ndu[2][Idu(i,j)];
     vb=(u_o*ndu[1][Idu(i,j)]+v_o*tdu[1][Idu(i,j)])/ndu[2][Idu(i,j)];

 //down, up wall boundary

	 if(j==2 || j==JT-2) {ub=0.; vb=0.;}
	 e=pb/(gb-1.)+0.5*db*(ub*ub+vb*vb);

// Eu_flux
			w1=0.;
			w2=0.;

     fldu[0][Idu(i,j)]=db*(ub-w1)*ndu[0][Idu(i,j)]+db*(vb-w2)*ndu[1][Idu(i,j)];
     fldu[1][Idu(i,j)]=(db*ub*(ub-w1)+pb)*ndu[0][Idu(i,j)]+db*(vb-w2)*ub*ndu[1][Idu(i,j)];
     fldu[2][Idu(i,j)]=db*vb*(ub-w1)*ndu[0][Idu(i,j)]+(db*vb*(vb-w2)+pb)*ndu[1][Idu(i,j)];
     fldu[3][Idu(i,j)]=(ub-w1)*(e+pb)*ndu[0][Idu(i,j)]+(vb-w2)*(e+pb)*ndu[1][Idu(i,j)];


// NS_flux
     if(Flag_turb == 1)
		{
		muef = mu_ref * (pow((pb/(db*Rair*Tem_ref)),0.76) + mu_turb_du[Idu(i,j)]);
		}
     else
		{
		muef=mu_ref*pow((pb/(db*Rair*Tem_ref)),0.76);
		}

	 t11=2.*(2.*dxfdu[1][Idu(i,j)]-dyfdu[2][Idu(i,j)])/3.;
	 t12=(dyfdu[1][Idu(i,j)]+dxfdu[2][Idu(i,j)]);   
	 t22=2.*(2.*dyfdu[2][Idu(i,j)]-dxfdu[1][Idu(i,j)])/3.;

     kapef=muef/0.72; 
     kapef*=3.5;
     if(mu_max < muef) mu_max=muef;

     t11*=muef; 
     t12*=muef; 
     t21=t12;
     t22*=muef; 

	 n=Idu(i,j);
	 if(coord > 0.)
		{
		t11-=2.*muef*vb/(3.*yc[IC(i,j)]);
		t22-=2.*muef*vb/(3.*yc[IC(i,j)]);
		}
     q1=kapef*dxfdu[4][Idu(i,j)];
     q2=kapef*dyfdu[4][Idu(i,j)];
     
     fldu[1][Idu(i,j)]-=t11*ndu[0][Idu(i,j)]+t12*ndu[1][Idu(i,j)];
     fldu[2][Idu(i,j)]-=t21*ndu[0][Idu(i,j)]+t22*ndu[1][Idu(i,j)];
     fldu[3][Idu(i,j)]-=(ub*t11+vb*t12+q1)*ndu[0][Idu(i,j)]+
                        (ub*t21+vb*t22+q2)*ndu[1][Idu(i,j)];

	 if( j > 2 && j < JT-2)
		{
		deta1=nlr[2][Ilr(i,j-1)];
		deta2=nlr[2][Ilr(i,j)];
		dt1=1.e17;
		if(s1 < 0.) dt1=deta1/(fabs(s1)+eforstep);
		dt2=1.e17;
		if(s3 > 0.) dt2=deta2/(fabs(s3)+eforstep);
		dt1=min(dt1,dt2);
		dt2=0.25*deta1*deta2/muef;
		dt_max=min(dt1,dt2);
		dt_du=min(dt_max,dt_du);
		}

                                                  
   }
return dt_du;
}


//*****************************************************

void new_gdf(void)
{
 int i,j,n,k;
 double p,e;
 double rhs[5],a0[5],a1[5],a_[5],fl[5],dxfc[5],dyfc[5];
 double q1,q2,mum,c,t11,t12,t21,t22,t33,muef,kapef;
 double fact;
 double fl5,rhs5,a05,a_5,al5;

 double gdftmp1,gdftmp2,gdftmp3,tmp1,tmp2,tmp3;


 for(k=0; k < 5; k++)
 {
	 a0[k]=0.; a1[k]=0.;  rhs[k]=0.; fl[k]=0.; dxfc[k]=0.; dyfc[k]=0.;
 }
 RES=0.;

 for(i=2; i < IT-2; i++)
  for(j=2; j < JT-2; j++)
   {
    n=IC(i,j);
    for(k=0; k < 5; k++)
     {
      dxfc[k]=0.25*(dxflr[k][Ilr(i,j)]+dxflr[k][Ilr(i+1,j)]+
                    dxfdu[k][Idu(i,j)]+dxfdu[k][Idu(i,j+1)]);
      dyfc[k]=0.25*(dyflr[k][Ilr(i,j)]+dyflr[k][Ilr(i+1,j)]+
                    dyfdu[k][Idu(i,j)]+dyfdu[k][Idu(i,j+1)]);
     }

	if(coord > 0.)
	{  
		/*  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  to check                                  
		if(Flag_turb == 1)
			{
				muef=mum=mu_ref*(pow((gdf[3][n]/gdf[0][n]),0.76) + 0.25 * (mu_turb_lr[0][IC(i,j)]+mu_turb_lr[1][IC(i,j)]+ 
                                                               mu_turb_du[0][IC(i,j)]+mu_turb_du[1][IC(i,j)]));
			}
		else {muef=mum=mu_ref*pow((gdf[3][n]/gdf[0][n]),0.76);}
		*/
		muef=mum=mu_ref*pow((gdf[3][n]/(gdf[0][n]*Rair*Tem_ref)),0.76);
		kapef=3.5*muef/0.72;
		t11=2.*muef*(2.*dxfc[1]-dyfc[2]-gdf[2][n]/yc[n])/3.;
		t12=muef*(dyfc[1]+dxfc[2]);
		t21=t12;
		t22=2.*muef*(2.*dyfc[2]-dxfc[1]-gdf[2][n]/yc[n])/3.;
		t33=2.*muef*(-dxfc[1]-dyfc[2]+2.*gdf[2][n]/yc[n])/3.;
		q1=kapef*dxfc[4];
		q2=kapef*dyfc[4];
		rhs[0]=-gdf[0][n]*gdf[2][n]/yc[n];
		rhs[1]=-(gdf[0][n]*gdf[1][n]*gdf[2][n]-t12)/yc[n];
		rhs[2]=-(gdf[0][n]*gdf[2][n]*gdf[2][n]-t22+t33)/yc[n];
		e=gdf[3][n]/(gdf[5][n]-1.)+0.5*gdf[0][n]*
			(gdf[1][n]*gdf[1][n]+gdf[2][n]*gdf[2][n]);
		rhs[3]=-(gdf[2][n]*(e+gdf[3][n])+q2-gdf[1][n]*t12-gdf[2][n]*t22)/yc[n];
	}

	muef=mum=mu_ref*pow((gdf[3][n]/(gdf[0][n]*Rair*Tem_ref)),0.76);
	c=sqrt(gdf[5][n]*gdf[3][n]/gdf[0][n]);
          
    a0[0]=gdf[0][n];
    a0[1]=gdf[0][n]*gdf[1][n];
    a0[2]=gdf[0][n]*gdf[2][n];
    e=gdf[3][n]/(gdf[5][n]-1.)+0.5*gdf[0][n]*(gdf[1][n]*gdf[1][n]+gdf[2][n]*gdf[2][n]);
    a0[3]=e;

    for(k=0; k < 4; k++)
    {
		fl[k]=(fllr[k][Ilr(i+1,j)]-fllr[k][Ilr(i,j)]+fldu[k][Idu(i,j+1)]-fldu[k][Idu(i,j)]);
		a1[k]=a0[k]+dt*(rhs[k]-fl[k]/vol[n]);
    }
    RES=max(RES,fabs(fl[0]));


    gdf[0][n]=a1[0];
    gdf[1][n]=a1[1]/a1[0];
    gdf[2][n]=a1[2]/a1[0];
    p=(gdf[5][n]-1.)*(a1[3]-0.5*(a1[1]*gdf[1][n]+a1[2]*gdf[2][n]));
    gdf[3][n]=p;
    gdf[4][n]=p/(gdf[0][n]*Rair);
    gdf[5][n]=GAM_cell(p,gdf[0][n]);
   
   }
}

//*****************************************************
void turb_lr(void) // Функция, заполняющая массив mu_turb_lr вязкостями;
{  
 int i, j;
 int Flag = 0;
 int JTH;
 double eps = 1.e-20;
 double mu_To, mu_Ti;
 double k = 0.4, k_caps = 0.0168, C_cp = 26., C_kleb = 0.3, C_wk = 0.25, A = 26.;
 double F_wake, Gamma, D;
 double F_max = 0., Y_max = 0, F_y, t_w, nu_;
 double U_max = -9.9e+9;
 double U_dif;
 double nu, omega, mu_m;
 double mu_To_max = 0., mu_Ti_max = 0., bug_min = 9.9e+9;
 double gdf_0, gdf_1, gdf_3, y_,YUP;

//JTH=(JT-3-2)/2;


 for(i=2; i < IT-2; i++)
 {

//  YUP=y[IN(i,JT-2)];

//DOWN => UP
	F_max = 0.;
	Y_max = 0.;
	U_max = 0.;
	for(j=2; j < JT-2; j++)
		{
		y_ = (yc[IC(i,j)] + yc[IC(i-1,j)]) / 2.;
		gdf_0 = (gdf[0][IC(i,j)] + gdf[0][IC(i-1,j)]) / 2.;
		gdf_1 = (gdf[1][IC(i,j)] + gdf[1][IC(i-1,j)]) / 2.;
		gdf_3 = (gdf[3][IC(i,j)] + gdf[3][IC(i-1,j)]) / 2.;

		mu_m=mu_ref*pow((gdf_3/(gdf_0*Rair*Tem_ref)),0.76);
		nu = mu_m / gdf_0;// динамическая вязкость
		if(j==2) //shouldnt this be i
			{
			t_w = nu * dyflr[1][Ilr(i,j)];// напряжение трения на стенке
			nu_ = sqrt(fabs(t_w / gdf_0));// динамическая скорость
			}
		omega = dxflr[2][Ilr(i,j)] - dyflr[1][Ilr(i,j)];
		F_y = y_ * fabs(omega) * (1 - exp((-y_ * nu_) / (A * nu))); // F(y)
		// Определяем Y_max, F_max
		if(F_max < F_y)
			{
			F_max = F_y;
			Y_max = y_;
			}
		// Определим U_dif = | U_max(x) - U_min(x) |
		U_max = max(U_max, gdf_1);
		}
	U_dif = fabs(U_max); //  - U_min
	F_wake = min(Y_max * F_max, C_wk * Y_max * pow(U_dif, 2.) / F_max );

	Flag=0;
	for(j=2; j < JT-2; j++)
		{
		y_ = (yc[IC(i,j)] + yc[IC(i-1,j)]) / 2.;
		gdf_0 = (gdf[0][IC(i,j)] + gdf[0][IC(i-1,j)]) / 2.;
		gdf_1 = (gdf[1][IC(i,j)] + gdf[1][IC(i-1,j)]) / 2.;
		gdf_3 = (gdf[3][IC(i,j)] + gdf[3][IC(i-1,j)]) / 2.;

		Gamma = pow((1 + 5.5 * pow((C_kleb * y_ / Y_max), 6.)), -1. );
		mu_m=mu_ref*pow((gdf_3/(gdf_0*Rair*Tem_ref)),0.76);
		nu = mu_m / gdf_0;// динамическая вязкость
		if(j==2) //shouldnt this be i
			{
			t_w = nu * dyflr[1][Ilr(i,j)];// напряжение трения на стенке
			nu_ = sqrt(fabs(t_w / gdf_0));// динамическая скорость
			}
		omega = dxflr[2][Ilr(i,j)] - dyflr[1][Ilr(i,j)];
		D = pow((1 - exp((-y_ * nu_) / (A * nu))), 2.);
 
		mu_Ti = pow((k * y_), 2) * D * fabs(omega); // Turb in
		mu_To = k_caps * C_cp * F_wake * Gamma; // Turb out

		//Отбор вязкости
		if(mu_To  > mu_Ti && Flag == 0)
			{
			mu_turb_lr[Ilr(i,j)] =  mu_Ti;
			}
		else
			{
			mu_turb_lr[Ilr(i,j)] = mu_To;
			Flag=1;
			}
		mu_turb_lr[Ilr(i,j)] = mu_turb_lr[Ilr(i,j)] * gdf_0;
		// printf("lr:  %d %d turb=%e \n",i,j,mu_turb_lr[Ilr(i,j)]);
		}

    }

}

//*****************************************************


void turb_du(void) // Функция, заполняющая массив mu_turb_du вязкостями;
{  
 int i, j;
 int Flag = 0;
 int JTH;

 double eps = 1.e-20;
 double mu_To, mu_Ti;
 double k = 0.4, k_caps = 0.0168, C_cp = 26., C_kleb = 0.3, C_wk = 0.25, A = 26.;
 double F_wake, Gamma, D;
 double F_max = 0, Y_max = 0, F_y, t_w, nu_;
 double U_max = -9.9e+9;
 double U_dif;
 double nu, omega, mu_m;
 double mu_To_max = 0., mu_Ti_max = 0., bug_min = 9.9e+9;
 double gdf_0, gdf_1, gdf_3, y_,YUP;

 //JTH=(JT-3-2+1)/2;

 for(i=2; i < IT-2; i++)
 {

//YUP=y[IN(i,JT-2)];

//DOWN => UP
	F_max = 0.;
    Y_max = 0.;
    U_max = 0.;
	for(j=2; j <=JT-2; j++)
		{
		y_ = y[IN(i,j)];
		gdf_0 = (gdf[0][IC(i,j)] + gdf[0][IC(i,j-1)]) / 2.;
		gdf_1 = (gdf[1][IC(i,j)] + gdf[1][IC(i,j-1)]) / 2.;
		gdf_3 = (gdf[3][IC(i,j)] + gdf[3][IC(i,j-1)]) / 2.;

		mu_m=mu_ref*pow((gdf_3/(gdf_0*Rair*Tem_ref)),0.76);
		nu = mu_m / gdf_0;// динамическая вязкость
		if(j==2)
			{
			t_w = nu * dyfdu[1][Idu(i,j)];// напряжение трения на стенке
			nu_ = sqrt(fabs(t_w / gdf_0));// динамическая скорость
			}
		omega = dxfdu[2][Idu(i,j)] - dyfdu[1][Idu(i,j)];
		F_y = y_ * fabs(omega) * (1. - exp((-y_ * nu_) / (A * nu))); // F(y)
		// Определяем Y_max, F_max 
		if(F_max < F_y)
			{
			F_max = F_y;
			Y_max = y_;
			}
		// Определим U_dif = | U_max(x) - U_min(x) |
		U_max = max(U_max, gdf_1);
		}

	U_dif = fabs(U_max); 
	if(Y_max < y[IN(i,3)])
		{
		Y_max = y[IN(i,3)];
		F_wake=0.;
		}
	else
		F_wake = min(Y_max * F_max, C_wk * Y_max * pow(U_dif, 2.) / F_max );

	Flag=0;

	for(j=2; j <= JT-2; j++)
		{
		y_ = y[IN(i,j)];
		gdf_0 = (gdf[0][IC(i,j)] + gdf[0][IC(i,j-1)]) / 2.;
		gdf_1 = (gdf[1][IC(i,j)] + gdf[1][IC(i,j-1)]) / 2.;
		gdf_3 = (gdf[3][IC(i,j)] + gdf[3][IC(i,j-1)]) / 2.;

		Gamma = pow((1 + 5.5 * pow((C_kleb * y_ / Y_max), 6.)), -1. );
		mu_m=mu_ref*pow((gdf_3/(gdf_0*Rair*Tem_ref)),0.76);
		nu = mu_m / gdf_0;// динамическая вязкость
		if(j==2)
			{
			t_w = nu * dyfdu[1][Idu(i,j)];// напряжение трения на стенке
			nu_ = sqrt(fabs(t_w / gdf_0));// динамическая скорость
			}
		omega = dxfdu[2][Idu(i,j)] - dyfdu[1][Idu(i,j)];
		D = pow((1 - exp((-y_ * nu_) / (A * nu))), 2);
   
		mu_Ti = pow((k * y_), 2) * D * fabs(omega); // Turb in
		mu_To = k_caps * C_cp * F_wake * Gamma; // Turb out

		//Отбор вязкости
		if(mu_Ti < mu_To && Flag == 0)
			{
			mu_turb_du[Idu(i,j)] =  mu_Ti;
			}
		else
			{
			mu_turb_du[Idu(i,j)] = mu_To;
			Flag=1;
			}                       

		mu_turb_du[Idu(i,j)] = mu_turb_du[Idu(i,j)] * gdf_0;
		}

  }
}

void TGAS4 (double P, double RHO, double *H)
{
    int JFLAG = 0, IFLAG = 0;
    double R0 = 1.292E00, P0 = 1.0133E05;
    double HLOW = 0, RSAVE = 0, YM = 0, YHIGH = 0, HHIGH = 0, YLOW = 0, Z1 = 0, GAMM = 0;
    double GAS1 = 0, GAS2 = 0, GAS3 = 0, GAS4 = 0, GAS5 = 0, GAS6 = 0, GAS7 = 0, GAS8 = 0, GAS9 = 0;
    double Y = log10(RHO/R0), X = log10(P/P0);
    if (fabs(Y+4.5E00) < 2.5E-02)
        goto p20;
    if (fabs(Y+0.5E00) < 5.0E-03)
        goto p50;
    IFLAG = -1;
    goto p90;
p10:	return;
p20:	IFLAG = 0;
    RSAVE = RHO;
    YM = Y;
    Y = -4.5E00+2.5E-02;
    YHIGH = Y;
    RHO = pow(10.,Y)*R0;
    JFLAG = -1;
    goto p90;
p30:	HHIGH = *H;
    Y = -4.5E00-2.5E-02;
    YLOW = Y;
    RHO = pow(10.,Y)*R0;
    JFLAG = 0;
    goto p90;
p40:    HLOW = *H;
    goto p80;
p50:	IFLAG = 1;
    RSAVE = RHO;
    YM = Y;
    Y = -0.5E00+0.5E-02;
    YHIGH = Y;
    RHO = pow(10.,Y)*R0;
    JFLAG = -1;
    goto p90;
p60:    HHIGH = *H;
    Y = -0.5E00-0.5E-02;
    YLOW = Y;
    RHO = pow(10.,Y)*R0;
    JFLAG = 0;
    goto p90;
p70:    HLOW = *H;
p80:	*H = HLOW+(HHIGH-HLOW)/(YHIGH-YLOW)*(YM-YLOW);
    RHO = RSAVE;
    goto p10;
p90:	Z1 = X-Y;
    if (Y > -0.5E00)
        goto p190;
    if (Y > -4.5E00)
        goto p140;
    if (Z1 > 0.1E00)
        goto p100;
    GAMM = 1.3986E00;
    goto p240;
p100:   if (Z1 > 0.85E00)
        goto p110;
    GAS1 = 2.53908E02+1.01491E02*Y;
    GAS2 = (-3.87199E02-1.54304E02*Y)*Z1;
    GAS3 = (7.28532E00-8.04378E00*Z1-1.82577E-03* Y)* Y*Y;
    GAS4 = (9.86233E01+4.63763E01*Y+2.18994E01*Z1)*Z1*Z1;
    GAS5 = -2.52423E02-1.01445E02*Y ;
    GAS6 = (3.87210E02+1.54298E02*Y)*Z1;
    GAS7 = (-7.2773E00+8.042277E00*Z1+2.28399E-03*Y)*Y*Y ;
    GAS8 = (-9.87576E01-4.63883E01*Y-2.19438E01*Z1)*Z1*Z1;
    GAS9 = exp(-11.E00+2.E00*Y+11.E00*Z1-2.E00*Y*Z1);
    GAMM=GAS1+GAS2+GAS3+GAS4+(GAS5+GAS6+GAS7+GAS8)/(1.-GAS9);
    goto p240;
p110:   if (Z1 > 1.30E00)
        goto p120;
    GAS1 = -1.05745E01-1.93693E00*Y;
    GAS2 = (3.07202E01+3.35578E00*Y)*Z1;
    GAS3 = (-7.79965E-02+6.68790E-02*Z1-9.86882E-04*Y)*Y*Y;
    GAS4 = (-2.60637E01-1.42391E00*Y+7.23223E00*Z1)*Z1*Z1;
    GAS5 = -1.86342E01+2.41997E-02*Y;
    GAS6 = (3.20880E01-7.46914E-01*Y)*Z1;
    GAS7 = (3.75161E-02-4.10125E-02*Z1+5.74637E-04*Y)*Y*Y;
    GAS8 = (-1.69985E01+5.39041E-01*Y+2.56253E00*Z1)*Z1*Z1;
    GAS9 = exp(2.768567E02+2.152383E01*Y-2.164837E02*Z1-1.394837E01*Y*Z1);
    goto p230;
p120:   if (Z1 > 1.95E00)
        goto p130;
    GAS1 = 6.17584E-01-2.40690E-01*Y;
    GAS2 = (1.95904E00+3.41644E-01*Y)*Z1;
    GAS3 = (-1.01073E-02+ 6.77631E-03*Z1-1.15922E-04*Y)*Y*Y;
    GAS4 = (-1.68951E00-1.10932E-01*Y+4.26058E-01*Z1)*Z1*Z1;
    GAS5 = -1.34222E01-5.43713E-01*Y;
    GAS6 = (1.81528E01+3.95928E-01*Y)*Z1;
    GAS7 = (-7.41105E-03+1.67768E-03*Z1-3.32714E-06*Y)*Y*Y;
    GAS8 = (-7.97425E00-5.80593E-02*Y+1.12448E00*Z1)*Z1*Z1;
    GAS9 = exp( 8.677803E01-8.370349E00*Y-4.074084E01*Z1+7.407405E00*Y*Z1);
    goto p230;
p130:   if (Z1 > 2.60E00)
        //fprintf (out, "RHO, P %f %f\n", RHO, P);
        //WRITE (*,*) RHO,P
        GAS1 = -8.32595E00-3.50219E-01*Y;
    GAS2 = (1.36455E01+3.59350E-01*Y)*Z1;
    GAS3 = (-3.70109E-03+3.30836E-03*Z1+1.10018E-04*Y)*Y*Y;
    GAS4 = (-6.49007E00-8.38594E-02*Y+1.02443E00*Z1)*Z1*Z1;
    GAS5 = -3.08441E01-1.49510E00*Y;
    GAS6 = (3.00585E01+9.19650E-01*Y)*Z1;
    GAS7 = (-3.60024E-02+1.02522E-02*Z1-4.68760E-04*Y)*Y*Y;
    GAS8 = (-9.33522E00-1.35228E-01*Y+8.92634E-01*Z1)*Z1*Z1;
    GAS9 = exp(8.800047E01-1.679356E01*Y-3.333353E01*Z1+8.465574E00*Y*Z1);
    goto p230;
p140:   if (Z1 > 0.1E00)
        goto p150;
    GAMM=1.4E00;
    goto p240;
p150:   if (Z1 > 0.95E00)
        goto p160;
    GAS1 = -1.33083E02-9.98707E00*Y;
    GAS2 = (3.94734E02+2.35810E01*Y)*Z1;
    GAS3 = (1.43957E00-1.43175E00*Z1+1.77068E-05*Y)*Y*Y;
    GAS4 = (-3.84712E02-1.36367E01*Y+1.24325E02*Z1)*Z1*Z1;
    GAS5 = 1.34486E02+9.99122E00*Y;
    GAS6 = (-3.94719E02-2.35853E01*Y)*Z1;
    GAS7 = (-1.43799E00+1.43039E00*Z1+1.44367E-04*Y)*Y*Y;
    GAS8 = (3.84616E02+1.36318E01*Y-1.24348E02*Z1)*Z1*Z1;
    GAS9 = exp(-2.141444E01+1.381584E00*Y+2.039473E01*Z1-1.315789E00*Y*Z1);
    GAMM = GAS1 + GAS2+GAS3 + GAS4+(GAS5 +GAS6 +GAS7+GAS8)/(1.-GAS9);
    goto p240;
p160:   if (Z1 > 1.50E00)
        goto p170;
    GAS1 = -7.36684E00-1.13247E00*Y;
    GAS2 = (2.47879E01 +1.99625E00*Y)*Z1;
    GAS3 = (-4.91630E-02+4.16673E-02*Z1-6.58149E-04*Y)*Y*Y;
    GAS4 = (-2.32990E01-8.59418E-01*Y+7.19016E00*Z1)*Z1*Z1;
    GAS5 = -2.42647E00+5.57912E-01*Y;
    GAS6 = (-2.03055E00-1.22031E00*Y)*Z1;
    GAS7 = (3.74866E-02-3.39278E-02*Z1+5.21042E-04*Y)*Y*Y;
    GAS8 = (7.75414E00+6.08488E-01*Y-3.68326E00*Z1)*Z1*Z1;
    GAS9 = exp(8.077385E01 -1.273807E01*Y - 6.547623E01*Z1+1.190475E01*Y*Z1);
    goto p230;
p170:   if (Z1 > 2.00E00)
        goto p180;
    GAS1 = 4.31520E-01-2.83857E-01*Y;
    GAS2 = (2.27791E00+3.99159E-01*Y)*Z1;
    GAS3 = (-1.29444E-02+8.78724E-03*Z1-1.60583E-04*Y)*Y*Y;
    GAS4 = (-1.84314E00-1.28136E-01*Y+4.45362E-01*Z1)*Z1*Z1;
    GAS5 = -1.03883E01-3.58718E-01*Y;
    GAS6 = (1.35068E01+1.87268E-01*Y)*Z1;
    GAS7 = (-4.28184E-03-9.52016E-04*Z1-4.10506E-05*Y)*Y*Y;
    GAS8 = (-5.63894E00-1.45626E-03*Y+7.39915E-01*Z1)*Z1*Z1;
    GAS9 = exp(2.949221E02+1.368660E01*Y-1.559335E02*Z1-3.787766E00*Y*Z1);
    goto p230;
p180:   if (Z1 > 2.5E00)
        //fprintf (out, "RHO, P %f %f\n", RHO, P);
        //WRITE (*,*) RHO,P
        GAS1 = -3.77766E00-5.53738E-01*Y;
    GAS2 = (6.60834E00+4.87181E-01*Y)*Z1;
    GAS3 = (-2.11045E-02+9.67277E-03*Z1-2.19420E-04*Y)*Y*Y;
    GAS4 = (-2.94754E00-1.02365E-01*Y+4.39620E-01*Z1)*Z1*Z1;
    GAS5 = 4.05813E01+3.25692E00*Y;
    GAS6 = (-4.79583E01-2.53660E00*Y)*Z1;
    GAS7 = (9.06436E-02-3.47578E-02*Z1+1.00077E-03*Y)*Y*Y;
    GAS8 = (1.89040E01+4.94114E-01*Y-2.48554E00*Z1)*Z1*Z1;
    GAS9 = exp(5.34718E02+7.495657E01*Y-2.219822E02*Z1-3.017229E01*Y*Z1);
    goto p230;
p190:   if (Z1 > 0.1E00)
        goto p200;
    GAMM = 1.4017E00;
    goto p240;
p200:   if (Z1 > 1.05E00)
        goto p210;
    GAS1 = -9.67488E01+2.05296E-01*Y;
    GAS2 = (2.69927E02-1.92887E00*Y)*Z1;
    GAS3 = (3.78392E-01-3.24965E-01*Z1-3.61036E-03*Y)*Y*Y;
    GAS4 = (-2.46711E02+1.54416E00*Y+7.48760E01*Z1)*Z1*Z1;
    GAS5 = 9.81502E01-2.05448E-01*Y;
    GAS6 = (-2.69913E02+1.93052E00*Y)*Z1;
    GAS7 = (-3.78527E-01+3.24832E-01*Z1+3.66182E-03*Y)*Y*Y;
    GAS8 = (2.46630E02-1.54646E00*Y-7.48980E01*Z1)*Z1*Z1;
    GAS9 = exp(-2.659865E01+1.564631E00*Y+2.312926E01*Z1-1.360543E00*Y*Z1);
    GAMM = GAS1+GAS2+GAS3+GAS4+(GAS5+GAS6+GAS7+GAS8)/(1.0-GAS9);
    goto p240;
p210:   if (Z1 > 1.60E00)
        goto p220;
    GAS1 = -2.67533E-01-1.87457E-01*Y;
    GAS2 = (5.07693E00+2.72286E-01*Y)*Z1;
    GAS3 = (1.04541E-02-1.42211E-02*Z1+6.38962E-4*Y)*Y*Y;
    GAS4 = (-5.08520E00-7.81935E-02*Y+1.58711E00*Z1)*Z1*Z1;
    GAS5 = 2.87969E00+ 3.9009E-01*Y;
    GAS6 = (-8.06179E00-5.51250E-01*Y)*Z1;
    GAS7 = (-1.01903E-02+1.35906E-02*Z1-8.97772E-04*Y)*Y*Y;
    GAS8 = (7.29592E00+1.83861E-01*Y-2.15153E00*Z1)*Z1*Z1;
    GAS9 = exp(1.828573E-02-3.428596E01*Y-1.51786E02*Z1+2.976212E01*Y*Z1);
    goto p230;
p220:   if (Z1 > 2.30E00)
        //printf ("RHO, P %f %f\n",  RHO, P);
        //WRITE (*,*) RHO,P
        GAS1 = 9.21537E-01-2.39670E-01*Y;
    GAS2 = (1.30714E00+3.42990E-01*Y)*Z1;
    GAS3 = (-2.18847E-02+1.36691E-02*Z1-4.90274E-04*Y)*Y*Y;
    GAS4 = (-1.20916E00-1.10206E-01*Y+3.087920E-01*Z1)*Z1*Z1;
    GAS5 = -6.77089E00-6.90476E-02*Y;
    GAS6 = (8.18168E00-9.52708E-02*Y)*Z1;
    GAS7 = (2.98487E-02-1.78706E-02*Z1+6.28419E-04*Y)*Y*Y;
    GAS8 = (-3.07662E00+6.60408E-02*Y+3.38590E-01*Z1)*Z1*Z1;
    GAS9 = exp(1.5916669E02+3.976192E01*Y-7.966199E01*Z1-1.66667E01*Y*Z1);
p230:	GAMM = GAS1+GAS2+GAS3+GAS4+(GAS5+GAS6+GAS7+GAS8)/(1.0+GAS9);
p240:	*H = GAMM/(GAMM-1.0E00)*P/RHO;
    if (IFLAG < 0)
        goto p10;
    else if (IFLAG == 0)
        goto p250;
    else
        goto p260;
p250:   if (JFLAG < 0)
        goto p30;
    else if (JFLAG == 0)
        goto p40;
    else
        goto p10;
p260:	if (JFLAG < 0)
        goto p60;
    else if (JFLAG == 0)
        goto p70;
    else
        goto p10;
}


double GAM_cell (double PIN, double DIN)
{
    //FILE *out = fopen ("output.txt", "w");
    double pa = P_ref, ra = D_ref;
    double H, GAMout = 0;

//!!!!!!!!!!!!!!!!!!!!!!!!!!For dim 11.11.23    
//    PIN = PIN * pa;
//    DIN = DIN * ra;
    TGAS4(PIN, DIN, &H);
    GAMout = H * DIN / (H * DIN - PIN);

    if (gam_ind != 0) return GAMout;
	else return 1.4;
}

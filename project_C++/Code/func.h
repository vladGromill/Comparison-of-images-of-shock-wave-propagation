void inp_mesh(void);
void inp_gdf(void);
void deriv_lr(void);
void deriv_du(void);
void step(void);
double flow_lr(void);
double flow_du(void);
void new_gdf(void);
void new_gdf_corr(void);

void process_fld(const char* name);
void turb_lr(void);
void turb_du(void);

void out(void);
void out_tec(char *name);
void memo_on(void);

void riemi(double *par,int *alarm, double *gam);

double GAM_cell (double PIN, double DIN);
void TGAS4 (double P, double RHO, double *H);


#define max(a,b)  (((a) > (b)) ? (a) : (b))
#define min(a,b)  (((a) < (b)) ? (a) : (b))

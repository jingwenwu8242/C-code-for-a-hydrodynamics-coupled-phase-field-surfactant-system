void initialization(double **psi, double **phi,double **u, double **v,double **p);
//*****************boundary condition****************
void augmenc(double **c, int nxt, int nyt);



//**************NS*******************//
void advection_step(double **u, double **v, double **c1, double **c2, double **adv_u, double **adv_v, double **adv_c1, double **adv_c2);
void augmenuv(double **u, double **v);
void augmentuv(double **u, double **v);
void sf_force(double **c1, double **mu, double **fx, double **fy);
//********** NS*****
void temp_uv(double **c1,double **c2,double **tu, double **tv, double **u, double **v, double **ou, double **ov, double **adv_u, double **adv_v, double **worku, double **workv, double **fx, double **fy);
void relax_uv(double **tu, double **tv, double **su, double **sv, int nx, int ny);
void Poisson(double **tu, double **tv, double **p, double **np);
void source_uv(double **tu, double **tv, double **p, double **divuv, int nxt, int nyt);
void div_uv(double **tu, double **tv, double **divuv, int nxt, int nyt);
void MG_Poisson(double **p, double **f);
void vcycle_uv(double **uf, double **ff, double wf, int nxf, int nyf, int ilevel);
void relax_p(double **p, double **f, double w, int ilevel, int nxt, int nyt);
void residual_den(double **r, double **u, double **f, double den, int nxt, int nyt);
void grad_p(double **p, double **dpdx, double **dpdy, int nxt, int nyt);
void prolong(double **u_coarse, double **u_fine, int nxt, int nyt);




/*****************function*****************************/
double functionMass(double **phi);
void functionWpsi(double **psi, double **phi,double **Wpsi);
void functionWphi(double **psi, double **phi, double **Wphi);
double functionQpsi(double **psi, double **dG, double **Wpsi);

void functionQphi(double **phi, double **psi ,double **dF, double **Wphi, double **Qphi);
void functiondF(double **phi,double **dF);
void functiondG(double **psi,double **dG);




/************  phase-field *******/
void cahnpsi(double **c_old, double **cc_old,double **adv_c,  double **c_new);

void sourcepsi(double **c_old, double **cc_old, double **adv_c, double **src_c, double **src_mu);
void cahnphi(double **c_old, double **cc_old,double **adv_c,   double **c_new);
void sourcephi(double **c_old, double **cc_old, double **adv_c, double **src_c, double **src_mu);

void laplace_ch(double **a, double **lap_a, int nxt, int nyt);

void vcycle(double **uf_new, double **wf_new, double **su, double **sw,double M,double S,double gam2, int nxf, int nyf, int ilevel);

void relax(double **c_new, double **mu_new, double **su, double **sw, double M,double S,double gam2,int ilevel, int nxt, int nyt);

void defect(double **duc, double **dwc, double **uf_new, double **wf_new,double **suf, double **swf,double M,double S, double gam2, int nxf, int nyf);

void nonL(double **ru, double **rw, double **c_new, double **mu_new,double M, double S, double gam2, int nxt, int nyt);

void restrict2(double **uf, double **uc, double **vf, double **vc, int nxt, int nyt);

void restrict1(double **vf, double **vc, int nxt, int nyt);

void prolong_ch(double **uc, double **uf, double **vc, double **vf, int nxt, int nyt);

double error(double **c_old, double **c_new, int nxt, int nyt);



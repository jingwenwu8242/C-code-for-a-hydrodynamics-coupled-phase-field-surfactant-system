/************************************************
cac psi,phi
shear flow
xy:periodic
************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "util.h"
#include "cacns.h"
#define NR_END 1

int it, max_it, nx, ny, n_level, c_relax, p_relax, count, Num;
double pi, xleft, xright, yleft, yright, volume, h, dt, gamphi, gampsi, Cahnphi, Cahnpsi, Mphi, Mpsi, alpha, theta, **Wpsi, **Wphi, Qpsi, **Qphi, **dFphi, **dGpsi,
    **mc, Sphi, Spsi, **ct, **sc, **smu, **mupsi, **muphi, **omupsi, **omuphi, **intmuphi, **intmupsi, **intpsi, **intphi, Masspsi, Massphi, **ou, **u, **nu, **ov, **v, **nv, **intu, **intv, **tu, **tv, **p, **np,
    **adv_u, **adv_v, **adv_phi, **adv_psi, **fx, **fy, **fxphi, **fyphi, **fxpsi, **fypsi, **worku, **workv, **workp, **nworku, **nworkv,
    Re, We, rho1, rho2, rho, gravity;
char bufferpsi[2000] = {0}, bufferphi[2000] = {0}, bufferu[2000] = {0}, bufferv[2000] = {0}, bufferp[2000] = {0};
int main()
{

    extern int it, max_it, nx, ny, n_level, c_relax, p_relax, count, Num;
    extern double pi, xleft, xright, yleft, yright, volume, h, dt, gamphi, gampsi, Cahnphi, Cahnpsi, Mphi, Mpsi, alpha, theta, **Wpsi, **Wphi, Qpsi, **Qphi, **dFphi, **dGpsi,
        **mc, Sphi, Spsi, **ct, **sc, **smu, **mupsi, **muphi, **omupsi, **omuphi, **intmuphi, **intmupsi, **intpsi, **intphi, Masspsi, Massphi, **ou, **u, **nu, **ov, **v, **nv, **intu, **intv, **tu, **tv, **p, **np,
        **adv_u, **adv_v, **adv_phi, **adv_psi, **fx, **fy, **fxphi, **fyphi, **fxpsi, **fypsi, **worku, **workv, **workp, **nworku, **nworkv,
        Re, We, rho1, rho2, rho, gravity;
    extern char bufferpsi[2000], bufferphi[2000], bufferu[2000], bufferv[2000], bufferp[2000];
    int i, j, ns;
    double **psi, **npsi, **opsi, **phi, **nphi, **ophi;

    FILE *fpsi, *fphi, *fmasspsi, *fmassphi, *fu, *fv, *fp;

    /*****/

    nx = gnx;
    ny = gny;
    n_level = (int)(log(nx) / log(2) - 0.9); 

    c_relax = 5; 
    p_relax = 7;

    pi = 4.0 * atan(1.0);

    xleft = 0.0, xright = 2.0 * pi;
    yleft = 0.0, yright = xright * (double)ny / (double)nx;
    volume = (xright - xleft) * (yright - yleft);

    count = 0;

    /***********************/

    h = xright / (double)nx;

    dt = 0.01;

   
    gamphi = 0.05; 
        gampsi = 0.01; 

    Cahnphi = pow(gamphi, 2); 
    Cahnpsi = pow(gampsi, 2);
    Mphi = 1;
    Mpsi = 1;
   alpha = 0.01;
    theta = 0.05;
    Sphi = 6;
    Spsi = 6;

    gravity = 0;
    rho1 = 1.0;
    rho2 = 1.0;
    rho = 1.0;
    Re = 1.0;
    We = 1.0;

    max_it = 1000;
    ns = max_it / 50;

    /***********************/

    printf("nx = %d, ny = %d\n", nx, ny);
    printf("n_level = %d\n", n_level);
    printf("dt      = %f\n", dt);
    printf("max_it  = %d\n", max_it);
    printf("ns      = %d\n", ns);

    intpsi = dmatrix(0, nx + 1, 0, ny + 1);
    intphi = dmatrix(0, nx + 1, 0, ny + 1);

    opsi = dmatrix(0, nx + 1, 0, ny + 1); 
    psi = dmatrix(0, nx + 1, 0, ny + 1);  
    npsi = dmatrix(0, nx + 1, 0, ny + 1); 

    ophi = dmatrix(0, nx + 1, 0, ny + 1); 
    phi = dmatrix(0, nx + 1, 0, ny + 1);  
    nphi = dmatrix(0, nx + 1, 0, ny + 1); 

    Wpsi = dmatrix(1, nx, 1, ny);
    Wphi = dmatrix(1, nx, 1, ny);
    Qphi = dmatrix(1, nx, 1, ny);
    dFphi = dmatrix(1, nx, 1, ny);
    dGpsi = dmatrix(1, nx, 1, ny);

    ct = dmatrix(1, nx, 1, ny); 
    sc = dmatrix(1, nx, 1, ny);
    smu = dmatrix(1, nx, 1, ny);            
    muphi = dmatrix(0, nx + 1, 0, ny + 1);  
    mupsi = dmatrix(0, nx + 1, 0, ny + 1);  
    omuphi = dmatrix(0, nx + 1, 0, ny + 1); 
    omupsi = dmatrix(0, nx + 1, 0, ny + 1); 

    intmuphi = dmatrix(0, nx + 1, 0, ny + 1);
    intmupsi = dmatrix(0, nx + 1, 0, ny + 1);

    
    
    ou = dmatrix(-1, nx + 1, 0, ny + 1); 
    u = dmatrix(-1, nx + 1, 0, ny + 1);  
    nu = dmatrix(-1, nx + 1, 0, ny + 1); 

    ov = dmatrix(0, nx + 1, -1, ny + 1); 
    v = dmatrix(0, nx + 1, -1, ny + 1);  
    nv = dmatrix(0, nx + 1, -1, ny + 1); 

    intu = dmatrix(-1, nx + 1, 0, ny + 1);
    intv = dmatrix(0, nx + 1, -1, ny + 1);

    p = dmatrix(0, nx + 1, 0, ny + 1);  
    np = dmatrix(0, nx + 1, 0, ny + 1); 

    adv_u = dmatrix(0, nx, 1, ny);
    adv_v = dmatrix(1, nx, 0, ny); 

    adv_psi = dmatrix(1, nx, 1, ny); 
    adv_phi = dmatrix(1, nx, 1, ny); 

    
    fx = dmatrix(0, nx, 1, ny);
    fy = dmatrix(1, nx, 0, ny);
    fxpsi = dmatrix(0, nx, 1, ny);
    fypsi = dmatrix(1, nx, 0, ny);
    fxphi = dmatrix(0, nx, 1, ny);
    fyphi = dmatrix(1, nx, 0, ny);

    
    tu = dmatrix(-1, nx + 1, 0, ny + 1);
    tv = dmatrix(0, nx + 1, -1, ny + 1);

    
    worku = dmatrix(0, nx, 1, ny);
    workv = dmatrix(1, nx, 0, ny);
    workp = dmatrix(1, nx, 1, ny); 
    nworku = dmatrix(0, nx, 1, ny);
    nworkv = dmatrix(1, nx, 0, ny);

    zero_matrix(mupsi, 1, nx, 1, ny);  
    zero_matrix(muphi, 1, nx, 1, ny);  
    zero_matrix(omupsi, 1, nx, 1, ny); 
    zero_matrix(omuphi, 1, nx, 1, ny); 

    sprintf(bufferpsi, "./data1/datapsi.m");
    sprintf(bufferphi, "./data1/dataphi.m");
    sprintf(bufferu, "./data1/datau.m");
    sprintf(bufferv, "./data1/datav.m");
    sprintf(bufferp, "./data1/datap.m");

    fpsi = fopen(bufferpsi, "w");
    fphi = fopen(bufferphi, "w");
    fu = fopen(bufferu, "w");
    fv = fopen(bufferv, "w");
    fp = fopen(bufferp, "w");

    fclose(fpsi);
    fclose(fphi);
    fclose(fu);
    fclose(fv);
    fclose(fp);

    initialization(psi, phi, u, v, p); 

    mat_copy(opsi, psi, 1, nx, 1, ny); 
    mat_copy(ophi, phi, 1, nx, 1, ny); 

    mat_copy(npsi, psi, 1, nx, 1, ny); 
    mat_copy(nphi, phi, 1, nx, 1, ny); 

    mat_copy(nu, u, 0, nx, 1, ny); 
    mat_copy(nv, v, 1, nx, 0, ny); 
    mat_copy(ou, u, 0, nx, 1, ny); 
    mat_copy(ov, v, 1, nx, 0, ny); 

    print_data(psi, phi, u, v, p); 
    Masspsi = functionMass(psi);
    Massphi = functionMass(phi);

    fmasspsi = fopen("./data1/masspsi.m", "w");
    fmassphi = fopen("./data1/massphi.m", "w");

    fprintf(fmasspsi, "%16.12f \n", Masspsi);
    fprintf(fmassphi, "%16.12f \n", Massphi);

    augmenc(psi, nx, ny);
    augmenc(phi, nx, ny);

    for (it =1; it <= max_it; it++)
    {

        i0jloop
        {
            intu[i][j] = 2.0 * u[i][j] - ou[i][j];
        }

        ij0loop
        {
            intv[i][j] = 2.0 * v[i][j] - ov[i][j];
        }

        ijloop
        {

            intpsi[i][j] = 2.0 * psi[i][j] - opsi[i][j];
            intphi[i][j] = 2.0 * phi[i][j] - ophi[i][j];
            intmuphi[i][j] = 2.0 * muphi[i][j] - omuphi[i][j];
            intmupsi[i][j] = 2.0 * mupsi[i][j] - omupsi[i][j];
        } 

        /*******************NSinitial********************/
        
       advection_step(intu, intv, intphi, intpsi, adv_u, adv_v, adv_phi, adv_psi);
        sf_force(intphi, intmuphi, fxphi, fyphi); 

        sf_force(intpsi,intmupsi, fxpsi, fypsi); 


        
        i0jloop
        {
            fx[i][j] = fxphi[i][j] / We + fxpsi[i][j] / We;
            
        }

        ij0loop
        {
            fy[i][j] = fyphi[i][j] / We + fypsi[i][j] / We;
            
        } 

        grad_p(p, worku, workv, nx, ny);                                                   
        temp_uv(intphi, intpsi, tu, tv, u, v, ou, ov, adv_u, adv_v, worku, workv, fx, fy); 
        augmentuv(tu, tv);
        Poisson(tu, tv, p, np); 

        grad_p(np, nworku, nworkv, nx, ny); 
        i0jloop
        {
            if (it==1)
            {
                nu[i][j] = tu[i][j] - dt * (nworku[i][j] - worku[i][j]) / rho;
            }
            else 
            {
                nu[i][j] = tu[i][j] - 2.0 * dt * (nworku[i][j] - worku[i][j]) / (3.0 * rho);
            }
        }

        ij0loop
        {
    if (it==1)
            {
                nv[i][j] = tv[i][j] - dt * (nworkv[i][j] - workv[i][j]) / rho;
            }
            else 
            {
                nv[i][j] = tv[i][j] - 2.0 * dt * (nworkv[i][j] - workv[i][j]) / (3.0 * rho);
            }
        } 
        // get u_{i+1/2,j},v_{i,j+1/2}//

        
  /*      ijloop{
              adv_phi[i][j]=0.0;
              adv_psi[i][j]=0.0;
          }   */

        functionWpsi(intpsi, intphi, Wpsi);
        functionWphi(intpsi, intphi, Wphi);
        functiondF(intphi, dFphi);
        functiondG(intpsi, dGpsi);
        Qpsi = functionQpsi(intpsi, dGpsi, Wpsi);

        functionQphi(intphi, intpsi, dFphi, Wphi, Qphi);

        
        cahnpsi(psi, opsi, adv_psi, npsi);

        
        cahnphi(phi, ophi, adv_phi, nphi);

        mat_copy(opsi, psi, 1, nx, 1, ny); 
        mat_copy(ophi, phi, 1, nx, 1, ny); 

        mat_copy(psi, npsi, 1, nx, 1, ny); 
        mat_copy(phi, nphi, 1, nx, 1, ny); 

        mat_copy(ou, u, 0, nx, 1, ny);
        mat_copy(u, nu, 0, nx, 1, ny);
        mat_copy(ov, v, 1, nx, 0, ny);
        mat_copy(v, nv, 1, nx, 0, ny);
        mat_copy(p, np, 1, nx, 1, ny);

        mat_copy(omuphi, muphi, 1, nx, 1, ny);
        mat_copy(omupsi, mupsi, 1, nx, 1, ny);
        printf("it = %d\n", it);

        if (it % ns == 0)
        {
            count++;
            Masspsi = functionMass(psi);
            Massphi = functionMass(phi);
            fprintf(fmasspsi, "%16.12f \n", Masspsi);
            fprintf(fmassphi, "%16.12f \n", Massphi);
            print_data(npsi, nphi, nu, nv, np); 
            printf("print out counts %d\n", count);
        }
    }

    return 0;
}
void initialization(double **psi, double **phi, double **u, double **v, double **p)
{
    extern double xright, yright, h, gamphi, pi;

    int i, j;
    double x, y, a, b;
    a = 2;
    b = 3;

    ijloop
    {
        x = ((double)i - 0.5) * h;
        y = ((double)j - 0.5) * h;

        psi[i][j] = 0.2 + 0.001 * (2.0 * rand() / (float)RAND_MAX - 1.0);
        phi[i][j] = 0.5 + 0.001 * (2.0 * rand() / (float)RAND_MAX - 1.0);
        //phi[i][j] = tanh(((0.4 * pi) - sqrt(pow(x - pi, 2) / a + pow(y - pi, 2) / b)) / (sqrt(2) * gamphi));
        p[i][j] = 0.0;
    }
    /*  i0jloop{
          x = ((double)i - 0.5) * h;
         y = ((double)j - 0.5) * h;
         u[i][j]=-0.7*(y-1);
     }
     ij0loop{
           x = ((double)i - 0.5) * h;
         y = ((double)j - 0.5) * h;
         v[i][j]=0;
     } */

    zero_matrix(u, 0, nx, 1, ny);
    zero_matrix(v, 1, nx, 0, ny);
}



void advection_step(double **u, double **v, double **c1, double **c2, double **adv_u, double **adv_v, double **adv_c1, double **adv_c2)
{
    extern int nx, ny;
    extern double h, **ou, **ov;
    double a, b, d, **ux, **uy, **vx, **vy;

    int i, j, k;
    

    ux = dmatrix(-1, nx, 1, ny);
    uy = dmatrix(0, nx, 0, ny);
    vx = dmatrix(0, nx, 0, ny);
    vy = dmatrix(1, nx, -1, ny);

    augmenuv(u, v);
    augmenc(c1, nx, ny);
    augmenc(c2, nx, ny);

    

    i0jloop
    {
        if (u[i][j] > 0.0)
            adv_u[i][j] = u[i][j] * (u[i][j] - u[i - 1][j]) / h;
        else
            adv_u[i][j] = u[i][j] * (u[i + 1][j] - u[i][j]) / h;

        if (v[i][j - 1] + v[i + 1][j - 1] + v[i][j] + v[i + 1][j] > 0.0)
            adv_u[i][j] += 0.25 * (v[i][j - 1] + v[i + 1][j - 1] + v[i][j] + v[i + 1][j]) * (u[i][j] - u[i][j - 1]) / h;
        else
            adv_u[i][j] += 0.25 * (v[i][j - 1] + v[i + 1][j - 1] + v[i][j] + v[i + 1][j]) * (u[i][j + 1] - u[i][j]) / h;
    }

    ij0loop
    {
        if (u[i - 1][j] + u[i][j] + u[i - 1][j + 1] + u[i][j + 1] > 0.0)
            adv_v[i][j] = 0.25 * (u[i - 1][j] + u[i][j] + u[i - 1][j + 1] + u[i][j + 1]) * (v[i][j] - v[i - 1][j]) / h;
        else
            adv_v[i][j] = 0.25 * (u[i - 1][j] + u[i][j] + u[i - 1][j + 1] + u[i][j + 1]) * (v[i + 1][j] - v[i][j]) / h;

        if (v[i][j] > 0.0)
            adv_v[i][j] += v[i][j] * (v[i][j] - v[i][j - 1]) / h;
        else
            adv_v[i][j] += v[i][j] * (v[i][j + 1] - v[i][j]) / h;
    }
    
    ijloop
    {
        adv_c1[i][j] = (u[i][j] * (c1[i + 1][j] + c1[i][j]) - u[i - 1][j] * (c1[i][j] + c1[i - 1][j])) / (2.0 * h) + (v[i][j] * (c1[i][j + 1] + c1[i][j]) - v[i][j - 1] * (c1[i][j] + c1[i][j - 1])) / (2.0 * h);

        adv_c2[i][j] = (u[i][j] * (c2[i + 1][j] + c2[i][j]) - u[i - 1][j] * (c2[i][j] + c2[i - 1][j])) / (2.0 * h) + (v[i][j] * (c2[i][j + 1] + c2[i][j]) - v[i][j - 1] * (c2[i][j] + c2[i][j - 1])) / (2.0 * h);
    }

    free_dmatrix(ux, -1, nx, 1, ny);
    free_dmatrix(uy, 0, nx, 0, ny);
    free_dmatrix(vx, 0, nx, 0, ny);
    free_dmatrix(vy, 1, nx, -1, ny);
}
void augmenuv(double **u, double **v)
{
    extern int nx, ny;

    int i, j;
    double aa;

    for (j = 1; j <= ny; j++)
    {
        aa = 0.5 * (u[0][j] + u[nx][j]);
        u[0][j] = u[nx][j] = aa;
        u[-1][j] = u[nx - 1][j];
        u[nx + 1][j] = u[1][j];
    }

    for (i = -1; i <= nx + 1; i++)
    {
        u[i][0] = -0.0 - u[i][1];
        u[i][ny + 1] = 0.0 - u[i][ny];
    }

    for (j = 0; j <= ny; j++)
    {
        v[0][j] = v[nx][j];
        v[nx + 1][j] = v[1][j];
    }

    for (i = 0; i <= nx + 1; i++)
    {
        v[i][0] = v[i][ny] = 0.0;
        v[i][-1] = -v[i][1];
        v[i][ny + 1] = -v[i][ny - 1];
    }
}

void augmentuv(double **u, double **v)
{
    extern int nx, ny;

    int i, j;
    double aa;

    for (j = 1; j <= ny; j++)
    {
        aa = 0.5 * (u[0][j] + u[nx][j]);
        u[0][j] = u[nx][j] = aa;
        u[-1][j] = u[nx - 1][j];
        u[nx + 1][j] = u[1][j];
    }

    for (i = -1; i <= nx + 1; i++)
    {
        u[i][0] = -0.0 - u[i][1];
        u[i][ny + 1] = 0.0 - u[i][ny];
    }

    for (j = 0; j <= ny; j++)
    {
        v[0][j] = v[nx][j];
        v[nx + 1][j] = v[1][j];
    }

    for (i = 0; i <= nx + 1; i++)
    {
        v[i][0] = v[i][ny] = 0.0;
        v[i][-1] = -v[i][1];
        v[i][ny + 1] = -v[i][ny - 1];
    }
}
void sf_force(double **c1, double **mu, double **fx, double **fy) 
{
    extern int nx, ny;
    extern double h;

    int i, j;

    augmenc(c1, nx, ny);
    augmenc(mu, nx, ny);

    i0jloop
    {
        fx[i][j] = 0.5 * (c1[i + 1][j] + c1[i][j]) * (mu[i + 1][j] - mu[i][j]) / (h);
    }

    ij0loop
    {
        fy[i][j] = 0.5 * (c1[i][j + 1] + c1[i][j]) * (mu[i][j + 1] - mu[i][j]) / (h);
    }
}


void temp_uv(double **c1, double **c2, double **tu, double **tv, double **u, double **v, double **ou, double **ov, double **adv_u, double **adv_v, double **worku, double **workv, double **fx, double **fy)
{
   int i, j, it_mg = 1, max_it = 500;
    extern double h, dt, rho, rho1, rho2, gravity,Re;
    extern int it;
    double residu = 1.0, residv = 1.0, tol = 1.0e-5, **soru, **sorv, **su, **sv;
  soru = dmatrix(0, nx, 1, ny);
    sorv = dmatrix(1, nx, 0, ny);

    su = dmatrix(0, nx, 1, ny);
    sv = dmatrix(1, nx, 0, ny);

    i0jloop
    {
        if (it==1){
            su[i][j] = u[i][j] / dt - adv_u[i][j] - worku[i][j] / rho - fx[i][j] / rho;
        }
        else{
            su[i][j] = (4.0 * u[i][j] - ou[i][j]) / (2.0 * dt) - adv_u[i][j] - worku[i][j] / rho - fx[i][j] / rho;
        }
        
    }  

    ij0loop
    {
        if (it==1){
            sv[i][j] = v[i][j] / dt - adv_v[i][j] - workv[i][j] / rho - fy[i][j] / rho - gravity * (rho1 * 0.5 * (c1[i][j + 1] + c1[i][j]) + rho2 * 0.5 * (c2[i][j + 1] + c2[i][j]) - rho) / rho;
        }
        else{
            sv[i][j] = (4.0 * v[i][j] - ov[i][j]) / (2.0 * dt) - adv_v[i][j] - workv[i][j] / rho - fy[i][j] / rho - gravity * (rho1 * 0.5 * (c1[i][j + 1] + c1[i][j]) + rho2 * 0.5 * (c2[i][j + 1] + c2[i][j]) - rho) / rho;
        }
        
    }
    while (it_mg <= max_it && residu >= tol && residv >= tol)
    {

        relax_uv(tu, tv, su, sv, nx, ny);

        mat_sub(soru, soru, tu, 0, nx, 1, ny); 
        mat_sub(sorv, sorv, tv, 1, nx, 0, ny);
        residu = mat_max(soru, 0, nx, 1, ny);
        residv = mat_max(sorv, 1, nx, 0, ny);
        mat_copy(soru, tu, 0, nx, 1, ny);
        mat_copy(sorv, tv, 1, nx, 0, ny);

        it_mg++;
    }

    printf("Velocity1 iteration = %d ", it_mg - 1);

    free_dmatrix(soru, 0, nx, 1, ny);
    free_dmatrix(sorv, 1, nx, 0, ny);

    free_dmatrix(su, 0, nx, 1, ny);
    free_dmatrix(sv, 1, nx, 0, ny);
}
void relax_uv(double **tu, double **tv, double **su, double **sv, int nx, int ny)
{
    extern double xright, dt, Re, h, rho;
    extern int it, max_it;

    int i, j, iter;
    double h2, sorc, coef;

    h2 = pow(h, 2);

    for (iter = 1; iter <= 5; iter++)
    {

        augmentuv(tu, tv);

        
        i0jloop
        {

            sorc = su[i][j] + (tu[i + 1][j] + tu[i - 1][j] + tu[i][j + 1] + tu[i][j - 1]) / (Re * h2 * rho);
            if (it==1)
            {
                coef = 1.0 / dt + 4.0 / (Re * h2 * rho);
            }
            else 
            {
                coef = 3.0 / (2.0 * dt) + 4.0 / (Re * h2 * rho);
            }

            tu[i][j] = sorc / coef;
        }

        ij0loop
        {

            sorc = sv[i][j] + (tv[i + 1][j] + tv[i - 1][j] + tv[i][j + 1] + tv[i][j - 1]) / (Re * h2 * rho);
             if (it==1)
            {
                coef = 1.0 / dt + 4.0 / (Re * h2 * rho);
            }
            else 
            {
                coef = 3.0 / (2.0 * dt) + 4.0 / (Re * h2 * rho);
            }

            tv[i][j] = sorc / coef;
        }
    }
}



void Poisson(double **tu, double **tv, double **p, double **np)
{
    extern int nx, ny;
    extern double **workp;

    source_uv(tu, tv, p, workp, nx, ny);

    MG_Poisson(np, workp);
}

void source_uv(double **tu, double **tv, double **p, double **divuv, int nxt, int nyt) 
{
    extern double dt, h, lam, rho;
    extern int it, max_it;

    int i, j;

    div_uv(tu, tv, divuv, nxt, nyt); 

    augmenc(p, nxt, nyt);

    ijloopt
    {
       if (it==1)
        {
            divuv[i][j] = divuv[i][j] / dt + (p[i + 1][j] + p[i - 1][j] + p[i][j + 1] + p[i][j - 1] - 4.0 * p[i][j]) / (h * h * rho);
        }
        else 
        {
           divuv[i][j] = 3.0 * divuv[i][j] / (2.0 * dt) + (p[i + 1][j] + p[i - 1][j] + p[i][j + 1] + p[i][j - 1] - 4.0 * p[i][j]) / (h * h * rho);
        }
    }
}

void div_uv(double **tu, double **tv, double **divuv, int nxt, int nyt) 
{
    extern double xright;

    int i, j;
    double ht;

    ht = xright / (double)nxt;
    
    

    ijloopt
    {
        divuv[i][j] = (tu[i][j] - tu[i - 1][j] + tv[i][j] - tv[i][j - 1]) / ht;
    }
}
void MG_Poisson(double **p, double **f)
{
    extern int nx, ny;
    extern double rho;

    int it_mg = 1, max_it = 100;
    double resid = 1.0, resid2 = 10.0, tol = 1.0e-5, **sor;

    sor = dmatrix(1, nx, 1, ny);

    mat_copy(sor, p, 1, nx, 1, ny); 

    while (it_mg <= max_it && resid >= tol)
    {

        vcycle_uv(p, f, rho, nx, ny, 1);

        pressure_update(p); 

        mat_sub(sor, sor, p, 1, nx, 1, ny);
        resid = mat_max(sor, 1, nx, 1, ny);
        mat_copy(sor, p, 1, nx, 1, ny);

        if (resid > resid2)
            it_mg = max_it;
        else
            resid2 = resid;

        it_mg++;
    }
    printf("Mac pressure iteration = %d   residual = %16.14f\n", it_mg - 1, resid);

    free_dmatrix(sor, 1, nx, 1, ny);

    return;
}

void vcycle_uv(double **uf, double **ff, double wf, int nxf, int nyf, int ilevel)
{
    extern int n_level;

    relax_p(uf, ff, wf, ilevel, nxf, nyf); 

    if (ilevel < n_level)
    {

        int nxc, nyc;
        double **rf, **fc, **uc;

        nxc = nxf / 2, nyc = nyf / 2; 

        rf = dmatrix(1, nxf, 1, nyf); 
        fc = dmatrix(1, nxc, 1, nyc); 
        uc = dmatrix(1, nxc, 1, nyc); 

        residual_den(rf, uf, ff, wf, nxf, nyf); 

        restrict1(rf, fc, nxc, nyc); 

        zero_matrix(uc, 1, nxc, 1, nyc);

        vcycle_uv(uc, fc, wf, nxc, nyc, ilevel + 1); 

        prolong(uc, rf, nxc, nyc); 

        mat_add(uf, uf, rf, 1, nxf, 1, nyf); 

        relax_p(uf, ff, wf, ilevel, nxf, nyf); 

        free_dmatrix(rf, 1, nxf, 1, nyf);
        free_dmatrix(fc, 1, nxc, 1, nyc);
        free_dmatrix(uc, 1, nxc, 1, nyc);
    }
}
void relax_p(double **p, double **f, double w, int ilevel, int nxt, int nyt)
{
    extern int ny, p_relax;
    extern double xright, Fr;

    int i, j, iter;
    double ht, ht2, a[4], sorc, coef;

    ht = xright / (double)nxt;
    ht2 = pow(ht, 2);

    for (iter = 1; iter <= p_relax; iter++)
    {

        ijloopt
        {
            a[0] = 1.0 / w;
            a[1] = 1.0 / w;
            a[2] = 1.0 / w;
            a[3] = 1.0 / w;

            sorc = f[i][j];
            coef = -(a[0] + a[1] + a[2] + a[3]) / ht2;
            
            if (i == 1)
            {
                sorc -= (a[0] * p[i + 1][j] + a[1] * p[nxt][j]) / ht2;
            }
            else if (i == nxt)
            {
                sorc -= (a[0] * p[1][j] + a[1] * p[i - 1][j]) / ht2;
            }
            else
            {
                sorc -= (a[0] * p[i + 1][j] + a[1] * p[i - 1][j]) / ht2;
            }

            /
            if (j == 1)
            {
                sorc -= (a[2] * p[i][j + 1] + a[3] * p[i][nyt]) / ht2;
            }
            else if (j == nyt)
            {
                sorc -= (a[2] * p[i][1] + a[3] * p[i][j - 1]) / ht2;
            }
            else
            {
                sorc -= (a[2] * p[i][j + 1] + a[3] * p[i][j - 1]) / ht2;
            }

            p[i][j] = sorc / coef;
        }
    }
}

void residual_den(double **r, double **u, double **f, double den, int nxt, int nyt)
{
    int i, j;
    double **dpdx, **dpdy;

    dpdx = dmatrix(0, nxt, 1, nyt);
    dpdy = dmatrix(1, nxt, 0, nyt);

    grad_p(u, dpdx, dpdy, nxt, nyt); 

    i0jloopt
    {
        dpdx[i][j] = dpdx[i][j] / den;
    }

    ij0loopt
    {
        dpdy[i][j] = dpdy[i][j] / den;
    } //(dpdx,dpdy)=1/rho*grad_p//

    div_uv(dpdx, dpdy, r, nxt, nyt);  
    mat_sub(r, f, r, 1, nxt, 1, nyt); 

    free_dmatrix(dpdx, 0, nxt, 1, nyt);
    free_dmatrix(dpdy, 1, nxt, 0, nyt);
}

void grad_p(double **p, double **dpdx, double **dpdy, int nxt, int nyt) 
{
    extern double xright;

    int i, j;
    double ht;

    ht = xright / (double)nxt;
    
    i0jloopt
    {
        if (i == 0)
            dpdx[0][j] = (p[1][j] - p[nxt][j]) / ht;
        else if (i == nxt)
            dpdx[nxt][j] = (p[1][j] - p[nxt][j]) / ht;
        else
            dpdx[i][j] = (p[i + 1][j] - p[i][j]) / ht;
    }

    /
    ij0loopt
    {
        if (j == 0)
        {
            dpdy[i][j] = (p[i][1] - p[i][nyt]) / ht;
        }
        else if (j == nyt)
        {
            dpdy[i][j] = (p[i][1] - p[i][nyt]) / ht;
        }

        else
            dpdy[i][j] = (p[i][j + 1] - p[i][j]) / ht;
    }
}

void prolong(double **u_coarse, double **u_fine, int nxt, int nyt)
{
    int i, j;

    ijloopt
    {
        u_fine[2 * i - 1][2 * j - 1] = u_fine[2 * i - 1][2 * j] =
            u_fine[2 * i][2 * j - 1] = u_fine[2 * i][2 * j] = u_coarse[i][j];
    }
}


void augmenc(double **c, int nxt, int nyt) 
{
    int i, j;

    for (j = 1; j <= nyt; j++)
    {
        c[0][j] = c[nxt][j];
        c[nxt + 1][j] = c[1][j];
    }

    for (i = 0; i <= nxt + 1; i++)
    {
        c[i][0] = c[i][nyt];
        c[i][nyt + 1] = c[i][1];
    }
}

/***************  phase field  *****************/

void cahnpsi(double **c_old, double **cc_old, double **adv_c, double **c_new)
{
    extern int nx, ny;
    extern double **ct, **sc, **smu, **mupsi, **dGpsi, **Wpsi, Qpsi, Mpsi, Spsi, Cahnpsi;

    int it_mg = 1, max_it_CH = 200;
    double resid = 1.0, tol = 1.0e-6;

    mat_copy(ct, c_old, 1, nx, 1, ny); 

    sourcepsi(c_old, cc_old, adv_c, sc, smu); 

    while (it_mg <= max_it_CH && resid > tol)
    {

        vcycle(c_new, mupsi, sc, smu, Mpsi, Spsi, Cahnpsi, nx, ny, 1);
        resid = error(ct, c_new, nx, ny);
        mat_copy(ct, c_new, 1, nx, 1, ny);

        it_mg++;
    }
    printf("cahnpsi %16.14f   %d\n", resid, it_mg - 1);
}

void sourcepsi(double **c_old, double **cc_old, double **adv_c, double **src_c, double **src_mu) 
{
    extern int nx, ny, it, max_it;
    extern double dt, h, Spsi, Cahnpsi, **dGpsi, **Wpsi, Qpsi, Mpsi;
    int i, j;

    ijloop
    {
        if (it==1)
        {
            src_c[i][j] = c_old[i][j] / dt + Mpsi * Qpsi - adv_c[i][j];
            src_mu[i][j] = 1 / Cahnpsi * dGpsi[i][j] + Wpsi[i][j] - Spsi / Cahnpsi * c_old[i][j];
        }
        else
        {
            src_c[i][j] = (4.0 * c_old[i][j] - cc_old[i][j]) / (2.0 * dt) + Mpsi * Qpsi - adv_c[i][j];
            src_mu[i][j] = 1 / Cahnpsi * dGpsi[i][j] + Wpsi[i][j] - Spsi / Cahnpsi * (2.0 * c_old[i][j] - cc_old[i][j]);
        }
    }
}
void cahnphi(double **c_old, double **cc_old, double **adv_c, double **c_new)
{
    extern int nx, ny;
    extern double **ct, **sc, **smu, **muphi, **dFphi, **Wphi, **Qphi, Mphi, Sphi, Cahnphi;

    int it_mg = 1, max_it_CH = 200;
    double resid = 1.0, tol = 1.0e-6;

    mat_copy(ct, c_old, 1, nx, 1, ny); 

    sourcephi(c_old, cc_old, adv_c, sc, smu); 

    while (it_mg <= max_it_CH && resid > tol)
    {

        vcycle(c_new, muphi, sc, smu, Mphi, Sphi, Cahnphi, nx, ny, 1);
        resid = error(ct, c_new, nx, ny);
        mat_copy(ct, c_new, 1, nx, 1, ny);

        it_mg++;
    }
    printf("cahnphi %16.14f   %d\n", resid, it_mg - 1);
}

void sourcephi(double **c_old, double **cc_old, double **adv_c, double **src_c, double **src_mu) 
{
    extern int nx, ny, it, max_it;
    extern double dt, h, Sphi, Cahnphi, **dFphi, **Wphi, **Qphi, Mphi;

    int i, j;

    ijloop
    {
        if (it==1)
        {
            src_c[i][j] = c_old[i][j] / dt + Mphi * Qphi[i][j] - adv_c[i][j];
            src_mu[i][j] = 1 / Cahnphi * dFphi[i][j] + Wphi[i][j] - Sphi / Cahnphi * c_old[i][j];
        }
        else
        {
            src_c[i][j] = (4.0 * c_old[i][j] - cc_old[i][j]) / (2.0 * dt) + Mphi * Qphi[i][j] - adv_c[i][j];
            src_mu[i][j] = 1 / Cahnphi * dFphi[i][j] + Wphi[i][j] - Sphi / Cahnphi * (2.0 * c_old[i][j] - cc_old[i][j]);
        }
    }
}

void vcycle(double **uf_new, double **wf_new, double **su, double **sw, double M, double S, double gam2, int nxf, int nyf, int ilevel)
{
    extern int n_level;

    relax(uf_new, wf_new, su, sw, M, S, gam2, ilevel, nxf, nyf); 

    if (ilevel < n_level)
    {

        int nxc, nyc;
        double **uc_new, **wc_new, **duc, **dwc,
            **uc_def, **wc_def, **uf_def, **wf_def;

        nxc = nxf / 2, nyc = nyf / 2;

        uc_new = dmatrix(1, nxc, 1, nyc);
        wc_new = dmatrix(1, nxc, 1, nyc);
        duc = dmatrix(1, nxc, 1, nyc);
        dwc = dmatrix(1, nxc, 1, nyc); 
        uc_def = dmatrix(1, nxc, 1, nyc);
        wc_def = dmatrix(1, nxc, 1, nyc);
        uf_def = dmatrix(1, nxf, 1, nyf);
        wf_def = dmatrix(1, nxf, 1, nyf);

        defect(duc, dwc, uf_new, wf_new, su, sw, M, S, gam2, nxf, nyf); 

        zero_matrix(uc_def, 1, nxc, 1, nyc);
        zero_matrix(wc_def, 1, nxc, 1, nyc);

        vcycle(uc_def, wc_def, duc, dwc, M, S, gam2, nxc, nyc, ilevel + 1); 

        prolong_ch(uc_def, uf_def, wc_def, wf_def, nxc, nyc); 

        mat_add2(uf_new, uf_new, uf_def, wf_new, wf_new, wf_def, 1, nxf, 1, nyf); 

        relax(uf_new, wf_new, su, sw, M, S, gam2, ilevel, nxf, nyf); 

        free_dmatrix(uc_new, 1, nxc, 1, nyc);
        free_dmatrix(wc_new, 1, nxc, 1, nyc);
        free_dmatrix(duc, 1, nxc, 1, nyc);
        free_dmatrix(dwc, 1, nxc, 1, nyc);
        free_dmatrix(uc_def, 1, nxc, 1, nyc);
        free_dmatrix(wc_def, 1, nxc, 1, nyc);
        free_dmatrix(uf_def, 1, nxf, 1, nyf);
        free_dmatrix(wf_def, 1, nxf, 1, nyf);
    }
}

void relax(double **c_new, double **mu_new, double **su, double **sw, double M, double S, double gam2, int ilevel, int nxt, int nyt)
{
    extern int c_relax, it, max_it;
    extern double xright, dt;

    int i, j, iter;
    double ht2, a[4], f[2], det;

    ht2 = pow(xright / (double)nxt, 2);

    for (iter = 1; iter <= c_relax; iter++)
    {

        ijloopt
        {
            if (it==1)
            {
                a[0] = 1.0 / dt;
            }
            else
            {
                a[0] = 3.0 / (2.0 * dt);
            }

            a[1] = M;

            a[2] = -4.0 / ht2 - S / gam2;

            a[3] = 1.0;

            f[0] = su[i][j];

            f[1] = sw[i][j];
            if (i > 1)
                f[1] -= c_new[i - 1][j] / ht2;
            else
                f[1] -= c_new[nxt][j] / ht2;

            if (i < nxt)
                f[1] -= c_new[i + 1][j] / ht2;
            else
                f[1] -= c_new[1][j] / ht2;

            if (j > 1)
                f[1] -= c_new[i][j - 1] / ht2;

            else
                f[1] -= c_new[i][nyt] / ht2;

            if (j < nyt)
                f[1] -= c_new[i][j + 1] / ht2;
            else
                f[1] -= c_new[i][1] / ht2;

            det = a[0] * a[3] - a[1] * a[2];

            c_new[i][j] = (a[3] * f[0] - a[1] * f[1]) / det;
            mu_new[i][j] = (-a[2] * f[0] + a[0] * f[1]) / det;
        }
    }
}

void defect(double **duc, double **dwc, double **uf_new, double **wf_new,
            double **suf, double **swf, double M, double S, double gam2, int nxf, int nyf) 
{
    double **ruf, **rwf, **ruc, **rwc, **rruf, **rrwf;

    ruf = dmatrix(1, nxf, 1, nyf);
    rwf = dmatrix(1, nxf, 1, nyf);
    rruf = dmatrix(1, nxf / 2, 1, nyf / 2);
    rrwf = dmatrix(1, nxf / 2, 1, nyf / 2);

    nonL(ruf, rwf, uf_new, wf_new, M, S, gam2, nxf, nyf); 

    mat_sub2(ruf, suf, ruf, rwf, swf, rwf, 1, nxf, 1, nyf); 

    restrict2(ruf, duc, rwf, dwc, nxf / 2, nyf / 2); 

    free_dmatrix(ruf, 1, nxf, 1, nyf);
    free_dmatrix(rwf, 1, nxf, 1, nyf);
    free_dmatrix(rruf, 1, nxf / 2, 1, nyf / 2);
    free_dmatrix(rrwf, 1, nxf / 2, 1, nyf / 2);
}

void nonL(double **ru, double **rw, double **c_new, double **mu_new, double M, double S, double gam2, int nxt, int nyt)
{
    extern double xright, dt;
    extern int it, max_it;

    int i, j;
    double ss, ht2, **lap_c, **lap_mu;

    lap_c = dmatrix(1, nxt, 1, nyt);

    laplace_ch(c_new, lap_c, nxt, nyt);
    ht2 = pow(xright / (double)nxt, 2);

    ijloopt
    {
        if (it==1)
        {
            ru[i][j] = c_new[i][j] / dt + M * mu_new[i][j];
        }
        else
        {
            ru[i][j] = 3.0 * c_new[i][j] / (2.0 * dt) + M * mu_new[i][j];
        }

        rw[i][j] = mu_new[i][j] + lap_c[i][j] - S / gam2 * c_new[i][j];
    }

    free_dmatrix(lap_c, 1, nxt, 1, nyt);
}

void laplace_ch(double **a, double **lap_a, int nxt, int nyt) 
{
    extern double xright;

    int i, j;
    double ht2, dadx_L, dadx_R, dady_B, dady_T;

    ht2 = pow(xright / (double)nxt, 2);

    ijloopt
    {

        if (i > 1)
            dadx_L = a[i][j] - a[i - 1][j];
        else
            dadx_L = a[i][j] - a[nxt][j];

        if (i < nxt)
            dadx_R = a[i + 1][j] - a[i][j];
        else
            dadx_R = a[1][j] - a[i][j];

        if (j > 1)
            dady_B = a[i][j] - a[i][j - 1];
        else
            dady_B = a[i][j] - a[i][nyt];

        if (j < nyt)
            dady_T = a[i][j + 1] - a[i][j];
        else
            dady_T = a[i][1] - a[i][j];

        lap_a[i][j] = (dadx_R - dadx_L + dady_T - dady_B) / ht2;
    }
}

void restrict2(double **uf, double **uc, double **vf, double **vc, int nxt, int nyt)
{
    int i, j;

    ijloopt
    {
        uc[i][j] = 0.25 * (uf[2 * i - 1][2 * j - 1] + uf[2 * i - 1][2 * j] + uf[2 * i][2 * j - 1] + uf[2 * i][2 * j]);
        vc[i][j] = 0.25 * (vf[2 * i - 1][2 * j - 1] + vf[2 * i - 1][2 * j] + vf[2 * i][2 * j - 1] + vf[2 * i][2 * j]);
    }
}

void restrict1(double **vf, double **vc, int nxt, int nyt)
{
    int i, j;

    ijloopt
    {
        vc[i][j] = 0.25 * (vf[2 * i - 1][2 * j - 1] + vf[2 * i - 1][2 * j] + vf[2 * i][2 * j - 1] + vf[2 * i][2 * j]);
    }
}

void prolong_ch(double **uc, double **uf, double **vc, double **vf, int nxt, int nyt)
{
    int i, j;

    ijloopt
    {
        uf[2 * i - 1][2 * j - 1] = uf[2 * i - 1][2 * j] = uf[2 * i][2 * j - 1] = uf[2 * i][2 * j] = uc[i][j];
        vf[2 * i - 1][2 * j - 1] = vf[2 * i - 1][2 * j] = vf[2 * i][2 * j - 1] = vf[2 * i][2 * j] = vc[i][j];
    }
}

double error(double **c_old, double **c_new, int nxt, int nyt)
{
    double **r, res;

    r = dmatrix(1, nxt, 1, nyt);

    mat_sub(r, c_new, c_old, 1, nxt, 1, nyt);
    res = mat_max(r, 1, nxt, 1, nyt);

    free_dmatrix(r, 1, nxt, 1, nyt);

    return res;
}

/*************** util ****************/
double **dmatrix(long nrl, long nrh, long ncl, long nch)
{
    double **m;
    long i, nrow = nrh - nrl + 1 + NR_END, ncol = nch - ncl + 1 + NR_END;

    m = (double **)malloc((nrow) * sizeof(double *));
    memset(m, 0, (nrow) * sizeof(double *));
    m += NR_END;
    m -= nrl;

    m[nrl] = (double *)malloc((nrow * ncol) * sizeof(double));
    memset(m[nrl], 0, (nrow * ncol) * sizeof(double));
    m[nrl] += NR_END;
    m[nrl] -= ncl;

    for (i = nrl + 1; i <= nrh; i++)
        m[i] = m[i - 1] + ncol;

    return m;
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
{
    free(m[nrl] + ncl - NR_END);
    free(m + nrl - NR_END);

    return;
}

void zero_matrix(double **a, int xl, int xr, int yl, int yr)
{
    int i, j;

    for (i = xl; i <= xr; i++)
    {
        for (j = yl; j <= yr; j++)
        {
            a[i][j] = 0.0;
        }
    }

    return;
}

void mat_copy(double **a, double **b, int xl, int xr, int yl, int yr)
{
    int i, j;

    for (i = xl; i <= xr; i++)
    {
        for (j = yl; j <= yr; j++)
        {
            a[i][j] = b[i][j];
        }
    }

    return;
}

void mat_copy2(double **a, double **b, double **a2, double **b2,
               int xl, int xr, int yl, int yr)
{
    int i, j;

    for (i = xl; i <= xr; i++)
    {
        for (j = yl; j <= yr; j++)
        {
            a[i][j] = b[i][j];
            a2[i][j] = b2[i][j];
        }
    }

    return;
}

void mat_add(double **a, double **b, double **psi,
             int xl, int xr, int yl, int yr)
{
    int i, j;

    for (i = xl; i <= xr; i++)
    {
        for (j = yl; j <= yr; j++)
        {
            a[i][j] = b[i][j] + psi[i][j];
        }
    }

    return;
}

void mat_add2(double **a, double **b, double **psi,
              double **a2, double **b2, double **phi,
              int xl, int xr, int yl, int yr)
{
    int i, j;

    for (i = xl; i <= xr; i++)
    {
        for (j = yl; j <= yr; j++)
        {
            a[i][j] = b[i][j] + psi[i][j];
            a2[i][j] = b2[i][j] + phi[i][j];
        }
    }

    return;
}

void mat_sub(double **a, double **b, double **psi,
             int nrl, int nrh, int ncl, int nch)
{
    int i, j;

    for (i = nrl; i <= nrh; i++)
    {
        for (j = ncl; j <= nch; j++)
        {
            a[i][j] = b[i][j] - psi[i][j];
        }
    }

    return;
}

void mat_sub2(double **a, double **b, double **psi,
              double **a2, double **b2, double **phi,
              int nrl, int nrh, int ncl, int nch)
{
    int i, j;

    for (i = nrl; i <= nrh; i++)
    {
        for (j = ncl; j <= nch; j++)
        {
            a[i][j] = b[i][j] - psi[i][j];
            a2[i][j] = b2[i][j] - phi[i][j];
        }
    }

    return;
}

double mat_max(double **a, int nrl, int nrh, int ncl, int nch)
{
    int i, j;
    double x = 0.0;

    for (i = nrl; i <= nrh; i++)
    {
        for (j = ncl; j <= nch; j++)
        {
            if (fabs(a[i][j]) > x)
                x = fabs(a[i][j]);
        }
    }

    return x;
}

void print_mat(FILE *fptr, double **a, int nrl, int nrh, int ncl, int nch)
{
    int i, j;

    for (i = nrl; i <= nrh; i++)
    {
        for (j = ncl; j <= nch; j++)
            fprintf(fptr, "   %16.14f", a[i][j]);

        fprintf(fptr, "\n");
    }

    return;
}

void print_data(double **psi, double **phi, double **u, double **v, double **p)
{
    extern char bufferpsi[2000], bufferphi[2000], bufferu[2000], bufferv[2000], bufferp[2000];
    int i, j;
    FILE *fpsi, *fphi, *fu, *fv, *fp;

    fpsi = fopen(bufferpsi, "a");
    fphi = fopen(bufferphi, "a");
    fu = fopen(bufferu, "a");
    fv = fopen(bufferv, "a");
    fp = fopen(bufferp, "a");

    iloop
    {
        jloop
        {

            fprintf(fpsi, "  %16.14f", psi[i][j]);
            fprintf(fphi, "  %16.14f", phi[i][j]);
            fprintf(fu, "  %16.14f", 0.5 * (u[i][j] + u[i - 1][j]));
            fprintf(fv, "  %16.14f", 0.5 * (v[i][j] + v[i][j - 1]));
            fprintf(fp, "  %16.14f", p[i][j]);
        }

        fprintf(fu, "\n");
        fprintf(fv, "\n");
        fprintf(fp, "\n");
        fprintf(fpsi, "\n");
        fprintf(fphi, "\n");
    }

    fclose(fpsi);
    fclose(fphi);
    fclose(fu);
    fclose(fv);
    fclose(fp);

    return;
}

void pressure_update(double **a)
{
    extern int nx, ny;

    int i, j;
    double ave = 0.0;

    ijloop
    {
        ave = ave + a[i][j];
    }
    ave /= (nx + 0.0) * (ny + 0.0);

    ijloop
    {
        a[i][j] -= ave;
    }

    return;
}

/*****************function*****************************/

double functionMass(double **phi)
{
    extern double h;
    int i, j;
    double res;
    res = 0;
    ijloop
    {
        res = res + phi[i][j];
    }
    res = pow(h, 2) * res;
    return res;
}

void functionWpsi(double **psi, double **phi, double **Wpsi)
{
    extern double theta, alpha, h;
    extern int nx, ny;
    int i, j;
    double **gadx, **gady;
    gadx = dmatrix(1, nx, 1, ny);
    gady = dmatrix(1, nx, 1, ny);

    augmenc(phi, nx, ny);
    ijloop { gadx[i][j] = (phi[i + 1][j] - phi[i - 1][j]) / (2.0 * h); }

    ijloop { gady[i][j] = (phi[i][j + 1] - phi[i][j - 1]) / (2.0 * h); }

    ijloop
    {

        Wpsi[i][j] = -theta * (pow(gadx[i][j], 2) + pow(gady[i][j], 2)) 
        + 0.5 * alpha * pow(phi[i][j], 2);
    }
    free_dmatrix(gadx, 1, nx, 1, ny);
    free_dmatrix(gady, 1, nx, 1, ny);
}

void functionWphi(double **psi, double **phi, double **Wphi)
{
    extern double theta, alpha, h;
    extern int nx, ny;
    int i, j;
    double a0, a1, a2, a3;
    augmenc(phi, nx, ny);
    augmenc(psi, nx, ny);
    ijloop
    {
        a0 = 0.5 * (psi[i + 1][j] + psi[i][j]);
        a1 = 0.5 * (psi[i - 1][j] + psi[i][j]);
        a2 = 0.5 * (psi[i][j + 1] + psi[i][j]);
        a3 = 0.5 * (psi[i][j - 1] + psi[i][j]);

        Wphi[i][j] = 2 * theta * ((a0 * phi[i + 1][j] + a1 * phi[i - 1][j] + a2 * phi[i][j + 1] + a3 * phi[i][j - 1]) / pow(h, 2) - (a0 + a1 + a2 + a3) * phi[i][j] / pow(h, 2)) 
        + alpha * psi[i][j] * phi[i][j];
    }
}
double functionQpsi(double **psi, double **dGpsi, double **Wpsi)
{
    extern double volume, Cahnpsi, h;
    int i, j;
    double res;

    res = 0.0;

    ijloop
    {
        res = res + 1.0 / Cahnpsi * dGpsi[i][j] + Wpsi[i][j];
    }

    res = res * pow(h, 2) / volume;

    return res;
}

void functionQphi(double **phi, double **psi, double **dFphi, double **Wphi, double **Qphi)
{
    extern double volume, Cahnphi, h, alpha;
    double res, beta, numerator, denominator;
    int i, j;

    res = 0;
    numerator = 0;
    denominator = 0;
    ijloop
    {
        numerator = numerator + dFphi[i][j];
        denominator = denominator + fabs(pow(phi[i][j], 2) - 1);

        

        res = res + alpha * psi[i][j] * phi[i][j];
    }
    beta = numerator * pow(h, 2) / (Cahnphi * denominator * pow(h, 2));
    res = res * pow(h, 2) / volume;
    ijloop
    {
        Qphi[i][j] = res + beta * fabs(pow(phi[i][j], 2) - 1);
    }
}
void functiondF(double **phi, double **dFphi)
{
    int i, j;
    ijloop
    {
        dFphi[i][j] = phi[i][j] * (pow(phi[i][j], 2) - 1);
    }
}
void functiondG(double **psi, double **dGpsi)
{
    int i, j;
    double kxi;

    ijloop
    {
        
        
        dGpsi[i][j] = log(psi[i][j] / (1 - psi[i][j]));
    }
    /* kxi = 0.00001;
           ijloop{

         if(psi[i][j] >= 1.0 - kxi)
         dGpsi[i][j] = log(psi[i][j]) - (1.0-psi[i][j])/kxi + 1.0 + log(kxi);
         else if (psi[i][j] > kxi && psi[i][j] < 1.0-kxi)
         dGpsi[i][j] = log(psi[i][j]) - log(1.0-psi[i][j]);
         else
         dGpsi[i][j] = -log(1.0-psi[i][j]) - 1.0 + psi[i][j]/kxi + log(kxi);

       }*/
}

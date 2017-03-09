//Created by Michael Grabchak and Lijuan Cao on 3/6/17.

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

double epsOne=.01; //Intervals near alpha=0 or alpha=1 that gets converted to 0 or 1
double tol = 1e-10;//tolerance
gsl_integration_qawo_table *tc;
gsl_integration_workspace *w;
gsl_integration_workspace *cw;
gsl_integration_workspace *wi;

void integrationSetup(){
    tc = gsl_integration_qawo_table_alloc(0, 1, GSL_INTEG_COSINE, 100);
    w = gsl_integration_workspace_alloc (10000);
    cw = gsl_integration_workspace_alloc (10000);
    gsl_set_error_handler_off();
}
void integrationCleanup(){
    gsl_integration_qawo_table_free(tc);
    gsl_integration_workspace_free (w);
    gsl_integration_workspace_free (cw);
}
//Evaluate the characteristic function
double cTSCharFunc(double t, void *params){
    
    double *par = (double *)params;
    
    double alpha = *par;
    double c = *(par+1);
    double ell = *(par+2);
    double s = t*ell;
    
    if(fabs(alpha-1)<epsOne)//if alpha is essentially 1
        return pow(1+s*s,c)*exp(-2*c*s*atan(s));
    else if(alpha < epsOne) //if alpha is less than 0
        return pow(1+s*s,-c);
    else //if alpha is greater than 0
        return exp( 2*c*gsl_sf_gamma(-alpha)*(pow(1+s*s,alpha/2) * cos(alpha*atan(s))-1));
}

//Evaluate the characteristic function divided by t, needed to evaluate the cdf
double  cTSCharFuncDiv(double t, void *params){
    
    return  cTSCharFunc(t, params)/t;
}

double integZero(double u, void *params)
{
    double *par = (double *)params;
    
    double omega = *par;
    double c = *(par+1);
    double ell = *(par+2);
    
    return exp(-2*u/ell)*pow(fabs(omega)*u+u*u, c-1);
}

double dCTSImp(double x, double mu, double *params){
    
    double result =0, result2 =0, error;
    double omega = x- mu;
    if(fabs(omega) < 1e-12)
        omega = 0;

    gsl_integration_qawo_table_set(tc,omega,1,GSL_INTEG_COSINE);
    
    gsl_function G;
    G.function = &cTSCharFunc;
    G.params = params;
    gsl_integration_qawo(&G, 0, tol, tol, 1000, w, tc, &result2, &error);
    gsl_integration_qawf (&G, 1, tol, 1000, w, cw, tc, &result, &error);
    
    return (result+result2)/M_PI;
}

//alpha = 0 second approach works for all c>0
double dCTSImp0(double x, double mu, double *params){
    
    double result =0, error;
    double omega = x-mu;
    if(fabs(omega) < 1e-12)
        omega = 0;
    
    double pars[3];
    double c = *(params+1);
    double ell = *(params+2);
    pars[0] = omega;
    pars[1] = c;
    pars[2] = ell;
    
    gsl_function G;
    G.function = &integZero;
    G.params = pars;
    gsl_integration_qagiu(&G, 0, tol, tol, 1000, wi, &result, &error);
    
    return result*exp(-fabs(omega)/ell)*pow(gsl_sf_gamma(c),-2)*pow(ell, -2*c);
}

//Evaluate pdf
void dCTS(double *x, int *len, double *mu, double *alpha, double *c, double *ell, double *result){
    
    double params[3];
    params[0] = *alpha;
    params[1] = *c;
    params[2] = *ell;
    
    integrationSetup();
    
    for(int i = 0; i < *len; i++) {
        
        result[i] = dCTSImp(x[i], *mu, params);
    }
    integrationCleanup();
}

void dCTS0(double *x, int *len, double *mu, double *alpha, double *c, double *ell, double *result){
    
    double params[3];
    params[0] = *alpha;
    params[1] = *c;
    params[2] = *ell;
    wi = gsl_integration_workspace_alloc (10000);
    for(int i = 0; i < *len; i++) {
        
        if(*c <= 0.5 && fabs(x[i] - *mu) < 0.00001)
            result[i] = 0;
        else
            result[i] = dCTSImp0(x[i], *mu, params);
    }
     gsl_integration_workspace_free (wi);
}

double pCTSImp(double x, double mu, double *params){
    
    double result =0, result2 = 0, error;
    double omega = x-mu;
    if(fabs(omega) < 1e-12)
        omega = 0;
    
    if(omega == 0)
        return 0.5;
    else {
            gsl_integration_qawo_table_set(tc,omega,1,GSL_INTEG_SINE);
            gsl_function G;
            G.function = & cTSCharFuncDiv;
            G.params = params;
        
            gsl_integration_qawo(&G, 0, tol, tol, 1000, w, tc, &result2, &error);
            gsl_integration_qawf (&G, 1, tol, 1000, w, cw, tc, &result, &error);
            return 0.5 + (result+result2)/M_PI;
        
    }
}

//Evaluate the cdf
void pCTS(double *x, int *len, double *mu, double *alpha, double *c, double *ell, double *result){
    
    double params[3];
    params[0] = *alpha;
    params[1] = *c;
    params[2] = *ell;
    
    integrationSetup();
    
    for(int i = 0; i < *len; i++) {
        
        result[i] = pCTSImp(x[i], *mu, params);
    }
    integrationCleanup();

}

double quantCTS(double y, double mu, double *params){
    double feps=1e-13;
    double peps=1e-13;
    double xLo = 0;
    double yLo = pCTSImp(xLo, mu, params);
    if(fabs(y-yLo)<feps)
        return xLo;
    double xHi = xLo, yHi = yLo;
    
    double alpha = *params;
    double c = *(params+1);
    double ell = *(params+2);
    
    while(y>yHi){
        yLo = yHi;
        xLo = xHi;
        xHi = xHi + 2*ell*sqrt(2*gsl_sf_gamma(2-alpha)*c);
        yHi = pCTSImp(xHi, mu, params);
    }
    
    double xMi = (xHi+xLo)/2;
    double yMi;
    int r = 0;
    while(fabs(yHi-yLo)>feps && fabs(xHi-xLo)>peps){
        r++;
        xMi = (xHi+xLo)/2;
        yMi = pCTSImp(xMi, mu, params);
        if(fabs(yMi-yLo)<feps){
            break;
        }
        else if(yMi < y){
            xLo = xMi;
            yLo = yMi;
        }
        else{
            xHi = xMi;
            yHi = yMi;
        }
        
    }
    return xMi;
}

double qCTSImp(double y,  double mu, double *params){
    
    if(y < 0.5){
        return -1*quantCTS(1-y, 0, params) + mu;
    }
    else
        return quantCTS(y, 0, params)+ mu;
    
}

//Evaluate the quantile function
void qCTS(double *y, int *len, double *mu, double *alpha, double *c, double *ell, double *result){
    
    double params[3];
    params[0] = *alpha;
    params[1] = *c;
    params[2] = *ell;
    
    integrationSetup();
    
    for(int i = 0; i < *len; i++) {
        result[i] = qCTSImp(y[i], *mu, params);
    }
    integrationCleanup();
}

//Evaluate the characteristic function
double saSCharFunc(double t, void *params){
    
    double *par = (double *)params;
    
    double alpha = *par;
    double c = *(par+1);
    
    return exp(-c*pow(t,alpha));
}

//Evaluate the characteristic function divided by t, needed to evaluate the cdf
double saSCharFuncDiv(double t, void *params){
    
    return saSCharFunc(t, params)/t;
}

double dSaSImp(double x, double mu, double *params){
    
    double result =0, result2 = 0,error;
    double omega = x-mu;
    if(fabs(omega) < 1e-12)
        omega = 0;
    
    gsl_integration_qawo_table_set(tc,omega,1,GSL_INTEG_COSINE);
    
    gsl_function G;
    G.function = &saSCharFunc;
    G.params = params;
    gsl_integration_qawo(&G, 0, tol, tol, 1000, w, tc, &result2, &error);
    gsl_integration_qawf (&G, 1, tol, 1000, w, cw, tc, &result, &error);
    
    return (result+result2)/M_PI;
}

//Evaluate pdf
void dSaS(double *x, int *len, double *mu, double *alpha, double *c,double *result){
    
    double params[2];
    params[0] = *alpha;
    params[1] = *c;
    integrationSetup();
    for(int i = 0; i < *len; i++) {
        
        result[i] = dSaSImp(x[i], *mu, params);
    }
    integrationCleanup();
}

double pSaSImp(double x, double mu, double *params){
    
    double result =0, result2 = 0, error;
    double omega = x-mu;
    if(fabs(omega) < 1e-12)
        omega = 0;
    
    if(omega == 0)
        return 0.5;
    else {
        gsl_integration_qawo_table_set(tc,omega,1,GSL_INTEG_SINE);
        
        gsl_function G;
        G.function = &saSCharFuncDiv;
        G.params = params;
        
        gsl_integration_qawo(&G, 0, tol, tol, 1000, w, tc, &result2, &error);
        gsl_integration_qawf (&G, 1, tol, 1000, w, cw, tc, &result, &error);
        return 0.5 + (result+result2)/M_PI;
    }
}

//Evaluate the cdf
void pSaS(double *x, int *len, double *mu, double *alpha, double *c,double *result){
    double params[2];
    params[0] = *alpha;
    params[1] = *c;
    integrationSetup();
    for(int i = 0; i < *len; i++) {
        
        result[i] = pSaSImp(x[i], *mu, params);
    }
    integrationCleanup();
}

double quantSaS(double y, double mu, double *params){
    double feps=1e-10;
    double peps=1e-10;
    double xLo = 0;
    double yLo = pSaSImp(xLo, mu, params);
    if(fabs(y-yLo)<feps)
        return xLo;
    double xHi = xLo, yHi = yLo;
    
    double newPar[2];
    newPar[0] = *params;
    newPar[1]=1;
    
    while(y>yHi){
        yLo = yHi;
        xLo = xHi;
        xHi = xHi + 2;
        yHi = pSaSImp(xHi, mu, newPar);
    }
    
    double xMi = (xHi+xLo)/2;
    double yMi;
    int r = 0;
    while(fabs(yHi-yLo)>feps && fabs(xHi-xLo)>peps){
        r++;
        xMi = (xHi+xLo)/2;
        yMi = pSaSImp(xMi, mu, newPar);
        if(fabs(yMi-yLo)<feps){
            break;
        }
        else if(yMi < y){
            xLo = xMi;
            yLo = yMi;
        }
        else{
            xHi = xMi;
            yHi = yMi;
        }
        
    }
    return xMi;
}

double qSaSImp(double y,  double mu, double *params){
    double alpha = *params;
    double c = *(params+1);
    
    if(y < 0.5){
        return (-pow(c,1/alpha))*quantSaS(1-y, 0, params) + mu;
    }
    else
        return pow(c,1/alpha)*quantSaS(y, 0, params) + mu;
    
}

//Evaluate the quantile function
void qSaS(double *y, int *len, double *mu, double *alpha, double *c,double *result){
    
    double params[2];
    params[0] = *alpha;
    params[1] = *c;
    integrationSetup();
    for(int i = 0; i < *len; i++) {
        result[i] = qSaSImp(y[i], *mu, params);
    }
    integrationCleanup();
}

double integ(double y, void *par){
    
    double *par2 = (double *)par;
    
    double alpha = *par2;
    double ell = *(par2+1);
    double t = *(par2+2);
    
    return pow(t*cos(y) + sin(y), -alpha)*pow(1+tan(y)/t, -ell)*cos((2-alpha)*y);
}

//Evaluate the characteristic function
double powCharFunc(double t, void *params){
    
    double *par = (double *)params;
    
    double alpha = *par;
    double c = *(par+1);
    double ell = *(par+2);
    double params2[3];
    params2[0] = alpha;
    params2[1] = ell;
    params2[2] = t;
    gsl_function G;
    G.function = &integ;
    G.params = params2;
    double result =0, error;
    
    gsl_integration_qag (&G, 0, M_PI/2, 0, tol, tol, 1000, wi, &result, &error);
    return exp (-result*2*c*gsl_sf_gamma(2-alpha)*pow(t,alpha+1));
}

//Evaluate the characteristic function divided by t, needed to evaluate the cdf
double powCharFuncDiv(double t, void *params){
    
    return powCharFunc(t, params)/t;
}

double dPowTSImp(double x, double mu, double *params){
    
    double result =0, result2 =0, error;
    double omega = x-mu;
    if(fabs(omega) < 1e-12)
        omega = 0;
    
    gsl_integration_qawo_table_set(tc,omega,1,GSL_INTEG_COSINE);
    
    gsl_function G;
    G.function = &powCharFunc;
    G.params = params;
    
    gsl_integration_qawo(&G, 0, tol, tol, 1000, w, tc, &result2, &error);
    gsl_integration_qawf (&G, 1, tol, 1000, w, cw, tc, &result, &error);
    
    return (result+result2)/M_PI;
}

//Evaluate pdf
void dPowTS(double *x, int *len, double *mu, double *alpha, double *c, double *ell, double *result){
    double params[3];
    params[0] = *alpha;
    params[1] = *c;
    params[2] = *ell;
    
    integrationSetup();
    wi = gsl_integration_workspace_alloc (10000);
    for(int i = 0; i < *len; i++) {
        
        result[i] = dPowTSImp(x[i], *mu, params);
    }
    integrationCleanup();
    gsl_integration_workspace_free (wi);
}

double pPowTSImp(double x, double mu, double *params){
    
    double result =0, result2 = 0, error;
    double omega = x-mu;
    if(fabs(omega) < 1e-12)
        omega = 0;
    
    if(omega == 0)
        return 0.5;
    else {
        gsl_integration_qawo_table_set(tc,omega,1,GSL_INTEG_SINE);
        
        gsl_function G;
        G.function = &powCharFuncDiv;
        G.params = params;
        
        gsl_integration_qawo(&G, 0, tol, tol, 1000, w, tc, &result2, &error);
        gsl_integration_qawf (&G, 1, tol, 1000, w, cw, tc, &result, &error);
        return 0.5 + (result+result2)/M_PI;
    }
}

//Evaluate the cdf
void pPowTS(double *x, int *len, double *mu, double *alpha, double *c, double *ell, double *result){
    double params[3];
    params[0] = *alpha;
    params[1] = *c;
    params[2] = *ell;
    
    integrationSetup();
    wi = gsl_integration_workspace_alloc (10000);
    for(int i = 0; i < *len; i++) {
        
        result[i] = pPowTSImp(x[i], *mu, params);
    }
    integrationCleanup();
    gsl_integration_workspace_free (wi);
}

double quantPowTS(double y, double mu, double *params){
    double feps=1e-13;
    double peps=1e-13;
    double xLo = 0;
    double yLo = pPowTSImp(xLo, mu, params);
    if(fabs(y-yLo)<feps)
        return xLo;
    double xHi = xLo, yHi = yLo;
    
    while(y>yHi){
        yLo = yHi;
        xLo = xHi;
        xHi = xHi + 2;
        yHi = pPowTSImp(xHi, mu, params);
    }
    
    double xMi = (xHi+xLo)/2;
    double yMi;
    int r = 0;
    while(fabs(yHi-yLo)>feps && fabs(xHi-xLo)>peps){
        r++;
        xMi = (xHi+xLo)/2;
        yMi = pPowTSImp(xMi, mu, params);
        if(fabs(yMi-yLo)<feps){
            break;
        }
        else if(yMi < y){
            xLo = xMi;
            yLo = yMi;
        }
        else{
            xHi = xMi;
            yHi = yMi;
        }
        
    }
    return xMi;
}

double qPowTSImp(double y,  double mu, double *params){
    
    if(y < 0.5){
        return -1*quantPowTS(1-y, 0, params) + mu;
    }
    else
        return quantPowTS(y, 0, params)+ mu;
    
}

//Evaluate the quantile function
void qPowTS(double *y, int *len, double *mu, double *alpha, double *c, double *ell, double *result){
    double params[3];
    params[0] = *alpha;
    params[1] = *c;
    params[2] = *ell;
    integrationSetup();
    wi = gsl_integration_workspace_alloc (10000);
    for(int i = 0; i < *len; i++) {
        result[i] = qPowTSImp(y[i], *mu, params);
    }
    integrationCleanup();
    gsl_integration_workspace_free (wi);
}
/* .C calls */

extern void dCTS(double *, int *, double *, double *, double *, double *, double *);
extern void dCTS0(double *, int *, double *, double *, double *, double *, double *);
extern void dPowTS(double *, int *, double *, double *, double *, double *, double *);
extern void dSaS(double *, int *, double *, double *, double *, double *);
extern void pCTS(double *, int *, double *, double *, double *, double *, double *);
extern void pPowTS(double *, int *, double *, double *, double *, double *, double *);
extern void pSaS(double *, int *, double *, double *, double *, double *);
extern void qCTS(double *, int *, double *, double *, double *, double *, double *);
extern void qPowTS(double *, int *, double *, double *, double *, double *, double *);
extern void qSaS(double *, int *, double *, double *, double *, double *);

static const R_CMethodDef CEntries[] = {
    {"dCTS",   (DL_FUNC) &dCTS,   7},
    {"dCTS0",  (DL_FUNC) &dCTS0,  7},
    {"dPowTS", (DL_FUNC) &dPowTS, 7},
    {"dSaS",   (DL_FUNC) &dSaS,   6},
    {"pCTS",   (DL_FUNC) &pCTS,   7},
    {"pPowTS", (DL_FUNC) &pPowTS, 7},
    {"pSaS",   (DL_FUNC) &pSaS,   6},
    {"qCTS",   (DL_FUNC) &qCTS,   7},
    {"qPowTS", (DL_FUNC) &qPowTS, 7},
    {"qSaS",   (DL_FUNC) &qSaS,   6},
    {NULL, NULL, 0}
};

void R_init_SymTS(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);

 }

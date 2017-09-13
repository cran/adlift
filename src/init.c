#include <stdlib.h>
#include <R_ext/Rdynload.h>

/* Symbol registration table */
#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

void amatdual(int *steps,int *po,int *re,int *nbrs,double *weights,double *alpha,int *lpo,int *lre,int *nn,double *adual);
void nbrsfromlca(double *lca,int *nc,int *rowno,int *nbrs);
void invtnp(double *X, double *coeff, double *lengths, double *lengthsremove, int *pointsin, double *lca,int *nadd,int *N, int *lr, int *nc,int *outpo,double *outlen);
void fwtnp(double *input,double *f,int *nkeep,int *intercept,int *initboundhandl,int *neighbours,int *closest,int *LocalPred,int *n,double *coeff,double *lengthsremove,double *lengths,double *lca,int *pointsin,int *nc,int *doW,double *W,int *varonly,double *v);
void pointsupdate(double *X,double *coeff,int *nn,int *index,int *remove,int *pointsin,double *wts,double *l,int *N,double *alpha,int *r);

static const R_CMethodDef R_CDef[] = {
   CALLDEF(amatdual, 10),
   CALLDEF(nbrsfromlca, 4),
   CALLDEF(invtnp, 12),
   CALLDEF(fwtnp, 19),
   CALLDEF(pointsupdate, 11),
   {NULL, NULL, 0}
};

void R_init_adlift(DllInfo *dll)
{
    R_registerRoutines(dll, R_CDef, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}

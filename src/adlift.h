
#ifndef ___ADLIFT_H___
#define ___ADLIFT_H___

void adaptneigh(int *pointsin, double *X, double *coeff, int *nbrs, int *remove, int *inter, int *nn, double *weights, int *scheme, int *clo, int *index, int *neighbours, int *N, int *docoeff);
void adaptpred(int *pointsin, double *X, double *coeff, int *nbrs, int *remove, int *inter, int *nn, double *weights, int *scheme, int *docoeff);
void amatdual(int *steps, int *po, int *re, int *nbrs, double *weights, double *alpha, int *lpo, int *lre, int *nn, double *adual);
void atimesb(double *a, double *b, int *n, double *prod, double *s);
void aug(double *A, int *ra, int *ca, double *Aaug);
void cubicpred(int *pointsin, double *X, double *coeff, int *nbrs, int *remove, int *inter, int *nn, double *weights, int *docoeff);
void fwtnp(double *input, double *f, int *nkeep, int *intercept, int *initboundhandl, int *neighbours, int *closest, int *LocalPred, int *n, double *coeff, double *lengthsremove, double *lengths, double *lca, int *pointsin, int *nc, int *doW, double *W, int *varonly, double *v);
void getnbrs2(double *X, int *remove, int *pointsin, int *lpo, int *neigh, int *closest, int *nbrs, int* index, int *nn);
void getridd(double *a, int *la, int *pos, double *b);
void getridi(int *a, int *la, int *pos, int *b);
void intervals(double *X, int *initboundhandl, int *n, double *inter);
void linearpred(int *pointsin, double *X, double *coeff, int *nbrs, int *remove, int *inter, int *nn, double *weights, int *docoeff);
void makelcaline(int *remove, int *nn, int *nbrs, double *alpha, double *weights, int *scheme, int *inter, int *closest, double *newline);
void mmult(double *x, double *y, int *nrx, int *ncx, int *ncy, double *ans);
void mycbind(double *a, double *b, int *ra, int *ca, int *cb, double *c);
void mycpyd(double *a, int *len, double *b);
void mycpyi(int *a, int *len, int *b);
void mydiag(double *d, int *n, int *ones, double *m);
void mymatchd(double *numa, double *numb, int *lnuma, int *lnumb, int *pos);
void mymatchi(int *numa, int *numb, int *lnuma, int *lnumb, int *pos);
void mymaxd(double *num, int *lnum, double *max, int *pos);
void mymaxi(int *num, int *lnum, int *max, int *pos);
void mymind(double *num, int *lnum, double *min, int *pos);
void mymini(int *num, int *lnum, int *min, int *pos);
void myrbind(double *a, double *b, int *ra, int *ca, int *rb, double *c);
void myrevd(double *dx, int *n, double *dy);
void myrevi(int *a, int *la, int *b);
void mysortd(double *a, int *la, double *sorted, int *order, int *inc);
void mysorti2(int *a, int *la, int* sorted, int *order, int *inc);
void mysvd(double *a, int *n, double *rvalues, double *rvectors, int *decreasing);
void myt(double *a, int *ra, int *ca, double *ta);
void mywhichd(double *num, int *lnum, double *a, int *pos);
void mywhichi(int *num, int *lnum, int *a, int *pos);
void pointsupdate(double *X, double *coeff, int *nn, int *index, int *remove, int *pointsin, double *wts, double *l, int *N, double *alpha, int *r);
void pts(double *input, double *start, int *n, double *X);
void quadpred(int *pointsin, double *X, double *coeff, int *nbrs, int *remove, int *inter, int *nn, double *weights, int *docoeff);
void rmatsolve(double *m, int *n, double *inv);
void transmatdual(double *lca, int *po, int *matno, int *lpo, int *nc, double *W, int *re);
void undopointsupdate(double *coeff, int *nbrs, int *index, int *remove, int *r, int *N, double *gamweights, double *l, double *lr, double *alpha, int *nn);
void updatelca(double *lca, int *nr, int *nc, double *newline, double *newlca);
void afromlca(double *lca, int *nc, int *rowno, double *alpha);
void clofromlca(double *lca, int *nr, int *nc, int *clo);
void interfromlca(double *lca, int *nr, int *nc, int *inter);
void nbrsfromlca(double *lca, int *nc, int *rowno, int *nbrs);
void schfromlca(double *lca, int *nr, int *nc, int *sch);
void wfromlca(double *lca, int *nc, int *rowno, double *weights);

void invtnp(double *X, double *coeff, double *lengths, double *lengthsremove, int *pointsin, double *lca, int *nadd, int *N, int *lr, int *nc, int *outpo, double *outlen);
void findadds(int *rem, int *l, double *lca, int *nc, int *index, int *li, int *a);
void delrow(double *M, int *nr, int *nc, int *i, double *Mnew);
void getnbrs(double *X, int *remove, int *pointsin, int *lpo, int *neigh, int *closest, int *nbrs, int* index, int *nn);
void mysorti(int *a, int *la, int* sorted, int *order, int *inc);
void mmult2(double *A, double *B, int *ra, int *ca, int *cb, double *C);
void mmult3(double *A, double *B, int *ra, int *ca, int *cb, double *C);
void fwtnpperm(double *input, double *f, int *nkeep, int *intercept, int *initboundhandl, int *neighbours, int *closest, int *LocalPred, int *n, double *coeff, double *lengthsremove, double *lengths, double *lca, int *pointsin, int *nc, int *traj, int *doW, double *W, int *varonly, double *v);

#endif

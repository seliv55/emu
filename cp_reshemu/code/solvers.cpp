//---------------------------------------------------------------------------
#include <iostream>
#include "nr.h"
#include "new.h"
#include "StiffIntegratorT.h"
#include "NonStiffIntegratorT.h"
#include "solvers.h"
//---------------------------------------------------------------------------
using namespace std;
   extern Ldistr emudyn;
DP dxsav;  
int kmax,kount;
Vec_DP *xp_p;
Mat_DP *yp_p;
Vec_INT *ija_p;
Vec_DP *sa_p;
int nrhs;   // counts function evaluations
Vec_DP *x_p;
Mat_DP *d_p;
double Vt;
string kin0;

void jac(const double& time,  double y[], const double yPrime[], double **PD, double& CJ, const double rPar[], const int iPar[]){}
void isores(const double& T, double y[],const double yprime[], double delta[], int& iRes, const double rPar[], const int iPar[]){
       int len=emudyn.glen();
       double dy[len];
	emudyn.mdistr(y, dy, T);
	for(unsigned i=0;i<len;i++) delta[i]=dy[i]-yprime[i];
}
void res(const double& T, double y[],const double yprime[], double delta[], int& iRes, const double rPar[], const int iPar[]){
       int len=emudyn.gnmet();
       double dy[len];
	emudyn.f(y, dy);
	for(unsigned i=0;i<len;i++) delta[i]=dy[i]-yprime[i];
}
/**/
double Ldistr::ddisolve() {
        int info[15],idid=0,lrw=800000,liw=1190,iwork[1190],ipar[2],ires=0,ikin=25;
double t=0.,xi=0., yprime[len],rtol=0.075,atol=1.0e-6,h0=0.1e-5, hmax=5.2, rpar[2], rwork[800000];
        for(int i=0;i<15;i++)   info[i]=0;
      info[6]=1; rwork[1]=hmax; rwork[2]=h0; info[10]=1;
 double tout,  tm;
//       shiso(xx);
    for(int i=0;i<len;i++) xx[i]=xinit[i];
        mdistr(xx, yprime,t);
	const int KMAX(2);
	for(int i=1;i<ntime;i++){
	tm=tex[i]/(double)ikin;
        for(int k=0;k<ikin;k++){ tout=t+tm;
ddassl_(isores,len,t,xx,yprime,tout,info,rtol,atol,idid,rwork,lrw,iwork, liw,  rpar, ipar, jac);
 if(idid<0) {  throw("dassl problem"); }
    t=tout;
    }
	 xi+=chisq();}
return xi;}

void derivsl(const DP x, Vec_IO_DP &y, Vec_O_DP &dydx){
	double *py=&y[0]; 
	DP *pdydt=&dydx[0];
        nrhs++; 
	emudyn.mdistr(py, pdydt, x);
}

double Ldistr::integrbs(){
  DP eps=1.0e-6,h1=0.00001,hmin=1.0e-11,x1=0.0, xfin, tm,xi=0.;
   const int KMAX(2);  Vec_DP yy(len); 
   for(int i=0;i<len;i++) {yy[i]=xinit[i]; cout<<yy[i]<<" ";} cout<<endl;
      xp_p=new Vec_DP(KMAX); yp_p=new Mat_DP(len,KMAX);
        Vec_DP &xp=*xp_p;  Mat_DP &yp=*yp_p;
    int nbad,nok,ikin=10; nrhs=0; kmax=KMAX;
//       shiso(pyinit);             //show isotopomers for all EMUs
       tm=1.;
	for(int i=1;i<ntime;i++){
      cout<<"ttime="<<x1<<endl;
	tm=tex[i]/ikin;
        dxsav = tm/((double)(KMAX-1));
        for(int k=0;k<ikin;k++){ xfin=x1+tm;cout<<"pass!!!***"<<endl;
    NR::odeint(yy,x1,xfin,eps,h1,hmin,nok,nbad,derivsl,NR::rkqs);
    x1=xfin;
    }//for(int i=0;i<len;i++) {xx[i]=pyinit[i]; cout<<xx[i]<<" "; } cout<<endl;
	}
      cout<<"ttime="<<x1<<endl;
   for(int i=1;i<len;i++) xx[i]=yy[i];
//        delete yp_p;
//        delete xp_p;
 return xi;}

void isT::Function(double x, double *y, double *dy){
	emudyn.f(y, dy);
}
void Jacobian(double x, double *y, double **J){}
void Mass(double **M){} // Mass

void isoT::Function(double x, double *y, double *dy){
	emudyn.mdistr(y, dy,x);}
	
void Ldistr::tsolve(const double tmax){
	// dimension of problem
    for(int i=0;i<len;i++) xx[i]=xinit[i];
	// initial value for x
	 int kmax=25;
   double xbeg(0.0), dx = tmax/((double)kmax), xend=dx;
	// rtoler and atoler are scalars
	int itoler(0);
	// relative tolerance
	double *rtoler = new double(1.0e-7);
	// absolute tolerance
	double *atoler = new double(1.0e-7);
	// use SolutionOutput routine
	const int iout(0);
	// initial step size
	double hinit(0.00001);
	// analytical Jacobian function provided
	 int ijac(0);
	// number of non-zero rows below main diagonal of Jacobian
	int mljac(len);
	// number of non-zero rows above main diagonal of Jacobian
	int mujac(len);
	// Mass matrix routine is identity
	const int imas(0);
	int mlmas(len);
	int mumas(0);
	
	// Use default values (see header files) for these parameters:
	double hmax(0.0);
	int nmax(0);
	double uround(0.0), safe(0.), facl(0.0), facr(0.0);
	int nit(0);
	bool startn(false);
	int nind1(0), nind2(0), nind3(0), npred(0), m1(0), m2(0);
	bool hess(false);
	double fnewt(0.0), quot1(0.0), quot2(0.0), thet(0.0);
	
ostringstream skin; 
skin<<"0 ";
//          for (int j=0;j<nmet;j++) cout<<xx[j]<<" "; cout<<endl;
//for(int i=0;i<kmax;i++){
//	isT stiffT(len, yy, xbeg, xend, dx, itoler, rtoler, atoler, iout, hinit, hmax, nmax, uround, safe, facl, facr, ijac, mljac,		mujac, imas, mlmas, mumas, nit, startn, nind1, nind2, nind3, npred, m1, m2, hess, fnewt, quot1, quot2, thet);
//  isT stiffT(len, xx, xbeg, xend, dx, ijac, mljac,mujac, imas, mlmas, mumas);
//	stiffT.Integrate();
//        skin<<xend<<" ";
//	for(int i=0;i<nmet; i++) skin<<xx[i]<<" "; skin<<endl;
//	xbeg=xend; xend += dx;
//	}	//for (int j=0;j<nmet;j++) cout<<xx[j]<<" "; cout<<endl;
//	kin0=skin.str();
//   int ni=nmet;
//   siso(xx,ni);
//       emudyn.sinit();
//    xbeg=0.0;  xend=dx;	 mljac=(len); mujac=(len); mlmas=(len);
	
//          shiso(yy);// kmax=1;
for(int i=0;i<kmax;i++){
  isoT stifT(len, xx, xbeg, xend, dx, itoler, rtoler, atoler, iout, hinit, hmax, nmax, uround, safe, facl, facr, ijac, mljac,   mujac, imas, mlmas, mumas, nit, startn, nind1, nind2, nind3, npred, m1, m2, hess, fnewt, quot1, quot2, thet);
//  isoT stifT(len, xx, xbeg, xend, dx, ijac, mljac,mujac, imas, mlmas, mumas);
	stifT.Integrate();
	xbeg=xend; xend += dx;
	}
	delete rtoler;
	delete atoler;
}


#include <iostream>
#include <sstream>
#include "nr.h"
#include "nums.hh"
#include "tk.hh"
#include "nv.hh"
#include "modlab.h"
#include "analis.h"
#include "solvers.h"
using namespace std;
double Analis::descent(double factor,int ip){ 
	double xi1, a, sens, slim=0.01, dp=factor-1.;
	const double xili(0.9998);
	int k, jfail(0), parcp[nrea];
  double tcal=tf;
  double xi0 = (chimin+suxx*suxx); cout<<"xi0="<<xi0<<endl;
  int npf = Problem.getListFit(parcp); ifn++;
  if(ip>=0) {for(k=0;k<npf;k++) if(parcp[k] = ip) break;
   cout<<"Descent: par="<<Problem.rea[parcp[k]].getname()<<endl;
                 for(int i=k;i<npf;i++) parcp[i] = parcp[i+1]; npf--;}
   cout<<"npf="<<npf<<"  "<<parcp[npf-1]<<endl;
  while (npf>0)  {int j(0), flag(0);
   int i = rand() % npf; npf--; 
    while (flag<2) { a = Problem.rea[parcp[i]].v();
     Problem.rea[parcp[i]].setVm(a*factor); cout << parcp[i] << ")";
      try{ xmin = solve(); xm=mader;  double tr=(double)tf/(double)tcal;
	xi1 = (xmin+suxx*suxx)*tr*mader/xm;
	 sens=(xi0-xi1)/xi0/dp;
	  cout<<Problem.rea[parcp[i]].getname()<<": sens="<<sens<<endl;
       } catch( char const* str )
	  {cout << "exception: "<< str <<endl; sens=0.;  jfail++; flag=2; Problem.rea[parcp[i]].setVm(a);
  if(!(jfail-5)){jfail=0; Problem.restoreVm(nrea,nv2); } }
 if (sens>slim) {
   xi0 = xi1; chimin=xmin; tcal=tf; xm=mader; jfail=0; j++; Problem.write(tf,ifn,xmin,suxx,0);
     if(!(j-3)) flag=2;   Problem.storeVms(nrea,nv1);
       	}
 else if(sens<-slim) {factor = 1./factor; ++flag; Problem.rea[parcp[i]].setVm(a);}
  else {flag=2; Problem.rea[parcp[i]].setVm(a);}
        }//end while flag
	for ( k=i;k<npf;k++) parcp[k] = parcp[k+1];}
 return xi0;}
 
 void Analis::rconfint(ifstream& fi, double ami[],double ama[]){
  string aaa;
   for(int i=0;i<nrea;i++) {fi>>aaa>>aaa>>aaa>>ami[i]>>ama[i]; getline(fi,aaa);}
 }
 
 void Analis::stepdown(double factor, int& i,double fdes){
   double oldp;
   cout << "; v["<<Problem.rea[i].getname()<<"]: "<<Problem.rea[i].v()<<endl;
  try { chimin = solve();
     for(int j=0;j<5;j++) {oldp = Problem.rea[i].chanVm(factor); chimin = solve(); 
       if((chimin-x00)>6.7) break;}   Problem.write(tf,ifn,chimin,suxx);   
        descent(fdes,i); Problem.write(tf,ifn,chimin,suxx); 
   if (chimin<x00) {x00=chimin; Problem.storeVms(nrea,nv1); Problem.write(tf,ifn,chimin,suxx);} }
  catch( char const* str ){ cout << "exception: "<< str <<endl; chimin=x00+1000.; Problem.rea[i].setVm(oldp);}
       }
 
void Analis::confidence(double factor,double fdes){ 
	double a0mi[nrea], a0ma[nrea], a1mi[nrea], a1ma[nrea];
	 ifstream fi("statfl"); rconfint(fi,a0mi,a0ma); fi.close();
	   fi.open("statfl1");  rconfint(fi,a1mi,a1ma); fi.close();
	int parcp[nrea], npf = Problem.getListFit(parcp);
	          Problem.storeVms(nrea,nv1);  chimin=x00;
  while (npf>0)  {   int i = rand() % npf;
   if((a0ma[parcp[i]]<a1mi[parcp[i]])&&(a0ma[parcp[i]]>7.e-7)){
       stepdown(factor,parcp[i],fdes);}
   else if((a0mi[parcp[i]]>a1ma[parcp[i]])&&(a1ma[parcp[i]>7.e-7])) {
       stepdown(1./factor,parcp[i],fdes);}
   npf--; for (int k=i;k<npf;k++) parcp[k]=parcp[k+1];
        Problem.restoreVm(nrea,nv1);  }
   }

void Fit::perturb(const double f1){
   int i=1; double sign,fact;
     while(par[i]>=0){// cout<<"par["<<i<<"]="<<par[i]<<endl;
      sign = (double)rand() / (double)RAND_MAX;
	fact = 1.- f1*(0.5 - sign);
	rea[par[i]].setVm(rea[par[i]].v()*fact);
	   i++;	}
	   }
	   
double Analis::dermax(){
  double dx[numx], amax, dmax(0.),dmin(0.),f=1.02, ff; int iax(0),iin(0),im;
   for(int i=0;i<numx;i++) xx[i]=xinit1[i];
    tsolve(3200.); Problem.f(xx,dx);
   for(int i=0;i<(numx-2);i++) if(dmax<dx[i]) {dmax=dx[i]; iax=i;}
                       else if(dmin>dx[i]) {dmin=dx[i]; iin=i;}
   dmin*=-1; if(dmax>dmin) {ff=f; amax=dmax; im=iax;}
        else {ff=1./f; amax=dmin; im=iin;}
   for(int i=0;i<numx;i++) xx[i]=xinit1[i];
         cout<<Problem.namex[im]<<", deriv="<<amax<<endl;
               return amax;}

               
void Fit::fitc(double dc,double dm,int iin,int iout){
   double f=0.95;
  if(dc>0) {if(dm<0)rea[iout].setVm(rea[iout].v()*f);   else rea[iin].setVm(rea[iin].v()/f);}
  else {if(dm>0) rea[iout].setVm(rea[iout].v()/f); else rea[iin].setVm(rea[iin].v()*f);}
 }

template<size_t l> void Ldistr::setmet(Metab<l>& met,data cmet[],data emet[][l+1],string sname,int vin,int vout){
   int tax=ntime-1;
 double s=met.gett(), m0=met.getm0(), dc=(s-cmet[tax].mean)/cmet[tax].sd, dm=m0-emet[tax][0].mean, lim=1., xi, f=0.95;
 cout<<sname; 
 for(int i=0;i<7;i++) {
  if(dc>lim){if(dm<0) Problem.rea[vout]*=f;   else Problem.rea[vin]/=f;}
  else if(dc<(-lim)){ if(dm>0) Problem.rea[vout]/=f; else Problem. rea[vin]*=f;}
   else break;
     xi=solve();
      s=met.gett(); m0=met.getm0(); dc=(s-cmet[tax].mean)/cmet[tax].sd; dm=m0-emet[tax][0].mean;
  cout<<"conc="<<s<<"; m0="<<m0<<endl;}
  if(xi>1.) Problem.write(tf,ifn,xi,suxx,0);
}
template<size_t l> void Ldistr::setm0(Metab<l>& met,data cmet[],data emet[][l+1],string sname,int vin,int vout){
   int tax=ntime-1;
 double s=met.gett(), m0=met.getm0(), dc=(s-cmet[tax].mean)/cmet[tax].sd,dm=(m0-emet[tax][0].mean)/emet[tax][0].sd, lim=1., xi, f=0.95; 
  cout<<sname;
   for(int i=0;i<7;i++){
  if(dm>lim) {if(dc>0.)  Problem.rea[vin]/=f; else Problem.rea[vout]/=f;}
  else if(dm<(-lim)) {if(dc<0.) Problem.rea[vin]*=f; else Problem.rea[vout]*=f;}
  else break;
    xi=solve();
 s=met.gett(); m0=met.getm0(); dc=(s-cmet[tax].mean)/cmet[tax].sd; dm=(m0-emet[tax][0].mean)/emet[tax][0].sd;
  cout<<"conc="<<s<<"; m0="<<m0<<endl;}
   if(xi>1.) Problem.write(tf,ifn,xi,suxx,0);
}
bool Parray::belong(int ip) const{
 bool b=false; int i=1; 
  while(par[i]+1) {if(!(ip-par[i])){b=true; break;} i++;}
  return b;}
   

void Analis::coord(const double f1,double fdes){
	double xi, xi0,dif0,limdx=20e-4,xlim=2.0;
	int  ifail(0);
cout<< "\nPerturbation+CoordinateDescent:\nPar#\txi2con\txi2iso\tParValue\tdif:" <<endl;
     for(;;){
        Problem.perturb(f1);
 for(int i=0;i<3;i++) {  try {
    xi=solve(); xm=mader; Problem.write(tf,ifn,xi,suxx);}
  catch( char const* str ){cout << "exception: "<< str <<endl; {Problem.restoreVm(nrea,nv2); }
 chimin=xi; dif0=dif;
//   cout<<"* reduce Lac *"<<endl; descent(fdes,-2);
     cout<<"* reduce xi total *"<<endl; descent(fdes);
 //      xi=horse.fitcon(); Problem.write(tf,ifn,xi,suxx); ifn++;
 //     horse.setcon();
       /*  while(xmax()>xlim) ; 
    double dm(1.);  while(dm>limdx) { dm=dermax(); if(dm>1.) {xi=solve();  Problem.write(tf,ifn,xi,suxx); break;}}*/
 chimin=solve(); xm=mader; Problem.write(tf,ifn,chimin,suxx); }
/*  while (chimin<xi) {
//  while (dif<dif0) {
    try {  while(xmax()>xlim) ;
    double dm(1.);  while(dm>limdx) { dm=dermax(); if(dm>1.) {xi=solve();  Problem.write(tf,ifn,xi,suxx); break;}}
//    double dm(1.);  while(dm>limdx) { dm=dermax(); if(dm>1.) {xi=solve();  Problem.write(tf,ifn,xi,suxx); break;}}
   for(int i=0;i<numx;i++) xinit1[i]=xx[i]; chimin=solve(); dif0=dif; Problem.write(tf,ifn,chimin,suxx);  }
  catch( char const* str ){cout << "exception: "<< str <<endl; get(nrea,Problem.nv,nv2); }
    xi=chimin; cout<<"* reduce Lac *"<<endl; descent(fdes,-2);
      cout<<"* reduce xi total *"<<endl; descent(fdes);}*/
	}
}}
/*
void Analis::sensitiv(const double tmax){ 
  double factor(1.1), xi1, xm1, a;
    get(nrea,nv2,Problem.nv); get(numx,xinit1,xx);//saves nv&xx
  chimin= solve();
  double xm0= getmax();
  int *par = Problem.getFitPar(), ip=0,jmi,mtb[15],nmet;
  string smtb[15];
   while (par[ip]+1)  ip++;
    jmi=horse.getmi(mtb,nmet,smtb); 
    cout<<"nmet="<<nmet<<endl;
    cout<<"mtb: "; for(int i=0;i<nmet;i++) cout<<mtb[i]<<" ";cout<<endl;
  double aa[ip+1][jmi];
   horse.miso(&aa[0][0]);
   for (int i=0;i<ip;i++) {	int flag(0); 
    while (flag<2) {cout<<"par["<<i<<"]="<<par[i]<<endl;
		a = Problem.getVal(par[i]);
     Problem.setVm(par[i], a*factor);
     try{ xmin = solve(); xm1 = getmax(); flag=2;
     } catch( char const* str ) {
         cout << "exception: "<< str <<endl; 
	  factor = 1./factor; ++flag;
	}
    }//end while flag
    horse.miso(&aa[i+1][0]);
    for(int j=0;j<jmi;j++) aa[i+1][j] -= aa[0][j];
     get(numx,xx,xinit2);//gets nv&xx
     Problem.setVm(par[i], a);
  }
  
    ofstream fi("sens.csv");
    fi<<"par# | "; for(int k=0;k<nmet;k++) fi<<smtb[k]<<"| "; fi<<endl;
      fi<<"init | ";
      for(int j=0;j<mtb[0];j++) fi<<aa[0][j]<<" "; fi<<"| ";
for(int k=1;k<nmet;k++){for(int j=mtb[k-1];j<mtb[k];j++) fi<<aa[0][j]<<" "; fi<<"| ";}
      fi<<endl;
    for(int i=1;i<ip;i++) {
      fi<<par[i]<<" | ";
      for(int j=0;j<mtb[0];j++) fi<<aa[i][j]<<" "; fi<<"| ";
for(int k=1;k<nmet;k++){for(int j=mtb[k-1];j<mtb[k];j++) fi<<aa[i][j]<<" "; fi<<"| ";}
      fi<<endl;
      }
}		

void Analis::cross(double p1[],double p2[]){
	int i,j1=rand() % nrea;
	for(i=0;i<j1;i++) p1[i]=p2[i];
	for(i=j1;i<nrea;i++) p2[i]=p1[i];
}
void Analis::mutate(int npf,int par[]){
	const double f1(0.3);
	       int sign = rand() % 100;
	       int imut=rand()%(npf-1);
		double fact = 1.- f1*(0.5 - sign*0.01);
		Problem.setVm(par[imut], Problem.getVal(par[imut])*fact); 
}
void Analis::genetic(const double tmax,const int ngen){
	char fn[11];
	double mean,xi, xi1,xm,xmm(100.);
	int itmp,par[nrea], ifail(0);
	get(nrea,nv2,Problem.nv); get(numx,xinit2,xx);//saves nv&xx
cout<<"evolution: crossing"<<endl;
		int npf = Problem.getListFit(par);
		int sign,imut,nfi;
 for(int k=0;k<9;k++){  nfi=ngen;
for(int j=1;j<ngen;j += 2){
	 sprintf(fn,"%i",j); cout<<j<<endl;
      Problem.read(itmp,mean,fn);
      mutate(npf,par);
	get(nrea,nv2,Problem.nv); get(numx,xinit1,xx);//saves nv&xx
	int j1= 1+ (rand() % (ngen-5));  cout<<"j1="<<j1<<endl;
	 sprintf(fn,"%i",j1); 
      Problem.read(itmp,mean,fn);
      mutate(npf,par);
      cross(Problem.nv,nv2);
	try { for(int i=0;i<2;i++){
      if(i) get(nrea,Problem.nv,nv2);
      xi=solve();	 xm=getmax();
if(xm<xmm){
        Problem.write(tf,nfi,xi,xm,0); nfi++;}	}
	} catch(char const * ierr){cout << "error " << endl; 
	 }
	} //NR::sort2a(a,b);
	nfi--;
	cout<<"nfi="<<nfi<<endl;
	Problem.stat(nfi);}
	for(int j=2;j<nfi;j+=5){
	 sprintf(fn,"%i",j);
	 cout<<"file:"<<fn<<endl;
      Problem.read(itmp,mean,fn);
	chimin = solve();
	descent(tmax);
	 xm=getmax();
                Problem.write(tmin,nfi,chimin,xm);}
}*/
void Analis::swarm(const double tmax,const int ngen){
	char fn[11];
double xi,r1,r2,f1,xmm(100.),c1(0.8),c2(1.),w(0.5),xm,v[nrea],f0[ngen];
double pbest[ngen][nrea],pcurr[ngen][nrea];
	int it,ifail(0),par[nrea];
cout<<"swarm:"<<endl;
		int sign,imut,nfi;
for(int i=1;i<ngen;i++){
 sprintf(fn,"%i",i);
  xi=Problem.read(it,xm,fn);
	Problem.storeVms(nrea, &pbest[i][0]);
	Problem.storeVms(nrea, &pcurr[i][0]);
  f0[i] = xi; 
  }
f0[0]=f0[1]; //get(nrea,&pbest[0][0],&pbest[1][0]);
		for(;;){
for(int i=2;i<ngen;i++){cout<<" i="<<i;
	Problem.restoreVm(nrea,&pcurr[i][0]);
	 r1 = (double)rand() / (double)RAND_MAX;
	 r2 = (double)rand() / (double)RAND_MAX;
   for(int j=0;j<nrea;j++){
 v[j]=w*v[j]+c1*r1*(pbest[i][j]-Problem.rea[j].v())+c2*r2*(pbest[0][j]-Problem.rea[j].v());
    Problem.rea[j].setVm(Problem.rea[j].v() + v[j]);
    }/**/
	try {
	//   get(numx,xx,xinit1); 
  xi=solve(); cout<<"; xi="<<xi<< endl;
    xm=suxx; Problem.storeVms(nrea,&pcurr[i][0]);
       f1 = xi;//+xm)*sqrt((dif/dif0)); 
	//   get(numx,xx,xinit1); 
    if(f1<f0[i])  {Problem.write(tf,i,xi,xm,0); 
       cout<<"f0="<<(f0[i]=f1)<<endl;
          Problem.storeVms(nrea,&pbest[i][0]); 
    if(f1<f0[0]){f0[0]=f1; Problem.storeVms(nrea,&pbest[0][0]);}}
	} catch(char const * ierr){cout << "error " << endl;}
		}}
}


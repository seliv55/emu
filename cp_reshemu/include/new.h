//---------------------------------------------------------------------------
#include <fstream> 
#include <string> 
#include <vector>
#ifndef newh
#define newh
const int tt=5;

class data {
public:
	 double mean, sd;
        data(double a=0.):mean(a),sd(a){}
        ~data(){}
};

class Reakcia {
 int npar;
 std::string naz;
 double flx[4], Vm, *par;
public:
 double *v(double x){flx[1]=par[0]/(par[1]+x); flx[0]=flx[1]*x; return flx;}
 
 double *v(double x,double y){flx[3]=par[0]/(par[1]+x)/(par[2]+y); flx[2]=flx[3]*x; flx[1]=flx[3]*y; flx[0]=flx[1]*x;  return flx;}
 
 double *gflx(){return flx;}
 
  std::string& gnaz(){return naz;}
  
 void rpar(std::ifstream& fi);

 void changeVm(double f){ par[0]*=f;}
 
Reakcia(){}
~Reakcia(){delete[] par;}
};

class Emu {
  std::string name;
  double *iso, *diso, sumiso, xi;
  int niso,nvec;
public:
 int set(std::string aaa){name=aaa; niso=name.length()+1; return niso; }
 
 void sinit(double c0){iso[0]=c0; for(int i=1;i<niso;i++) iso[i]=0.;}
 
 void volume(double vol){ for(int i=0;i<niso;i++) diso[i] /= vol;}

 double sumt(){ sumiso=iso[0]; for(int i=1;i<niso;i++) sumiso += iso[i]; return sumiso;}
 
 double sfrac(data *eda);
    
 double chisq(data *eda);
    
 int gniso(){return niso;}
 
 int gnvec(){return nvec;}
 
 std::string gname(){return name;}
 
 bool cmname(std::string sss){return (sss==name);}
 
 void markinit(std::string mark, double mval);
           
 int siso(double vec[],int ni){nvec=ni; iso=&vec[ni]; return ni+niso;}
 
 int sdiso(double vec[],int ni){diso=&vec[ni]; for(int i=0;i<niso;i++) diso[i]=0.; return ni+niso;}
 
  double zero(double flx) {double sum(0); for(int i=0;i<niso;i++) {
              double x = flx*iso[i]; diso[i] -= x; sum += x;} return sum;}
              
  double uni(double flx,Emu* pr) {double sum(0); for(int i=0;i<niso;i++) {
              double x = flx*iso[i]; diso[i] -= x; pr->diso[i] += x; sum += x;} return sum;}
              
  double uni1(double flx,Emu* pr) {double sum(0); for(int i=0;i<niso;i++) {
              double x = flx*iso[i]; pr->diso[i] += x; sum += x;} return sum;}
              
  double bi(double flx,Emu* su2,Emu* pr) {double sum(0); for(int i=0;i<niso;i++) {
             double x = flx*iso[i]; for(int j=0;j<su2->gniso();j++){ double x1=x*su2->iso[j];
             diso[i] -= x1; su2->diso[j] -= x1;
              pr->diso[i+j] += x1; sum += x1;}}  return sum;}
              
Emu(){nvec=0;}
~Emu(){}
};

class Edata {
  data econc[tt], *emiso[tt];
  double xicon[tt], xi[tt];
  int ntime, niso;
  std::string sfrg;
public:
 void setex0();
 
 double getmi0(int imi){return emiso[0][imi].mean;}
 
 data *geda(int nt){return emiso[nt];}
 
 int gniso(){return this->niso;}

 void setsfrg(std::string sss){ sfrg=sss; }
 
 std::string getsfrg(){return sfrg;}
 
 void setrav(int nt1,int nt2);
 
 void read(std::ifstream& fi,int nt);
 
 void readc(std::ifstream& fi,  int nt);
 
 Edata():niso(0){}
 ~Edata(){for(int i=0;i<ntime;i++) if(niso) delete[] emiso[i];}
};

class Metab {
  double conc;
  std::string name;
  int nemus,len, ieref;
  Emu *emu;
  Edata edata;
public:
  Edata* getex(){return &edata;}
  int read(std::ifstream& fi);
     
  int glen(){return len;} //length of array of isotopomers
  
  double shiso();
        
  double chisq();
        
  std::string& gname(){return name;}
 
  void scon(double *vec,int ni){ vec[ni]=conc;}
  
  int markinit(std::string mname,std::string mark,double mval);
                           
  double gcon(){return emu[0].sumt();}
  
  Emu* gemu(int iem){return &emu[iem];}
  
  int siso(double *vec,int ni){for(int i=0;i<nemus;i++) ni=emu[i].siso(vec,ni); return ni;}
  
  int sdiso(double *vec,int ni){for(int i=0;i<nemus;i++) ni=emu[i].sdiso(vec,ni); return ni;}
  
  void volume(double vol){  if(name.at(0)!='E') for(int i=0;i<nemus;i++) emu[i].volume(vol);}
  
  int fndemu(std::string sss);
       
  int findemu(std::string exemu);
       
  void sinit(){for(int i=0;i<nemus;i++) emu[i].sinit(conc);}
  
  void concor(double& a){ a= emu[0].sumt();}
  
  void sumt(){for(int i=0;i<nemus;i++) emu[i].sumt();}
  
  void readc0(std::ifstream& fi);
  
Metab(){ieref=-1;}
~Metab(){  delete[] emu;}
};

class Ldistr {
 int hkf, hkr, pfk, fbase, aldf, aldr, t3pep, pept3, pk, pyrlac, lacpyr, citdr, citdf, csyn, akgcit, citakg, akgdr, akgdf, citoxc, liase, ppp, pdh, pc, malic, malicc, pyrdf, pyrdr, maloa, oamal, oadr, oadf, akgfum, glnin, glnout, gluin, gluout, tkp5k, tkh6k, tks7k, tkt3a, tke4a, tkp5a, tas7k, tah6k, tae4a, tat3a, nre;
 int Eglc, h6p, fbp, t3p, pep, pyrc, Elac, citc, cit, oa, accoa, akgc, co2, akg, oac, p5, pyr, fumal, Egln, Eglu, gae, e4p, s7p, dhe, nmet;
 int ntime, len; //number of timepoints, accounted mass isotopomers
 Metab *met;
 Reakcia *rr;
 std::vector<int> vpar;
 double  *xx, *xinit, tex[tt], Vi, Vt, mu, mval;
 std::string mname, mark;
public:
 void setmetrea();
 
 double* getxx(){return xx;};
 
 double* getxinit(){return xinit;};
 
 void f(const double *y,double *dydx);
 
 void ff(const double *y,double *dydx);
 
 void mdistr(double *py,double *pdydt,double t);
 
 void markinit();
 
 void rpar(std::string fn);
 
 void read(std::string fn);
 
 double readex(std::string fn,int itp);
 
 int glen() { return len;}
 
 int gnmet() { return nmet;}
 
 int scon(double *vec){for(int i=0;i<nmet;i++) met[i].scon(vec,i); return nmet; }
 
 void gcon();
 
 void concor(double *vec){for(int i=0;i<nmet;i++) met[i].concor(vec[i]);  }
 
 void siso(double *vec,int& ni){for(int i=0;i<nmet;i++) ni=met[i].siso(vec,ni); }
 
 void sdiso(double *vec,int& ni){for(int i=0;i<nmet;i++) ni=met[i].sdiso(vec,ni); }
 
 void volume(double vol){for(int i=0;i<nmet;i++) met[i].volume(vol); }
 
 void sinit(){for(int i=0;i<nmet;i++) met[i].sinit(); }
 
 void sumt(){for(int i=0;i<nmet;i++) met[i].sumt(); }
 
 void shiso();
 
 double chisq();
 
 void tsolve(const double tmax);
 
 void tisolve(const double tmax);
 
 double ddisolve();
 
 double integrbs();
 
  void setinit(std::string fimod="metemu",std::string fipar="1",std::string fiex="mglc",int ntp=1);
  
Ldistr(){}
~Ldistr(){delete[] met; delete[] rr;  delete[] xx; delete[] xinit;}
};

//---------------------------------------------------------------------------
#endif


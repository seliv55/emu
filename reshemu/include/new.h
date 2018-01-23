//---------------------------------------------------------------------------
#include <fstream> 
#include <string> 

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
 double flx[4], *par;
public:
 double *v(double x){flx[1]=par[0]/(par[1]+x); flx[0]=flx[1]*x; return flx;}
 double *v(double x,double y){flx[3]=par[0]/(par[1]+x)/(par[2]+y); flx[2]=flx[3]*x; flx[1]=flx[3]*y; flx[0]=flx[1]*x;  return flx;}
 double *gflx(){return flx;}
  std::string& gnaz(){return naz;}
 void rpar(std::ifstream& fi){ std::string aaa;
    fi>>aaa>>naz>>npar; par=new double[npar];
      for(int i=0;i<npar;i++) fi>>par[i];
        }
Reakcia(){}
~Reakcia(){delete[] par;}
};

class Emu {
  std::string name;
  double *iso, *diso, sumiso;
  int niso,nvec;
public:
 int set(std::string aaa){name=aaa; niso=name.length()+1; return niso; }
 void sinit(double c0){iso[0]=c0; for(int i=1;i<niso;i++) iso[i]=0.;}
 void volume(double vol){ for(int i=0;i<niso;i++) diso[i] /= vol;}
// void shiso(){std::cout<<"\t"<<name<<": ";for(int i=0;i<niso;i++) std::cout<<iso[i]<<" ";std::cout<<"\n";}
 double sumt(){ sumiso=iso[0]; for(int i=1;i<niso;i++) sumiso += iso[i]; return sumiso;}
 
 double sfrac(data *eda){ double a, xi(0);  sumt();
     std::cout<<name<<": "; 
     for(int i=0;i<niso;i++) { double fr=iso[i]/sumiso;
         if(eda[i].sd>0.009) {a=(fr-eda[i].mean)/eda[i].sd; xi += a*a;}
         std::cout<<eda[i].mean<<"-"<<fr<<", ";}//<<" xi:"<<xi
     std::cout<<"; sum="<<sumiso<<"; xi="<<xi<<"\n";   
    return xi;}
    
 int gniso(){return niso;}
 int gnvec(){return nvec;}
 std::string gname(){return name;}
 bool cmname(std::string sss){return (sss==name);}
 
 void markinit(std::string mark, double mval){  // mark: positions containing ¹³C
   int l=mark.length(), nmark=0;                // nmark: number of ¹³C in a givem EMU
      for(int i=0;i<l;i++) {  int npos=name.find(mark.at(i)); // check if the EMU contains the labeled position
               if(npos+1) nmark++;      }
      if(nmark){ iso[nmark]=iso[0]*mval; iso[0] -= iso[nmark];}
        std::cout<<name<<" nmark: "<<nmark<<" val: "<<iso[nmark]<<'\n';
           }
           
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
//             diso[i] -= x1; su2->diso[j] -= x1;
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
 void setex0(){ emiso[0]=new data[niso];
     emiso[0][0].mean=1.; emiso[0][0].sd=0.01;
      for(int i=1;i<niso;i++){ emiso[0][i].mean=0.; emiso[0][i].sd=0.01;  }
                         }
 double getmi0(int imi){return emiso[0][imi].mean;}
 
 data *geda(int nt){return emiso[nt];}
 
 int gniso(){return this->niso;}

// double getxi(double iso[],int liso){ double xi(0);  for(int i=0;i<liso;i++) {
//                  double a=(iso[i]-emiso[i]->mean)/emiso[i]->sd; xi +=a*a; }
//                   return xi; }
                   
 void setsfrg(std::string sss){ sfrg=sss; }
 std::string getsfrg(){return sfrg;}
 void setrav(int nt1,int nt2){// sets edata at t=0
      for(int i=0;i<=niso;i++){ emiso[nt1][i].mean=emiso[nt2][i].mean; emiso[nt1][i].sd=emiso[nt2][i].sd;  }
                         }
 void read(std::ifstream& fi,int nt){
	std::string aaa;//i: isotopomer; nt: time; imt: a number assigned to isotopomer
	fi>>sfrg; niso=sfrg.length()+1;
       emiso[nt]=new data[niso];
   for(int i=0;i<niso;i++) fi>>emiso[nt][i].mean;  fi>>aaa; fi>>aaa; // std::cout<<aaa<<" "; 
   for(int i=0;i<niso;i++) {fi>>emiso[nt][i].sd; if(emiso[nt][i].sd<0.01) emiso[nt][i].sd=0.01;}
	}
 void readc(std::ifstream& fi,  int nt){ std::string aaa;
     fi>>aaa;
     for(int i=0;i<nt;i++) fi>>econc[i].mean;
     for(int i=0;i<nt;i++) {fi>>econc[i].sd; if(econc[i].sd<0.01) econc[i].sd=0.01;}
	}
 int shex(int nt){if(niso) {std::cout<<"Edata "<<sfrg<<": ";
   for(int i=0;i<niso;i++) std::cout<<emiso[nt][i].mean<<" "; std::cout<<"\n\n";}
     return niso; }
 Edata(){}
 ~Edata(){for(int i=0;i<ntime;i++) delete[] emiso[i];}
};

class Metab {
  double conc;
  std::string name;
  int nemus,len, ieref;
  Emu *emu;
  Edata edata;
public:
  Edata* getex(){return &edata;}
  int read(std::ifstream& fi){ std::string aaa; len=0;
   fi>>name>>nemus; emu=new Emu[nemus]; // std::cout<<name<<" "<<nemus<<std::endl;
    for(int i=0;i<nemus;i++) {fi>>aaa; len+=emu[i].set(aaa);} 
     return len; }
     
  int glen(){return len;} //length of array of isotopomers
  
  double shiso(double *vec){ double xi(0); int emunum;//show  EMUs
      if(edata.gniso()){std::cout<<name<<":\t"<<edata.getsfrg()<<"-";
        emunum=findemu(edata.getsfrg()); if(emunum>=0) xi=emu[emunum].sfrac(edata.geda(1));}
       if(emunum<0) std::cout<<'\n';
//        edata.shex(1);
        return xi;}
        
  std::string& gname(){return name;}
 
  void scon(double *vec,int ni){ vec[ni]=conc;}
  
  int markinit(std::string mname,std::string mark,double mval) {
     if(mname==name){  for(int i=0;i<nemus;i++) emu[i].markinit(mark,mval);
      std::cout<<name<<" mark="<<mark<<" %"<<mval<<'\n'; return 1; } 
           return 0;}
                           
//  double getxi(){ if(ieref == -1) return 0.;
//                   else return edata.getxi(gemu(ieref)->sfrac(),gemu(ieref)->gniso());}
  
  double gcon(){return conc;}
  Emu* gemu(int iem){return &emu[iem];}
  int siso(double *vec,int ni){for(int i=0;i<nemus;i++) ni=emu[i].siso(vec,ni); return ni;}
  int sdiso(double *vec,int ni){for(int i=0;i<nemus;i++) ni=emu[i].sdiso(vec,ni); return ni;}
  void volume(double vol){  if(name.at(0)!='E') for(int i=0;i<nemus;i++) emu[i].volume(vol);}
  
  int fndemu(std::string sss){
       for(int i=0;i<nemus;i++) if(emu[i].cmname(sss)) {ieref=i; break;}  return ieref;}
       
  int findemu(std::string exemu){  int emunum(-1);
       for(int i=0;i<nemus;i++) if(emu[i].gname()==exemu)  {emunum=i; break;}
         return emunum;}
       
  void sinit(){for(int i=0;i<nemus;i++) emu[i].sinit(conc);}
  void concor(double& a){ a=conc= emu[0].sumt();}
  void sumt(){for(int i=0;i<nemus;i++) emu[i].sumt();}
  void readc0(std::ifstream& fi) {std::string aaa;
    fi>>aaa>>conc; if(aaa!=name) std::cout<<aaa<<"name conflict!"<<name<<std::endl;}
Metab(){ieref=-1;}
~Metab(){  delete[] emu;}
};

class Ldistr {
 int hkf, hkr, pfk, fbase, aldf, aldr, t3pep, pept3, pk, pyrlac, lacpyr, citdr, citdf, csyn, akgcit, citakg, akgdr, akgdf, citoxc, liase, ppp, pdh, pc, malic, malicc, pyrdf, pyrdr, maloa, oamal, oadr, oadf, akgfum, glnin, glnout, gluin, gluout, tkp5k, tkh6k, tks7k, tkt3a, tke4a, tkp5a, tas7k, tah6k, tae4a, tat3a, nre;
 int Eglc, h6p, fbp, t3p, pep, pyrc, Elac, citc, cit, oa, accoa, akgc, co2, akg, oac, p5, pyr, fumal, Egln, Eglu, gae, e4p, s7p, dhe, nmet;
 int ntime, len; //number of timepoints, accounted mass isotopomers
 Metab *met;
 Reakcia *rr;
 double  tex[tt], Vi, Vt, mu, mval;
 std::string mname, mark;
public:
 void setmetrea();
 void f(const double *y,double *dydx);
 void ff(const double *y,double *dydx);
 void mdistr(double *py,double *pdydt,double t);
 void markinit(){
        for(int i=0;i<nmet;i++) if(met[i].markinit(mname,mark,mval)) break;
           }
 void rpar(std::string fn){std::ifstream fi(fn.c_str());
        for(int i=0;i<nre;i++)     rr[i].rpar(fi);      // read parameters
            for(int i=0;i<nmet;i++) met[i].readc0(fi);  // read concentrations
       fi>>mname>>mname>>mark>>mval;                    // labeled substrate
             }
 void read(std::string fn){ std::ifstream fi(fn.c_str());
     len=nmet; for(int i=0;i<nmet;i++) len+=met[i].read(fi); std::cout<<"; len "<<len<<std::endl;}
 double readex(std::string fn,int itp);
 int glen() { return len;}
 int gnmet() { return nmet;}
 int scon(double *vec){for(int i=0;i<nmet;i++) met[i].scon(vec,i); return nmet; }
 void gcon(){for(int i=0;i<nmet;i++) std::cout<<met[i].gname()<<"="<<met[i].gcon()<<" "; std::cout<<std::endl; }
 void concor(double *vec){for(int i=0;i<nmet;i++) met[i].concor(vec[i]);  }
 void siso(double *vec,int& ni){for(int i=0;i<nmet;i++) ni=met[i].siso(vec,ni); }
 void sdiso(double *vec,int& ni){for(int i=0;i<nmet;i++) ni=met[i].sdiso(vec,ni); }
 void volume(double vol){for(int i=0;i<nmet;i++) met[i].volume(vol); }
//  double getxi(){ double sum=0.; for(int i=0;i<nmet;i++) sum += met[i].getxi();  return sum;}
 void sinit(){for(int i=0;i<nmet;i++) met[i].sinit(); }
 void sumt(){for(int i=0;i<nmet;i++) met[i].sumt(); }
 void shiso(double *vec){double xi(0.);std::cout<<"** Data experiment-calculation (emu: m0,m1,...): **\n";
    for(int i=0;i<nmet;i++) xi += met[i].shiso(vec); 
    std::cout<<"xi="<<xi<<std::endl;}
 void showflx(){for(int i=0;i<nre;i++) std::cout<<rr[i].gnaz()<<" "<<rr[i].gflx()[0]<<"; "; std::cout << std::endl;}
 void tsolve(const double tmax);
 void tisolve(const double tmax);
 double ddisolve();
 double integrbs();
Ldistr(){}
~Ldistr(){delete[] met; delete[] rr; }
};

//---------------------------------------------------------------------------
#endif


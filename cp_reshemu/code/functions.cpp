#include <iostream>
#include <cmath>
#include "new.h"
using namespace std;
double xribi, xasp=1.;
extern double dt;
double Ldistr::readex (string fn,int itp) {
  string aaa, scar;
  ifstream fi(fn.c_str()); double Ti,ts1;  mu=0.;
//  start reading
   fi>> dt; getline(fi,aaa);      tex[0]=0.;  
      for(ntime=1;;ntime++) {fi>>tex[ntime]; tex[ntime] *= 60.; if(tex[ntime]<0) break;}//time points
        fi>> Vi; getline(fi,aaa); //third line
 double Nc[ntime]; for(int i=0;i<ntime;i++) fi>>Nc[i]; getline(fi,aaa);//cells number
        for(int i=1;i<ntime;i++) mu += log(Nc[i]/Nc[0])/tex[i];
    mu /= ((double)(ntime-1)*1.0); getline(fi,aaa); // exponential growth of cell number
	 cout<<"Edata: ";
   while(!fi.eof()){ fi>>aaa;  //reading metabolite name, positions of carbons
     if(fi.eof()) break;
   int i; for(i=0;i<nmet;i++){ // find metabolite in model corresponding to data
    if(aaa.find(met[i].gname())+1) { cout<<" "<<met[i].gname();
       met[i].getex()->read(fi,itp); // set edata for time point itp
        met[i].getex()->setex0(); //set initial m0=1, and m1,...=0
         break;}
                             }
      if(i==nmet){getline(fi,aaa);getline(fi,aaa);} //if no correspondence to model metabolites, go further
    } cout<<endl;
    if(itp==1) met[Eglc].getex()->setrav(0,1);   //set initially labeled substrate
return ts1;}

 void Ldistr::rpar(std::string fn){std::ifstream fi(fn.c_str());
        for(int i=0;i<nre;i++) rr[i].rpar(fi);      // read parameters
                              int ipar; vpar.clear(); 
        for(;;) {fi>>ipar; if(ipar<0) break; vpar.push_back(ipar);} 
                              string aaa; getline(fi,aaa);
        for(int i=0;i<nmet;i++) met[i].readc0(fi);  // read concentrations
        fi>>mname>>mname>>mark>>mval;               // labeled substrate
             }
             
 void Ldistr::read(std::string fn){ std::ifstream fi(fn.c_str()); //read "metemu" 
     len=nmet; for(int i=0;i<nmet;i++) len+=met[i].read(fi);
      std::cout<<"; len "<<len<<std::endl;
       xinit=new double[len]; xx=new double[len];}
       
 int Metab::read(std::ifstream& fi){ std::string aaa; len=0;
   fi>>name>>nemus; emu=new Emu[nemus]; // std::cout<<name<<" "<<nemus<<std::endl;
    for(int i=0;i<nemus;i++) {fi>>aaa; len+=emu[i].set(aaa);} 
     return len; }
     
 void Edata::read(std::ifstream& fi,int nt){
	std::string aaa;//i: isotopomer; nt: time; imt: a number assigned to isotopomer
	fi>>sfrg; niso=sfrg.length()+1;
       emiso[nt]=new data[niso];
   for(int i=0;i<niso;i++) fi>>emiso[nt][i].mean;  fi>>aaa; fi>>aaa; // std::cout<<aaa<<" "; 
   for(int i=0;i<niso;i++) {fi>>emiso[nt][i].sd; if(emiso[nt][i].sd<0.01) emiso[nt][i].sd=0.01;}
	}
	
 void Metab::readc0(std::ifstream& fi) {std::string aaa;
    fi>>aaa>>conc; if(aaa!=name) std::cout<<aaa<<"name conflict!"<<name<<std::endl;}

 void Edata::readc(std::ifstream& fi,  int nt){ std::string aaa;
     fi>>aaa;
     for(int i=0;i<nt;i++) fi>>econc[i].mean;
     for(int i=0;i<nt;i++) {fi>>econc[i].sd; if(econc[i].sd<0.01) econc[i].sd=0.01;}
	}
	
 void Edata::setrav(int nt1,int nt2){// sets edata at t=0
      for(int i=0;i<=niso;i++){ emiso[nt1][i].mean=emiso[nt2][i].mean; emiso[nt1][i].sd=emiso[nt2][i].sd;  }
                         }

void Edata::setex0(){ emiso[0]=new data[niso];
     emiso[0][0].mean=1.; emiso[0][0].sd=0.01;
      for(int i=1;i<niso;i++){ emiso[0][i].mean=0.; emiso[0][i].sd=0.01;  }
                         }

 double Emu::sfrac(data *eda){
     std::cout<<name<<": "; 
     for(int i=0;i<niso;i++) { double fr=iso[i]/sumiso;
         std::cout<<eda[i].mean<<"-"<<fr<<", ";}//<<" xi:"<<xi
     std::cout<<" sum="<<sumiso<<"; xi="<<xi<<"\n";   
  return xi; }

 double Emu::chisq(data *eda){ double a; xi=0.;  sumt();
     for(int i=0;i<niso;i++) { double fr=iso[i]/sumiso;
         if(eda[i].sd>0.009) {a=(fr-eda[i].mean)/eda[i].sd; xi += a*a;}
         }//<<" xi:"<<xi
    return xi;}

 void Emu::markinit(std::string mark, double mval){  // mark: positions containing ¹³C
   int l=mark.length(), nmark=0;                // nmark: number of ¹³C in a givem EMU
      for(int i=0;i<l;i++) {  int npos=name.find(mark.at(i)); // check if the EMU contains the labeled position
               if(npos+1) nmark++;      }
      if(nmark){ iso[nmark]=iso[0]*mval; iso[0] -= iso[nmark];}
           }

 void Reakcia::rpar(std::ifstream& fi){ std::string aaa;
    fi>>aaa>>naz>>npar; par=new double[npar];
      for(int i=0;i<npar;i++) fi>>par[i]; Vm=par[0];
        }

  double Metab::shiso(){ double xi(0); int emunum;//show  EMUs
      if(edata.gniso()){std::cout<<name<<":\t"<<edata.getsfrg()<<"-";
        emunum=findemu(edata.getsfrg()); if(emunum>=0) xi=emu[emunum].sfrac(edata.geda(1));}
       if(emunum<0) std::cout<<'\n';
        return xi;}
        
  double Metab::chisq(){ double xi(0); int emunum;//show  EMUs
      if(edata.gniso()){
        emunum=findemu(edata.getsfrg()); if(emunum>=0) xi=emu[emunum].chisq(edata.geda(1));}
        return xi;}
        
  int Metab::markinit(std::string mname,std::string mark,double mval) {
     if(mname==name){  for(int i=0;i<nemus;i++) emu[i].markinit(mark,mval);
      std::cout<<name<<" mark="<<mark<<" %"<<mval<<'\n'; return 1; } 
           return 0;}
                           
  int Metab::fndemu(std::string sss){
       for(int i=0;i<nemus;i++) if(emu[i].cmname(sss)) {ieref=i; break;}  return ieref;}
       
  int Metab::findemu(std::string exemu){  int emunum(-1);
       for(int i=0;i<nemus;i++) if(emu[i].gname()==exemu)  {emunum=i; break;}
         return emunum;}
       
 void Ldistr::markinit(){
        for(int i=0;i<nmet;i++) if(met[i].markinit(mname,mark,mval)) break;
           }

 void Ldistr::shiso(){double xi(0.);std::cout<<"** Data experiment-calculation (emu: m0,m1,...): **\n";
    for(int i=0;i<nmet;i++) xi += met[i].shiso(); 
    }
 double Ldistr::chisq(){double xi(0.);
    for(int i=0;i<nmet;i++) xi += met[i].chisq(); 
    return xi;}

 void Ldistr::gcon(){std::cout<<"*** Metabolite concentrations: ***\n";
      for(int i=0;i<nmet;i++) std::cout<<met[i].gname()<<"="<<met[i].gcon()<<" "; std::cout<<std::endl; }



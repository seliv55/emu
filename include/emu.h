//---------------------------------------------------------------------------
#include <fstream> 
#include <vector>

#ifndef emuh
#define emuh
 extern std::ostringstream ostr;
 extern int dlinemu;
 
class Metab {
 std::string imya, emu[22];
  int nat, nemus;// # of C atoms, # of emus 
public:
  void read(std::ifstream& fi,std::string& aa){fi>>imya>>nat>>aa;}
//
  void winit(std::ostringstream& os) {os<<imya<<" 1.0\n";}
//   
  std::string& getimya(){return imya;}
//      
  int getnum(std::vector<Metab*>& mlst ){ int i(0), nmax=mlst.size();
   for(i=0;i<nmax;i++) if(mlst[i]->getimya()==this->getimya()) return i;
      return i; }
//      
  void admet(std::vector<Metab*>& met){int nmet=met.size();
   int j=getnum(met); if(!(j-nmet)) {Metab *a=this; met.push_back(a);} }
//
  int getnat(){return nat;}
  void showmet(std::ostringstream& sfo){ sfo<<imya<<" "<<nemus;
    for(int j=0;j<nemus;j++) sfo<<" "<<emu[j]; sfo<<std::endl;}
//
  int ademu(std::string chemu);
//
  std::string symm(int je);
  void symal() {for(int i=0;i<nemus;i++){std::string a=symm(i); ademu(a); }}
//
  void setic(std::string eobj,std::vector<Metab*>& oblst){
   admet(oblst); ademu(eobj); std::cout<<" obj"<<nemus<<" "<<emu[nemus-1]<<std::endl;}
//
std::string stemu(Metab *su){return "met["+su->getimya()+"].gemu(";}
std::string strr(std::string& rimya){return "(rr["+rimya+"].gflx()[";}
  void equal(Metab *su) { for(int j=0;j<this->nemus;j++) {int i=su->ademu(this->emu[j]);}}
  void equal_rr(std::string& rimya,Metab *su) { for(int j=0;j<nemus;j++) {
          int i=su->ademu(emu[j]);
      ostr<<stemu(su)<<i<<")->uni"<<strr(rimya)<<"1],"<<stemu(this)<<j<<("));// e: "+su->emu[i]+" -> "+emu[j]+";\n");}}
//
  std::string split(int je,int pat[],int sat[], int star);
  void splital(int *pat,int *sat,Metab* su,int star=0){  for(int j=0;j<nemus;j++){
         std::string a=this->split(j,pat,sat,star); int i=su->ademu(a);}}

  void split_rr(std::string& rimya,int *pat,int *sat,Metab* su,int star=0){
    for(int j=0;j<su->nemus;j++){int i; for(i=0;i<nemus;i++){std::string a=split(i,pat,sat,star); if(su->emu[j]==a) break;}
      if(i<nemus) ostr<<stemu(su)<<j<<")->uni"<<strr(rimya)<<"1],"<<stemu(this)<<i<<("));// s: "+su->emu[j]+" -> "+emu[i]+";\n");
             else ostr<<stemu(su)<<j<<")->zero"<<strr(rimya)<<"1]);// s: "+su->emu[j]+" -> ;\n";}}
//
  void condens(int je,std::string a[],int pat[],int sat[],int nas1);
  void condensal(int *pat,int *sat,Metab *mets1,Metab *mets2,std::string a[]){
   int nas1=mets1->getnat();
    for(int j=0;j<nemus;j++){ this->condens(j,a,pat,sat,nas1);
   mets1->ademu(a[0]); mets2->ademu(a[1]); }}
  void condens_rr(std::string& rimya,int *pat,int *sat,Metab *mets1,Metab *mets2,std::string a[]){
   int nas1=mets1->getnat(),i1st[5],i2st[5];
    int i=0; for(int j=0;j<nemus;j++){ this->condens(j,a,pat,sat,nas1);
   int i1=mets1->ademu(a[0]); int i2=mets2->ademu(a[1]);
    if((i1>=0)&&(i2>=0)) { if(mets2->getimya()=="co2a") ostr<<stemu(mets1)<<i1<<")->uni"<<strr(rimya)<<"1],"<<stemu(this)<<j<<"));// c: "<<mets1->emu[i1]<<" + "<<mets2->emu[i2]<<" -> "<<emu[j]<<";\n";
                      else ostr<<stemu(mets1)<<i1<<")->bi"<<strr(rimya)<<"3],"<<stemu(mets2)<<i2<<"),"<<stemu(this)<<j<<"));// c: "<<mets1->emu[i1]<<" + "<<mets2->emu[i2]<<" -> "<<emu[j]<<";\n"; i1st[i]=i1; i2st[i]=i2; i++;}
    else if(i1>=0) {int k; for(k=0;k<i;k++) if(!(i1-i1st[k])) break;
      if(i-k) ostr<<stemu(mets1)<<i1<<")->uni1"<<strr(rimya)<<"1],"<<stemu(this)<<j<<"));// c: "<<mets1->emu[i1]<<" -> "<<emu[j]<<";\n";
      else  ostr<<stemu(mets1)<<i1<<")->uni"<<strr(rimya)<<"1],"<<stemu(this)<<j<<"));// c: "<<mets1->emu[i1]<<" -> "<<emu[j]<<";\n";}
    else if(i2>=0) {int k; for(k=0;k<i;k++) if(!(i2-i2st[k])) break;
      if(i-k) ostr<<stemu(mets2)<<i2<<")->uni1"<<strr(rimya)<<"2],"<<stemu(this)<<j<<"));// c: "<<mets2->emu[i2]<<" -> "<<emu[j]<<";\n";
        else  ostr<<stemu(mets2)<<i2<<")->uni"<<strr(rimya)<<"2],"<<stemu(this)<<j<<"));// c: "<<mets2->emu[i2]<<" -> "<<emu[j]<<";\n";}
    }}
//
  Metab(){imya=""; nat=0; nemus=0;}
  ~Metab(){}
};
/***/
class React {
 int *sat, *pat,npar,ns;
  std::string rname,rtip;
   Metab s[4],*pr;
public:
 bool chekrob(Metab* a,int ip=0){ return (pr[ip].getimya()==a->getimya()); }
//
 Metab& getsu(int ns=0){return s[ns];}
//
 int* getpat(){return pat;}
 int* getsat(){return sat;}
 std::string& getrname(){return rname;}
//
 void condensal(Metab* pr,Metab* su1,Metab* su2,std::string a[]){ pr->condensal(pat,sat,su1,su2,a);}
//
 void admetrr(std::vector<Metab*>& met){
       s[0].admet(met);
       pr[0].admet(met);
   if(s[1].getnat()) s[1].admet(met);
   else if(pr[1].getnat()) pr[1].admet(met);
 }
//
  void read(std::ifstream& fi);
  void wreq(std::ostringstream& seq,std::ostringstream& sff);
   void wpar(std::ostringstream& seq,int i);
  React(){pr=&s[2];}
  ~React(){}
};
//---------------------------------------------------------------------------
#endif


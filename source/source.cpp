#include <iostream> 
#include <sstream> 
#include <string>
#include "emu.h"
using namespace std;

inline int chatoin(char a){return (int)a-(int)'0';}
//
inline void sort(int len,int a[]){int aa; for(int i=1;i<len;i++)
  for(int j=0;j<i;j++) if(a[i]<a[j]) {aa=a[j]; a[j]=a[i]; a[i]=aa;}}
//
 std::string Metab::symm(int je) {
  int a[nat], len=emu[je].length();
   for(int i=0;i<len;i++) for(int j=0;j<nat;j++) if(chatoin(emu[je].at(i))==j) {
     a[i]=nat-j-1; break;} sort(len,a);
    std::ostringstream lst; for(int i=0;i<len;i++) lst<<a[i];
   return lst.str();}
//            
 int Metab::ademu(std::string chemu){ int i;
  for(i=0;i<this->nemus;i++) if(chemu==this->emu[i]) break;
   if(!chemu.length()) i=-1;
    else if(!(this->nemus-i)){
        this->emu[this->nemus]=chemu; this->nemus++;}
          return i;}
//
 string Metab::split(int je,int pat[],int sat[], int star) {
  int len=this->emu[je].length(), a[this->nat];
    for(int i=0;i<len;i++) 
     for(int k=0;k<this->nat;k++)
      if(chatoin(this->emu[je].at(i))==pat[k+star]) {a[i]=sat[k+star]; break;}
                     sort(len,a);
   ostringstream lst; for(int l=0;l<len;l++) lst<<a[l];
    return lst.str();}
//
 void Metab::condens(int je,std::string a[],int pat[],int sat[],int nas1) {
     int len=emu[je].length(), as1[len], as2[len],j(0),l(0);
      for(int i=0;i<len;i++)
       for(int k=0;k<nat;k++) if(chatoin(emu[je].at(i))==pat[k]) {
              if(k<nas1) {as1[j]=sat[k]; j++; break;} else {as2[l]=sat[k]; l++; break;}}
                     sort(j,as1); sort(l,as2); //sorting
                 ostringstream lst,lst1; 
   for(int k=0;k<j;k++) lst<<as1[k]; for(int k=0;k<l;k++) lst1<<as2[k];
      a[0]=lst.str(); a[1]=lst1.str(); }
//
  void React::read(ifstream& fi) { string aaa;
  fi>>rname; s[0].read(fi,aaa);
    if(aaa=="+") { s[1].read(fi,aaa); pr[0].read(fi,aaa); rtip="cond";//condensation
     int nat=pr[0].getnat();
      sat=new int[nat]; pat=new int[nat];
        for(int i=0;i<nat;i++) fi>>sat[i]>>aaa>>pat[i];}
    else {pr[0].read(fi,aaa); rtip=aaa;
          if(aaa=="+") {  pr[1].read(fi,aaa); rtip="spli";//split
            int nat=s[0].getnat();
             sat=new int[nat]; pat=new int[nat];
               for(int i=0;i<nat;i++) fi>>sat[i]>>aaa>>pat[i]; }
          }
      }
   void React::wreq(ostringstream& seq,ostringstream& sff){ string dydx="", su=s[0].getimya(), pp=pr[0].getimya(), dff="";
    if(su.at(0)=='E') {dydx="dydx["+pp+"] += rr["+rname+"].v(y["+su+"])[0];"; dff="dydx["+su+"] -= rr["+rname+"].v(y["+su+"])[0];\n";}
     else {dydx="dydx["+su+"] -= rr["+rname+"].v(y["+su+"])[0];";
    if(pp.at(0)=='E') dff="dydx["+pp+"] += rr["+rname+"].v(y["+su+"])[0];\n";   else dydx +=" dydx["+pp+"] += rr["+rname+"].gflx()[0];";}
     npar=2; ns=1;
   if(rtip=="spli") {dydx += "dydx["+pr[1].getimya()+"] += rr["+rname+"].gflx()[0];";}
   else if(rtip=="cond") { npar++; ns++;
      dydx="dydx["+s[0].getimya()+"] -= rr["+rname+"].v(y["+s[0].getimya()+"],y["+s[1].getimya()+"])[0]; dydx["+pr[0].getimya()+"] += rr["+rname+"].gflx()[0]; dydx["+s[1].getimya()+"] -= rr["+rname+"].gflx()[0];";}
    seq<<dydx<<"\n"; sff<<dff;}
    
   void React::wpar(ostringstream& seq,int nr){
    seq<<nr<<"   "<<rname<<" par "<<npar; for(int i=0;i<npar;i++) seq<<" 0.1"; seq<<"\n";}
   


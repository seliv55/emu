#include <iostream> 
#include <sstream> 
#include <string>
#include "emu.h"
using namespace std;
        ostringstream ostr;
        
 int main () {
       vector<Metab*> met, oblst;
    ifstream fi("model"); string aaa; 
        vector<React> rr;
  // read and form lists of reactions & metabolites:
  ostringstream kin,kin0;
   kin0<<"void Ldistr::ff(const double *y,double *dydx) {\n";
   kin<<"void Ldistr::f(const double *y,double *dydx) { for(int i=0;i<nmet;i++) dydx[i]=0.;\n";
   while(true) {React* r0=new React; r0->read(fi); if(fi.eof()) break; r0->admetrr(met); r0->wreq(kin,kin0); rr.push_back(*r0);}
   kin0<<"}\n\n";
   kin<<"for(int i=0;i<nmet;i++) if(met[i].gname().at(0) != 'E') dydx[i] /= Vt;}\n\n";

  int rnum=rr.size();
   cout<<"metabolites: "<<met.size()<<"; reactions"<<rnum<<endl;
   string rlist="\nvoid Ldistr::setmetrea(){ "+rr[0].getrname()+"=0, "; 
  for(int i=1;i<rnum;i++) rlist += rr[i].getrname()+"="+rr[i-1].getrname()+"+1, ";
   rlist+="nre="+rr[rnum-1].getrname()+"+1;\n "+met[0]->getimya()+"=0, ";
  for(int i=1;i<met.size();i++) rlist += met[i]->getimya()+"="+met[i-1]->getimya()+"+1, ";
   rlist+="nmet="+met[met.size()-1]->getimya()+"+1;\nrr=new Reakcia[nre];  met=new Metab[nmet]; cout<<\"reactions \"<<nre<<\"; metabolites \"<<nmet;}\n\n//  file 'include/new.h':\n// int "; 
  for(int i=0;i<rnum;i++) rlist += rr[i].getrname()+", "; rlist+="nre;\n// int "+met[0]->getimya()+", ";
  for(int i=1;i<met.size();i++) rlist += met[i]->getimya()+", ";rlist+="nmet;\n\n";
// save the formed lists:
    ofstream fo("reshemu/code/nv.cpp");
    fo<<"#include <iostream>\n#include \"new.h\"\nusing namespace std;\n\tdouble dt;\n";
     fo<<rlist<<kin0.str()<<kin.str(); fo.close();

   for(int i=0;i<met.size();i++) { aaa=met[i]->getimya(); // define the emu that are measured:
     cout<<i<<" "<<aaa<<"; ";
      if(aaa=="Elac")  met[i]->setic("012",oblst); 
      if(aaa=="Eglu")  met[i]->setic("1234",oblst); 
      if(aaa=="p5")  met[i]->setic("01234",oblst); 
     if(aaa=="Egln")   met[i]->setic("01234",oblst); }
      int nc(7);// number of passages through the reaction scheme
   // emu to follow
  for(int ij=0;ij<nc;ij++){
   int k=0, ipr, isu;
   while (k<oblst.size()) {  //adding new metabolites to the list of objects to follow
 for(int i=0;i<rnum;i++) { //search through all the products of all the reactions
  if(rr[i].chekrob(oblst[k])) { //search through first products of the reactions
          ipr=oblst[k]->getnum(met); // number of the product in the list "met"
          isu=rr[i].getsu().getnum(met);  // number of the substrate in the list "met"
     met[isu]->admet(oblst); //add substrate to the list of objects to follow
          
       if(met[ipr]->getimya()=="fumal") met[ipr]->symal(); // accounting fot symmetry
       
 if(!(rr[i].getsu(1).getnat()+rr[i].getsu(3).getnat())) {       // equal
  met[ipr]->equal(met[isu]);   //makes sure that all emus of product is presented in substrate
   if(!(ij-nc+1)) met[ipr]->equal_rr(rr[i].getrname(),met[isu]);//construct ODEs corresponding to the reaction.
    }
       
 if(rr[i].getsu(3).getnat()) {                                  // split
   met[ipr]->splital(rr[i].getpat(),rr[i].getsat(),met[isu]); // construct emu of substrate corresponding to that of product.
    if(!(ij-nc+1)) met[ipr]->split_rr(rr[i].getrname(), rr[i].getpat(),rr[i].getsat(),met[isu]);} 
     
 if(rr[i].getsu(1).getnat()){ int isu2=rr[i].getsu(1).getnum(met);  // condence
     met[isu2]->admet(oblst);      string a[2];
   met[ipr]->condensal(rr[i].getpat(),rr[i].getsat(),met[isu],met[isu2],a);// construct emus of su corresponding to that of pr
    if(!(ij-nc+1)) met[ipr]->condens_rr(rr[i].getrname(),rr[i].getpat(),rr[i].getsat(),met[isu],met[isu2],a);}
            }
             if(rr[i].chekrob(oblst[k],1)) {//if the second product exists, add it to the list and construct emus
          ipr=oblst[k]->getnum(met);
          isu=rr[i].getsu().getnum(met);
       if(rr[i].getsu(3).getnat()) {
        met[ipr]->splital(rr[i].getpat(),rr[i].getsat(),met[isu],rr[i].getsu(2).getnat()); }
            }
       }
                   k++;}
  }
      fo.open("reshemu/asdf.cpp"); fo<<"#include <iostream>\n#include <cmath>\n#include \"new.h\"\nusing namespace std;\nvoid Ldistr::mdistr(double *py,double *pdydt,double t) {\n\tint nx=nmet, ni=nmet;\n//	siso(py,ni);\n\tsdiso(pdydt,nx);\n\tVt=Vi*exp(mu*t);\n\tconcor(py);\n\tf(py,pdydt);\n\tff(py,pdydt);\n";
      fo<<ostr.str()<<"\tvolume(Vt); }\n"; fo.close();
//       for(int i=0;i<obnum;i++) { oblst[i]->showmet(); }/**/
    ostringstream sfo;
       for(int i=0;i<met.size();i++) { met[i]->showmet(sfo); }/**/
    ostringstream parol;
       for(int i=0;i<rnum;i++) { rr[i].wpar(parol,i); }/**/
       for(int i=0;i<met.size();i++) { met[i]->winit(parol); }/**/
       fo.open("par"); fo<<parol.str(); fo.close();
       fo.open("reshemu/metemu"); fo<<sfo.str(); fo.close(); cout<<'\n';
  return 0;
}

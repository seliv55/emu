#include <iostream> 
#include <sstream> 
#include <fstream> 
#include <string>
#include "nr.h"
#include "solvers.h"
#include "new.h"
using namespace std;
  Ldistr emudyn;
  extern string foc, kin0, kin, kinflx;

void Ldistr::setinit(char* fimod,char *fipar,char* fiex,int ntp){
     setmetrea();//code in nv.cpp: assign numbers for metabolites and reactions
  read(fimod);//define emus
   rpar(fipar);// read reaction constants, initial concentrations of metabolites, labeled substrate
   readex(fiex,ntp);
   int ni=scon(xinit); 
   siso(xinit,ni);
       sinit();
       markinit();          //set initial substrate labeling
  }

int main () {
 emudyn.setinit();
//cout<<"files="<<emudyn.chekifn()<<'\n';
  double xi=emudyn.ddisolve();
//   double xi=emudyn.integrbs();
//   emudyn.tsolve(100.);
   emudyn.shiso();
	emudyn.gcon();
	emudyn.wpar(1);
	emudyn.coord(0.1,1.05);
return 0;}

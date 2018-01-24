#include <iostream> 
#include <sstream> 
#include <string>
#include "nr.h"
#include "solvers.h"
#include "new.h"
using namespace std;
  Ldistr emudyn; int  ifn=0;
    time_t ts,tf; string fex1="../files/mglc", fex2="../files/mglc";
  extern string foc, kin0, kin, kinflx;
  void Ldistr::setinit(){
     setmetrea();//code in nv.cpp: assign numbers for metabolites and reactions
  read("metemu");//define emus
   rpar("1");// read reaction constants, initial concentrations of metabolites, labeled substrate
   readex("mglc",1);
   int ni=scon(xinit); 
   siso(xinit,ni);
       sinit();
       markinit();          //set initial substrate labeling
  }

int main () {
   char fn[15];
   for(int i=1;;i++) { sprintf(fn,"%i",i);
	   ifstream checkfi(fn);
	   if(!checkfi.good()) { ifn=i-1; break;}
	   checkfi.close();
   }
 emudyn.setinit();
  time_t ts=clock(); 
   emudyn.ddisolve();
//   emudyn.integrbs();
//   emudyn.tsolve(100.);
   emudyn.shiso(emudyn.getxx());
   tf= (clock()-ts);
	emudyn.gcon();
	 cout<<"time="<< tf / (double) CLOCKS_PER_SEC<<" s\n"; 
return 0;}

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

int main () {
   char fn[15];
   for(int i=1;;i++) { sprintf(fn,"%i",i);
	   ifstream checkfi(fn);
	   if(!checkfi.good()) { ifn=i-1; break;}
	   checkfi.close();
   }

     emudyn.setmetrea();//code in nv.cpp: assign numbers for metabolites and reactions
  emudyn.read("metemu");//define emus
   emudyn.rpar("1");// read reaction constants, initial concentrations of metabolites, labeled substrate
   emudyn.readex("mglc",1);
  time_t ts=clock(); 
   emudyn.ddisolve();
//   emudyn.tsolve(100.);
//   emudyn.integrbs();
   tf= (clock()-ts);
	emudyn.gcon();
	 cout<<"time="<< tf / (double) CLOCKS_PER_SEC<<" s\n"; 
return 0;}

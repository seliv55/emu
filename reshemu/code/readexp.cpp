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
	 bool a=0; cout<<"Edata: ";
   while(!a){ fi>>aaa;  //reading metabolite name, positions of carbons
     if(a=fi.eof()) break;
   int i; for(i=0;i<nmet;i++){ // find metabolite in model corresponding to data
      int imatch=aaa.find(met[i].gname());
    if(imatch+1) { cout<<" "<<met[i].gname();
       met[i].getex()->read(fi,itp); // set edata for time point itp
        met[i].getex()->setex0(); //set initial m0=1, and m1,...=0
         break;}
                             }
      if(i==nmet){getline(fi,aaa);getline(fi,aaa);} //if no correspondence to model metabolites, go further
    } cout<<endl;
    if(itp==1) met[Eglc].getex()->setrav(0,1);   //set initially labeled substrate
return ts1;}

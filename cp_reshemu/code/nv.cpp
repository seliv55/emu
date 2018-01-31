#include <iostream>
#include "new.h"
using namespace std;
	double dt;

void Ldistr::setmetrea(){ hkf=0, hkr=hkf+1, pfk=hkr+1, fbase=pfk+1, aldf=fbase+1, aldr=aldf+1, t3pep=aldr+1, pept3=t3pep+1, pk=pept3+1, pyrlac=pk+1, lacpyr=pyrlac+1, citdr=lacpyr+1, citdf=citdr+1, csyn=citdf+1, akgcit=csyn+1, citakg=akgcit+1, akgdr=citakg+1, akgdf=akgdr+1, citoxc=akgdf+1, liase=citoxc+1, ppp=liase+1, pdh=ppp+1, pc=pdh+1, malic=pc+1, malicc=malic+1, pyrdf=malicc+1, pyrdr=pyrdf+1, maloa=pyrdr+1, oamal=maloa+1, oadr=oamal+1, oadf=oadr+1, akgfum=oadf+1, glnin=akgfum+1, glnout=glnin+1, gluin=glnout+1, gluout=gluin+1, tkp5k=gluout+1, tkh6k=tkp5k+1, tks7k=tkh6k+1, tkt3a=tks7k+1, tke4a=tkt3a+1, tkp5a=tke4a+1, tas7k=tkp5a+1, tah6k=tas7k+1, tae4a=tah6k+1, tat3a=tae4a+1, nre=tat3a+1;
 Eglc=0, h6p=Eglc+1, fbp=h6p+1, t3p=fbp+1, pep=t3p+1, pyrc=pep+1, Elac=pyrc+1, citc=Elac+1, cit=citc+1, oa=cit+1, accoa=oa+1, akgc=accoa+1, co2=akgc+1, akg=co2+1, oac=akg+1, p5=oac+1, pyr=p5+1, fumal=pyr+1, Egln=fumal+1, Eglu=Egln+1, gae=Eglu+1, e4p=gae+1, s7p=e4p+1, dhe=s7p+1, nmet=dhe+1;
rr=new Reakcia[nre];  met=new Metab[nmet]; cout<<"reactions "<<nre<<"; metabolites "<<nmet;}

//  file 'include/new.h':
// int hkf, hkr, pfk, fbase, aldf, aldr, t3pep, pept3, pk, pyrlac, lacpyr, citdr, citdf, csyn, akgcit, citakg, akgdr, akgdf, citoxc, liase, ppp, pdh, pc, malic, malicc, pyrdf, pyrdr, maloa, oamal, oadr, oadf, akgfum, glnin, glnout, gluin, gluout, tkp5k, tkh6k, tks7k, tkt3a, tke4a, tkp5a, tas7k, tah6k, tae4a, tat3a, nre;
// int Eglc, h6p, fbp, t3p, pep, pyrc, Elac, citc, cit, oa, accoa, akgc, co2, akg, oac, p5, pyr, fumal, Egln, Eglu, gae, e4p, s7p, dhe, nmet;

void Ldistr::ff(const double *y,double *dydx) {
dydx[Eglc] -= rr[hkf].v(y[Eglc])[0];
dydx[Eglc] += rr[hkr].v(y[h6p])[0];
dydx[Elac] += rr[pyrlac].v(y[pyrc])[0];
dydx[Elac] -= rr[lacpyr].v(y[Elac])[0];
dydx[Egln] -= rr[glnin].v(y[Egln])[0];
dydx[Egln] += rr[glnout].v(y[akg])[0];
dydx[Eglu] -= rr[gluin].v(y[Eglu])[0];
dydx[Eglu] += rr[gluout].v(y[akg])[0];
}

void Ldistr::f(const double *y,double *dydx) { for(int i=0;i<nmet;i++) dydx[i]=0.;
dydx[h6p] += rr[hkf].v(y[Eglc])[0];
dydx[h6p] -= rr[hkr].v(y[h6p])[0];
dydx[h6p] -= rr[pfk].v(y[h6p])[0]; dydx[fbp] += rr[pfk].gflx()[0];
dydx[fbp] -= rr[fbase].v(y[fbp])[0]; dydx[h6p] += rr[fbase].gflx()[0];
dydx[fbp] -= rr[aldf].v(y[fbp])[0]; dydx[t3p] += rr[aldf].gflx()[0];dydx[t3p] += rr[aldf].gflx()[0];
dydx[t3p] -= rr[aldr].v(y[t3p],y[t3p])[0]; dydx[fbp] += rr[aldr].gflx()[0]; dydx[t3p] -= rr[aldr].gflx()[0];
dydx[t3p] -= rr[t3pep].v(y[t3p])[0]; dydx[pep] += rr[t3pep].gflx()[0];
dydx[pep] -= rr[pept3].v(y[pep])[0]; dydx[t3p] += rr[pept3].gflx()[0];
dydx[pep] -= rr[pk].v(y[pep])[0]; dydx[pyrc] += rr[pk].gflx()[0];
dydx[pyrc] -= rr[pyrlac].v(y[pyrc])[0];
dydx[pyrc] += rr[lacpyr].v(y[Elac])[0];
dydx[citc] -= rr[citdr].v(y[citc])[0]; dydx[cit] += rr[citdr].gflx()[0];
dydx[cit] -= rr[citdf].v(y[cit])[0]; dydx[citc] += rr[citdf].gflx()[0];
dydx[oa] -= rr[csyn].v(y[oa],y[accoa])[0]; dydx[cit] += rr[csyn].gflx()[0]; dydx[accoa] -= rr[csyn].gflx()[0];
dydx[akgc] -= rr[akgcit].v(y[akgc],y[co2])[0]; dydx[citc] += rr[akgcit].gflx()[0]; dydx[co2] -= rr[akgcit].gflx()[0];
dydx[cit] -= rr[citakg].v(y[cit])[0]; dydx[akg] += rr[citakg].gflx()[0];dydx[co2] += rr[citakg].gflx()[0];
dydx[akgc] -= rr[akgdr].v(y[akgc])[0]; dydx[akg] += rr[akgdr].gflx()[0];
dydx[akg] -= rr[akgdf].v(y[akg])[0]; dydx[akgc] += rr[akgdf].gflx()[0];
dydx[citc] -= rr[citoxc].v(y[citc])[0]; dydx[akgc] += rr[citoxc].gflx()[0];dydx[co2] += rr[citoxc].gflx()[0];
dydx[citc] -= rr[liase].v(y[citc])[0]; dydx[oac] += rr[liase].gflx()[0];dydx[accoa] += rr[liase].gflx()[0];
dydx[h6p] -= rr[ppp].v(y[h6p])[0]; dydx[p5] += rr[ppp].gflx()[0];dydx[co2] += rr[ppp].gflx()[0];
dydx[pyr] -= rr[pdh].v(y[pyr])[0]; dydx[accoa] += rr[pdh].gflx()[0];dydx[co2] += rr[pdh].gflx()[0];
dydx[pyr] -= rr[pc].v(y[pyr],y[co2])[0]; dydx[oa] += rr[pc].gflx()[0]; dydx[co2] -= rr[pc].gflx()[0];
dydx[fumal] -= rr[malic].v(y[fumal])[0]; dydx[pyr] += rr[malic].gflx()[0];dydx[co2] += rr[malic].gflx()[0];
dydx[oac] -= rr[malicc].v(y[oac])[0]; dydx[pyrc] += rr[malicc].gflx()[0];dydx[co2] += rr[malicc].gflx()[0];
dydx[pyrc] -= rr[pyrdf].v(y[pyrc])[0]; dydx[pyr] += rr[pyrdf].gflx()[0];
dydx[pyr] -= rr[pyrdr].v(y[pyr])[0]; dydx[pyrc] += rr[pyrdr].gflx()[0];
dydx[fumal] -= rr[maloa].v(y[fumal])[0]; dydx[oa] += rr[maloa].gflx()[0];
dydx[oa] -= rr[oamal].v(y[oa])[0]; dydx[fumal] += rr[oamal].gflx()[0];
dydx[oac] -= rr[oadr].v(y[oac])[0]; dydx[fumal] += rr[oadr].gflx()[0];
dydx[fumal] -= rr[oadf].v(y[fumal])[0]; dydx[oac] += rr[oadf].gflx()[0];
dydx[akg] -= rr[akgfum].v(y[akg])[0]; dydx[fumal] += rr[akgfum].gflx()[0];dydx[co2] += rr[akgfum].gflx()[0];
dydx[akg] += rr[glnin].v(y[Egln])[0];
dydx[akg] -= rr[glnout].v(y[akg])[0];
dydx[akg] += rr[gluin].v(y[Eglu])[0];
dydx[akg] -= rr[gluout].v(y[akg])[0];
dydx[p5] -= rr[tkp5k].v(y[p5])[0]; dydx[t3p] += rr[tkp5k].gflx()[0];dydx[gae] += rr[tkp5k].gflx()[0];
dydx[h6p] -= rr[tkh6k].v(y[h6p])[0]; dydx[e4p] += rr[tkh6k].gflx()[0];dydx[gae] += rr[tkh6k].gflx()[0];
dydx[s7p] -= rr[tks7k].v(y[s7p])[0]; dydx[p5] += rr[tks7k].gflx()[0];dydx[gae] += rr[tks7k].gflx()[0];
dydx[t3p] -= rr[tkt3a].v(y[t3p],y[gae])[0]; dydx[p5] += rr[tkt3a].gflx()[0]; dydx[gae] -= rr[tkt3a].gflx()[0];
dydx[e4p] -= rr[tke4a].v(y[e4p],y[gae])[0]; dydx[h6p] += rr[tke4a].gflx()[0]; dydx[gae] -= rr[tke4a].gflx()[0];
dydx[p5] -= rr[tkp5a].v(y[p5],y[gae])[0]; dydx[s7p] += rr[tkp5a].gflx()[0]; dydx[gae] -= rr[tkp5a].gflx()[0];
dydx[s7p] -= rr[tas7k].v(y[s7p])[0]; dydx[e4p] += rr[tas7k].gflx()[0];dydx[dhe] += rr[tas7k].gflx()[0];
dydx[h6p] -= rr[tah6k].v(y[h6p])[0]; dydx[t3p] += rr[tah6k].gflx()[0];dydx[dhe] += rr[tah6k].gflx()[0];
dydx[e4p] -= rr[tae4a].v(y[e4p],y[dhe])[0]; dydx[s7p] += rr[tae4a].gflx()[0]; dydx[dhe] -= rr[tae4a].gflx()[0];
dydx[t3p] -= rr[tat3a].v(y[t3p],y[dhe])[0]; dydx[h6p] += rr[tat3a].gflx()[0]; dydx[dhe] -= rr[tat3a].gflx()[0];
for(int i=0;i<nmet;i++) if(met[i].gname().at(0) != 'E') dydx[i] /= Vt;}


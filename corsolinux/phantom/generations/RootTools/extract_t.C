/* In VV scattering 6f events extract_t computes the Mandelstam t between
   vectors
   ptot: total lab final state momentum
   pVf1,2 : momentum of reconstructed final state vector bosons
   ptag1,2 : momentum of tag jets
*/

#if !defined(_CINT_) || defined(_MAKECINT_)
#include <Riostream.h>
#include <stdio.h>
#include "TMath.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#endif

double extract_t(TLorentzVector &ptot, TLorentzVector &pVf1,
                 TLorentzVector &pVf2, TLorentzVector &ptag1,
                 TLorentzVector &ptag2)
{
  double x1TimesEcoll,x2TimesEcoll,t,aux1,aux2;
  TLorentzVector p1,p2,pVi1;
//  TLorentzVector pVi2; 
  
  x1TimesEcoll = 0.5*(ptot.E()+ptot.Pz());
  x2TimesEcoll = 0.5*(ptot.E()-ptot.Pz());

  p1 = TLorentzVector(0.0,0.0,x1TimesEcoll,x1TimesEcoll);
  p2 = TLorentzVector(0.0,0.0,-x2TimesEcoll,x2TimesEcoll);


  aux1=TMath::Abs((p1-ptag1)*(p1-ptag1)+(p2-ptag2)*(p2-ptag2));
  aux2=TMath::Abs((p1-ptag2)*(p1-ptag2)+(p2-ptag1)*(p2-ptag1));
  if(aux1 <= aux2){
     pVi1 = p1-ptag1;
//     pVi2 = p2-ptag2;
  }else{
     pVi1 = p1-ptag2;
//     pVi2 = p2-ptag1;
  }

//  t = TMath::Min(
//        TMath::Abs((pVi1-pVf1)*(pVi1-pVf1)),
//        TMath::Abs((pVi1-pVf2)*(pVi1-pVf2))
                );
  t = TMath::Max(
        (pVi1-pVf1)*(pVi1-pVf1),
        (pVi1-pVf2)*(pVi1-pVf2)
                );
/*
cout<<"p1 "<<endl;
printf("px = %.18g, py= %.18g, pz = %.18g, pe = %.18g\n",
        p1.X(),p1.Y(),p1.Z(),p1.E());
cout<<"p2 "<<endl;
printf("px = %.18g, py= %.18g, pz = %.18g, pe = %.18g\n",
        p2.X(),p2.Y(),p2.Z(),p2.E());
cout<<"pVi1 "<<endl;
printf("px = %.18g, py= %.18g, pz = %.18g, pe = %.18g\n",
        pVi1.X(),pVi1.Y(),pVi1.Z(),pVi1.E());
cout<<"pVi2 "<<endl;
printf("px = %.18g, py= %.18g, pz = %.18g, pe = %.18g\n",
        pVi2.X(),pVi2.Y(),pVi2.Z(),pVi2.E());

*/
  return t;
}

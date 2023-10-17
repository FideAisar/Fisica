// Adapted from integro.c (by S.Bolognesi)

double integro(TH1D *histo)
{
  int nbins;
  double width,integral;
  width = histo->GetBinWidth(1);
  //integral = (histo->Integral())*width;
  nbins=histo->GetNbinsX();
  integral = (histo->Integral(1,nbins,"width"));
  return integral;
}

double integro(TH1D *histo,double low,double high)
{
  double width,min,max,integral;
  int highBin,lowBin;
  width = histo->GetBinWidth(1);
  
  TAxis *xaxis = histo->GetXaxis();
  min=xaxis->GetXmin();
  max=xaxis->GetXmax();
  
  if(low>=min) lowBin=(int)TMath::Ceil((low-min)/width);
  else
    {
      cout<<"*** integro.c ERROR: lower integ. limit .LT.< histogram range"<<endl;
//      return 0;
      exit(0);
    }
  
  if(high<=max) highBin=(int)TMath::Ceil((high-min)/width);
  else 
    {
//      cout<<"limite sup. > del range"<<endl;
      highBin=(int)TMath::Ceil((max-min)/width);
    }

  low = xaxis->GetBinLowEdge(lowBin);
  high = ((xaxis->GetBinLowEdge(highBin))+width);

  cout<<"integrale"<<endl<<"da "<<low<<" (bin "<<lowBin<<")"<<endl
      <<"a "<<high<<" (bin "<<highBin<<")"<<endl; 
  integral= (histo->Integral(lowBin,highBin))*width;
  return integral;
}

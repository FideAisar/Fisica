double integro_nowidth(TH1D *histo,double low,double high)
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
  integral= (histo->Integral(lowBin,highBin));
  return integral;
}


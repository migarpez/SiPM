#include "TSpectrum.h"

const int NMAXFILES = 3;
const int WFLENGTH  = 1250;
bool plot_results   = true;

//**********************************************************
int getMultipleEvents(TH1F* hxt,
		      TLine* tl){
//**********************************************************

  //create tspectrum for finding peak
  TSpectrum *s = new TSpectrum(1);

  //look for single pe value and compute 1.5 pe;
  s->Search(hxt,0.0001,"",0.1);
  double* pe = s->GetPositionX();
  double oneandahalfpe =  3*pe[0]/2;
  tl->SetY1(oneandahalfpe);
  tl->SetY2(oneandahalfpe);

  double multiple = 0;

  for(int i = 0; i < hxt->GetNbinsX(); i++)if(hxt->GetBinCenter(i) > oneandahalfpe)multiple = multiple + hxt->GetBinContent(i);

  return multiple;
}

//**********************************************************
double getBaseline(double* time, double* V){
//**********************************************************

  //fit all wf to an horizontal line and get mean value
  TGraph* tg = new TGraph(WFLENGTH,time,V);
  tg->Fit("pol0","Q");
  //tg->Draw("ap");gPad->WaitPrimitive();
  double baseline = static_cast<TF1*>(tg->GetFunction("pol0"))->GetParameter(0);
  tg->Delete();
  return baseline;
}

//**********************************************************
double getMax(double* time, double* V){
//**********************************************************

  double Vmax = 0;
  for(int i = 0; i < WFLENGTH; i++)if(V[i] > Vmax)Vmax = V[i];

  if(isnan(Vmax) || isinf(Vmax))return -1;
  else return Vmax;
}

//**********************************************************
bool isNoise(double* time, double* V, double baseline){
//**********************************************************

  bool isNoise = false;

  for(int i = 0; i < WFLENGTH; i++)if(time[i] > 0.05e-6 && time[i] < 0.2e-6 && V[i] < baseline){
      isNoise = true;
      break;
    }

  /*if(!isNoise){
    TGraph* tg = new TGraph(WFLENGTH,time,V);
    TLine* tl = new TLine(-0.4e-7,baseline,2e-6,baseline);
    tg->Draw("ap");tl->Draw();gPad->WaitPrimitive();
    }*/
  
  return isNoise;
}

//**********************************************************
void getVoltageFromName(const std::string& name,
			double &V){
//**********************************************************

  //the input file name
  std::string filename(name);

  //get Voltage from name
  V = stod(filename.substr(filename.find("V")+1,5)) - 41.69;
}

//**********************************************************
void analizeTree(const std::string& name,
		 double &DCR , double &eDCR ,
		 double &fDCR, double &feDCR,
		 double &XT  , double &eXT  ,
		 const double OV,
		 const int NFILES){
//**********************************************************

  //the input file name
  std::string filename(name);
  std::cout << "running on file '" << filename << "'" << endl;

  //get OV into a string to print plots with appropiate name
  std::stringstream sOV;
  sOV << OV;
  
  //histogram for dcr. Prepare logarithmic binning
  double binmin      = 1e-7;
  double binmax      = 1e2;
  const int nbins    = 100;
  double bins[nbins] = {0};
  for(int i = 0; i < nbins + 1 ; i++){
    bins[i] = TMath::Power(10,log10(binmin) + i*(log10(binmax/binmin)/nbins));
  }
  TH1F* hdcr = new TH1F("hdcr","hdcr",nbins,bins);

  //histogram for dcr (mHz)
  TH1F* hmhz = new TH1F("hmhz","hmhz",80,0,400);

  //histogram for xtalk
  TH1F* hxt  = new TH1F("hxt" ,"hxt" ,400,0,0.04);
  
  //TGraph for 2D plot
  TGraph* tg = new TGraph();

  //TLine for XT level
  TLine* tl = new TLine();

  //open root file
  TFile* rfile = TFile::Open(filename.c_str());
  TTree* wf    = (TTree*)rfile->Get("wf");

  //set addres for variables
  double fwftime   = 0;
  double ftime[WFLENGTH] = {0};
  double fV[WFLENGTH]    = {0};
  wf->SetBranchAddress("wftime", &fwftime);
  wf->SetBranchAddress("time"  , &ftime  );
  wf->SetBranchAddress("V"     , &fV     );
  
  //set variables for dcr
  double triggertime = 1;
  double fulltime    = 1;
  double baseline    = 0;
  int nevents        = 0;
  
  //Get number of wf
  int nentries = wf->GetEntries();

  //loop over wf
  for(int i = 0; i < nentries; i++){
    wf->GetEntry(i);
    triggertime = triggertime + fwftime;
    fulltime    = fulltime + fwftime;
    baseline    = getBaseline(ftime,fV);
    
    //if waveform is not noise, add one more event and fill graphs
    if(!isNoise(ftime,fV,baseline)){
      nevents++;
      hdcr->Fill(triggertime);
      hmhz->Fill(1000/(triggertime*36));
      hxt ->Fill(getMax(ftime,fV)-baseline);
      tg->SetPoint(nevents-1,triggertime,getMax(ftime,fV)-baseline);
      triggertime = 0;
    }
  }
  
  //compute DCR
  DCR  = nevents/fulltime;
  DCR  = DCR/36;
  DCR  = DCR*1000;
  eDCR = 1000*sqrt(nevents)/(fulltime*36);

  //compute xt
  //get events with more that one pe
  int multiple = getMultipleEvents(hxt,tl);
  XT  = 100*(double)multiple/nevents;
  eXT = 100*sqrt(pow(sqrt(multiple)/nevents,2)+pow(multiple*sqrt(nevents)/(nevents*nevents),2));

  //compute the fit rate
  hmhz->GetXaxis()->SetTitle("#it{Rate} (mHz/mm^{2})");
  hmhz->GetYaxis()->SetTitle("#it{Events}");
  hmhz->Draw("e");
  gPad->SetLogx(0);
  gPad->SetLogy();
  TF1* f = new TF1("f","landau(0)+pol0(3)");
  f->SetParameters(300,10,4,1);
  hmhz->Fit("f");
  fDCR  = f->GetParameter(1);
  feDCR = f->GetParError(1);
  if(plot_results)gPad->Print(("rate_"+sOV.str()+"OV_2.pdf").c_str());
  
  //plot 1D dcr log scale
  hdcr->GetXaxis()->SetTitle("#it{Delay time} (s)");
  hdcr->GetYaxis()->SetTitle("#it{Events}");
  hdcr->Draw();
  gPad->SetLogx();
  if(plot_results)gPad->Print(("delaytime_"+sOV.str()+"OV_2.pdf").c_str());
  
  //plot 2D plot with separation line
  tg->GetXaxis()->SetTitle("#it{Delay time} (s)");
  tg->GetYaxis()->SetTitle("#it{Amplitude} (V)");
  tg->SetMarkerStyle(20);
  gPad->SetLogy(0);
  tg->Draw("ap");
  tl->SetX1(tg->GetXaxis()->GetXmin());
  tl->SetX2(tg->GetXaxis()->GetXmax());
  tl->SetLineWidth(3);
  tl->SetLineColor(2);
  tl->Draw();
  if(plot_results)gPad->Print(("2dplot_"+sOV.str()+"OV_2.pdf").c_str());
    
  //delete histograms and other stuff
  gPad->SetLogx(0);
  gPad->SetLogy(0);
  gPad->Clear();
  hxt->Delete();
  hmhz->Delete();
  hdcr->Delete();
  tl->Delete();
  tg->Delete();
  //fit->Delete();
  
}

//**********************************************************
void dcr_from_tree(){
//**********************************************************

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  
  std::string listname = "/home/miguel/Documents/PhD/PD/macros/first_test_split/lists/dcr.list";
  
  // Open the file
  FILE *pFile = fopen (listname.c_str(), "r");
  if( pFile == NULL ){
    std::cout << "Cannot open File '" <<  listname << "'" << std::endl;
    exit(1);
  }

  //variables for main plot
  double DCR[NMAXFILES]   = {0};
  double eDCR[NMAXFILES]  = {0};
  double fDCR[NMAXFILES]  = {0};
  double feDCR[NMAXFILES] = {0};
  double XT[NMAXFILES]    = {0};
  double eXT[NMAXFILES]   = {0};
  double V[NMAXFILES]     = {0};
  
  //Read the list file
  char filename[200];
  int NFILES = 0;
  
  while(fscanf(pFile,"%s",filename) == 1 && NFILES < NMAXFILES){  
    if(filename[0] == '/' && filename[1] == '/')continue;
    getVoltageFromName(filename,V[NFILES]);
    analizeTree(filename,DCR[NFILES],eDCR[NFILES],fDCR[NFILES],feDCR[NFILES],XT[NFILES],eXT[NFILES],V[NFILES],NFILES);
    NFILES++;
  }

  std::cout << "OV(V) \t" << "DCR(mHz/mm²) \t Error" << std::endl;
  for(int i = 0; i < NFILES; i++)std::cout << DCR[i] << " " << eDCR[i] << std::endl;

  std::cout << "OV(V) \t" << "DCR(mHz/mm²) fit Error" << std::endl;
  for(int i = 0; i < NFILES; i++)std::cout << fDCR[i] << " " << feDCR[i] << std::endl;

  std::cout << "OV(V) \t" << "XT(%) \t \t Error" << std::endl;
  for(int i = 0; i < NFILES; i++)std::cout << XT[i]   << " " << eXT[i] << std::endl;
  

  TCanvas* c1 = new TCanvas("c1","c1",600,600);
  TGraphErrors* tg = new TGraphErrors(NFILES,V,DCR,0,eDCR);
  tg->SetMarkerStyle(20);
  tg->GetXaxis()->SetTitle("#it{OV} (V)");
  tg->GetYaxis()->SetTitle("#it{DCR} (mHz)");
  tg->Draw("ap");
  if(plot_results)c1->Print("DCR_2.pdf");

  TCanvas* c2 = new TCanvas("c2","c2",600,600);
  TGraphErrors* tg2 = new TGraphErrors(NFILES,V,fDCR,0,feDCR);
  tg2->SetMarkerStyle(20);
  tg2->GetXaxis()->SetTitle("#it{OV} (V)");
  tg2->GetYaxis()->SetTitle("#it{DCR from fit} (mHz)");
  tg2->Draw("ap");
  if(plot_results)c2->Print("DCRfromfit_2.pdf");

  TCanvas* c3 = new TCanvas("c3","c3",600,600);
  TGraphErrors* tg3 = new TGraphErrors(NFILES,V,XT,0,eXT);
  tg3->SetMarkerStyle(20);
  tg3->GetXaxis()->SetTitle("#it{OV} (V)");
  tg3->GetYaxis()->SetTitle("#it{X-Talk} (%)");
  tg3->Draw("ap");
  if(plot_results)c3->Print("XT_2.pdf");
}

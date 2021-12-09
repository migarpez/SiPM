#include "TSpectrum.h"
#include "Utils.C"

//general parameters
const int NMAXFILES = 100;
const int NMAXPEAKS = 5;
const int WFLENGTH  = 1500;

//parameters for ploting/printing results
bool plot_intermediate = false;
bool plot_is_noise     = false;
bool plot_clean_peaks  = false;
bool plot_results      = false;
bool use_abs_voltage   = false;
bool print_results     = true;

//**********************************************************
double getPE(TH1F* hxt){
//**********************************************************

  //create tspectrum for finding peak
  TSpectrum *s = new TSpectrum(1);

  //look for single pe value;
  s->Search(hxt,0.0001,"",0.15);
  double* pe = s->GetPositionX();

  std::cout << 3*pe[0]/2 << std::endl;
  if(plot_intermediate){hxt->Draw();gPad->Update();gPad->WaitPrimitive();}
  
  return pe[0];
}


//**********************************************************
int getMultipleEvents(TH1F* hxt,
		      TLine* tl,
		      double pe){
//**********************************************************

  tl->SetY1(3*pe/2);
  tl->SetY2(3*pe/2);

  double multiple = 0;

  for(int i = 0; i < hxt->GetNbinsX(); i++)if(hxt->GetBinCenter(i) > 3*pe/2)multiple = multiple + hxt->GetBinContent(i);

  return multiple;
}

//**********************************************************
double getBaseline(double*time, double* V){
//**********************************************************

  //fit all wf to an horizontal line and get mean value
  TGraph* tg = new TGraph(WFLENGTH,time,V);
  tg->Fit("pol0","Q","",2e-6,9e-6);
  tg->Draw("ap");
  gPad->Update();gPad->WaitPrimitive();
  if(!static_cast<TF1*>(tg->GetFunction("pol0"))){
    double baseline = 0;
    tg->Delete();
    return baseline;
  }
  else{
    double baseline = static_cast<TF1*>(tg->GetFunction("pol0"))->GetParameter(0);
    tg->Delete();
    return baseline;
  }
}

//**********************************************************
void getVoltageFromName(const std::string& name,
			double &V){
//**********************************************************

  //the input file name
  std::string filename(name);

  //get Voltage from name
  V = stod(filename.substr(filename.find("V")+1,5));
}

//**********************************************************
bool isNoise(double* time, double* V){
//**********************************************************

  bool isNoise = false;
  
  for(int i = 0; i < WFLENGTH; i++){
    if(time[i] > 0.15e-6 && time[i] < 0.2e-6 && V[i]<0.03){
      isNoise = true;
      break;
    }
  }

  if(isNoise && plot_is_noise){
    TGraph* tg = new TGraph(WFLENGTH,time,V);
    tg->Draw("ap");gPad->Update();gPad->WaitPrimitive();
  }
  
  return isNoise;
}

//**********************************************************
void analizeWf(double* ftime, double* fV,
	       const double absolute_time,
	       std::vector<std::pair<int,std::pair<double,double>>> &pulses,
	       int iwf){
//**********************************************************

  //create a histogram for the wf and another for the baseline
  TH1F* h   = new TH1F("h"  ,"h"  ,WFLENGTH,ftime[0],ftime[WFLENGTH-1]);

  //fill wf and baseline
  for(int i = 0; i < WFLENGTH; i++)h->SetBinContent(i,fV[i]);

  h->Rebin(16);//h->Draw();gPad->WaitPrimitive();

  //create tspectrum for finding peaks
  TSpectrum *s = new TSpectrum(NMAXPEAKS);

  //look for peaks
  int npeaks = s->Search(h,1,"",0.4);
  double* tpeak = s->GetPositionX();
  double* Vpeak = s->GetPositionY();
  //h->Draw();gPad->WaitPrimitive();

  //order peaks temporaly
  double aux1,aux2;
  for(int i = 0 ; i < npeaks ; i++){
    for(int j = i+1 ; j < npeaks ; j++){
      if(tpeak[i] > tpeak[j] && tpeak[j] != 0){
	aux1    = tpeak[i];
	tpeak[i] = tpeak[j];
	tpeak[j] = aux1;
	aux2    = Vpeak[i];
	Vpeak[i]  = Vpeak[j];
	Vpeak[j]  = aux2;
      }
    }
  }

  if(npeaks>1){
    h->Draw();
    if(plot_clean_peaks){gPad->Update();gPad->WaitPrimitive();}
  }

  //store info in the tgraph
  for(int i = 0; i < npeaks; i++){
    pulses.push_back(std::make_pair(iwf,std::make_pair(tpeak[i]+absolute_time,Vpeak[i])));
  }
  
  h->Delete();
}

//**********************************************************
void analizeTree(const std::string& name,
		 double &DCR , double &eDCR ,
		 double &fDCR, double &feDCR,
		 double &XT  , double &eXT  ,
		 double &AP  , double &eAP  ,
		 const int NFILES){
//**********************************************************

  //the input file name
  std::string filename(name);
  std::cout << "running on file '" << filename << "'" << endl;

  //histogram for dcr. Prepare logarithmic binning
  double binmin      = 1e-7;
  double binmax      = 1e2;
  const int nbins    = 50;
  double bins[nbins] = {0};
  for(int i = 0; i < nbins + 1 ; i++){
    bins[i] = TMath::Power(10,log10(binmin) + i*(log10(binmax/binmin)/nbins));
  }
  TH1F* hdcr = new TH1F("","",nbins,bins);

  //histogram for dcr (mHz)
  TH1F* hmhz = new TH1F("hmhz","hmhz",80,0,150);

  //histogram for xtalk
  TH1F* hxt  = new TH1F("hxt" ,"hxt" ,1000,0,5);
  
  //vector for storing peak information
  std::vector<std::pair<int,std::pair<double,double>>> pulses;
  pulses.clear();

  //TLine for XT level
  TLine* tlxt = new TLine();
  TLine* tlap = new TLine();
  
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

  //set variables for cnoise
  double absolute_time = 1;
  int ntriggers        = 0;
  int nevents          = 0;
  
  //Get number of wf
  int nentries = wf->GetEntries();
  
  //loop over the wfs
  for(int i = 0; i < nentries; i++){

    //get info
    wf->GetEntry(i);
    absolute_time = absolute_time + fwftime;

    //if waveform is not noise, add one more event and fill graphs
    if(!isNoise(ftime,fV)){
      //get peaks position and amplitude from wf and fill histograms
      analizeWf(ftime,fV,absolute_time,pulses,i);
      ntriggers++;
      if(i%500==0)std::cout << i+1 << "/" << nentries << std::endl;
    }
  }

  //fill amplitude histogram
  for(int i = 0; i < pulses.size(); i++){
      hxt->Fill(pulses[i].second.second);
  }

  //get photoelectron level to set boundaries
  double pe = getPE(hxt);
  
  //create new tgraph to do the main plot. 
  TGraph* mtg = new TGraph();

  //run over the stored peaks in and get info
  double V0,t0,Vnext,tnext,triggertime;
  int point = 0;
  //get point "zero"
  int zero = 0;
  for(int i = 0; i < pulses.size(); i++){
    if(pulses[i].second.second<pe/2)continue;
    t0 = pulses[i].second.first;
    V0 = pulses[i].second.second;
    zero = i;
    nevents++;
    break;
  }

  for(int i = zero+1; i < pulses.size(); i++){
    //get "next" point
    tnext = pulses[i].second.first;
    Vnext = pulses[i].second.second;
    //skip it if below half PE
    if(Vnext<pe/2)continue;

    //compute time between points and increase AP if needed
    triggertime = tnext-t0;
    if(triggertime < 0)continue;
    if(triggertime > 0 && triggertime<5e-6)AP++;
    
    //fill histograms and main plot
    mtg->SetPoint(nevents-1,triggertime,Vnext);
    hdcr->Fill(triggertime);
    hmhz->Fill(1000/((triggertime)*36));

    V0 = Vnext;
    t0 = tnext;
    nevents++;
  }

  //compute AP
  eAP = 100*sqrt(pow(sqrt(AP)/nevents,2)+pow(AP*sqrt(nevents)/(nevents*nevents),2));
  AP  = 100*AP/nevents;

  //compute xt
  //get events with more that one pe
  int multiple = getMultipleEvents(hxt,tlxt,pe);
  XT  = 100*(double)multiple/nevents;
  eXT = 100*sqrt(pow(sqrt(multiple)/nevents,2)+pow(multiple*sqrt(nevents)/(nevents*nevents),2));
  
  //draw and print 1D rate plot. Compute DCR
  gPad->SetLogx();
  hdcr->GetXaxis()->SetTitle("#it{Delay time} (s)");
  hdcr->GetYaxis()->SetTitle("#it{Events}");
  hdcr->Draw();if(plot_intermediate){gPad->Update();gPad->WaitPrimitive();}
  DCR  = nevents/absolute_time;//ntriggers/absolute_time;
  DCR  = DCR/36;
  DCR  = DCR*1000;
  eDCR = 1000*sqrt(nevents)/(absolute_time*36);

  //draw hz distribution, fit and compute fit DCR
  //gPad->SetLogx(0);
  gPad->SetLogy(0);
  gPad->SetLogx(0);
  hmhz->GetXaxis()->SetTitle("#it{Rate} (mHz/mm^{2})");
  hmhz->GetYaxis()->SetTitle("#it{Events}");
  hmhz->Draw("e");
  TF1* f = new TF1("f","landau(0)+pol0(3)");
  f->SetParameters(60,10,5,0);
  f->SetParLimits(1,0,30);
  hmhz->Fit("f","Q");if(plot_intermediate){gPad->Update();gPad->WaitPrimitive();}
  fDCR  = f->GetParameter(1);
  feDCR = f->GetParError(1);
  
  //draw and print main plot
  gPad->SetLogx();gPad->SetLogy(0);
  mtg->SetMarkerStyle(20);
  mtg->GetXaxis()->SetTitle("#it{Delay time} (s)");
  mtg->GetYaxis()->SetTitle("#it{Amplitude} (V)");
  mtg->Draw("ap");
  tlxt->SetX1(mtg->GetXaxis()->GetXmin());
  tlxt->SetX2(mtg->GetXaxis()->GetXmax());
  tlxt->SetLineWidth(3);
  tlxt->SetLineColor(2);
  tlxt->Draw();
  tlap->SetY1(mtg->GetYaxis()->GetXmin());
  tlap->SetY2(mtg->GetYaxis()->GetXmax());
  tlap->SetX1(5e-7);
  tlap->SetX2(5e-7);
  tlap->SetLineWidth(3);
  tlap->SetLineColor(2);
  tlap->Draw();
  if(plot_intermediate){gPad->Update();gPad->WaitPrimitive();}
  
  gPad->SetLogx(0);gPad->SetLogy(0);
  
  //delete stuff
  gPad->Clear();
  hxt->Delete();
  hmhz->Delete();
  hdcr->Delete();
  tlxt->Delete();
  tlap->Delete();
  mtg->Delete();
}

//**********************************************************
void cnoise_from_tree(){
//**********************************************************
  
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  std::string listname = "lists/cnoise.list";
  
  // Open the file
  FILE *pFile = fopen (listname.c_str(), "r");
  if( pFile == NULL ){
    std::cout << "Cannot open File '" <<  listname << "'" << std::endl;
    exit(1);
  }

  //variables for main plots
  double DCR[NMAXFILES]     = {0};
  double eDCR[NMAXFILES]    = {0};
  double fDCR[NMAXFILES]    = {0};
  double feDCR[NMAXFILES]   = {0};
  double XT[NMAXFILES]      = {0};
  double eXT[NMAXFILES]     = {0};
  double AP[NMAXFILES]      = {0};
  double eAP[NMAXFILES]     = {0};
  
  //Read the list file
  char filename[200];
  int NFILES = 0;
  int sipm   = 0;
  
  while(fscanf(pFile,"%s",filename) == 1 && NFILES < NMAXFILES){  
    if(filename[0] == '/' && filename[1] == '/')continue;
    analizeTree(filename,DCR[NFILES],eDCR[NFILES],fDCR[NFILES],feDCR[NFILES],XT[NFILES],eXT[NFILES],AP[NFILES],eAP[NFILES],NFILES);
    NFILES++;
  }
  
  //print results
  if(print_results){
    std::cout << "------------------------------------------------------------------" << std::endl;
    std::cout << "RESULTS" << std::endl;
    std::cout << "------------------------------------------------------------------" << std::endl;
    std::cout << "DCR(mHz/mm²) \t Error" << std::endl;
    for(int i = 0; i < NFILES; i++)std::cout << DCR[i] << " " << eDCR[i] << std::endl;

    std::cout << "DCR(mHz/mm²) fit Error" << std::endl;
    for(int i = 0; i < NFILES; i++)std::cout << fDCR[i] << " " << feDCR[i] << std::endl;

    std::cout << "XT(%) \t \t Error" << std::endl;
    for(int i = 0; i < NFILES; i++)std::cout << XT[i] << " " << eXT[i] << std::endl;

    std::cout << "AP(%) \t \t Error" << std::endl;
    for(int i = 0; i < NFILES; i++)std::cout << AP[i] << " " << eAP[i] << std::endl;
  }

  TH1F* hdcr = new TH1F("hdcr","hdcr",20,0,50);
  TH1F* hxt = new TH1F("hxt","hxt",20,0,30);
  TH1F* hap = new TH1F("hap","hap",20,0,5);
  for(int i = 0; i < NFILES; i++){
    hdcr->Fill(fDCR[i]);
    hxt->Fill(XT[i]);
    hap->Fill(AP[i]);
  }

  hdcr->Draw();gPad->Update();gPad->WaitPrimitive();
  hxt->Draw();gPad->Update();gPad->WaitPrimitive();
  hap->Draw();gPad->Update();gPad->WaitPrimitive();
}


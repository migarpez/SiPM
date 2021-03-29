#include "TSpectrum.h"

//general parameters
const int NMAXFILES = 4;
const int NMAXPEAKS = 5;
const int WFLENGTH  = 2500;

//sipm pitch related
const int SIPM75[]           = {11,15,16,17,19,21,12}; //12 is from Madrid
const double OV50[NMAXFILES] = {3,4  ,5};
const double OV75[NMAXFILES] = {2,2.5,3};

//
const double cut_level = 0.004; //zero if not wanted a lower limit

//parameters for ploting/printing results
bool plot_intermediate = true;
bool plot_results      = true;
bool use_abs_voltage   = false;
bool print_results     = true;

//**********************************************************
int getSiPMFromName(char* name){
//**********************************************************
  //the input file name
  std::string filename(name);

  //get sipm number from filename
  int sipm = stoi(filename.substr(filename.find("sipm")+4,2));
  return sipm;
}

//**********************************************************
bool is75Pitch(const int sipm){
//**********************************************************

  bool isIt = false;
  for(int i = 0; i < sizeof(SIPM75); i++){
    if(sipm == SIPM75[i]){
      isIt = true;
      break;
    }
  }
  return isIt;
}

//**********************************************************
void initializeByPitch(const int sipm, double* OV){
//**********************************************************
  if(is75Pitch(sipm)){
    for(int i = 0; i < NMAXFILES; i++)OV[i] = OV75[i];
    std::cout << "------------------------------------------------------------------" << std::endl;
    std::cout << "Analyzing Dark Noise from 75 um pitch sipm " << sipm << std::endl;
    std::cout << "------------------------------------------------------------------" << std::endl;
  }
  else{
    for(int i = 0; i < NMAXFILES; i++)OV[i] = OV50[i];
    std::cout << "------------------------------------------------------------------" << std::endl;
    std::cout << "Analyzing Dark Noise from 50 um pitch sipm " << sipm << std::endl;
    std::cout << "------------------------------------------------------------------" << std::endl;
  }
}

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
  std::cout << oneandahalfpe << std::endl;
  hxt->Draw();gPad->WaitPrimitive();
  double multiple = 0;

  for(int i = 0; i < hxt->GetNbinsX(); i++)if(hxt->GetBinCenter(i) > oneandahalfpe)multiple = multiple + hxt->GetBinContent(i);

  return multiple;
}

//**********************************************************
double getBaseline(double*time, double* V){
//**********************************************************

  //fit all wf to an horizontal line and get mean value
  TGraph* tg = new TGraph(WFLENGTH,time,V);
  tg->Fit("pol0","Q","",4e-6,9e-6);
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
bool isNoise(double* time, double* V, const double baseline){
//**********************************************************

  bool isNoise = false;

  for(int i = 0; i < WFLENGTH; i++)if(time[i] > 0.05e-6 && time[i] < 0.3e-6 && V[i] < baseline){
      isNoise = true;
      //Graph* tg = new TGraph(WFLENGTH,time,V);tg->Draw("ap");gPad->WaitPrimitive();
      break;
    }
  
  return isNoise;
}

//**********************************************************
void analizeWf(double* ftime, double* fV,
	       const double absolute_time,
	       const double baseline,
	       int &nevents,
	       TGraph* tg){
//**********************************************************

  //create a histogram for the wf and another for the baseline
  TH1F* h   = new TH1F("h"  ,"h"  ,WFLENGTH,ftime[0],ftime[WFLENGTH-1]);
  TH1F* hbl = new TH1F("hbl","hbl",200     ,-0.05  ,0.05);

  //fill wf and baseline
  for(int i = 0; i < WFLENGTH; i++){
    h->SetBinContent(i,fV[i]-baseline);
    hbl->Fill(fV[i]-baseline);
  }

  h->Rebin(2);

  //create tspectrum for finding peaks
  TSpectrum *s = new TSpectrum(NMAXPEAKS);

  //look for peaks
  int npeaks = s->Search(h,0.00001,"",0.4);
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

  //Get baseline rms
  double blrms    = hbl->GetRMS();
  
  //check that found peaks are not noise
  int clean_peaks = npeaks;
  for(int i = 0; i < npeaks; i++){
    //set absolute time for peaks
    tpeak[i] = tpeak[i] + absolute_time;
    if(Vpeak[i] < 3*blrms){
      Vpeak[i] = -999;
      clean_peaks = clean_peaks-1;
    }
  }

  if(clean_peaks>1){
    h->Draw();
    std::cout << baseline << "+/- " << hbl->GetRMS() << " clean = " << clean_peaks << std::endl;
    //gPad->WaitPrimitive();
  }

  //store info in the tgraph
  for(int i = 0; i < npeaks; i++){
    if(Vpeak[i] == -999 || tpeak[i] < 0 || Vpeak[i] < cut_level)continue;
    tg->SetPoint(nevents,tpeak[i],Vpeak[i]);
    //tg->Draw("ap");gPad->WaitPrimitive();
    nevents = nevents+1;
  }

  hbl->Delete();
  h->Delete();

}

//**********************************************************
void analizeTree(const std::string& name,
		 double &DCR , double &eDCR ,
		 double &fDCR, double &feDCR,
		 double &XT  , double &eXT  ,
		 double &AP  , double &eAP  ,
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
  TH1F* hmhz = new TH1F("hmhz","hmhz",200,0,1000);

  //histogram for xtalk
  TH1F* hxt  = new TH1F("hxt" ,"hxt" ,1000,0,0.5);
  
  //TGraph for 2D plot
  TGraph* tg = new TGraph();

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
  double baseline      = 0;
  int ntriggers        = 0;
  int nevents          = 0;
  
  //Get number of wf
  int nentries = wf->GetEntries();
  
  //loop over the wfs
  for(int i = 0; i < nentries; i++){

    //get info
    wf->GetEntry(i);
    absolute_time = absolute_time + fwftime;
    baseline      = getBaseline(ftime,fV);

    //if waveform is not noise, add one more event and fill graphs
    if(!isNoise(ftime,fV,baseline)){
      //get peaks position and amplitude from wf and fill histograms
      analizeWf(ftime,fV,absolute_time,baseline,nevents,tg);
      ntriggers++;
      if(i%5000==0)std::cout << i+1 << "/" << nentries << std::endl;
    }
  }

  //create new tgraph to do the main plot. 
  TGraph* mtg = new TGraph();

  //run over the stored peaks in the histogram and get info
  double V0,t0,Vnext,tnext,triggertime;
  for(int i = 1; i < nevents; i++){

    //get point "zero"
    tg->GetPoint(i-1,t0,V0);
    //get "next" point
    tg->GetPoint(i,tnext,Vnext);

    //compute time between points and increase AP if needed
    triggertime = tnext-t0;
    if(triggertime < 0 || triggertime > 100)continue;
    if(triggertime>0 && triggertime<5e-6)AP++;
    
    //fill histograms and main plot
    mtg->SetPoint(i-1,triggertime,Vnext);
    hdcr->Fill(triggertime);
    hmhz->Fill(1000/((triggertime)*36));
    hxt ->Fill(Vnext);
  }

  //compute AP
  eAP = 100*sqrt(pow(sqrt(AP)/nevents,2)+pow(AP*sqrt(nevents)/(nevents*nevents),2));
  AP  = 100*AP/nevents;
  
  //compute xt
  //get events with more that one pe
  int multiple = getMultipleEvents(hxt,tlxt);
  XT  = 100*(double)multiple/nevents;
  eXT = 100*sqrt(pow(sqrt(multiple)/nevents,2)+pow(multiple*sqrt(nevents)/(nevents*nevents),2));
  
  //draw and print 1D rate plot. Compute DCR
  gPad->SetLogx();
  hdcr->GetXaxis()->SetTitle("#it{Delay time} (s)");
  hdcr->GetYaxis()->SetTitle("#it{Events}");
  hdcr->Draw();if(plot_intermediate){gPad->Update();gPad->WaitPrimitive();}
  if(plot_results)gPad->Print(("delaytime_"+sOV.str()+"OV.pdf").c_str());
  DCR  = ntriggers/absolute_time;
  DCR  = DCR/36;
  DCR  = DCR*1000;
  eDCR = 1000*sqrt(nevents)/(absolute_time*36);

  //draw hz distribution, fit and compute fit DCR
  gPad->SetLogx(0);gPad->SetLogy();
  hmhz->GetXaxis()->SetTitle("#it{Rate} (mHz/mm^{2})");
  hmhz->GetYaxis()->SetTitle("#it{Events}");
  hmhz->Draw("e");
  TF1* f = new TF1("f","landau(0)+pol0(3)");
  f->SetParameters(400,100,40,1);
  hmhz->Fit("f","Q");if(plot_intermediate){gPad->Update();gPad->WaitPrimitive();}
  fDCR  = f->GetParameter(1);
  feDCR = f->GetParError(1);
  if(plot_results)gPad->Print(("rate_"+sOV.str()+"OV.pdf").c_str());

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
  tlap->Draw();if(plot_intermediate){gPad->Update();gPad->WaitPrimitive();}
  if(plot_results)gPad->Print(("2dplot_"+sOV.str()+"OV.pdf").c_str());

  //delete stuff
  gPad->Clear();
  hxt->Delete();
  hmhz->Delete();
  hdcr->Delete();
  tlxt->Delete();
  tlap->Delete();
  mtg->Delete();
  tg->Delete();
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
  double V[NMAXFILES]       = {0};
  double OV[NMAXFILES]      = {0};
  double voltage[NMAXFILES] = {0};
  
  //Read the list file
  char filename[200];
  int NFILES = 0;
  int sipm   = 0;
  
  while(fscanf(pFile,"%s",filename) == 1 && NFILES < NMAXFILES){  
    if(filename[0] == '/' && filename[1] == '/')continue;
    if(NFILES == 0){
      sipm = getSiPMFromName(filename);
      initializeByPitch(sipm,OV);
    }
    //getx operating voltage and decide wether to use absolute or relative
    getVoltageFromName(filename,V[NFILES]);
    if(use_abs_voltage)voltage[NFILES] = V[NFILES];
    else voltage[NFILES] = OV[NFILES];
    analizeTree(filename,DCR[NFILES],eDCR[NFILES],fDCR[NFILES],feDCR[NFILES],XT[NFILES],eXT[NFILES],AP[NFILES],eAP[NFILES],voltage[NFILES],NFILES);
    NFILES++;
  }

  //print results
  if(print_results){
    std::cout << "------------------------------------------------------------------" << std::endl;
    std::cout << "RESULTS" << std::endl;
    std::cout << "------------------------------------------------------------------" << std::endl;
    std::cout << "OV(V) \t" << "DCR(mHz/mm²) \t Error" << std::endl;
    for(int i = 0; i < NFILES; i++)std::cout << DCR[i] << " " << eDCR[i] << std::endl;

    std::cout << "OV(V) \t" << "DCR(mHz/mm²) fit Error" << std::endl;
    for(int i = 0; i < NFILES; i++)std::cout << fDCR[i] << " " << feDCR[i] << std::endl;

    std::cout << "OV(V) \t" << "XT(%) \t \t Error" << std::endl;
    for(int i = 0; i < NFILES; i++)std::cout << XT[i] << " " << eXT[i] << std::endl;

    std::cout << "OV(V) \t" << "AP(%) \t \t Error" << std::endl;
    for(int i = 0; i < NFILES; i++)std::cout << AP[i] << " " << eAP[i] << std::endl;
  }

  //plot results
  TCanvas* c1 = new TCanvas("c1","c1",600,600);
  TGraphErrors* tg = new TGraphErrors(NFILES,V,DCR,0,eDCR);
  tg->SetMarkerStyle(20);
  tg->GetXaxis()->SetTitle("#it{OV} (V)");
  tg->GetYaxis()->SetTitle("#it{DCR} (mHz)");
  tg->Draw("ap");
  if(plot_results)c1->Print("DCR.pdf");

  TCanvas* c2 = new TCanvas("c2","c2",600,600);
  TGraphErrors* tg2 = new TGraphErrors(NFILES,V,fDCR,0,feDCR);
  tg2->SetMarkerStyle(20);
  tg2->GetXaxis()->SetTitle("#it{OV} (V)");
  tg2->GetYaxis()->SetTitle("#it{DCR from fit} (mHz)");
  tg2->Draw("ap");
  if(plot_results)c2->Print("DCRfromfit.pdf");

  TCanvas* c3 = new TCanvas("c3","c3",600,600);
  TGraphErrors* tg3 = new TGraphErrors(NFILES,V,XT,0,eXT);
  tg3->SetMarkerStyle(20);
  tg3->GetXaxis()->SetTitle("#it{OV} (V)");
  tg3->GetYaxis()->SetTitle("#it{X-Talk} (%)");
  tg3->Draw("ap");
  if(plot_results)c3->Print("XT.pdf");

  TCanvas* c4 = new TCanvas("c4","c4",600,600);
  TGraphErrors* tg4 = new TGraphErrors(NFILES,V,AP,0,eAP);
  tg4->SetMarkerStyle(20);
  tg4->GetXaxis()->SetTitle("#it{OV} (V)");
  tg4->GetYaxis()->SetTitle("#it{AP} (%)");
  tg4->Draw("ap");
  if(plot_results)c4->Print("AP.pdf");
}


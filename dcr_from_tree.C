#include "TSpectrum.h"
#include "Utils.C"

//general parameters
const int NMAXFILES = 3;
const int WFLENGTH  = 3125;

//sipm pitch related
const int SIPM75[]           = {11,15,16,17,19,21,12}; //12 is from Madrid
const double OV50[NMAXFILES] = {3,4  ,5};
const double OV75[NMAXFILES] = {2,2.5,3};

//
const double cut_level = 0.005;//0.0035; //zero if not wanted a lower limit
const double noise_level = 0.000; //zero if not wanted a lower limit

//parameters for ploting/printing results

bool plot_intermediate = true;
bool plot_is_noise     = false;
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

  std::cout << pe << std::endl;
  
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
  tg->Fit("pol0","Q","");
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
  double vmax  = 0;
  
  for(int i = 0; i < WFLENGTH; i++){
    if(time[i] > 0.05e-6 && time[i] < 0.1e-6 && V[i]-baseline<0){
      isNoise = true;
      break;
    }
  }

  if(plot_is_noise && isNoise){
    std::cout << "noise: " << isNoise << std::endl;
    TGraph* tg = new TGraph(WFLENGTH,time,V);tg->Draw("ap");gPad->Update();gPad->WaitPrimitive();
  }

  return isNoise;
}

//**********************************************************
double getMaximum(double* time, double* V){
//**********************************************************

  double max = 0;

  for(int i = 0; i < WFLENGTH; i++){
    if(time[i]<2e-6)if(V[i]>max)max=V[i];
  }

  return max;
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
  const int nbins    = 50;
  double bins[nbins] = {0};
  for(int i = 0; i < nbins + 1 ; i++){
    bins[i] = TMath::Power(10,log10(binmin) + i*(log10(binmax/binmin)/nbins));
  }
  TH1F* hdcr = new TH1F("","",nbins,bins);

  //histogram for dcr (mHz)
  TH1F* hmhz = new TH1F("hmhz","hmhz",250,0.1,250);

  //histogram for xtalk
  TH1F* hxt  = new TH1F("hxt" ,"hxt" ,1000,0,0.1);
  
  //TLine for XT level
  TLine* tlxt = new TLine();
  
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
  double vmax          = 0;
  std::vector<double> vmaxs;
  std::vector<double> times;
  vmaxs.clear();
  times.clear();
  
  //Get number of wf
  int nentries = wf->GetEntries();
  
  //loop over the wfs
  //for(int i = 0; i < 1000; i++){
  for(int i = 0; i < nentries; i++){

    //get info
    wf->GetEntry(i);
    absolute_time = absolute_time + fwftime;
    baseline      = getBaseline(ftime,fV);

    //if waveform is not noise, add one more event and fill graphs
    if(!isNoise(ftime,fV,baseline)){
      //get maximum value of the wf
      vmax = getMaximum(ftime,fV)-baseline;
      vmaxs.push_back(vmax);
      times.push_back(fwftime);
      hxt->Fill(vmax);
      if(i%5000==0)std::cout << i+1 << "/" << nentries << std::endl;
    }
  }
  
  std::cout << absolute_time << std::endl;
  
  //get photoelectron level to set boundaries
  double pe = getPE(hxt);
  std::cout << pe << std::endl;
  
  //run over the stored peaks
  double triggertime = 0;
  int point = 0;
  for(int i = 0; i < vmaxs.size(); i++){

    triggertime = triggertime+times[i];
    //skip it if below half PE
    if(vmaxs[i]<pe/2 || vmaxs[i] < cut_level)continue;

    if(times[i] < 0 || times[i] > 100)continue;
    
    //fill histograms and main plot
    hdcr->Fill(triggertime);
    hmhz->Fill(1000/((triggertime)*36));
    triggertime = 0;
    point++;
  }

  //compute xt
  //get events with more that one pe
  int multiple = getMultipleEvents(hxt,tlxt,pe);
  XT  = 100*(double)multiple/point;
  eXT = 100*sqrt(pow(sqrt(multiple)/point,2)+pow(multiple*sqrt(point)/(point*point),2));
  
  //draw and print 1D rate plot. Compute DCR
  gPad->SetLogx();
  hdcr->GetXaxis()->SetTitle("#it{Delay time} (s)");
  hdcr->GetYaxis()->SetTitle("#it{Events}");
  hdcr->Draw();if(plot_intermediate){gPad->Update();gPad->WaitPrimitive();}
  if(plot_results)gPad->Print(("delaytime_"+sOV.str()+"OV.pdf").c_str());
  DCR  = point/absolute_time;
  DCR  = DCR/36;
  DCR  = DCR*1000;
  eDCR = 1000*sqrt(point)/(absolute_time*36);

  //draw hz distribution, fit and compute fit DCR
  //gPad->SetLogx(0);
  gPad->SetLogy();
  hmhz->GetXaxis()->SetTitle("#it{Rate} (mHz/mm^{2})");
  hmhz->GetYaxis()->SetTitle("#it{Events}");
  hmhz->Draw("e");
  TF1* f = new TF1("f","landau(0)+pol0(3)");
  f->SetParameters(60,1,1,2);
  hmhz->Fit("f","Q");if(plot_intermediate){gPad->Update();gPad->WaitPrimitive();}
  fDCR  = f->GetParameter(1);
  feDCR = f->GetParError(1);
  if(plot_results)gPad->Print(("rate_"+sOV.str()+"OV.pdf").c_str());

  gPad->SetLogy(0);
  gPad->SetLogx(0);
  
  //delete stuff
  gPad->Clear();
  hxt->Delete();
  hmhz->Delete();
  hdcr->Delete();
  tlxt->Delete();
}

//**********************************************************
void dcr_from_tree(){
//**********************************************************
  
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  std::string listname = "lists/dcr.list";
  
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
  double V[NMAXFILES]       = {0};
  double OV[NMAXFILES]      = {3,4,5};
  double voltage[NMAXFILES] = {0};
  
  //Read the list file
  char filename[200];
  int NFILES = 0;
  int sipm   = 0;
  
  while(fscanf(pFile,"%s",filename) == 1 && NFILES < NMAXFILES){  
    if(filename[0] == '/' && filename[1] == '/')continue;
    if(NFILES == 0){
      sipm = SiPMUtils::getSiPMFromName(filename);
    }
    analizeTree(filename,DCR[NFILES],eDCR[NFILES],fDCR[NFILES],feDCR[NFILES],XT[NFILES],eXT[NFILES],voltage[NFILES],NFILES);
    NFILES++;
  }

  //print results
  if(print_results){
    std::cout << "------------------------------------------------------------------" << std::endl;
    std::cout << "RESULTS" << std::endl;
    std::cout << "------------------------------------------------------------------" << std::endl;
    //std::cout << "OV(V) \t" << "DCR(mHz/mm²) \t Error" << std::endl;
    //for(int i = 0; i < NFILES; i++)std::cout << DCR[i] << " " << eDCR[i] << std::endl;

    //std::cout << "OV(V) \t" << "DCR(mHz/mm²) fit Error" << std::endl;
    //for(int i = 0; i < NFILES; i++)std::cout << fDCR[i] << " " << feDCR[i] << std::endl;

    //std::cout << "OV(V) \t" << "XT(%) \t \t Error" << std::endl;
    //for(int i = 0; i < NFILES; i++)std::cout << XT[i] << " " << eXT[i] << std::endl;

    for(int i = 0; i < NFILES; i++)std::cout << DCR[i] << " " << eDCR[i] << " " << fDCR[i] << " " << feDCR[i]<< " " << XT[i] << " " << eXT[i]<< std::endl;
  }

  //plot results
  TCanvas* c1 = new TCanvas("c1","c1",600,600);
  TGraphErrors* tg = new TGraphErrors(NFILES,OV,DCR,0,eDCR);
  tg->SetMarkerStyle(20);
  tg->GetXaxis()->SetTitle("#it{OV} (V)");
  tg->GetYaxis()->SetTitle("#it{DCR} (mHz)");
  tg->Draw("ap");
  if(plot_results)c1->Print("DCR.pdf");

  TCanvas* c2 = new TCanvas("c2","c2",600,600);
  TGraphErrors* tg2 = new TGraphErrors(NFILES,OV,fDCR,0,feDCR);
  tg2->SetMarkerStyle(20);
  tg2->GetXaxis()->SetTitle("#it{OV} (V)");
  tg2->GetYaxis()->SetTitle("#it{DCR from fit} (mHz)");
  tg2->Draw("ap");
  if(plot_results)c2->Print("DCRfromfit.pdf");

  TCanvas* c3 = new TCanvas("c3","c3",600,600);
  TGraphErrors* tg3 = new TGraphErrors(NFILES,OV,XT,0,eXT);
  tg3->SetMarkerStyle(20);
  tg3->GetXaxis()->SetTitle("#it{OV} (V)");
  tg3->GetYaxis()->SetTitle("#it{X-Talk} (%)");
  tg3->Draw("ap");
  if(plot_results)c3->Print("XT.pdf");
}


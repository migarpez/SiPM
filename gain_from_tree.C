//IMPORTANT!
//this script needs an input list to work, list/gain.list
//it assumes files are provided from lower OV to higher OV and that only
//three files are being analized. It also assumes a single SiPM is being
//analized

#include "Utils.C"

//**********************************************************
//******PREVIOUS DEFINITIONS********************************
//**********************************************************

//general parameters
const int NMAXFILES    = 3;
const int NSETFILES    = 3;
const int NMAXPEAKS    = 15;
const int NMAXPEAKSFIT = 3;
const int WFLENGTH     = 1300;

//constants for calculating
const double AMPFACTOR_LN2_50O = 2940;
const double AMPFACTOR_LN2_1MO = 6740;
const double AMPFACTOR         = 1;//AMPFACTOR_LN2_1MO;
const double ECHARGE           = 1.60217662e-19;

//sipm pitch related

//parameters for defining QvsPE histogram
double hqmin1;
double hqmin2;
double hqmax1;
double hqmax2;

//parameters for cleaning data sample
const bool lookForDCR = false; //if integrating whole wf no need to use it

//parameters for numerical integration. If integrateUpTo90 is
//true, it will use time_limit_up between min and max. If false,
//it will integrate till waveform end time_limit_up_fix
const bool integrateUpToPer    = false;
const double PER               = 0.90;
const double area_factor       = 1;//1.233; //1.233 for 75%
const double time_limit_down   = 0.00000005;
//const double time_limit_up_fix[3] = {0.0000012,0.0000013,0.0000016};
const double time_limit_up_fix = 0.000005;
const double time_limit_up_min = 0.0000004;
const double time_limit_up_max = 0.000005;
const bool use_fit             = false;  //if true it fits the wf to a gausºexp
                                        //and integrates it

//parameters for ploting/printing results
bool plot_intermediate = true;
bool plot_results      = false;
bool use_abs_voltage   = false;
bool print_results     = true;


//**********************************************************
//******FIT CHARGE HISTOGRAM********************************
//**********************************************************

//**********************************************************
Double_t fpeaks(Double_t *x, Double_t *par){
//**********************************************************
  Double_t result = 0;
  for (int i = 0 ; i < NMAXPEAKS ; i++){
    Double_t norm  = par[3*i+0];
    Double_t mean  = par[3*i+1];
    Double_t sigma = par[3*i+2];
    result += norm*TMath::Gaus(x[0],mean,sigma);
  }
  return result;
}

//**********************************************************
void computeGain(TH1F* hq,
		 double* SNR,
		 double &GAIN, double &eGAIN,
		 const double OV,
		 const int NFILES){
//**********************************************************

  std::cout << std::endl;
  
  //get OV into a string to print plots with appropiate name
  std::stringstream sOV;
  sOV << OV;
  
  Double_t par[NMAXPEAKS*3] = {0};

  //create an instance of TSpectrum
  TSpectrum *s = new TSpectrum(NMAXPEAKS);
  hq->GetXaxis()->SetTitle("#it{Charge} (nV s)");
  hq->GetYaxis()->SetTitle("#it{Counts}");

  hq->Draw();
  if(plot_intermediate){
    gPad->Update();
    gPad->WaitPrimitive();
  }
  
  //substract background
  TH1 *hb = s->Background(hq,20,"same");
  hq->Add(hb,-1);

  //rebin
  hq->Rebin(4);
  hq->Draw();
  
  //go for peaks
  int npeaks = s->Search(hq,1,"",0.05);
  Double_t *peakx = s->GetPositionX();
  Double_t *peaky = s->GetPositionY();
  for(int i = 0 ; i < npeaks ; i++){
    par[i*3+0] = peaky[i];
    par[i*3+1] = peakx[i];
    par[i*3+2] = 0.00000001;
  }

  if(plot_intermediate){
    gPad->Update();
    gPad->WaitPrimitive();
  }
  
  //second gaussian sumatory
  TF1* fgaus = new TF1("fgaus",fpeaks,hqmin1+NFILES*(hqmin2-hqmin1)/NMAXFILES,hqmax1+NFILES*(hqmax2-hqmax1)/NMAXFILES,3*npeaks);
  fgaus->SetParameters(par);

  //for(int i = 0 ; i < npeaks ; i++){
  //  fgaus->SetParLimits(i*3+0,1,1000);
  //  fgaus->SetParLimits(i*3+2,0,0.1);
  //}

  
  //fit histogram without bg
  hq->GetXaxis()->SetTitle("#it{Charge} (Vs)");
  hq->GetYaxis()->SetTitle("#it{Entries}");
  hq->Fit("fgaus","Q");
  if(plot_intermediate){
    gPad->Update();
    gPad->WaitPrimitive();
  }
  if(plot_results)gPad->Print(("spectra_"+sOV.str()+"OV.pdf").c_str());

  //store peaks position and rms
  double mean[NMAXPEAKS] = {0}, emean[NMAXPEAKS] = {0}, rms[NMAXPEAKS] = {0}, number[NMAXPEAKS] = {0};
  for(int i = 0; i < npeaks; i++){
    mean[i]   = fgaus->GetParameter(3*i+1);
    emean[i]  = fgaus->GetParError(3*i+1);
    rms[i]    = fgaus->GetParameter(3*i+2);
    number[i] = i;
  }
  
  //order peaks for lower to higher
  double aux1,aux2,aux3;
  for(int i = 0 ; i < npeaks ; i++){
    for(int j = i+1 ; j < npeaks ; j++){
      if(mean[i] > mean[j] && mean[j] != 0){
	aux1     = mean[i];
	mean[i]  = mean[j];
	mean[j]  = aux1;
	aux2     = emean[i];
	emean[i] = emean[j];
	emean[j] = aux2;
	aux3     = rms[i];
	rms[i]   = rms[j];
	rms[j]   = aux3;
      }
    }
  }

  //compute SNR
  SNR[0] = (mean[1]-mean[0])/rms[1];
  SNR[1] = (mean[1]-mean[0])/rms[0];
  SNR[2] = (mean[1]-mean[0])/sqrt(pow(rms[0],2)+pow(rms[1],2));
  
  //draw mean vs peak number for getting gaining
  TGraphErrors* tg = new TGraphErrors(NMAXPEAKSFIT,number,mean,0,emean);
  tg->Draw("ap");
  tg->GetXaxis()->SetTitle("#it{pe}");
  tg->GetYaxis()->SetTitle("#it{Charge} (Vs)");
  tg->SetMarkerStyle(20);
  tg->Fit("pol1","Q");
  if(plot_intermediate){
    gPad->Update();
    gPad->WaitPrimitive();
  }
  if(plot_results)gPad->Print(("gainfit_"+sOV.str()+"OV.pdf").c_str());
  
  //get gain and error
  TF1* fgain = (TF1*)tg->GetListOfFunctions()->FindObject("pol1");
  GAIN  = fgain->GetParameter(1)/(AMPFACTOR*ECHARGE);
  eGAIN = fgain->GetParError(1)/(AMPFACTOR*ECHARGE);

  //sanity chech for errors
  if(eGAIN>GAIN){
    std::cout << "problem computing error, assuming zero as error" << std::endl;
    eGAIN = 0;
  }

  //tg->Delete();
  //fgain->Delete();
}

//**********************************************************
//******ANALIZE DATAFILE, FILL CHARGE HISTOGRAM*************
//**********************************************************

//**********************************************************
Double_t integrateWf(double* time, double* V,
		     const int NFILES){
//**********************************************************
 
  //create wf tgraph and histogram for baseline
  TGraph* tg = new TGraph();
  TH1F* h    = new TH1F("h","h",1000,-0.1,0.1);

  double time_limit_up = -999;
  double max           = 0;

  //fill histogram and look for max
  for(int i = 0; i < WFLENGTH; i++){
    if(time[i] < time_limit_down)h->Fill(V[i]);
    if(V[i]>max)max = V[i];
  }

  //get baseline
  double baseline = h->GetMean();
  
  //fill tg substracting baseline and look for integration limit
  for(int i = 0; i < WFLENGTH; i++){
    tg->SetPoint(i,time[i],V[i]-baseline);
    if(time[i] > time_limit_up_min && (V[i]-baseline) < (max-baseline)*(1-PER) && time_limit_up == -999)time_limit_up = time[i];
  }

  //decide time limit
  if(integrateUpToPer){
    if(time_limit_up == -999)time_limit_up = time_limit_up_max;
    if(time_limit_up > time_limit_up_max)time_limit_up = time_limit_up_max;
  }
  else time_limit_up = time_limit_up_fix;//[NFILES%NSETFILES];

  //create the function to be integrated
  TF1* f;
  if(use_fit){
    if(SiPMUtils::isNoise(time,V))f = new TF1("f",[&](double*x, double *par){return tg->Eval(x[0]);},time_limit_down,time_limit_up,0);
    else {
      f = new TF1("f","[3]*([1]/2)*exp(([1]/2)*(2*[0]+[1]*[2]*[2]-2*x))*TMath::Erfc(([0]+[1]*[2]*[2]-x)/([2]*sqrt(2)))+[4]",time[0],time[WFLENGTH-1]);
      f->SetParameters(2e-7,1000000,5e-8,0.000000002,0.001);
      tg->Fit("f","Q","",time[0],1e-6);//tg->Draw("ap");gPad->Update();gPad->WaitPrimitive();
    }
  }
  else{
    //it is just an interpolation of the waveform
    f = new TF1("f",[&](double*x, double *par){return tg->Eval(x[0]);},time_limit_down,time_limit_up,0);
  }
  
  //numerically integrate it
  double val = f->Integral(time_limit_down,time_limit_up);
  if(integrateUpToPer && !SiPMUtils::isNoise(time,V))val = val*area_factor;
  
  //detele stuff
  h->Delete();
  f->Delete();
  tg->Delete();

  return val;
}
  
//**********************************************************
void analizeTree(const std::string& name,
		 TH1F* hq,
		 double* SNR,
		 double &GAIN, double &eGAIN,
		 const double OV,
		 const int NFILES){
//**********************************************************
  
  //the input file name
  std::string filename(name);
  std::cout << "------------------------------------------------------------------" << std::endl;
  std::cout << "running on file '" << filename << "'" << endl;

  //open root file
  TFile* rfile = TFile::Open(filename.c_str());
  TTree* wf    = (TTree*)rfile->Get("wf");

  //set addres for variables
  double fwftime         = 0;
  double ftime[WFLENGTH] = {0};
  double fV[WFLENGTH]    = {0};
  wf->SetBranchAddress("wftime", &fwftime);
  wf->SetBranchAddress("time"  , &ftime  );
  wf->SetBranchAddress("V"     , &fV     );
  
  //Get number of wf
  int nentries = wf->GetEntries();
  
  //loop over wf
  for(int i = 0; i < nentries; i++){
    wf->GetEntry(i);

    //check if there has been any DC in the waveform
    //if(lookForDCR){if(hasDCR(ftime,fV))continue;}
        
    //if not integrate each waveform and fill histogram
    hq->Fill(integrateWf(ftime,fV,NFILES));
    if(i%5000==0)std::cout << "looped over " << i+1 << "/" << nentries << " entries" << std::endl;
  }

  //compute gain from histogram
  computeGain(hq,SNR,GAIN,eGAIN,OV,NFILES);
  
}

//**********************************************************
//******MAIN FUNCTION***************************************
//**********************************************************

//**********************************************************
void gain_from_tree(){
//**********************************************************

  gStyle->SetOptStat(1);
  gStyle->SetOptFit(0);

  //listname
  std::string listname = "lists/gain.list";
  
  // Open the file
  FILE *pFile = fopen (listname.c_str(), "r");
  if( pFile == NULL ){
    std::cout << "Cannot open File '" <<  listname << "'" << std::endl;
    exit(1);
  }

  //variables for main plot
  double GAIN[NMAXFILES]             = {0};
  double eGAIN[NMAXFILES]            = {0};
  double V[NMAXFILES]                = {0};
  double OV[NMAXFILES]               = {0};
  double voltage[NMAXFILES]          = {0};
  double SNR[NMAXFILES][3]           = {0};
  
  //Read the list file
  char filename[200];
  int NFILES     = 0;
  int sipm       = 0;

  //main loop
  while(fscanf(pFile,"%s",filename) == 1){  
    if(filename[0] == '/' && filename[1] == '/')continue;

    //create histogram with variable limits, depending on SiPM pitch
    //select overvoltages depending on SiPM pitch
    if(NFILES == 0){
      sipm = SiPMUtils::getSiPMFromName(filename);
      SiPMUtils::initializeChargeHistogramByPitch(sipm,OV,hqmin1,hqmin2,hqmax1,hqmax2);
    }
    TH1F* hq = new TH1F("hq","hq",500,
			hqmin1+(NFILES%NSETFILES)*(hqmin2-hqmin1)/NSETFILES,
			hqmax1+(NFILES%NSETFILES)*(hqmax2-hqmax1)/NSETFILES);

    //getx operating voltage and decide wether to use absolute or relative
    SiPMUtils::getVoltageFromName(filename,V[NFILES]);
    if(use_abs_voltage)voltage[NFILES] = V[NFILES];
    else voltage[NFILES] = OV[NFILES];
   
    //analize datafile
    analizeTree(filename,hq,SNR[NFILES],GAIN[NFILES],eGAIN[NFILES],voltage[NFILES],NFILES);
    
    //clean memory
    hq->Delete();

    //increase file number
    NFILES++;
    if(NFILES >= NMAXFILES)break;
  }

  fclose(pFile);

  //print results
  if(print_results){
    std::cout << "------------------------------------------------------------------" << std::endl;
    std::cout << "RESULTS" << std::endl;
    std::cout << "------------------------------------------------------------------" << std::endl;
    std::cout << "GAIN \t\t ERROR" << std::endl;
    for(int i = 0; i < NFILES; i++){
      std::cout << GAIN[i] << " " << eGAIN[i] << std::endl;
    }
    std::cout << "SNR" << std::endl;
    for(int i = 0; i < NFILES; i++){
      std::cout << SNR[i][0] << " " << SNR[i][1] << " " << SNR[i][2] << std::endl;
    }
  }

  /*std::cout << "OV \t SNR1 \t\t SNR2 \t\t SNR3" << std::endl;
  for(int i = 0; i < NFILES; i++){
    std::cout << V[i] << " \t " << SNR[i][0] << " \t " << SNR[i][1] << " \t " << SNR[i][2] << std::endl;
    }*/

  //draw gain vs ov
  TCanvas* c1 = new TCanvas("c1","c1",900,600);
  TGraphErrors* tg = new TGraphErrors(NFILES,voltage,GAIN,0,eGAIN);
  tg->SetMarkerStyle(20);
  tg->GetXaxis()->SetTitle("#it{OV} (V)");
  tg->GetYaxis()->SetTitle("#it{GAIN}");
  tg->Draw("ap");
  if(plot_results)gPad->Print("GAINvsOV.pdf");

  /*//draw SNR
  TCanvas* c2 = new TCanvas("c2","c2",600,600);
  TGraph* tgsnr[NFILES];
  double peaks[NMAXPEAKS-1] = {0};
  for(int i = 0; i < NMAXPEAKS-1; i++)peaks[i] = i + 1;
  for(int i = 0; i < 3; i++){
    tgsnr[i] = new TGraph(NMAXPEAKS-1,peaks,SNR[i]);
    tgsnr[i]->SetLineColor(i+1);
    if(i == 0){
      tgsnr[i]->GetXaxis()->SetTitle("#it{PE}");
      tgsnr[i]->GetYaxis()->SetTitle("#it{SNR}");
      tgsnr[i]->GetYaxis()->SetRangeUser(5,40);
      tgsnr[i]->SetLineWidth(3);
      tgsnr[i]->Draw("al");
    }
    else{
      tgsnr[i]->SetLineWidth(3);
      tgsnr[i]->Draw("lsame");
    }
    }*/

}

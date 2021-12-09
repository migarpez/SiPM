#include "Utils.C"

//macro for ploting voltage vs intensity data and getting breakdown voltage
const int NMAXFILES     = 500;
const int NBOARDS       = 12;
const int NSIPMPERBOARD = 6;

const double NPIXELS50  = 14331;
const double NPIXELS75  = 6364;
const double NPIXELSFBK = 11188;
const double IMIN       = 5.0;
const double IMAX       = 6.2; 

//parameters for ploting/printing results
bool plot_intermediate = false;
//bool plot_results      = false;
//bool use_abs_voltage   = false;
//bool print_results     = true;

//**********************************************************
void getRQCell(char* name,
	       double RQ, double eRQ,
	       double& RQ_cell, double& eRQ_cell){
//**********************************************************

  if(SiPMUtils::isBoard(name)){
    //RQ_cell  = NPIXELS50*RQ;
    //eRQ_cell = NPIXELS50*eRQ;
    RQ_cell  = NPIXELSFBK*RQ;
    eRQ_cell = NPIXELSFBK*eRQ;
  }
  else if(SiPMUtils::is75Pitch(name) && !SiPMUtils::isFBK(name)){
    RQ_cell  = NPIXELS75*RQ;
    eRQ_cell = NPIXELS75*eRQ;
  }
  else if(!SiPMUtils::is75Pitch(name) && !SiPMUtils::isFBK(name)){
    RQ_cell  = NPIXELS50*RQ;
    eRQ_cell = NPIXELS50*eRQ;
  }
  else if(SiPMUtils::isFBK(name)){
    RQ_cell  = NPIXELSFBK*RQ;
    eRQ_cell = NPIXELSFBK*eRQ;
  }
}


//**********************************************************
void getFitLimits(double& VMIN, double& VMAX,
		  const TGraph* tg){
//**********************************************************

  //set funcion boundaries
  double fmin = TMath::MinElement(tg->GetN(),tg->GetX());
  double fmax = TMath::MaxElement(tg->GetN(),tg->GetX());

  //define function for minimum voltage
  TF1* tfmin = new TF1("tfmin",[&](double*x, double *par){return abs(tg->Eval(x[0])-IMIN);},fmin,fmax,0);
  VMIN = tfmin->GetMinimumX();

  //define function for maximum voltage
  TF1* tfmax = new TF1("tfmax",[&](double*x, double *par){return abs(tg->Eval(x[0])-IMAX);},fmin,fmax,0);
  VMAX = tfmax->GetMinimumX();

  tfmin->Delete();
  tfmax->Delete();
}

//**********************************************************
void analizeCurve(char* name,
		  TCanvas* c1,
		  double& RQ, double& eRQ){
//**********************************************************

  TGraph* tg = new TGraph();

  //the input file name for indirect curva
  std::string filename(name);
  
  // Open the file
  FILE *iFile = fopen (filename.c_str(),"r");
  if(iFile == NULL){
    cout << "Cannot open File '" << filename << "'" << endl;
    exit(1);
  }

  //std::cout << "running on file '" << filename << "'" << std::endl;
  
  //Read the input file and fill the i:v graph
  int i = 0;
  double v,c;
  char line[100];

  //skip header
  fscanf(iFile,"%[^\n]s",line);

  //read data
  while(fscanf(iFile,"%lf %lf",&v,&c) == 2){
    tg->SetPoint(i,v,c);
    i++;
  }

  fclose(iFile);

  //draw the graph
  tg->GetXaxis()->SetTitle("#it{V} (V)");
  tg->GetYaxis()->SetTitle("#it{I} (mA)");
  tg->GetYaxis()->SetRangeUser(0,20);
  //tg->SetMarkerStyle(20);
  //tg->SetMarkerSize(0.25);
  tg->Draw("apc");
  //gPad->WaitPrimitive();

  //get fit limits
  double fitmin,fitmax;
  getFitLimits(fitmin,fitmax,tg);

  //fit data
  tg->Fit("pol1","Q","",fitmin,fitmax);
  gPad->SetGridy();
  gPad->SetGridx();
  if(plot_intermediate){gPad->Update();gPad->WaitPrimitive();}
  
  //compute RQ and error
  TF1* fit = (TF1*)tg->GetListOfFunctions()->FindObject("pol1");

  //in KOhms
  RQ  = 1000000000/(fit->GetParameter(1))*pow(10,-6);
  eRQ = 1000000000/(pow(fit->GetParameter(1),2))*(fit->GetParError(1))*pow(10,-6);
}

//**********************************************************
void forward(){
//**********************************************************
  //the input file name
  std::string listname= "lists/forward.list";
  
  //define canvas and other stuff
  TCanvas* c1 = new TCanvas("c1","c1",600,600);
  TCanvas* c2 = new TCanvas("c2","c2",600,600);
  c2->SetGrid();

  //define variables
  double RQ[NMAXFILES]       = {0};
  double eRQ[NMAXFILES]      = {0};
  double RQ_cell[NMAXFILES]  = {0};
  double eRQ_cell[NMAXFILES] = {0};
  double sipm[NMAXFILES]     = {0};
  double board[NMAXFILES]    = {0};

  // Open the file
  FILE *pFile = fopen (listname.c_str(),"r");
  if(pFile == NULL){
    std::cout << "Cannot open File '" <<  listname << "'" << std::endl;;
    exit(1);
  }

  //Read the input file
  c1->cd();
  char file_name[200];
  int counter = 0;

  std::cout << "Board \t \t SiPM \t \t" << "RQ \t \t \t" << "RQ/cell" << std::endl;
  while(fscanf(pFile,"%s",file_name) == 1 && counter < NMAXFILES){

    if(file_name[0] == '/' && file_name[1] == '/')continue;

    sipm[counter]  = SiPMUtils::getSiPMFromName(file_name);
    board[counter] = SiPMUtils::getBoardFromName(file_name);

    analizeCurve(file_name,c1,RQ[counter],eRQ[counter]);
    //getRQCell(file_name,RQ[counter],eRQ[counter],RQ_cell[counter],eRQ_cell[counter]);
    RQ_cell[counter] = RQ[counter]/11188;
    eRQ_cell[counter] = eRQ[counter]/11188;
    
    
    std::cout << RQ[counter] << std::endl;
    
    counter++;
  }  
  
  /*c2->cd();
  TGraphErrors* tg = new TGraphErrors(counter,sipm,RQ_cell,0,eRQ_cell);
  tg->GetXaxis()->SetTitle("#it{SiPM}");
  tg->GetYaxis()->SetTitle("#it{R_{q}} (M#Omega)");
  tg->GetYaxis()->SetMaxDigits(2);
  tg->SetMarkerStyle(20);
  tg->Draw("ap");
  gStyle->SetOptFit();
  //tg->Fit("pol1");0

  //fill a common histogram
  TH1F* ht = new TH1F("ht","ht",50,30,80);
  for(int i = 0; i < counter; i++)ht->Fill(RQ[i]);//hr[i]->Fill(RQ_cell[i*6+j]);
  
  //draw stacked histogram for boards
  THStack* hs = new THStack("hs",""); 
  TH1F* hr[NBOARDS];
  TLegend* lg = new TLegend (0.52,0.72,0.9,0.9);

  ht->GetXaxis()->SetTitle("#it{R_{Q}} (#Omega)");
  ht->Draw();

  std::string names[NBOARDS] = {"M3_7","P11_2","P11_3","P11_4","P11_5","P11_6","P11_8","P11_9","P11_10","P12_1","P12_2","P12_3"};
  
  for(int i = 0; i < NBOARDS; i++){
    std::stringstream ssi, ssii;
    ssi << i;
    ssii << i+1;
    //hr[i] = new TH1F(("hb"+ssi.str()+"").c_str(),("hb"+ssi.str()+"").c_str(),50,380,450);
    hr[i] = new TH1F(("hb"+ssi.str()+"").c_str(),("hb"+ssi.str()+"").c_str(),50,30,80);
    for(int j = 0; j < NSIPMPERBOARD; j++)hr[i]->Fill(RQ[i*6+j]);//hr[i]->Fill(RQ_cell[i*6+j]);
    hr[i]->SetFillColor(i+1);
    lg->AddEntry(hr[i],(names[i]).c_str(),"f");
    hs->Add(hr[i]);
  }
  hs->Draw("same");
  hs->GetXaxis()->SetTitle("#it{R_{Q}} (#Omega)");
  lg->Draw();*/
}

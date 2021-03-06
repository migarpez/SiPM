//macro for ploting voltage vs intensity data and getting breakdown voltage

#include "Utils.C"

const int NMAXFILES     = 500;
const int NBOARDS       = 12;
const int NSIPMPERBOARD = 6;

const double VMIN = 26.5;
const double VMAX = 27.5;
const double RANGE = 0.08 ;

bool plot_intermediate = true;

//**********************************************************
void analizeCurve(char* name,
		  TCanvas* c1,TCanvas* c2,
		  double& BKV, double& eBKV){
//**********************************************************
  
  //the input file name for indirect curva
  std::string filename(name);

  //tgraphs
  TGraph* tg  = new TGraph();
  TGraph* tgi = new TGraph();
  TGraph* tgf = new TGraph();
  
  // Open the file
  FILE *iFile = fopen (filename.c_str(),"r");
  if(iFile == NULL){
    cout << "Cannot open File '" << filename << "'" << endl;
    exit(1);
  }

  //skip header
  char line[100];
  fscanf(iFile,"%[^\n]s",line);
  
  // Read the input file and fill the I and the 1/I graph
  int i = 0;
  double v,c;
  
  while(fscanf(iFile,"%lf %lf",&v,&c) == 2){
    tg ->SetPoint(i,v,c);
    tgi->SetPoint(i,v,1/c);
    i++;
  }

  //draw IvsV graph
  c1->cd();
  tg->GetXaxis()->SetTitle("#it{Voltage} (V)");
  tg->GetYaxis()->SetTitle("#it{I} (mA)");
  tg->GetYaxis()->SetTitleOffset(1.4);
  tg->Draw("ac");
  //if(plot_intermediate){gPad->Update();gPad->WaitPrimitive();}

  //close file
  fclose(iFile);
  
  //define IvsV function to compute derivative
  TF1* IvsV = new TF1("IvsV",[&](double*x, double *par){return tg->Eval(x[0]);},VMIN,VMAX,0);

  //compute final function
  double inverse, voltage, value, maxvalue = 0, maxvoltage;
  int counter = 0;
  for(int j = 0; j < i; j++){
    tgi->GetPoint(j,voltage,inverse);
    value = inverse*IvsV->Derivative(voltage);
    if(voltage < VMIN || voltage > VMAX)continue;
    tgf->SetPoint(counter,voltage,value);
    if(value > maxvalue && voltage>VMIN){
      maxvalue = value;
      maxvoltage = voltage;
    }
    counter++;
  }

  c2->cd();
  tgf->GetXaxis()->SetTitle("#it{Voltage} (V)");
  tgf->GetXaxis()->SetRangeUser(VMIN,VMAX);
  tgf->GetYaxis()->SetTitle("#it{#frac{dI}{IdV}} (V^{-1})");
  tgf->Draw("al");
  tgf->Fit("landau","Q","",maxvoltage-RANGE,maxvoltage+RANGE*1.5);

  BKV = static_cast<TF1*>(tgf->GetFunction("landau"))->GetParameter(1);
  eBKV = static_cast<TF1*>(tgf->GetFunction("landau"))->GetParameter(2);

  if(plot_intermediate){gPad->Update();gPad->WaitPrimitive();}

  //delete created things
  tg->Delete();
  tgi->Delete();
  tgf->Delete();
  IvsV->Delete();
}

//**********************************************************
void inverse(){
//**********************************************************

  //gStyle->SetOptStat(1111);
  
  //the input file name
  std::string listname= "lists/inverse.list";
  
  //define canvas
  TCanvas* c1 = new TCanvas("c1","c1",600,600);
  TCanvas* c2 = new TCanvas("c2","c2",600,600);
  TCanvas* c3 = new TCanvas("c3","c3",600,600);
  c3->SetGrid();

  //define variables
  double BKV[NMAXFILES]   = {0};
  double eBKV[NMAXFILES]  = {0};
  double sipm[NMAXFILES]  = {0};
  double board[NMAXFILES]  = {0};

  //Open the file
  FILE *pFile = fopen (listname.c_str(),"r");
  if(pFile == NULL){
    std::cout << "Cannot open File '" <<  listname << "'" << std::endl;;
    exit(1);
  }

  //Read the input file
  c1->cd();
  char file_name[200];
  int counter = 0;

  std::cout << "SiPM \t BV(V) \t error" << std::endl;
  
  while(fscanf(pFile,"%s",file_name) == 1 && counter < NMAXFILES){

    if(file_name[0] == '/' && file_name[1] == '/')continue;
    sipm[counter] = SiPMUtils::getSiPMFromName(file_name);
    board[counter] = SiPMUtils::getBoardFromName(file_name);

    analizeCurve(file_name,c1,c2,BKV[counter],eBKV[counter]);
    
    //std::cout << board[counter] << " " << sipm[counter] << " " << BKV[counter] << " " << eBKV[counter] << std::endl;
    std::cout << BKV[counter] << std::endl;
    
    counter++;
  }
  
  /*c3->cd();
  TGraphErrors* tg = new TGraphErrors(counter,sipm,BKV,0,eBKV);
  tg->GetXaxis()->SetTitle("#it{SiPM}");
  tg->GetYaxis()->SetTitle("#it{V_{breakdown}} (V)");
  tg->GetYaxis()->SetMaxDigits(2);
  tg->SetMarkerStyle(20);
  tg->Draw("ap");
  gStyle->SetOptFit();
  //tg->Fit("pol1");
  
  TH1F* hr = new TH1F("hr","hr",50,32.5,33);
  for(int i = 0; i < counter; i++)hr->Fill(BKV[i]);
  hr->GetXaxis()->SetTitle("#it{V_{B}} (V)");
  hr->Draw();

  //draw stacked histogram for boards
  THStack* hs = new THStack("hs",""); 
  TH1F* hb[NBOARDS];
  TLegend* lg = new TLegend (0.52,0.72,0.9,0.9);

  std::string names[NBOARDS] = {"M3_7","P11_2","P11_3","P11_4","P11_5","P11_6","P11_8","P11_9","P11_10","P12_1","P12_2","P12_3"};
  
  for(int i = 0; i < NBOARDS; i++){
    std::stringstream ssi, ssii;
    ssi << i;
    ssii << i+1;
    hb[i] = new TH1F(("hb"+ssi.str()+"").c_str(),("hb"+ssi.str()+"").c_str(),50,32.6,33.6);
    for(int j = 0; j < NSIPMPERBOARD; j++)hb[i]->Fill(BKV[i*6+j]);
    hb[i]->SetFillColor(i+1);
    lg->AddEntry(hb[i],(names[i]).c_str(),"f");
    hs->Add(hb[i]);
  }
  hs->Draw("same");
  hs->GetXaxis()->SetTitle("#it{V_{B}} (V)");
  lg->Draw();

  hr->Draw();*/
}

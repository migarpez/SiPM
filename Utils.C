namespace SiPMUtils{

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
  int getSiPMFromName(char* name){
  //**********************************************************

    //the input file name
    std::string filename(name);
    
    //get sipm from name
    if(filename.find("sipm") != -1){
      return stod(filename.substr(filename.find("sipm")+4,2));
    }
    else if(filename.find("fbk") != -1){
      return stod(filename.substr(filename.find("fbk")+3,2));
    }
    else return -1;
  }
  
  //**********************************************************
  bool is75Pitch(const int sipm){
  //**********************************************************

    const int SIPM75[] = {11,15,16,17,19,21,12}; //12 is from Madrid
    
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
  bool is75Pitch(char* name){
  //**********************************************************

    return is75Pitch(getSiPMFromName(name));
  }

  //**********************************************************
  bool isBoard(char* name){
  //**********************************************************

    //the input file name
    std::string filename(name);

    bool isIt = false;
    
    //look for board in name
    if(filename.find("board") != -1){
      isIt = true;
    }

    return isIt;
  }

  //**********************************************************
  int getBoardFromName(char* name){
  //**********************************************************

    //the input file name
    std::string filename(name);
    
    //get sipm from name
    if(filename.find("board") != -1){
      return stod(filename.substr(filename.find("board")+5,2));
    }
    else return -1;
  }

  //**********************************************************
  bool isNoise(double* time, double* V){
  //**********************************************************

    bool isNoise = false;
    TH1F* hy = new TH1F("hy","hy",100,-0.05,0.05);
    for(int i = 0; i < sizeof(time); i++)hy->Fill(V[i]);
    double baseline = hy->GetMean();
    
    for(int i = 0; i < sizeof(time); i++)if(time[i] > 0.2e-6 && time[i] < 0.4e-6 && V[i] < baseline){
	isNoise = true;
	break;
      }

    hy->Delete();
    return isNoise;
  }

  //**********************************************************
  void initializeChargeHistogramByPitch(const int board, const int sipm,
					double* OV,
					double &hqmin1, double &hqmin2,
					double &hqmax1, double &hqmax2){
  //**********************************************************

    const double hqmin1_50 = -10e-9;
    const double hqmin2_50 = -20e-9;
    const double hqmax1_50 = 100e-9;
    const double hqmax2_50 = 140e-9;
    const double hqmin1_75 = -40e-9;
    const double hqmin2_75 = -40e-9;
    const double hqmax1_75 = 500e-9;
    const double hqmax2_75 = 500e-9;

    const double OV50[] = {3,4  ,5};
    const double OV75[] = {2,2.5,3};
    
    if(SiPMUtils::is75Pitch(sipm) && board==-1){
      hqmin1 = hqmin1_75;
      hqmin2 = hqmin2_75;
      hqmax1 = hqmax1_75;
      hqmax2 = hqmax2_75;
      for(int i = 0; i < sizeof(OV); i++)OV[i] = OV75[i];
      std::cout << "------------------------------------------------------------------" << std::endl;
      std::cout << "Analyzing GAIN from 75 um pitch sipm " << sipm << std::endl;
      std::cout << "------------------------------------------------------------------" << std::endl;
    }
    else{
      hqmin1 = hqmin1_50;
      hqmin2 = hqmin2_50;
      hqmax1 = hqmax1_50;
      hqmax2 = hqmax2_50;
      for(int i = 0; i < sizeof(OV); i++)OV[i] = OV50[i];
      std::cout << "------------------------------------------------------------------" << std::endl;
      std::cout << "Analyzing GAIN from 50 um pitch sipm " << sipm << std::endl;
    std::cout << "------------------------------------------------------------------" << std::endl;
    }
  }

  //**********************************************************
  bool isColdFromName(const std::string& name){
  //**********************************************************

  //the input file name
    std::string filename(name);
    
    //look for 'cold' char
    if(filename.find("cold")!=-1)return true;
    else return false;
  }

  //**********************************************************
  bool isFBK(char* name){
  //**********************************************************

    //the input file name
    std::string filename(name);
    
    if(filename.find("fbk") != -1)return true;
    else return false;
  }

  //**********************************************************
  //bool is75Pitch(const int sipm){
  //**********************************************************
  /*
    bool isIt = false;
    for(int i = 0; i < sizeof(SIPM75); i++){
      if(sipm == SIPM75[i]){
	isIt = true;
	break;
      }
    }
    return isIt;
    }*/
  
}


//**********************************************************
//bool hasDCR(double* time, double* V){
//**********************************************************
/*
  TH1F* h   = new TH1F("h","h",WFLENGTH,time[0],time[WFLENGTH-1]);
  TH1F* hbl = new TH1F("hbl","hbl",200     ,-0.05  ,0.05);

  //fill wf and baseline
  for(int i = 0; i < WFLENGTH; i++){
    h->SetBinContent(i,V[i]);
    hbl->Fill(V[i]);
  }

  //create an instance of TSpectrum
  TSpectrum *s = new TSpectrum(NMAXPEAKS);
  
  //go for peaks
  //int npeaks = s->Search(h,0.1,"",0.5);
  int npeaks = s->Search(h,0.000001,"",0.2);
  double* Vpeak = s->GetPositionY();

  h->Draw();gPad->WaitPrimitive();
  //Get baseline rms
  double blmean = hbl->GetMean();
  double blrms  = hbl->GetRMS();
  double noise_limit = blmean + 3*blrms;
  
  if(npeaks > 1 && Vpeak[0] > noise_limit && Vpeak[1] > noise_limit){
    std::cout << noise_limit << std::endl;
    h->Draw();gPad->WaitPrimitive();
    h->Delete();
    hbl->Delete();
    return true;
  }

  else{
    h->Delete();
    hbl->Delete();
    return false;
  }
  
}*/

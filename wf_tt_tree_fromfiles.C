//macro for storing waveform data into a tree
//It creates a single tree from different input files

const int NMAXFILES         = 10;
const int WFLENGHT          = 1250;
const double default_wftime = 0.000002;
bool use_default            = false;

//**********************************************************
void filltree(TTree* wf,
	      const std::string& filename,
	      double &fwftime, double *ftime, double *fV){
//**********************************************************  

  std::cout << "Processing file '" << filename << "'" << std::endl;

  // Open the file
  FILE *pFile = fopen (filename.c_str(), "r");
  if( pFile == NULL ){
    std::cout << "Cannot open File '" <<  filename << "'" << std::endl;
    exit(1);
  }

  //read header and get relevant info
  int NFRAMES,NPOINTS;
  char line[100];
  while(fscanf(pFile,"%s",line) == 1){
    std::string sline(line);

    if(sline.find( "Length," ) != -1){
      NPOINTS = stoi(sline.substr(filename.find("Length,")+8,filename.size()-(filename.find("Length,")+8)));
      if(WFLENGHT != NPOINTS){
	std::cout << "WARNING! Waveform length in file does not match with WFLENGHT value! Please, change it" <<  std::endl;
	exit(1);
      }
    }
    
    if(sline.find( "Count," ) != -1){
      NFRAMES = stoi(sline.substr(filename.find("Count,")+7,filename.size()-(filename.find("Count,")+7)));
      std::cout << "This file has " << NFRAMES << " waveforms"<< std::endl;
    }

    if(sline.find( "TIME" ) != -1)break;
  }

  //prepare time trend file
  std::string ttname = filename;
  ttname.replace(ttname.find("wf"),2,"tt");

  //check wheter it exists or not
  FILE *tFile = fopen (ttname.c_str(), "r");
  if(tFile == NULL){
    use_default = true;
    fwftime = default_wftime;
    std::cout << "Not time trend associated file found, using default time between waveform" << std::endl;
    std::cout << "default time = " << default_wftime << " s" << std::endl;
  }

  //skip time trend file header
  if(!use_default)fgets(line,100,tFile);
 
  //fill tree
  double wft,t,v;

  //read waveforms. Each one is a new entry in the ttree
  for(int i = 0; i < NFRAMES; i++){
    for(int j = 0; j < WFLENGHT; j++){
      fscanf(pFile,"%lf,%lf",&t,&v);
      ftime[j] = t;
      fV[j] = v;
    }
    if(!use_default){
      fscanf(tFile,"%*lf,%lf",&wft);
      fwftime = wft;
    }
    if(i%100==0)std::cout << i+1 << "/" << NFRAMES << " frames stored" << std::endl;
    wf->Fill();
  }

  if(tFile)fclose(tFile);
  fclose(pFile);
}

//**********************************************************
void wf_tt_tree_fromfiles(){
//**********************************************************
  
  std::string listname = "lists/wf_tt_fromfiles.list";
  
  // Open the file
  FILE *pFile = fopen (listname.c_str(), "r");
  if( pFile == NULL ){
    std::cout << "Cannot open File '" <<  listname << "'" << std::endl;
    exit(1);
  }

  //create output tree
  const std::string rootfilename = "/home/miguel/Documents/PhD/PD/first_split_tests/HQR50/sipm24/dcr/after/2020_08_05_sipm26_cold_dcr_V45.65.root";
  
  //define tree variables
  double fwftime,ftime[WFLENGHT],fV[WFLENGHT];

  //Create and open the output root file
  TFile* rfile = new TFile(rootfilename.c_str(),"NEW");
  
  //Define the output tree
  TTree* wf = new TTree("wf","wf");
  std::stringstream sWFLENGHT;
  sWFLENGHT << WFLENGHT;
  wf->Branch("wftime", &fwftime, "wftime/D"    );
  wf->Branch("time"  , &ftime  , ("time["+sWFLENGHT.str()+"]/D").c_str());
  wf->Branch("V"     , &fV     , ("V["+sWFLENGHT.str()+"]/D").c_str()   );
  
  wf->SetBranchAddress("wftime", &fwftime);
  wf->SetBranchAddress("time"  , &ftime  );
  wf->SetBranchAddress("V"     , &fV     );

  char filename[200];

  //Read the list file  
  while(fscanf(pFile,"%s",filename) == 1){
    filltree(wf,filename,fwftime,ftime,fV);
  }    
    
  wf->Write();
  rfile->Close();
  fclose(pFile);

  std::cout << "ROOT file stored as '" << rootfilename << "'" << std::endl;
}

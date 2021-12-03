//macro for storing waveform data into a tree
//it creates a different tree for each input file
//IMPORTANT! the WFLENGTH parameter should match the input file!!

const int NMAXFILES         = 100;
const int WFLENGHT          = 1500; //short files, gain-like, 1250. Long files, cnoise-like, 3125.
const double default_wftime = 0.000002;
bool use_default            = false; //if false, look for a tt-file. if not found, default will be used

//**********************************************************
void readHeader(const std::string& filename,
		FILE* pFile,
		int &NFRAMES, int &NPOINTS,
		bool &isLengthCorrect){
//**********************************************************  

  char line[100];
  while(fscanf(pFile,"%s",line) == 1){
    std::string sline(line);

    if(sline.find( "Length," ) != -1){
      NPOINTS = stoi(sline.substr(filename.find("Length,")+8,filename.size()-(filename.find("Length,")+8)));
      if(WFLENGHT != NPOINTS){
	isLengthCorrect = false;
	break;
      }
    }
    
    if(sline.find( "Count," ) != -1){
      NFRAMES = stoi(sline.substr(filename.find("Count,")+7,filename.size()-(filename.find("Count,")+7)));
      std::cout << "This file has " << NFRAMES << " waveforms"<< std::endl;
    }

    if(sline.find( "TIME" ) != -1)break;
  }
}

//**********************************************************
void fillTree(const std::string& filename){	      
//**********************************************************  

  std::cout << "------------------------------------------------------------------" << std::endl;
  std::cout << "Processing file '" << filename << "'" << std::endl;

  //Open the file
  FILE *pFile = fopen (filename.c_str(), "r");
  if( pFile == NULL ){
    std::cout << "Cannot open File '" <<  filename << "'" << std::endl;
    exit(1);
  }

  //read header
  bool isLenghtCorrect = true;
  int NFRAMES,NPOINTS;
  readHeader(filename,pFile,NFRAMES,NPOINTS,isLenghtCorrect);
  if(!isLenghtCorrect){
    std::cout << "WARNING! Waveform length in file does not match with WFLENGHT value! Please, change it" <<  std::endl;
    std::cout << "Waveform length in file is " << NPOINTS << " and WFLENGHT in macro is " << WFLENGHT << std::endl;
    exit(1);
  }
  
  //create output tree
  std::string rootfilename = filename;
  rootfilename.replace(rootfilename.find("wf"),3,"");
  rootfilename.replace(rootfilename.find("csv"),3,"root");
  std::cout << "Creating rootfile '" << rootfilename << "'" << std::endl;

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

  //prepare time trend file
  std::string ttname = filename;
  ttname.replace(ttname.find("wf"),2,"tt");

  //check wheter it exists or not
  FILE *tFile = fopen (ttname.c_str(), "r");
  if(tFile == NULL){
    use_default = true;
    fwftime = default_wftime;
    std::cout << "Not time trend associated file found, using default time between waveforms" << std::endl;
    std::cout << "default time = " << default_wftime << " s" << std::endl;
  }

  //skip time trend file header
  char line[100];
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
    if(i%1000==0)std::cout << i+1 << "/" << NFRAMES << " frames stored" << std::endl;
    wf->Fill();
  }

  //close files
  if(tFile)fclose(tFile);
  fclose(pFile);

  //write tree and close root file
  wf->Write();
  rfile->Close();

  std::cout << "ROOT file stored" << std::endl;
  std::cout << "-----------------------------------" << std::endl;
}

//**********************************************************
void wf_tt_tree(){
//**********************************************************
  
  std::string listname = "lists/wf_tt.list";
  
  // Open the file
  FILE *pFile = fopen (listname.c_str(), "r");
  if( pFile == NULL ){
    std::cout << "Cannot open File '" <<  listname << "'" << std::endl;
    exit(1);
  }

  char filename[200];

  //Read the list file  
  while(fscanf(pFile,"%s",filename) == 1){
    if(filename[0] == '/' && filename[1] == '/')continue;
    fillTree(filename);
  }    
    
  fclose(pFile);
}

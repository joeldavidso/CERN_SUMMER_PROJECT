#include <cmath>
#include <algorithm>


#include <sstream>
#include <string>
#include <iostream>

template <typename Type> std::string to_str(const Type & t) // Used for converting floats to strings without the annoying .00000000
{
  std::ostringstream os;
  os << t;
  return os.str ();
}


TString datadirmain = "/eos/atlas/atlascerngroupdisk/phys-exotics/jdm/svjets-schannel/"; // Strings for accessing the ntuples
TString datadirsub = "v1.0/";
TString filename1 = "user.ebusch.5085";
TString filename2 = ".MGPy8EG_SVJSChan_";
TString filevarM[2] = {"1500", "750"};
TString filevarr[2] = {"8", "3"};
TString filename3 = "_tree.root";
TString mcs[3] = {".mc16a", ".mc16d", ".mc16e"};
TString fileN[4] = {"47", "48", "49", "50"};
TString ndN[12] = {"3733", "3907", "3996", "3741", "3916", "4007", "3753", "3924", "4015", "3764", "3929", "4025"};

int COLOR[4] = {810,860,800,870}; // The Colours of the 4 signal points

TString multijetdir = "user.ebusch.364702.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ2WithSW."; // More location stuff but for multijet bkgrd
TString multijetndN[3] = {"198","165","147"};
TString multiNO[3] = {"22","04","10"};
TString fileN1[3] = {"98","65","47"};
TString fileN2[3] = {"22","04","10"};


TString tree_entry = "outTree/nominal";

    /////////////////////////////////////////////////////////////////
    ////                                                         ////
    ////           Different Variables For Later Use             ////
    ////                                                         ////
    /////////////////////////////////////////////////////////////////

int nsamples = 4;

double matchcut = 0.4; // The matching criteria for jets with dark quarks

int ptc = 20;     // The small R jet minimum pt
int ptlc = 200;   // The large R jet minimum pt

int nintcuts[2] = {20, 40};      // The region edges for the interactions per bunch crossing plots
int ptcuts[3] = {ptc,60,100};  // The region edges for the average jets close to dark quark at diff pt

int nbins = 5040; // number of bins for each histogram

    /////////////////////////////////////////////////////////////////
    ////                                                         ////
    ////                Start Main Body Of Code                  ////
    ////                                                         ////
    /////////////////////////////////////////////////////////////////

void project_Pv3(){  


  gInterpreter->GenerateDictionary("vector<vector<int> >", "vector");

    /////////////////////////////////////////////////////////////////
    ////                                                         ////
    ////                    Defines Histograms                   ////
    ////                                                         ////
    /////////////////////////////////////////////////////////////////


    TH1D* JETPT[nsamples][2]; // 0 -> Leading | 1 -> Subleading
    TH2D* JETPT2D[nsamples];
    TH1D* DQPT[nsamples][2]; // 0 -> Leading | 1 -> Subleading
    TH2D* DQPT2D[nsamples];
    TH1D* DRJJ[nsamples];
    TH1D* DRQQ[nsamples];
    TH2D* JETMATCHDR[nsamples][2];     // 0 -> Leading | 1 -> Subleading
    TH2D* JETMATCHINT[nsamples][2][3]; // 0 -> 1st int bin | 1 -> 2nd int bin | 2 -> 3rd int bin
    TH2D* JETMATCHPT[nsamples][2][3];  // 0 -> jet pt > ptc | 1 -> jet pt > 60 | 2 -> jet pt > 100

    TH1D* MAXDR[nsamples][2]; // 0 -> Leading Quark | 1 -> Subleading Quark
    TH1D* NVIS[nsamples];
    TH1D* NINV[nsamples];
    TH1D* NSTAB[nsamples];
    TH1D* NUNST[nsamples];
    TH1D* VISUNSTRAT[nsamples];
    TH1D* UNSTINVRAT[nsamples];
    TH1D* STABTOTRAT[nsamples];
    TH1D* INVTOTRAT[nsamples];
    
    TH1D* VDHPT[nsamples][2];        // 0 -> Leading Quark | 1 -> Subleading Quark
    TH1D* DPHIMATCHMET[nsamples][2]; // 0 -> Leading Quark | 1 -> Subleading Quark
    TH2D* DPHIMATCHMET2D[nsamples];   
    TH1D* DPHIJETMET[nsamples][3];   // 0 -> Leading Jet | 1 -> Subleading Jet | 2 -> Sum jet
    TH1D* MINDPHIMET[nsamples];
    TH1D* MAXDPHIMET[nsamples];

    TH1D* MTJETSMET[nsamples][2]; // 0 -> Aligned Jet + Anti-Aligned Jet + MET | 1 -> Leading Jet + Subleading Jet
    TH1D* SUMOFPT[nsamples];

    TH1D* DPHIJJ[nsamples];
    TH1D* DPHIQQ[nsamples];

    /////////////////////////////////////////////////////////////////
    ////                                                         ////
    ////                 Signal Grid Loop &                      ////
    ////                   ntuple access                         ////
    /////////////////////////////////////////////////////////////////

    for (int j = 0; j<nsamples; j++){
      int M = 0;
    int r = 0;
    if(j == 1) r = 1;
    if(j == 2) M = 1;
    if(j == 3){
      M = 1;
      r = 1;
    }
    


    int pass1 = 0;

    int success = 0;

    TChain *chain = new TChain(tree_entry);

    cout << "Z" << filevarM[M] << "r" << filevarr[r] << ":" << endl;
    for (int i = 0;i<3;i++){
      TString ndfilenam = "user.ebusch.2909"+ndN[3*j+i]+"._000001.tree.root";
      TString filenam = datadirmain+datadirsub+filename1+fileN[j]+filename2+filevarM[M]+"_"+filevarr[r]+mcs[i]+filename3;
      chain-> Add(filenam+"/"+ndfilenam);
      cout << mcs[i] << endl;
    }
    
    /////////////////////////////////////////////////////////////////
    ////                                                         ////
    ////                 Gets Variable Adresses                  ////
    ////                                                         ////
    /////////////////////////////////////////////////////////////////


    int ntruthBSM[22];
    chain->SetBranchAddress("ntruthBSM",&ntruthBSM);

    vector<int> *truthBSM_status = nullptr;
    chain->SetBranchAddress("truthBSM_status",&truthBSM_status);

    vector<int> *truthBSM_pdgId = nullptr;
    chain->SetBranchAddress("truthBSM_pdgId",&truthBSM_pdgId);

    vector<float> *truthBSM_pt = nullptr;
    chain->SetBranchAddress("truthBSM_pt",&truthBSM_pt);

    vector<float> *truthBSM_eta = nullptr;
    chain->SetBranchAddress("truthBSM_eta",&truthBSM_eta);

    vector<float> *truthBSM_phi = nullptr;
    chain->SetBranchAddress("truthBSM_phi",&truthBSM_phi);

    vector<float> *truthBSM_e = nullptr;
    chain->SetBranchAddress("truthBSM_e",&truthBSM_e);

    vector<float> *truthBSM_m = nullptr;
    chain->SetBranchAddress("truthBSM_m",&truthBSM_m);

    vector<int> *truthBSM_barcode = nullptr;
    chain->SetBranchAddress("truthBSM_barcode",&truthBSM_barcode);

    vector<int> *truthBSM_nChildren = nullptr;
    chain->SetBranchAddress("truthBSM_nChildren",&truthBSM_nChildren);

    vector<vector<int>> *truthBSM_child_status = nullptr;
    chain->SetBranchAddress("truthBSM_child_status",&truthBSM_child_status);

    vector<vector<int>> *truthBSM_child_pdgId = nullptr;
    chain->SetBranchAddress("truthBSM_child_pdgId",&truthBSM_child_pdgId);

    vector<vector<int>> *truthBSM_child_barcode = nullptr;
    chain->SetBranchAddress("truthBSM_child_barcode",&truthBSM_child_barcode);

    vector<int> *truthBSM_nParents = nullptr;
    chain->SetBranchAddress("truthBSM_nParents",&truthBSM_nParents);

    vector<vector<int>> *truthBSM_parent_status = nullptr;
    chain->SetBranchAddress("truthBSM_parent_status",&truthBSM_parent_status);

    vector<vector<int>> *truthBSM_parent_pdgId = nullptr;
    chain->SetBranchAddress("truthBSM_parent_pdgId",&truthBSM_parent_pdgId);

    int na4_pflowjets[22];
    chain->SetBranchAddress("na4_pflowjets",&na4_pflowjets);

    vector<float> *a4_pflowjets_pt = nullptr;
    chain->SetBranchAddress("a4_pflowjets_pt",&a4_pflowjets_pt);

    vector<float> *a4_pflowjets_eta = nullptr;
    chain->SetBranchAddress("a4_pflowjets_eta",&a4_pflowjets_eta);

    vector<float> *a4_pflowjets_phi = nullptr;
    chain->SetBranchAddress("a4_pflowjets_phi",&a4_pflowjets_phi);

    vector<float> *a4_pflowjets_E = nullptr;
    chain->SetBranchAddress("a4_pflowjets_E",&a4_pflowjets_E);

    int na10_lctopojets[22];
    chain->SetBranchAddress("na10_lctopojets",&na10_lctopojets);

    vector<float> *a10_lctopojets_pt = nullptr;
    chain->SetBranchAddress("a10_lctopojets_pt",&a10_lctopojets_pt);

    vector<float> *a10_lctopojets_eta = nullptr;
    chain->SetBranchAddress("a10_lctopojets_eta",&a10_lctopojets_eta);

    vector<float> *a10_lctopojets_phi = nullptr;
    chain->SetBranchAddress("a10_lctopojets_phi",&a10_lctopojets_phi);
    
    float metFinalTrk;
    chain->SetBranchAddress("metFinalTrk",&metFinalTrk);

    float metFinalTrkPx;
    chain->SetBranchAddress("metFinalTrkPx",&metFinalTrkPx);
    
    float metFinalTrkPy;
    chain->SetBranchAddress("metFinalTrkPy",&metFinalTrkPy);

    float metFinalTrkPhi;
    chain->SetBranchAddress("metFinalTrkPhi",&metFinalTrkPhi);

    float metFinalTrkSumEt;
    chain->SetBranchAddress("metFinalTrkSumEt",&metFinalTrkSumEt);

    float actualInteractionsPerCrossing;
    chain->SetBranchAddress("actualInteractionsPerCrossing",&actualInteractionsPerCrossing);

    /////////////////////////////////////////////////////////////////
    ////                                                         ////
    ////                    Create Histograms                    ////
    ////                                                         ////
    /////////////////////////////////////////////////////////////////



    JETPT[j][0] = new TH1D((const TString)("JETPT0"+to_str(j)),"Leading Jet p_{T} [GeV]",nbins,0,1000);
    JETPT[j][1] = new TH1D((const TString)("JETPT1"+to_str(j)),"Subleading Jet p_{T} [GeV]",nbins,0,700);
    JETPT2D[j] = new TH2D((const TString)("JETPT2D"+to_str(j)),"Leading Jet p_{T} VS Subleading Jet p_{T}",nbins,0,1000,nbins,0,700);
    DQPT[j][0] = new TH1D((const TString)("DQPT0"+to_str(j)),"Leading Dark Quark p_{T} [GeV]",nbins,0,1200);
    DQPT[j][1] = new TH1D((const TString)("DQPT1"+to_str(j)),"Subleading Dark Quark p_{T} [GeV]",nbins,0,1000);
    DQPT2D[j] = new TH2D((const TString)("DQPT2D"+to_str(j)),"Leading Dark Quark p_{T} VS SUbleading Dark Quaark p_{T}",nbins,0,1200,nbins,0,1000);
    DRJJ[j] = new TH1D((const TString)("DRJJ"+to_str(j)),"DR(Leading Jet, Subleading Jet)",nbins,0,6);   
    DRQQ[j] = new TH1D((const TString)("DRQQ"+to_str(j)),"DR(Leading Dark Quark, Subleading Dark Quark)",nbins,0,6);
    JETMATCHDR[j][0] = new TH2D((const TString)("JETMATCHDR0"+to_str(j)),"Average NO. Jets With DR(J,\\chi) < X",30,0,2,10,0,10);
    JETMATCHDR[j][1] = new TH2D((const TString)("JETMATCHDR1"+to_str(j)),"Average NO. Jets With DR(J,\\chi) < X",30,0,2,10,0,10);
    JETMATCHINT[j][0][0] = new TH2D((const TString)("JETMATCHINT00"+to_str(j)),"",30,0,2,10,0,10);
    JETMATCHINT[j][0][1] = new TH2D((const TString)("JETMATCHINT01"+to_str(j)),"",30,0,2,10,0,10);
    JETMATCHINT[j][0][2] = new TH2D((const TString)("JETMATCHINT02"+to_str(j)),"",30,0,2,10,0,10);
    JETMATCHINT[j][1][0] = new TH2D((const TString)("JETMATCHINT10"+to_str(j)),"",30,0,2,10,0,10);
    JETMATCHINT[j][1][1] = new TH2D((const TString)("JETMATCHINT11"+to_str(j)),"",30,0,2,10,0,10);
    JETMATCHINT[j][1][2] = new TH2D((const TString)("JETMATCHINT12"+to_str(j)),"",30,0,2,10,0,10);
    JETMATCHPT[j][0][0] = new TH2D((const TString)("JETMATCHPT00"+to_str(j)),"",30,0,2,10,0,10);
    JETMATCHPT[j][0][1] = new TH2D((const TString)("JETMATCHPT01"+to_str(j)),"",30,0,2,10,0,10);
    JETMATCHPT[j][0][2] = new TH2D((const TString)("JETMATCHPT02"+to_str(j)),"",30,0,2,10,0,10);
    JETMATCHPT[j][1][0] = new TH2D((const TString)("JETMATCHPT10"+to_str(j)),"",30,0,2,10,0,10);
    JETMATCHPT[j][1][1] = new TH2D((const TString)("JETMATCHPT11"+to_str(j)),"",30,0,2,10,0,10);
    JETMATCHPT[j][1][2] = new TH2D((const TString)("JETMATCHPT12"+to_str(j)),"",30,0,2,10,0,10);

    MAXDR[j][0] = new TH1D((const TString)("MAXDR0"+to_str(j)),"Max DR(VDH_{i}, VDH_{j}) From Leading Dark Quark",nbins,0,3);
    MAXDR[j][1] = new TH1D((const TString)("MAXDR1"+to_str(j)),"Max DR(VDH_{i}, VDH_{j}) From Subleading Dark Quark",nbins,0,3);
    NVIS[j] = new TH1D((const TString)("NVIS"+to_str(j)),"NO. Visible Dark Hadrons",nbins,0,10);
    NINV[j] = new TH1D((const TString)("NINV"+to_str(j)),"NO. Invisible Dark Hadrons",nbins,0,20);
    NSTAB[j] = new TH1D((const TString)("NSTAB"+to_str(j)),"NO. Stable Dark Hadrons",nbins,0,18);
    NUNST[j] = new TH1D((const TString)("NUNST"+to_str(j)),"NO. Unstable Dark Hadrons",nbins,0,18);
    VISUNSTRAT[j] = new TH1D((const TString)("VISUNSTRAT"+to_str(j)),"Fraction Of Unstable Dark AHdrons Which Are Visible",nbins,0,1);
    UNSTINVRAT[j] = new TH1D((const TString)("UNSTINVRAT"+to_str(j)),"Fraction Of Invisible Dark Hadrons Which Are Unstable",nbins,0,1);
    STABTOTRAT[j] = new TH1D((const TString)("STABTOTRAT"+to_str(j)),"Fraction Of Total Dark Hadrons Which Are Stable",nbins,0,1);
    INVTOTRAT[j] = new TH1D((const TString)("INVTOTRAT"+to_str(j)),"Fraction Of Total Dark Hadrons Which Are Invisible",nbins,0,1);

    VDHPT[j][0] = new TH1D((const TString)("VDHPT0"+to_str(j)),"Visible Dark Hadron p_{T} From The Leading Dark Quark [GeV]",nbins,0,1000);
    VDHPT[j][1] = new TH1D((const TString)("VDHPT1"+to_str(j)),"Visible Dark Hadron p_{T} From The Subleading Dark Quark [GeV]",nbins,0,1000);
    DPHIMATCHMET[j][0] = new TH1D((const TString)("DPHIMATCHMET0"+to_str(j)),"DPhi(Matched Jet, MET) From Leading Dark Quark",nbins,0,TMath::Pi());
    DPHIMATCHMET[j][1] = new TH1D((const TString)("DPHIMATCHMET1"+to_str(j)),"DPhi(Matched Jet, MET) From Subleading Dark Quark",nbins,0,TMath::Pi());
    DPHIMATCHMET2D[j] = new TH2D((const TString)("DPHIMATCHMET2D"+to_str(j)),"DPhi(Matched Jet, MET), Leading Dark Quark VS Subleading Dark Quark",nbins,0,TMath::Pi(),nbins,0,TMath::Pi());
    DPHIJETMET[j][0] = new TH1D((const TString)("DPHIJETMET0"+to_str(j)),"DPhi(Leading Jet, MET)",nbins,0,TMath::Pi());
    DPHIJETMET[j][1] = new TH1D((const TString)("DPHIJETMET1"+to_str(j)),"DPhi(Subleading Jet, MET)",nbins,0,TMath::Pi());
    MINDPHIMET[j] = new TH1D((const TString)("MINDPHIMET"+to_str(j)),"Minimum DPhi(Jet, MET) For All Jets",nbins,0,TMath::Pi());
    MAXDPHIMET[j] = new TH1D((const TString)("MAXDPHIMET"+to_str(j)),"Maximum Dphi(Jet, MET) For All Jets",nbins,0,TMath::Pi());
    
    MTJETSMET[j][0] = new TH1D((const TString)("MTJETSMET0"+to_str(j)),"Transverse Mass Of Aligned Jet + Anti-Aligned Jet + MET",nbins,0,2000);
    MTJETSMET[j][1] = new TH1D((const TString)("MTJETSMET1"+to_str(j)),"Transverse Mass Of Leading Jet + Subleading Jet + MET",nbins,0,2000);
    SUMOFPT[j] = new TH1D((const TString)("SUMOFPT"+to_str(j)),"h_{T} (Scalar Sum Of Jet p_{T}) [GeV]",nbins,0,2000);

    DPHIJJ[j] = new TH1D((const TString)("DPHIJJ"+to_str(j)),"DPhi(Leading Jet, Subleading Jet)",nbins,0,TMath::Pi());
    DPHIQQ[j] = new TH1D((const TString)("DPHIQQ"+to_str(j)),"DPhi(Leading Dark Quark, Subleading Dark Quark)",nbins,0,TMath::Pi());
    
    /////////////////////////////////////////////////////////////////
    ////                                                         ////
    ////                        Event Loop                       ////
    ////                                                         ////
    /////////////////////////////////////////////////////////////////


    int passes = 0;
    int nout = 0;



    TH1D* matchescut = new TH1D((const TString)("matchescut"+std::to_string(j)),"",4,0,4);  // BIN1 -> 0JETS | BIN2 -> JET1 | BIN3 -> JET2 | BIN4 -> BOTHJETS

    TH1D* matchesalignedcut = new TH1D((const TString)("matchesalignedcut"+std::to_string(j)),"",5,0,5); // BIN1 -> 0JETS | BIN2 -> J||MET | BIN3 -> J!|MET | BIN4 -> BOTHJETS | BIN5 -> J1 == J!|MET

    TH1D* matchesmetcut = new TH1D((const TString)("matchesmetcut"+std::to_string(j)),"",6,0,6);  // BIN1,2 -> NONE / BIN3,4 -> J1||MET / BIN5,6 -> J2||MET

    TH1D* matchesmetaacut = new TH1D((const TString)("matchesmetaacut"+std::to_string(j)),"",6,0,6);  // BIN1,2 -> NONE / BIN3,4 -> J1!|MET / BIN5,6 -> J2!|MET 

    TH1D* matchesfinalecut = new TH1D((const TString)("matchesfinalecut"+std::to_string(j)),"",8,0,8); //  |            Leading Jet              |              Other Jet              |                       
                                                                                                         //  |       Aligned    |        Anti-A    |       Aligned    |        Anti-A    |
                                                                                                         //  | No Match | Match | No Match | Match | No Match | Match | No Match | Match |
                                                                                                         //  |  BIN  1  | BIN 2 |  BIN  3  | BIN 4 |  BIN  5  | BIN 6 |  BIN  7  | BIN 8 |

    TH1D* matchesfinaleqcut = new TH1D((const TString)("matcehsfinaleqcut"+std::to_string(j)),"",4,0,4); // BIN1 -> Both Leading Jet & Other Jet Are Matched To Seperate Dark Quarks | BIN2 -> Same Dark Quarks
                                                                                                         // BIN3 -> Only Leading Jet Matched | BIN4 -> Only Other Jet Matched

    for(int e = 0; e < chain->GetEntries(); e++){

      chain->GetEntry(e); // Gets the event

      int nbsm = ntruthBSM[0];       // Number of each type of object
      int njets = na4_pflowjets[0];

      TLorentzVector TEMP(0,0,0,0); // Temp lorentz vector for use during object loops


                                                             /////////////////////////////////////////////////////
      ///////////////////////////////////////////////////////////                MET Vector                  ///////////////////////////////////////////////////////////
                                                             /////////////////////////////////////////////////////

      TLorentzVector MET(0,0,0,0);
      MET.SetPtEtaPhiM(metFinalTrk,0,metFinalTrkPhi,0);
      MET.SetPz(0);
      MET.SetE(MET.Pt());

                                                             /////////////////////////////////////////////////////
      ///////////////////////////////////////////////////////////                  BSM Loop                   ///////////////////////////////////////////////////////////
                                                             /////////////////////////////////////////////////////
 
      vector <TLorentzVector> Quarks;
      vector <TLorentzVector> VisHadrons;
      vector <TLorentzVector> InvHadrons;
      vector <TLorentzVector> StabHadrons;
      vector <TLorentzVector> UnstHadrons;
      Quarks.clear();
      VisHadrons.clear();
      StabHadrons.clear();
      UnstHadrons.clear();

      for(int B = 0; B < nbsm; B++){ // BSM main loop
	
	TEMP.SetPtEtaPhiE(truthBSM_pt->at(B),truthBSM_eta->at(B),truthBSM_phi->at(B),truthBSM_e->at(B));

	int absid = std::abs(truthBSM_pdgId->at(B));

	if(absid == 5000001){;}                               // Z'
	if(absid == 4900101 && truthBSM_status->at(B) == 23){ // Dark Quark
	  Quarks.push_back(TEMP);
	} 
	if(absid == 4900111 || absid == 4900113){             // Unstable Dark Hadron
	  bool vis = false;
	  for(int C = 0; C < truthBSM_nChildren->at(B); C++){
	    int childID = std::abs(truthBSM_child_pdgId->at(B).at(C));
	    if(childID < 50 || (childID > 53 && childID < 4900000)){vis = true;} // Chekcs if any of the children are SM -> Visible Hadron
	  }
	  if(vis){VisHadrons.push_back(TEMP);}
	  else{InvHadrons.push_back(TEMP);}
	  UnstHadrons.push_back(TEMP);
	} 
	if(absid == 4900211 || absid == 4900213){             // Stable Dark Hadron
	  InvHadrons.push_back(TEMP);
	  StabHadrons.push_back(TEMP);
	} 

      }
 
      if(Quarks.size() == 2){
	if(Quarks.at(0).Pt() < Quarks.at(1).Pt()){ // Orders the dark quarks in decreasing pt order
	  Quarks.push_back(Quarks.at(0));
	  Quarks.erase(Quarks.begin());
	}
      }

      double maxdrvdh[2] = {0,0}; // Calculating the max DR between the visible dark hadrons from each dark quark
      
      if(Quarks.size() == 2){
	for(TLorentzVector v1 : VisHadrons){
	  if(Quarks.at(0).DeltaR(v1) < Quarks.at(1).DeltaR(v1)){
	    for(TLorentzVector v2 : VisHadrons){
	      if(Quarks.at(0).DeltaR(v2) < Quarks.at(1).DeltaR(v2) && v1.DeltaR(v2) > maxdrvdh[0]){
		maxdrvdh[0] = v1.DeltaR(v2);
	      }
	    }	  
	  }
	  else{
	    for(TLorentzVector v2 : VisHadrons){
	      if(Quarks.at(0).DeltaR(v2) > Quarks.at(1).DeltaR(v2) && v1.DeltaR(v2) > maxdrvdh[1]){
		maxdrvdh[1] = v1.DeltaR(v2);
	      }
	    }	  
	  }
	}
      }
                                                             /////////////////////////////////////////////////////
      ///////////////////////////////////////////////////////////             Small R Jet Loop                ///////////////////////////////////////////////////////////
                                                             /////////////////////////////////////////////////////
      vector <TLorentzVector> Jets;
      vector <TLorentzVector> LeadQuarkMatches;
      vector <TLorentzVector> SubQuarkMatches;

      Jets.clear();
      LeadQuarkMatches.clear();
      SubQuarkMatches.clear();

      TLorentzVector MinPhiV(0,0,0,0);
      TLorentzVector MaxPhiV(0,0,0,0);
	
      int matchjetpt[2] = {-1,-1};
      int matchjetaa[2] = {-1,-1};

      double mindphi = TMath::Pi();
      double maxdphi = 0;
      double hT = 0;
      for(int J = 0; J < njets; J++){
        
	TEMP.SetPtEtaPhiE(a4_pflowjets_pt->at(J),a4_pflowjets_eta->at(J),a4_pflowjets_phi->at(J),a4_pflowjets_E->at(J));
	
	if(TEMP.Pt() > ptc && std::abs(TEMP.Phi()) < TMath::Pi() && std::abs(TEMP.Eta()) < 4.5){
	  if(std::abs(TEMP.DeltaPhi(MET)) < mindphi){mindphi = std::abs(TEMP.DeltaPhi(MET)); MinPhiV = TEMP; matchjetaa[0] = -1;}
	  if(std::abs(TEMP.DeltaPhi(MET)) > maxdphi){maxdphi = std::abs(TEMP.DeltaPhi(MET)); MaxPhiV = TEMP; matchjetaa[1] = -1;}
	  
	  if(Quarks.size() == 2){
	    if(Quarks.at(0).DeltaR(TEMP) < matchcut && Quarks.at(0).DeltaR(TEMP) < Quarks.at(1).DeltaR(TEMP)){ // Checks if jet is matching 
	      LeadQuarkMatches.push_back(TEMP);
	      if(J == 0){matchjetpt[0] = 1;}
	      else if(J == 1){matchjetpt[1] = 1;}
	      if(TEMP == MinPhiV){matchjetaa[0] = 1;}
	      if(TEMP == MaxPhiV){matchjetaa[1] = 1;}
	    }
	    else if(Quarks.at(1).DeltaR(TEMP) < matchcut){
	      SubQuarkMatches.push_back(TEMP);
	      if(J == 0){matchjetpt[0] = 2;}
	      else if(J == 1){matchjetpt[1] = 2;}
	      if(TEMP == MinPhiV){matchjetaa[0] = 2;}
	      if(TEMP == MaxPhiV){matchjetaa[1] = 2;}
	    }
	    
	  }
	  hT += TEMP.Pt();
	  Jets.push_back(TEMP);
	}

      }

      double mt[3] = {0,0,0};

      if(Jets.size() >= 2 && Quarks.size() == 2){

	TEMP.SetPtEtaPhiE(0,0,0,0); // Calculates MT JJ+MET (A + AA)
	TEMP += MaxPhiV;
	TEMP += MinPhiV;
	//TEMP += MET;
	//mt[0] = TEMP.Mt();
	mt[0] = std::sqrt(std::pow(std::sqrt(std::pow(TEMP.M(),2)+std::pow(MaxPhiV.Pt()+MinPhiV.Pt(),2))+MET.Pt(),2) - std::pow(MaxPhiV.Pt()+MinPhiV.Pt()+MET.Pt(),2));

        TEMP.SetPtEtaPhiE(0,0,0,0); // Calculates MT JJ+MET (A + AA)
	TEMP += Jets.at(0);
	TEMP += Jets.at(1);
	//TEMP += MET;
	//mt[1] = TEMP.Mt();
	mt[1] = std::sqrt(std::pow(std::sqrt(std::pow(TEMP.M(),2)+std::pow(Jets.at(0).Pt()+Jets.at(1).Pt(),2))+MET.Pt(),2) - std::pow(Jets.at(0).Pt()+Jets.at(1).Pt()+MET.Pt(),2));

      }
    
    /////////////////////////////////////////////////////////////////
    ////                                                         ////
    ////                         Cut Tests                       ////
    ////                                                         ////
    /////////////////////////////////////////////////////////////////

      int TEST = -1; // If this stays as -1 then the cuts are passed


      //if(!() && TEST == -1){TEST = 1;}    <-- Template cut for copy pasting

      if(!(Quarks.size() == 2) && TEST == -1){TEST = 1;} // New test value can be used to produce a cutflow histogram
      if(!(Jets.size() >= 2) && TEST == -1){TEST = 1;}
      if(!(metFinalTrk > 200) && TEST == -1){TEST = 1;}

      //if(!(Jets.size() <= 4) && TEST == -1){TEST = 1;}

      if(TEST == -1){
	//if(!(Jets.at(0).Pt() > 150) && TEST == -1){TEST = 1;}
	//if(!(Jets.at(1).Pt() > 33) && TEST == -1){TEST = 1;}

	//if(!(Jets.at(0).Pt() < 66) && TEST == -1){TEST = 1;}
	//if(!(Jets.at(1).Pt() < 66) && TEST == -1){TEST = 1;}

	//if(!(MinPhiV != Jets.at(0)) && TEST == -1){TEST = 1;}
	//if(!(MinPhiV == Jets.at(1)) && TEST == -1){TEST = 1;}

	//if(!(MaxPhiV == Jets.at(0)) && TEST == -1){TEST = 1;}
	//if(!(MaxPhiV != Jets.at(1)) && TEST == -1){TEST = 1;}

	//if(!(hT > 200) && TEST == -1){TEST = 1;}
      }

      if(TEST == -1){    //Succesfully passed cuts
	passes++;
	if(passes%17 == 0 && nout < 5){
	  nout++;
	  cout << "----->> Event " << e << " Passed" << endl;
	}

    /////////////////////////////////////////////////////////////////
    ////                                                         ////
    ////                        Fill Numbers                     ////
    ////                                                         ////
    /////////////////////////////////////////////////////////////////

	if(matchjetpt[0] != -1 && matchjetpt[1] == -1){matchescut->AddBinContent(2);}               // Match Eff for Leading / Subleading Jets
	else if(matchjetpt[0] == -1 && matchjetpt[1] != -1){matchescut->AddBinContent(3);}
	else if(matchjetpt[0] != -1 && matchjetpt[1] != -1){matchescut->AddBinContent(4);}
	else{matchescut->AddBinContent(1);}

	if(matchjetaa[0] != -1 && matchjetaa[1] == -1){matchesalignedcut->AddBinContent(2);}        // Match Eff for Most Aligned / Most Anti-Aligned
	else if(matchjetaa[0] == -1 && matchjetaa[1] != -1){matchesalignedcut->AddBinContent(3);}
	else if(matchjetaa[0] != -1 && matchjetaa[1] != -1){matchesalignedcut->AddBinContent(4);}
	else{matchesalignedcut->AddBinContent(1);}
	if(MaxPhiV == Jets.at(0)){matchesalignedcut->AddBinContent(5);}

	if(Jets.at(0) == MinPhiV){                                      // Matching of Leading/Subleading when Most aligned
	  if(matchjetpt[0] == -1){matchesmetcut->AddBinContent(3);}
	  else{matchesmetcut->AddBinContent(4);}
	}
	else if(Jets.at(1) == MinPhiV){
	  if(matchjetpt[1] == -1){matchesmetcut->AddBinContent(5);}
	  else{matchesmetcut->AddBinContent(6);}
	}
	else{
	  if(Quarks.at(0).DeltaR(MinPhiV) < 0.4 || Quarks.at(1).DeltaR(MinPhiV) < 0.4){matchesmetcut->AddBinContent(1);}
	  else{matchesmetcut->AddBinContent(2);}
	}
	
	if(Jets.at(0) == MaxPhiV){                                      // Matching of Leading/Subleading when Most aligned
	  if(matchjetpt[0] == -1){matchesmetaacut->AddBinContent(3);}
	  else{matchesmetaacut->AddBinContent(4);}
	}
	else if(Jets.at(1) == MaxPhiV){
	  if(matchjetpt[1] == -1){matchesmetaacut->AddBinContent(5);}
	  else{matchesmetaacut->AddBinContent(6);}
	}
	else{
	  if(Quarks.at(0).DeltaR(MaxPhiV) < 0.4 || Quarks.at(1).DeltaR(MaxPhiV) < 0.4){matchesmetaacut->AddBinContent(1);}
	  else{matchesmetaacut->AddBinContent(2);}
	}

	bool matchtracker[2] = {false,false}; 

	if(
	   Jets.at(0) == MaxPhiV
	   //std::abs(Jets.at(0).DeltaPhi(MET)) >= (double)TMath::Pi()/(double)2 && Jets.at(0) != MinPhiV
	   ){  // Matching and alignment of leading jet and other jet
	  if(Jets.at(0).DeltaR(Quarks.at(0)) < matchcut || Jets.at(0).DeltaR(Quarks.at(1)) < matchcut){matchesfinalecut->AddBinContent(4); matchtracker[0] = true;}
	  else{matchesfinalecut->AddBinContent(3);}
	  
	  if(MinPhiV.DeltaR(Quarks.at(0)) < matchcut || MinPhiV.DeltaR(Quarks.at(1)) < matchcut){matchesfinalecut->AddBinContent(6); matchtracker[1] = true;}
	  else{matchesfinalecut->AddBinContent(5);}

	  if(matchtracker[0] && matchtracker[1]){ // checks if matching to same dark quark if both matched
	    if(std::count(LeadQuarkMatches.begin(),LeadQuarkMatches.end(),Jets.at(0)) == std::count(LeadQuarkMatches.begin(),LeadQuarkMatches.end(),MinPhiV)){matchesfinaleqcut->AddBinContent(2);}
	    else{matchesfinaleqcut->AddBinContent(1);}
	  }
	  else if(matchtracker[0]){matchesfinaleqcut->AddBinContent(3);}
	  else if(matchtracker[1]){matchesfinaleqcut->AddBinContent(4);}
	}
	else if(
		Jets.at(0) == MinPhiV
		//std::abs(Jets.at(0).DeltaPhi(MET)) < (double)TMath::Pi()/(double)2 && Jets.at(0) != MaxPhiV
		){
	  if(Jets.at(0).DeltaR(Quarks.at(0)) < matchcut || Jets.at(0).DeltaR(Quarks.at(1)) < matchcut){matchesfinalecut->AddBinContent(2); matchtracker[0] = true;}
	  else{matchesfinalecut->AddBinContent(1);}
	  
	  if(MaxPhiV.DeltaR(Quarks.at(0)) < matchcut || MaxPhiV.DeltaR(Quarks.at(1)) < matchcut){matchesfinalecut->AddBinContent(8); matchtracker[1] = true;}
	  else{matchesfinalecut->AddBinContent(7);}
	
	  if(matchtracker[0] && matchtracker[1]){ // checks if matching to same dark quark if both matched
	    if(std::count(LeadQuarkMatches.begin(),LeadQuarkMatches.end(),Jets.at(0)) == std::count(LeadQuarkMatches.begin(),LeadQuarkMatches.end(),MaxPhiV)){matchesfinaleqcut->AddBinContent(2);}
	    else{matchesfinaleqcut->AddBinContent(1);}
	  }
	  else if(matchtracker[0]){matchesfinaleqcut->AddBinContent(3);}
	  else if(matchtracker[1]){matchesfinaleqcut->AddBinContent(4);}
	}

    /////////////////////////////////////////////////////////////////
    ////                                                         ////
    ////                   Fill varying DRCut                    ////
    ////                                                         ////
    /////////////////////////////////////////////////////////////////
  
	if(true){                        // This is for filling the varying drcut hists but increases runtime by a large amount so turn off in not using it
	  for(int B = 0; B < 2; B++){
	    int N1 = JETMATCHDR[j][B]->GetXaxis()->GetNbins();
	    double binA1[N1+1];
	    for(int BIN = 0; BIN < N1; BIN++){binA1[BIN] = JETMATCHDR[j][B]->GetXaxis()->GetBinLowEdge(BIN+1);}
	    binA1[N1] = 2;
	    double wbinA1 = (double)std::abs(binA1[N1] - binA1[N1-1])/(double)2;

	    for(int bn = 0; bn < N1; bn++){
	      int bincount = 0;
	      int bincountP[3] = {0,0,0};
	      for(TLorentzVector J: Jets){
		if(J.DeltaR(Quarks.at(B)) < binA1[bn+1]){
		  bincount++;
		  if(J.Pt() > ptcuts[0]){bincountP[0]++;}
		  if(J.Pt() > ptcuts[1]){bincountP[1]++;}
		  if(J.Pt() > ptcuts[2]){bincountP[2]++;}
		}
	      }
	      JETMATCHDR[j][B]->Fill(binA1[bn],bincount);
	      if(actualInteractionsPerCrossing < nintcuts[0]){JETMATCHINT[j][B][0]->Fill(binA1[bn],bincount);}
	      else if(actualInteractionsPerCrossing > nintcuts[0] && actualInteractionsPerCrossing < nintcuts[1]){JETMATCHINT[j][B][1]->Fill(binA1[bn],bincount);}
	      else if(actualInteractionsPerCrossing > nintcuts[1]){JETMATCHINT[j][B][2]->Fill(binA1[bn],bincount);}
	      JETMATCHPT[j][B][0]->Fill(binA1[bn],bincountP[0]);
	      JETMATCHPT[j][B][1]->Fill(binA1[bn],bincountP[1]);
	      JETMATCHPT[j][B][2]->Fill(binA1[bn],bincountP[2]);
	    }
	  }
	}

     
    /////////////////////////////////////////////////////////////////
    ////                                                         ////
    ////                      Fill Histograms                    ////
    ////                                                         ////
    /////////////////////////////////////////////////////////////////
 
        JETPT[j][0]->Fill(Jets.at(0).Pt());
	JETPT[j][1]->Fill(Jets.at(1).Pt());
	JETPT2D[j]->Fill(Jets.at(0).Pt(),Jets.at(1).Pt());
	DQPT[j][0]->Fill(Quarks.at(0).Pt());
	DQPT[j][1]->Fill(Quarks.at(1).Pt());
	DQPT2D[j]->Fill(Quarks.at(0).Pt(),Quarks.at(1).Pt());
	DRJJ[j]->Fill(Jets.at(0).DeltaR(Jets.at(1)));
	DRQQ[j]->Fill(Quarks.at(0).DeltaR(Quarks.at(1)));
	//	JETMATCHDR[j]->Fill();
        if(maxdrvdh[0] > 0){MAXDR[j][0]->Fill(maxdrvdh[0]);}
	if(maxdrvdh[1] > 0){MAXDR[j][1]->Fill(maxdrvdh[1]);}
	NVIS[j]->Fill(VisHadrons.size());
	NINV[j]->Fill(InvHadrons.size());
	NSTAB[j]->Fill(StabHadrons.size());
	NUNST[j]->Fill(UnstHadrons.size());
	VISUNSTRAT[j]->Fill((double)VisHadrons.size()/(double)UnstHadrons.size());
	UNSTINVRAT[j]->Fill((double)UnstHadrons.size()/(double)InvHadrons.size());
	STABTOTRAT[j]->Fill((double)StabHadrons.size()/(double)(StabHadrons.size()+UnstHadrons.size()));
	INVTOTRAT[j]->Fill((double)InvHadrons.size()/(double)(StabHadrons.size()+UnstHadrons.size()));

	for(TLorentzVector v : VisHadrons){
	  if(Quarks.at(0).DeltaR(v) < Quarks.at(1).DeltaR(v)){VDHPT[j][0]->Fill(v.Pt());}
	  else{VDHPT[j][1]->Fill(v.Pt());}
	}
	for(TLorentzVector m : LeadQuarkMatches){
	  DPHIMATCHMET[j][0]->Fill(std::abs(m.DeltaPhi(MET)));
	}
	for(TLorentzVector m : SubQuarkMatches){
	  DPHIMATCHMET[j][1]->Fill(std::abs(m.DeltaPhi(MET)));
	}
	if(LeadQuarkMatches.size() == 1 && SubQuarkMatches.size() == 1){
	  DPHIMATCHMET2D[j]->Fill(std::abs(LeadQuarkMatches.at(0).DeltaPhi(MET)), std::abs(SubQuarkMatches.at(0).DeltaPhi(MET)));
	}
	DPHIJETMET[j][0]->Fill(std::abs(Jets.at(0).DeltaPhi(MET)));
	DPHIJETMET[j][1]->Fill(std::abs(Jets.at(1).DeltaPhi(MET)));
	MINDPHIMET[j]->Fill(mindphi);
	MAXDPHIMET[j]->Fill(maxdphi);

	MTJETSMET[j][0]->Fill(mt[0]);
	MTJETSMET[j][1]->Fill(mt[1]);
	SUMOFPT[j]->Fill(hT);

	DPHIJJ[j]->Fill(std::abs(Jets.at(0).DeltaPhi(Jets.at(1))));
	DPHIQQ[j]->Fill(std::abs(Quarks.at(0).DeltaPhi(Quarks.at(1))));

      }
    }
    /////////////////////////////////////////////////////////////////
    ////                                                         ////
    ////                     Write Histograms                    ////
    ////                                                         ////
    /////////////////////////////////////////////////////////////////


    TString opname="Z"+filevarM[M]+"r"+filevarr[r]+".root";
    TFile *f=(TFile*)gROOT->GetListOfFiles()->FindObject(opname);
    if(!f){
      f = new TFile(opname, "RECREATE");
    }    


    JETPT[j][0]->Write();
    JETPT[j][1]->Write();
    JETPT2D[j]->Write();
    DQPT[j][0]->Write();
    DQPT[j][1]->Write();
    DQPT2D[j]->Write();
    DRJJ[j]->Write();
    DRQQ[j]->Write();

    MAXDR[j][0]->Write();
    MAXDR[j][1]->Write();
    NVIS[j]->Write();
    NINV[j]->Write();
    NSTAB[j]->Write();
    NUNST[j]->Write();
    VISUNSTRAT[j]->Write();
    UNSTINVRAT[j]->Write();
    STABTOTRAT[j]->Write();
    INVTOTRAT[j]->Write();

    VDHPT[j][0]->Write();
    VDHPT[j][1]->Write();
    DPHIMATCHMET[j][0]->Write();
    DPHIMATCHMET[j][1]->Write();
    DPHIMATCHMET2D[j]->Write();
    DPHIJETMET[j][0]->Write();
    DPHIJETMET[j][1]->Write();
    MINDPHIMET[j]->Write();
    MAXDPHIMET[j]->Write();

    MTJETSMET[j][0]->Write();
    MTJETSMET[j][1]->Write();
    SUMOFPT[j]->Write();

    DPHIJJ[j]->Write();
    DPHIQQ[j]->Write();

    for(int B = 0; B < 2; B++){
      JETMATCHDR[j][B]->Write();
      for(int b = 0; b < 3; b++){
	JETMATCHINT[j][B][b]->Write();
	JETMATCHPT[j][B][b]->Write();
      }
    }

    /////////////////////////////////////////////////////////////////
    ////                                                         ////
    ////                    Output Numbers                       ////
    ////                                                         ////
    /////////////////////////////////////////////////////////////////

    double efficiencies[3][5];


    cout << endl;
    cout << endl;
    cout << "NO. Events: " << chain->GetEntries() << endl;
    cout << "NO. Passes: " << passes << endl;
    cout << endl;
    cout << endl;

    float normval = (float)100/(float)matchescut->GetSumOfWeights();
    normval = 1;
    cout << "Overall Jet Matching After Applied Cuts:" << endl;
    cout << endl;
    cout << std::setw(20) << std::right << "Matched Jets |";
    cout << std::setw(20) << std::right << "No Jets";
    cout << std::setw(20) << std::right << "Leading Jet";
    cout << std::setw(20) << std::right << "Subleading Jet";
    cout << std::setw(20) << std::right << "Both Jets";
    cout << endl;
    cout << std::setw(20) << std::right << "NO. Matches  |";
    cout << std::setw(20) << std::right << matchescut->GetBinContent(1)*normval;
    cout << std::setw(20) << std::right << matchescut->GetBinContent(2)*normval;
    cout << std::setw(20) << std::right << matchescut->GetBinContent(3)*normval;
    cout << std::setw(20) << std::right << matchescut->GetBinContent(4)*normval;
    cout << endl;
    cout << endl;

    normval = (float)100/(float)matchescut->GetSumOfWeights();
    efficiencies[0][0] = matchescut->GetBinContent(2)*normval + matchescut->GetBinContent(4)*normval;
    efficiencies[0][1] = matchescut->GetBinContent(3)*normval + matchescut->GetBinContent(4)*normval;
    efficiencies[0][2] = matchescut->GetBinContent(4)*normval;
    efficiencies[0][3] = 100 - matchescut->GetBinContent(4)*normval - matchescut->GetBinContent(1)*normval;
    efficiencies[0][4] = matchescut->GetBinContent(1)*normval;
    

    normval = (float)100/(float)(matchesalignedcut->GetBinContent(1) + matchesalignedcut->GetBinContent(2) + matchesalignedcut->GetBinContent(3) + matchesalignedcut->GetBinContent(4));
    normval = 1;
    cout << "Aligned/Anti-Aligned Jet Matching After Applied Cuts:" << endl;
    cout << endl;
    cout << std::setw(20) << std::right << "Matched Jets |";
    cout << std::setw(20) << std::right << "No Jets";
    cout << std::setw(20) << std::right << "Aligned Jet";
    cout << std::setw(20) << std::right << "Anti-Aligned Jet";
    cout << std::setw(20) << std::right << "Both Jets";
    cout << endl;
    cout << std::setw(20) << std::right << "NO. Matches  |";
    cout << std::setw(20) << std::right << matchesalignedcut->GetBinContent(1)*normval;
    cout << std::setw(20) << std::right << matchesalignedcut->GetBinContent(2)*normval;
    cout << std::setw(20) << std::right << matchesalignedcut->GetBinContent(3)*normval;
    cout << std::setw(20) << std::right << matchesalignedcut->GetBinContent(4)*normval;
    cout << endl;
    cout << endl;
    cout << "NO. Events where the Anti-Aligned Jet Is The Leading Jet:   " << matchesalignedcut->GetBinContent(5) << endl;
    cout << endl;
    cout << endl;

    normval = (float)100/(float)(matchesalignedcut->GetBinContent(1) + matchesalignedcut->GetBinContent(2) + matchesalignedcut->GetBinContent(3) + matchesalignedcut->GetBinContent(4));
    efficiencies[1][0] = matchesalignedcut->GetBinContent(2)*normval + matchesalignedcut->GetBinContent(4)*normval;
    efficiencies[1][1] = matchesalignedcut->GetBinContent(3)*normval + matchesalignedcut->GetBinContent(4)*normval;
    efficiencies[1][2] = matchesalignedcut->GetBinContent(4)*normval;
    efficiencies[1][3] = 100 - matchesalignedcut->GetBinContent(4)*normval - matchesalignedcut->GetBinContent(1)*normval;
    efficiencies[1][4] = matchesalignedcut->GetBinContent(1)*normval;

    normval = (float)100/(float)matchesmetcut->GetSumOfWeights();
    cout << "Overall Jet||MET After Applied Cuts:" << endl;
    cout << endl;
    cout << std::setw(20) << std::right << "Jet||MET |";
    cout << std::setw(30) << std::right << "Other Jet";
    cout << std::setw(30) << std::right << "Leading Jet";
    cout << std::setw(30) << std::right << "Subleading Jet";
    cout << endl;
    cout << std::setw(20) << std::right << "";
    cout << std::setw(15) << std::right << "No Match";
    cout << std::setw(15) << std::right << "Match";
    cout << std::setw(15) << std::right << "No Match";
    cout << std::setw(15) << std::right << "Match";
    cout << std::setw(15) << std::right << "No Match";
    cout << std::setw(15) << std::right << "Match";
    cout << endl;
    cout << std::setw(20) << std::right << "NO. ||MET  |";
    for(int b = 0; b < matchesmetcut->GetNbinsX(); b++){
      cout << std::setw(15) << std::right << matchesmetcut->GetBinContent(b+1)*normval;
    }
    cout << endl;
    cout << endl;

    normval = (float)100/(float)matchesmetaacut->GetSumOfWeights();
    cout << "Overall Jet!|MET After Applied Cuts:" << endl;
    cout << endl;
    cout << std::setw(20) << std::right << "Jet||MET |";
    cout << std::setw(30) << std::right << "Other Jet";
    cout << std::setw(30) << std::right << "Leading Jet";
    cout << std::setw(30) << std::right << "Subleading Jet";
    cout << endl;
    cout << std::setw(20) << std::right << "";
    cout << std::setw(15) << std::right << "No Match";
    cout << std::setw(15) << std::right << "Match";
    cout << std::setw(15) << std::right << "No Match";
    cout << std::setw(15) << std::right << "Match";
    cout << std::setw(15) << std::right << "No Match";
    cout << std::setw(15) << std::right << "Match";
    cout << endl;
    cout << std::setw(20) << std::right << "NO. ||MET  |";
    for(int b = 0; b < matchesmetaacut->GetNbinsX(); b++){
      cout << std::setw(15) << std::right << matchesmetaacut->GetBinContent(b+1)*normval;
    }
    cout << endl;
    cout << endl;

    normval = (float)200/(float)matchesfinalecut->GetSumOfWeights();
    cout << "Leading Jet + Other Jet MET Alignment & Dark Quark Matching After Applied cuts:" << endl;
    cout << endl;
    cout << std::setw(20) << std::right << "Jet   |";
    cout << std::setw(60) << std::right << "Leading Jet            ";
    cout << std::setw(60) << std::right << "Other Jet             |";
    cout << endl;
    cout << std::setw(20) << std::right << "Alignment |";
    cout << std::setw(30) << std::right << "Aligned     ";
    cout << std::setw(30) << std::right << "Anti-Al     ";
    cout << std::setw(30) << std::right << "Aligned     ";
    cout << std::setw(30) << std::right << "Anti-Al    |";
    cout << endl;
    cout << std::setw(20) << std::right << "Matching |";
    cout << std::setw(15) << std::right << "No Match ";
    cout << std::setw(15) << std::right << "Match ";
    cout << std::setw(15) << std::right << "No Match ";
    cout << std::setw(15) << std::right << "Match ";
    cout << std::setw(15) << std::right << "No Match ";
    cout << std::setw(15) << std::right << "Match ";
    cout << std::setw(15) << std::right << "No Match ";
    cout << std::setw(15) << std::right << "Match |";
    cout << endl;
    cout << std::setw(20) << std::right << "Percentage |";
    for(int b = 0; b < matchesfinalecut->GetNbinsX(); b++){
      cout << std::setw(15) << std::right << matchesfinalecut->GetBinContent(b+1)*normval;
    }
    cout << endl;
    cout << endl;
    // normval = (double)100/(double)matchesfinaleqcut->GetSumOfWeights();
    cout << "NO. events for this table = " << (double)matchesfinalecut->GetSumOfWeights()/(double)2 << endl;
    cout << "Fraction of these events where both the leading jet & other jet are matched to seperate dark quarks = " << matchesfinaleqcut->GetBinContent(1)*normval << endl;
    cout << "Fraction of these events where both the leading jet & other jet are matched to the same dark quark = " << matchesfinaleqcut->GetBinContent(2)*normval << endl;
    cout << "Fraction of these events where only the leading jet is matched to a dark quark = " << matchesfinaleqcut->GetBinContent(3)*normval << endl;
    cout << "Fraction of these events where only the other jet is matched to a dark quark = " << matchesfinaleqcut->GetBinContent(4)*normval << endl;
    cout << endl;
    cout << endl;

    efficiencies[2][0] = matchesfinalecut->GetBinContent(2)*normval + matchesfinalecut->GetBinContent(4)*normval;
    efficiencies[2][1] = matchesfinalecut->GetBinContent(6)*normval + matchesfinalecut->GetBinContent(8)*normval;
    efficiencies[2][2] = matchesfinaleqcut->GetBinContent(1)*normval + matchesfinaleqcut->GetBinContent(2)*normval;
    efficiencies[2][3] = matchesfinaleqcut->GetBinContent(3)*normval + matchesfinaleqcut->GetBinContent(4)*normval;
    efficiencies[2][4] = 100 - efficiencies[2][2] - efficiencies[2][3]; 

    cout << "Pt Eff Table:" << endl;
    cout << std::setw(20) << std::right << "Jet(s) |";
    cout << std::setw(20) << std::right << "Leading Jet Match";
    cout << std::setw(20) << std::right << "Subleading Jet Match";
    cout << std::setw(20) << std::right << "Both Match";
    cout << std::setw(20) << std::right << "one Match";
    cout << std::setw(20) << std::right << "No Match";
    cout << endl;   
    cout << std::setw(20) << std::right << "Efficiency (%) |";
    cout << std::setw(20) << std::right << efficiencies[0][0];
    cout << std::setw(20) << std::right << efficiencies[0][1];
    cout << std::setw(20) << std::right << efficiencies[0][2];
    cout << std::setw(20) << std::right << efficiencies[0][3];
    cout << std::setw(20) << std::right << efficiencies[0][4];
    cout << endl;

    cout << "Alignment Eff Table:" << endl;
    cout << std::setw(20) << std::right << "Jet(s) |";
    cout << std::setw(25) << std::right << "Most Aligned Jet Match";
    cout << std::setw(25) << std::right << "Most Anti-Aligned Jet Match";
    cout << std::setw(20) << std::right << "Both Match";
    cout << std::setw(20) << std::right << "one Match";
    cout << std::setw(20) << std::right << "No Match";
    cout << endl;   
    cout << std::setw(20) << std::right << "Efficiency (%) |";
    cout << std::setw(25) << std::right << efficiencies[1][0];
    cout << std::setw(25) << std::right << efficiencies[1][1];
    cout << std::setw(20) << std::right << efficiencies[1][2];
    cout << std::setw(20) << std::right << efficiencies[1][3];
    cout << std::setw(20) << std::right << efficiencies[1][4];
    cout << endl;

    cout << "Leading + Other Eff Table:" << endl;
    cout << std::setw(20) << std::right << "Jet(s) |";
    cout << std::setw(20) << std::right << "Leading Jet Match";
    cout << std::setw(20) << std::right << "Other Jet Match";
    cout << std::setw(20) << std::right << "Both Match";
    cout << std::setw(20) << std::right << "one Match";
    cout << std::setw(20) << std::right << "No Match";
    cout << endl;   
    cout << std::setw(20) << std::right << "Efficiency (%) |";
    cout << std::setw(20) << std::right << efficiencies[2][0];
    cout << std::setw(20) << std::right << efficiencies[2][1];
    cout << std::setw(20) << std::right << efficiencies[2][2];
    cout << std::setw(20) << std::right << efficiencies[2][3];
    cout << std::setw(20) << std::right << efficiencies[2][4];
    cout << endl;
    cout << endl;

  

}
  

}

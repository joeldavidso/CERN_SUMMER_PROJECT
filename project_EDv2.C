
#include "ATLASSTYLE/AtlasStyle.C"

#include "ATLASSTYLE/AtlasUtils.C"
#include "ATLASSTYLE/AtlasLabels.C"

#include <sstream>
#include <string>
#include <iostream>

template <typename Type> std::string to_str(const Type & t) // used to turn floats to strings without the pesky .000000000 showing
{
  std::ostringstream os;
  os << t;
  return os.str ();
}

const TString datadirmain = "/eos/atlas/atlascerngroupdisk/phys-exotics/jdm/svjets-schannel/"; // Strings for locating the nTuples (Change this if using different files)
const TString datadirsub = "v1.0/";
const TString filename1 = "user.ebusch.5085";
const TString filename2 = ".MGPy8EG_SVJSChan_";
const TString filename3 = "_tree.root";
const TString mcs[3] = {".mc16a", ".mc16d", ".mc16e"};

const TString tree_entry = "outTree/nominal";


const int ptc = 20; // Pt Minimum for small R jets [GeV]
const int lptc = 200; // Pt Minimum for large R jets [GeV]



/////////////////////////////////////////////////////////////////                  //////
////                                                         ////              ///////
////          Input Variables For Event Displays             ////           //////////////////////////////////////////////////////
////                                                         ////              ///////
/////////////////////////////////////////////////////////////////                  //////


// Only need to change these four constants for running the code!!!

const int nevents = 3; // The Number of events per signal point you want displayed & saved
const int nsignals = 1; // The Number of signal points you want 

const TString filevar[nsignals][2] = {{"1500","8"}}; // The Signal points you want to look at {MZ',r_inv}

const int events[nevents][nsignals] = {{108},{185},{295}}; // The Event number(s) you want displayed and saved for each signal point {event# for signal point 1, ... , event# for signal point nsignals}



//Numbers specific to the signal location
//The Three Commented out elements correspond to 1500|3, 750|8, 750|3 (Mz'|r_inv)

const TString fileN[nsignals] = {"47"
				 // , "48", "49", "50"
};                                                                        
const TString ndN[nsignals][3] = {{"3733", "3907", "3996"}, 
				  //{"3741", "3916", "4007"}, {"3753", "3924", "4015"}, {"3764", "3929", "4025"}
};


void project_EDv2(){ // start program

  gInterpreter->GenerateDictionary("vector<vector<int> >", "vector"); // Loads Vector<Vector<int>> dictionary for child pdgId's

  SetAtlasStyle(); // Sets Style For ATLAS

  TCanvas* C[nsignals][nevents];
  TPad* P[nsignals][nevents];
  TPad* P1[nsignals][nevents];
  TPad* P2[nsignals][nevents];
  TLegend* L[nsignals][nevents];
   
  /////////////////////////////////////////////////////////////////
  ////                                                         ////
  ////        Run Over The Selected Signal Points &            ////
  ////                    Access nTuples                       ////
  /////////////////////////////////////////////////////////////////

  for (int j = 0; j<nsignals; j++){ // Loops over the selected signal points

    TChain *chain = new TChain(tree_entry);

    cout << "Z" << filevar[j][0] << "r" << filevar[j][1] << ":" << endl;
    for (int i = 0;i<3;i++){
      TString ndfilenam = "user.ebusch.2909"+ndN[j][i]+"._000001.tree.root";
      TString filenam = datadirmain+datadirsub+filename1+fileN[j]+filename2+filevar[j][0]+"_"+filevar[j][1]+mcs[i]+filename3;
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

    vector<int> *truthBSM_nChildren = nullptr;
    chain->SetBranchAddress("truthBSM_nChildren",&truthBSM_nChildren);

    vector<vector<int>> *truthBSM_child_pdgId = nullptr;
    chain->SetBranchAddress("truthBSM_child_pdgId",&truthBSM_child_pdgId);
    
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

    float metFinalTrkPhi;
    chain->SetBranchAddress("metFinalTrkPhi",&metFinalTrkPhi);

    int eventNumber[22];
    chain->SetBranchAddress("eventNumber",&eventNumber);

    /////////////////////////////////////////////////////////////////
    ////                                                         ////
    ////                 Selected Events Loop                    ////
    ////                                                         ////
    /////////////////////////////////////////////////////////////////

    for(int e = 0; e < nevents; e++){

      C[j][e] = new TCanvas((const TString)(std::to_string(e)+"_"+std::to_string(j)),(const TString)("E"+std::to_string(events[j][e])+"P"+std::to_string(e)),800,600); // Sets up the canvas

      P1[j][e] = new TPad((const TString)("P1"+std::to_string(e)+"_"+std::to_string(j)),(const TString)("P1"+std::to_string(e)),0,0,0.75,1); // The first main pad where the markers are drawn
      P1[j][e]->Draw();

      C[j][e]->cd(0);
      P1[j][e]->cd(0);
      auto frame = P1[j][e]->DrawFrame(-TMath::Pi(),-4.5,TMath::Pi(),4.5,";\\phi;\\eta"); // The frame in the first pad where the drawing happens


      chain->GetEntry(events[j][e]); // Gets Event

      TLorentzVector TEMP(0,0,0,0);
     

                                                             /////////////////////////////////////////////////////
      ///////////////////////////////////////////////////////////                Draw MET Line                ///////////////////////////////////////////////////////////
                                                             /////////////////////////////////////////////////////
    
      
      TLine *METL = new TLine(metFinalTrkPhi,-4.5,metFinalTrkPhi,4.5);
      METL->SetLineColor(kRed);
      METL->SetLineStyle(2);
      METL->DrawClone();
    

                                                             /////////////////////////////////////////////////////
      ///////////////////////////////////////////////////////////            Run over Large R Jets            ///////////////////////////////////////////////////////////
                                                             /////////////////////////////////////////////////////

      vector <TLorentzVector> LJets;  // LorentzVectors of (Large R) jets (Used for jet cone drawing later)
      int LLeadPos = 0;               // Leading jet's position in LJets
      int LSubPos = 0;                // Subleading jet's position in LJets
      int LLeadPt = 0;                // Leading jet's Pt
      int LSubPt = 0;                 // Subleading jet's Pt 
    
      for(int J = 0; J < na10_lctopojets[0]; J++){ // Loop over all large R jets to add to the vectors
	
	TEMP.SetPtEtaPhiE(a10_lctopojets_pt->at(J),a10_lctopojets_eta->at(J),a10_lctopojets_phi->at(J),0);
	
	if(TEMP.Pt() > LLeadPt){LSubPos = LLeadPos; LSubPt = TEMP.Pt(); LLeadPos = J; LLeadPt = TEMP.Pt();}
	else if(TEMP.Pt() > LSubPt){LSubPos = J; LSubPt = TEMP.Pt();}
	if(TEMP.Pt() > lptc && std::abs(TEMP.Phi()) < TMath::Pi() && std::abs(TEMP.Eta()) < 4.5){LJets.push_back(TEMP);}
      }


      /////////////////////////////////////////////////////
      ////           Large R Jets Main Loop            ////
      /////////////////////////////////////////////////////

      
      float jetr = 1;

      // These are all the variables which determine how to draw the Jet radii in the display

      vector<float> theta;       // This is the opening angle with the |Phi| = Pi edges
      vector<float> beta;        // This is the opening angle with the |Eta| = 4.5 edges
      vector<float> rotate;      // This is the rotation angle for theta > 0
      vector<float> rotateb;     // This is the rotation angle for beta > 0
      vector<int> cross;         // This is the index of which |Phi| = Pi edge is hit
      vector<int> crossb;        // This is the index of which |Eta| = 4.5  edge is hit
      
      vector<float> cornerdist;  // This is the distance to the nearest corner from the jet centre position

      float jetphi = 0;
      float jeteta = 0;

      float thet, bet, rott, rotb, ph, et; // renamed variables so code takes up less space

      int LS = 2; int LC = 40; int FS = 0; int FC = 38; int MST = 26; double MSI = 1; // Numbers for drawing the jet radii in the correct style

      for(int J = 0; J < LJets.size(); J++){ // Runs over all large R jets

	/////////////////////////////////////////////////////
	////        Large R Jets Edge Detection          ////
	/////////////////////////////////////////////////////

	jetphi = LJets.at(J).Phi(); // Initialises the values
        jeteta = LJets.at(J).Eta();
        theta.push_back(0);
	beta.push_back(0);
	rotate.push_back(0);
	rotateb.push_back(0);
	cross.push_back(-1);
	crossb.push_back(-1);
	cornerdist.push_back(100);

	for(int s = 0;  s < 4; s++){ // Loop of all 4 edges
	  int A = s %2;
	  float root = 0;
	  float rootb = 0;
	  float cdist = std::sqrt((jeteta - (4.5*(2*(s/2)-1)))*(jeteta - (4.5*(2*(s/2)-1)))+(jetphi - TMath::Pi()*(1-2*(((s+1)/2)%2)))*(jetphi - TMath::Pi()*(1-2*(((s+1)/2)%2)))); // Distance to a corner
	  if(cdist < cornerdist.at(J)){cornerdist.at(J) = cdist;}  // Saves corner dist if smallest so far for that jet
	  
	  if(A == 0){ // Checks if jet radius interescts with |Phi| = Pi
	    root = jetr*jetr - (jetphi-(1-s)*TMath::Pi())*(jetphi-(1-s)*TMath::Pi()); 
	    if(root > 0){
	      cross.at(J) = s;
	      theta.at(J) = asin((double)((double)std::sqrt(root)/(double)jetr))* (double)((double)360/(double)(2*TMath::Pi()));
	      rotate.at(J) = theta.at(J);
	      if(s == 2){rotate.at(J) = 180+theta.at(J);}
	    }
	  }
	  if(A == 1){ // Checks if jet radius interescts with |Eta| = 4.5
	    rootb = jetr*jetr - (jeteta+(2-s)*4.5)*(jeteta+(2-s)*4.5);
	    if(rootb > 0){
	      crossb.at(J) = s;
	      beta.at(J) = asin((double)((double)std::sqrt(rootb)/(double)jetr))* (double)((double)360/(double)(2*TMath::Pi()));
	      rotateb.at(J) = 270 + beta.at(J);
	      if(s == 3){rotateb.at(J) = 90 + beta.at(J);}
	    }
	  }
	}

	/////////////////////////////////////////////////////
	////         Drawing Large R Jets  Radii         ////
	/////////////////////////////////////////////////////

      
	thet = theta.at(J); bet = beta.at(J); 
	rott = rotate.at(J); rotb = rotateb.at(J);
	ph = LJets.at(J).Phi(); et = LJets.at(J).Eta();
	
	TEllipse* JARC;
	JARC = new TEllipse(ph,et,jetr); // creates arc as if no edge crossing
	TGraph* FRAMEJ = new TGraph();
	FRAMEJ->AddPoint(ph,et);
	FRAMEJ->GetXaxis()->SetLimits(-TMath::Pi(),TMath::Pi());
	FRAMEJ->GetYaxis()->SetLimits(-4.5,4.5);

	
	if(cross.at(J) != -1 || crossb.at(J) != -1){ // Checks if it crosses any edge

	  int narc = 0; // The number of extra arks to draw (e.g it loops over the |Phi| = Pi edge)
	  TEllipse* CJARC[2];
	  
	  if(cross.at(J) != -1 && crossb.at(J) == -1){ // Crosses right or left only
	    JARC = new TEllipse(ph,et,jetr,jetr,0,360-2*thet,rott);
	    CJARC[0] = new TEllipse(ph-(2-cross.at(J)*2)*TMath::Pi(),et,jetr,jetr,0,2*thet,360-rott); 
	    narc = 1;
	  }
	  
	  if(cross.at(J) == -1 && crossb.at(J) != -1){ // Crosses top or bottom only
	    JARC = new TEllipse(ph,et,jetr,jetr,0,360-2*bet,rotb);
	  }

	  if(cross.at(J) != -1 && crossb.at(J) != -1){ // Crosses a corner
	    float rotatey;
	    float rotateyr;
	    float rotatey2;
	    int cn = 10*cross.at(J)+crossb.at(J); // 1 -> BR| 20+1 -> BL| 23 -> TL| 3 -> TR (8 Different total corner cases)
	    if(cn == 1){
	      rotatey = thet; 
	      rotateyr = 360-thet;
	      rotatey2 = 270 + bet;
	      if(90-thet-bet <= 0){ // If corner is witin jet radius
		rotateyr = 270 + bet;
		rotatey2 = -1;
	      }
	    }
	    if(cn == (20+1)){
	      rotatey = 270+bet+1.5; 
	      rotateyr = 180-thet;
	      rotatey2 = 180+thet;
	      if(90-thet-bet <= 0){ // If corner is witin jet radius
		rotateyr = 180-thet-1;
		rotatey2 = -1;
	      }
	    }
	    if(cn == 23){
	      rotatey = 180+thet+1.5; 
	      rotateyr = 180-thet;
	      rotatey2 = 90 + bet;
	      if(90-thet-bet <= 0){ // If corner is witin jet radius
		rotateyr = 90+bet-1;
		rotatey2 = -1;
	      }
	    }
	    if(cn == 3){
	      rotatey = 90+bet; 
	      rotateyr = 360-thet;
	      rotatey2 = thet;
	      if(90-thet-bet <= 0){ // If corner is witin jet radius
		rotateyr = 360-thet;
		rotatey2 = -1;
	      }
	    }
	    JARC = new TEllipse(ph,et,jetr,jetr,0,270-bet-thet,rotatey); // Overrides original jet radius

	    if(rotatey2 == -1){ // Corner is within jet radius
	      CJARC[0] = new TEllipse(ph-(2-cross.at(J)*2)*TMath::Pi(),et,jetr,jetr,0,90+thet-bet,rotateyr); // Additional arc to draw
	      narc = 1;
	    }
	    else{ // Corner is not within jet radius
	      CJARC[0] = new TEllipse(ph-(2-cross.at(J)*2)*TMath::Pi(),et,jetr,jetr,0,90-bet+thet,rotateyr); // Additional arcs to draw
	      CJARC[1] = new TEllipse(ph,et,jetr,jetr,0,90-thet-bet,rotatey2);
	      narc = 2;
	    }
	  }
	  

	  for(int p = 0; p < narc; p++){ // Draws the additional jet radii
	    CJARC[p]->SetNoEdges();
	    CJARC[p]->SetLineColor(LC);
	    CJARC[p]->SetLineStyle(2);
	    CJARC[p]->SetFillStyle(FS);
	    CJARC[p]->SetFillColor(FC);
	    CJARC[p]->SetLineWidth(2);
	    CJARC[p]->DrawClone("SAME");
	  }

	} 

	JARC->SetNoEdges();     // Draws the main jet radius
	JARC->SetLineColor(LC);
	JARC->SetFillStyle(FS);
	JARC->SetFillColor(FC);
	JARC->SetLineWidth(2);
	JARC->SetLineStyle(2);
	JARC->DrawClone("SAME"); 
      }

                                                             /////////////////////////////////////////////////////
      ///////////////////////////////////////////////////////////            Run over Small R Jets            ///////////////////////////////////////////////////////////
                                                             /////////////////////////////////////////////////////

      vector <TLorentzVector> Jets;  // LorentzVectors of (Small R) jets (Used for jet cone drawing later)
      int LeadPos = 0;               // Leading jet's position in Jets
      int SubPos = 0;                // Subleading jet's position in Jets
      int LeadPt = 0;                // Leading jet's Pt
      int SubPt = 0;                 // Subleading jet's Pt 

      for(int J = 0; J < na4_pflowjets[0]; J++){ // Loop over all small R jets to add to the vectors
	
	TEMP.SetPtEtaPhiE(a4_pflowjets_pt->at(J),a4_pflowjets_eta->at(J),a4_pflowjets_phi->at(J),a4_pflowjets_E->at(J));
	
	if(TEMP.Pt() > LeadPt){SubPos = LeadPos; SubPt = TEMP.Pt(); LeadPos = J; LeadPt = TEMP.Pt();}
	else if(TEMP.Pt() > SubPt){SubPos = J; SubPt = TEMP.Pt();}
	if(TEMP.Pt() > ptc && std::abs(TEMP.Phi()) < TMath::Pi() && std::abs(TEMP.Eta()) < 4.5){Jets.push_back(TEMP);}
      }


      /////////////////////////////////////////////////////
      ////           Small R Jets Main Loop            ////
      /////////////////////////////////////////////////////

      
      jetr = 0.4;

      // These are all the variables which determine how to draw the Jet radii in the display (Reused from large R jets)

      theta.clear();       // This is the opening angle with the |Phi| = Pi edges
      beta.clear();        // This is the opening angle with the |Eta| = 4.5 edges
      rotate.clear();      // This is the rotation angle for theta > 0
      rotateb.clear();     // This is the rotation angle for beta > 0
      cross.clear();         // This is the index of which |Phi| = Pi edge is hit
      crossb.clear();        // This is the index of which |Eta| = 4.5  edge is hit

      cornerdist.clear();  // This is the distance to the nearest corner from the jet centre position

      for(int J = 0; J < Jets.size(); J++){ // Runs over all small R jets

      /////////////////////////////////////////////////////
      ////        Small R Jets Edge Detection          ////
      /////////////////////////////////////////////////////

	jetphi = Jets.at(J).Phi(); // Initialises the values
        jeteta = Jets.at(J).Eta();
        theta.push_back(0);
	beta.push_back(0);
	rotate.push_back(0);
	rotateb.push_back(0);
	cross.push_back(-1);
	crossb.push_back(-1);
	cornerdist.push_back(100);

	for(int s = 0;  s < 4; s++){ // Loop of all 4 edges
	  int A = s %2;
	  float root = 0;
	  float rootb = 0;
	  float cdist = std::sqrt((jeteta - (4.5*(2*(s/2)-1)))*(jeteta - (4.5*(2*(s/2)-1)))+(jetphi - TMath::Pi()*(1-2*(((s+1)/2)%2)))*(jetphi - TMath::Pi()*(1-2*(((s+1)/2)%2)))); // Distance to a corner
	  if(cdist < cornerdist.at(J)){cornerdist.at(J) = cdist;}  // Saves corner dist if smallest so far for that jet
	  
	  if(A == 0){ // Checks if jet radius interescts with |Phi| = Pi
	    root = jetr*jetr - (jetphi-(1-s)*TMath::Pi())*(jetphi-(1-s)*TMath::Pi()); 
	    if(root > 0){
	      cross.at(J) = s;
	      theta.at(J) = asin((double)((double)std::sqrt(root)/(double)jetr))* (double)((double)360/(double)(2*TMath::Pi()));
	      rotate.at(J) = theta.at(J);
	      if(s == 2){rotate.at(J) = 180+theta.at(J);}
	    }
	  }
	  if(A == 1){ // Checks if jet radius interescts with |Eta| = 4.5
	    rootb = jetr*jetr - (jeteta+(2-s)*4.5)*(jeteta+(2-s)*4.5);
	    if(rootb > 0){
	      crossb.at(J) = s;
	      beta.at(J) = asin((double)((double)std::sqrt(rootb)/(double)jetr))* (double)((double)360/(double)(2*TMath::Pi()));
	      rotateb.at(J) = 270 + beta.at(J);
	      if(s == 3){rotateb.at(J) = 90 + beta.at(J);}
	    }
	  }
	}

      /////////////////////////////////////////////////////
      ////    Drawing Small R Jets Markers & Radii     ////
      /////////////////////////////////////////////////////

      
	thet = theta.at(J); bet = beta.at(J); 
	rott = rotate.at(J); rotb = rotateb.at(J);
	ph = Jets.at(J).Phi(); et = Jets.at(J).Eta();
	
	TEllipse* JARC;
	JARC = new TEllipse(ph,et,jetr); // creates arc as if no edge crossing
	TGraph* FRAMEJ = new TGraph();
	FRAMEJ->AddPoint(ph,et);
	FRAMEJ->GetXaxis()->SetLimits(-TMath::Pi(),TMath::Pi());
	FRAMEJ->GetYaxis()->SetLimits(-4.5,4.5);

	
	if(cross.at(J) != -1 || crossb.at(J) != -1){ // Checks if it crosses any edge

	  int narc = 0; // The number of extra arks to draw (e.g it loops over the |Phi| = Pi edge)
	  TEllipse* CJARC[2];
	  
	  if(cross.at(J) != -1 && crossb.at(J) == -1){ // Crosses right or left only
	    JARC = new TEllipse(ph,et,jetr,jetr,0,360-2*thet,rott);
	    CJARC[0] = new TEllipse(ph-(2-cross.at(J)*2)*TMath::Pi(),et,jetr,jetr,0,2*thet,360-rott); 
	    narc = 1;
	  }
	  
	  if(cross.at(J) == -1 && crossb.at(J) != -1){ // Crosses top or bottom only
	    JARC = new TEllipse(ph,et,jetr,jetr,0,360-2*bet,rotb);
	  }

	  if(cross.at(J) != -1 && crossb.at(J) != -1){ // Crosses a corner
	    float rotatey;
	    float rotateyr;
	    float rotatey2;
	    int cn = 10*cross.at(J)+crossb.at(J); // 1 -> BR| 20+1 -> BL| 23 -> TL| 3 -> TR (8 Different total corner cases)
	    if(cn == 1){
	      rotatey = thet; 
	      rotateyr = 360-thet;
	      rotatey2 = 270 + bet;
	      if(90-thet-bet <= 0){ // If corner is witin jet radius
		rotateyr = 270 + bet;
		rotatey2 = -1;
	      }
	    }
	    if(cn == (20+1)){
	      rotatey = 270+bet+1.5; 
	      rotateyr = 180-thet;
	      rotatey2 = 180+thet;
	      if(90-thet-bet <= 0){ // If corner is witin jet radius
		rotateyr = 180-thet-1;
		rotatey2 = -1;
	      }
	    }
	    if(cn == 23){
	      rotatey = 180+thet+1.5; 
	      rotateyr = 180-thet;
	      rotatey2 = 90 + bet;
	      if(90-thet-bet <= 0){ // If corner is witin jet radius
		rotateyr = 90+bet-1;
		rotatey2 = -1;
	      }
	    }
	    if(cn == 3){
	      rotatey = 90+bet; 
	      rotateyr = 360-thet;
	      rotatey2 = thet;
	      if(90-thet-bet <= 0){ // If corner is witin jet radius
		rotateyr = 360-thet;
		rotatey2 = -1;
	      }
	    }
	    JARC = new TEllipse(ph,et,jetr,jetr,0,270-bet-thet,rotatey); // Overrides original jet radius

	    if(rotatey2 == -1){ // Corner is within jet radius
	      CJARC[0] = new TEllipse(ph-(2-cross.at(J)*2)*TMath::Pi(),et,jetr,jetr,0,90+thet-bet,rotateyr); // Additional arc to draw
	      narc = 1;
	    }
	    else{ // Corner is not within jet radius
	      CJARC[0] = new TEllipse(ph-(2-cross.at(J)*2)*TMath::Pi(),et,jetr,jetr,0,90-bet+thet,rotateyr); // Additional arcs to draw
	      CJARC[1] = new TEllipse(ph,et,jetr,jetr,0,90-thet-bet,rotatey2);
	      narc = 2;
	    }
	  }
	  

	  LC = 1; FS = 0; FC = 1;LS = 1; // Sets the draw style and options for the radii
	  MST = 23;
	  MSI = 1;
	  if(J == 0 || J == 1){LS = 2;LC = 1; FS = 0; FC = 38; MST = 32; MSI = 1;}

	  for(int p = 0; p < narc; p++){ // Draws the radii
	    CJARC[p]->SetNoEdges();
	    CJARC[p]->SetLineColor(LC);
	    CJARC[p]->SetLineStyle(LS);
	    CJARC[p]->SetFillStyle(FS);
	    CJARC[p]->SetFillColor(FC);
	    CJARC[p]->SetLineWidth(2);
	    CJARC[p]->DrawClone("SAME");
	  }

	}

	LC = 1; FS = 0; FC = 1; LS = 1; // Sets the draw style and options for the radii
	MST = 23;
	MSI = 1;
	if(J == 0 || J == 1){LS = 2;LC = 1; FS = 0; FC = 38; MST = 32; MSI = 1;}

	JARC->SetNoEdges();     // Draws the radii
	JARC->SetLineColor(LC);
	JARC->SetFillStyle(FS);
	JARC->SetFillColor(FC);
	JARC->SetLineWidth(2);
	JARC->SetLineStyle(LS);
	JARC->DrawClone("SAME");
	
	
	FRAMEJ->SetMarkerStyle(MST);    // Sets the draw style and options for the jet marker & draws it
	FRAMEJ->SetMarkerSize(1.3*MSI*((float)std::sqrt(Jets.at(J).Pt())/(float)10));
	FRAMEJ->SetMarkerColor(1);
	FRAMEJ->DrawClone("SAME P");	  
     

      }

                                                             /////////////////////////////////////////////////////
      ///////////////////////////////////////////////////////////             Run over All BSM                ///////////////////////////////////////////////////////////
                                                             /////////////////////////////////////////////////////

      vector <TLorentzVector> BSMs;    // vector of TLorentzVectors of all BSM objets in the event which require markers int the Event Display 
      vector <int> IDS;                // vector of BSM marker id number for later identification
                                       // IDS == 0 -> Z', 1 -> Dark Quark, 2 -> Vis Dark Hadron, 3 -> Inv Dark Hadrons

      float LeadQPt = 0;  //Stores the leading & subleading pt of the dark quarks for later use
      float SubQPt = 0;
      for(int B = 0; B < ntruthBSM[0]; B++){ // Loop over all BSM objets to find the dark quarks, Z' & Dark Hadrons 

	TEMP.SetPtEtaPhiE(truthBSM_pt->at(B),truthBSM_eta->at(B),truthBSM_phi->at(B),truthBSM_e->at(B));

	if(std::abs(truthBSM_pdgId->at(B)) == 5000001){ // Z'
	  BSMs.push_back(TEMP);
	  IDS.push_back(0);
	}
	if(std::abs(truthBSM_pdgId->at(B)) == 4900101 && truthBSM_status->at(B) == 23){ // Dark Quarks
	  if(truthBSM_pt->at(B) > LeadQPt){SubQPt = LeadQPt; LeadQPt = truthBSM_pt->at(B);}
	  else if(truthBSM_pt->at(B) >  SubQPt){SubQPt = truthBSM_pt->at(B);}
	  BSMs.push_back(TEMP);
	  IDS.push_back(1);
	}
	if(std::abs(truthBSM_pdgId->at(B)) == 4900211 || std::abs(truthBSM_pdgId->at(B)) == 4900213){ // Stable Dark Hadrons
	  BSMs.push_back(TEMP);
	  IDS.push_back(3);
	}
	if(std::abs(truthBSM_pdgId->at(B)) == 4900111 || std::abs(truthBSM_pdgId->at(B)) == 4900113){ // Unstable Dark Hadrons
	  BSMs.push_back(TEMP);
	  bool vis = false;
	  for(int C = 0; C < truthBSM_nChildren->at(B); C++){
	    int childID = std::abs(truthBSM_child_pdgId->at(B).at(C));
	    if(childID < 50 || (53 < childID && childID < 4900000)){vis = true;} // Checks if Hadron is Visible
	  }
	  if(vis){IDS.push_back(2);}
	  else{IDS.push_back(3);}
	}
      }


	/////////////////////////////////////////////////////
	////    Ordering BSM Objetcs In Descending Pt    ////
	/////////////////////////////////////////////////////


      float temppt = 0;
      int tempn = 0;
      for(int i = 0; i < BSMs.size();i++){  
	temppt = 0;
	tempn = 0;
	for(int l = 0; l < BSMs.size() - i; l++){
	  if(BSMs.at(l).Pt() > temppt){tempn = l; temppt = BSMs.at(i).Pt();}
	}
	BSMs.push_back(BSMs.at(tempn));
	BSMs.erase(BSMs.begin()+tempn);
	IDS.push_back(IDS.at(tempn));
	IDS.erase(IDS.begin()+tempn);
      }


	/////////////////////////////////////////////////////
	////             Drawing BSM Objetcs             ////
	/////////////////////////////////////////////////////


      int style[4] = {42,24,20,20};  // Style for BSM markers
      double size[4] = {1,1,1,1};    // Size for BSM markers
      int color[4] = {1,1,3,2};      // Colour for BSM markers

      
      for(int B = 0; B < BSMs.size(); B++){ // LOOPS OVER ALL BSM AND DRAWS INTO THE CANVAS
        
	TGraph* FRAME = new TGraph();
	FRAME->AddPoint(BSMs.at(B).Phi(),BSMs.at(B).Eta());
	FRAME->GetXaxis()->SetLimits(-TMath::Pi(),TMath::Pi());
	FRAME->GetYaxis()->SetLimits(-4.5,4.5);

	float fst = style[IDS.at(B)];
	float fsi = size[IDS.at(B)];
	float fcl = color[IDS.at(B)];
	
	FRAME->SetMarkerStyle(fst);
	FRAME->SetMarkerSize(1.3*fsi*((float)std::sqrt(BSMs.at(B).Pt())/(float)10));	  
	FRAME->SetMarkerColor(fcl);
	FRAME->DrawClone("SAME P");
      }
      
                                                             /////////////////////////////////////////////////////
      ///////////////////////////////////////////////////////////             Legend & Numbers                ///////////////////////////////////////////////////////////
                                                             /////////////////////////////////////////////////////




      C[j][e]->cd(0);   // Creates a new pad for the legend (top right) and the numbers (bottom right)

      P2[j][e] = new TPad((const TString)("P2"+std::to_string(e)+"_"+std::to_string(j)),(const TString)("P2"+std::to_string(e)),0.75,0,1,1); 
      P2[j][e]->Divide(1,2);
      P2[j][e]->Draw();
      

      P2[j][e]->cd(1); // Focus on top of new pad
      

	/////////////////////////////////////////////////////
	////      Creating & Drawing The Legend          //// (Could Be Improved By Just Making The Legend At The Beginning Instead Of Remaking It Each Time)
	/////////////////////////////////////////////////////


      // All the names, styles sizes and colours of the merkers and lines implemented in the display
      const int LN = 8;
      TString LNAMES[LN] = {"Z'","Dark Quarks","Visible Dark Hadron","Invisible Dark Hadron", "Jet","Leading/Subleading Jets","Large Jet Cones","MET"}; 
      int LCL[LN] = {1,1,3,2,1,1,40,46};
      int LST[LN] = {42,24,20,20,23,32,2,2};
      int LSI[LN] = {1,1,1,1,1,1,2,2};
 
      L[j][e] = new TLegend(0,0,1,1); 

      for(int l = 0; l < LN; l++){ // Adds each marker/line to the legend
	
	if(l != 7 && l != 6){ // For non large jet & non MET entries to the legend
	  TH2D* TEMPH = new TH2D((const TString)("TEMPH"+LNAMES[l]+std::to_string(e)),"",1,0,1,1,0,1); // makes temp TGraph with desiered marker for addition to legend
	  TEMPH->SetMarkerColor(LCL[l]);
	  TEMPH->SetMarkerStyle(LST[l]);
	  TEMPH->SetMarkerSize(1.3*LSI[l]);
	  L[j][e]->AddEntry(TEMPH,LNAMES[l],"P");
	  if(l == 5){  // For the jet radii
	    TArc* TEMPARC[2];
	    for(int t = 0; t < 2; t++){
	      TEMPARC[t] = new TArc(0,0,0.4); // Makes a temp arc for addition to legend
	      TEMPARC[t]->SetLineColor(1);
	      TEMPARC[t]->SetLineStyle(t+1);
	      TEMPARC[t]->SetFillStyle(0);
	      TEMPARC[t]->SetLineWidth(2);
	      if(t == 1){L[j][e]->AddEntry(TEMPARC[t],"Leading/Subleading Jet Cones","L");}
	      else{L[j][e]->AddEntry(TEMPARC[t],"Jet Cones","L");}
	    }
	  }
	}
	else{ // For large jet cone & MET line
	  TLine *TEMPL = new TLine(0,0,1,1);
	  TEMPL->SetLineColor(LCL[l]);
	  TEMPL->SetLineStyle(LST[l]);
	  TEMPL->SetLineWidth(LSI[l]);
	  if(l!=1){L[j][e]->AddEntry(TEMPL,LNAMES[l],"L");}
	}
      }

      L[j][e]->SetBorderSize(0); // Sets the style of the legend and draws it
      L[j][e]->SetTextFont(42);
      L[j][e]->SetTextSize(0.05);
      L[j][e]->SetFillStyle(0);
      L[j][e]->SetHeader((const TString)("Z"+filevar[j][0]+"r"+filevar[j][1]+"   Event:"+to_str(eventNumber[0])));
     
      L[j][e]->Draw();


	/////////////////////////////////////////////////////
	////              Draws The Numbers              ////
	/////////////////////////////////////////////////////

      P2[j][e]->cd(2);
      
      const int ntext = 5;
      TString TNAMES[ntext] = {"MET","Leading Jet Pt","Subleading Jet Pt","Leading Dark Quark Pt","Subleading Dark Quark Pt"}; // The number's names
      float TVALS[ntext] = {metFinalTrk, (float)Jets.at(0).Pt(), (float)Jets.at(1).Pt(), LeadQPt, SubQPt};                     // The actual numbers
      string TVNAMES[ntext];
      for(int t = 0; t < ntext; t++){
	TVNAMES[t] = to_str(TVALS[t]);
      }
      TText * TEXT[ntext];
      TText * TEXT2[ntext];
      for(int t = 0; t < ntext; t++){
	TEXT[t] = new TText(0, 0.9-0.06*(2*t), TNAMES[t]+" [GeV] :"); // Positions & draws the name
	TEXT[t]->SetTextAlign(11);
	TEXT[t]->SetTextSize(0.06);
	TEXT[t]->Draw("SAME");	
	TEXT2[t] = new TText(0, 0.9-0.06*(2*t+1)+0.015,(const TString)("  " + TVNAMES[t])); // Positions & draws the number
	TEXT2[t]->SetTextAlign(11);
	TEXT2[t]->SetTextSize(0.06);
	TEXT2[t]->Draw("SAME");
      }

	/////////////////////////////////////////////////////
	////             Saves The Display               ////
	/////////////////////////////////////////////////////



      TString filenames = "EVENTDISPLAYS/EVENT"+std::to_string(e)+"Z"+filevar[j][0]+"r"+filevar[j][1]; // Saves to the EVENTDISPLAYS directory


      C[j][e]->Print(filenames+".pdf");
      
    
    }
  }


}

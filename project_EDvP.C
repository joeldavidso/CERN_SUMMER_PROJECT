
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

const int nevents = 1; // The Number of events per signal point you want displayed & saved
const int nsignals = 1; // The Number of signal points you want 

const TString filevar[nsignals][2] = {{"1500","8"}}; // The Signal points you want to look at {MZ',r_inv}

const int events[nevents][nsignals] = {{3042}}; // The Event number(s) you want displayed and saved for each signal point {event# for signal point 1, ... , event# for signal point nsignals}



//Numbers specific to the signal location
//The Three Commented out elements correspond to 1500|3, 750|8, 750|3 (Mz'|r_inv)

const TString fileN[nsignals] = {"47"
				 // , "48", "49", "50"
};                                                                        
const TString ndN[nsignals][3] = {{"3733", "3907", "3996"}, 
				  //{"3741", "3916", "4007"}, {"3753", "3924", "4015"}, {"3764", "3929", "4025"}
};


double D(double phi){
  return phi*(double)360/(double)(2*TMath::Pi());
}


double R(double phi){
  return phi*(double)(2*TMath::Pi())/(double)360;
}

void project_EDvP(){ // start program

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

      chain->GetEntry(events[j][e]); // Gets Event


      double rad;

      TEllipse* circ;
      for(int R = 2; R < 5; R++){
	rad = 0.4*((double)R/(double)4)+0.00000001;
	circ = new TEllipse(0.5,0.5,rad,rad,0,360,0);
	circ->SetLineWidth(2);
	circ->SetFillStyle(0);
	circ->DrawClone("SAME");
      }

      TLine* marks;
      
      double scale[4] = {0.02,0.005,0.01,0.005};

      for(int q = 0; q < 8; q++){
	for(int t = 0; t < 4; t++){

	  double minx = 0.5 + (0.4 - scale[t])*sin(q*TMath::Pi()/(double)(4) + t*TMath::Pi()/(double)(16));
	  double maxx = 0.5 + (0.4 + scale[t])*sin(q*TMath::Pi()/(double)(4) + t*TMath::Pi()/(double)(16));

	  double miny = 0.5 + (0.4 - scale[t])*cos(q*TMath::Pi()/(double)(4) + t*TMath::Pi()/(double)(16));
	  double maxy = 0.5 + (0.4 + scale[t])*cos(q*TMath::Pi()/(double)(4) + t*TMath::Pi()/(double)(16));

	  marks = new TLine(minx,miny,maxx,maxy);

	  marks->SetLineWidth(2);
	  marks->DrawClone("SAME");
	}
      }
  

      TLatex *LABL;

      LABL = new TLatex(0.462,0.925,"\\phi = 0");
      LABL->SetTextSize(0.04);
      LABL->DrawClone("SAME");

      LABL = new TLatex(0.93,0.49,"#frac{\\pi}{2}");
      LABL->SetTextSize(0.04);
      LABL->DrawClone("SAME");

      LABL = new TLatex(0.478,0.05,"\\pm \\pi");
      LABL->SetTextSize(0.04);
      LABL->DrawClone("SAME");

      LABL = new TLatex(0.035,0.49,"- #frac{\\pi}{2}");
      LABL->SetTextSize(0.04);
      LABL->DrawClone("SAME");


      TLorentzVector TEMP(0,0,0,0);
     

                                                             /////////////////////////////////////////////////////
      ///////////////////////////////////////////////////////////                 Run over Jets               ///////////////////////////////////////////////////////////
                                                             /////////////////////////////////////////////////////

      double jetr;
      
      double difrad = 0.4*(double)(360)/(double)(2*TMath::Pi());

      TLine *JETL;
      TEllipse *JETE;

      /////////////////////////////////////////////////////
      ////           Small R Jets Main Loop            ////
      /////////////////////////////////////////////////////


      jetr = 0.4;

      vector <TLorentzVector> Jets;  // LorentzVectors of (Small R) jets (Used for jet cone drawing later)
      int LeadPt = 0;                // Leading jet's Pt
      int SubPt = 0;                 // Subleading jet's Pt 

      for(int J = 0; J < na4_pflowjets[0]; J++){ // Loop over all small R jets to add to the vectors
	
	TEMP.SetPtEtaPhiE(a4_pflowjets_pt->at(J),a4_pflowjets_eta->at(J),a4_pflowjets_phi->at(J),a4_pflowjets_E->at(J));
	
        if(J == 0){LeadPt = TEMP.Pt();}
	if(J == 1){SubPt = TEMP.Pt();}
	
	if(TEMP.Pt() > ptc && std::abs(TEMP.Phi()) < TMath::Pi() && std::abs(TEMP.Eta()) < 4.5){

	  Jets.push_back(TEMP);


	  double arcr = 0.1 + 0.27*(double)(TEMP.Pt() - 20)/(double)(Jets.at(0).Pt() - 20);

	  JETE = new TEllipse(0.5,0.5,arcr,arcr,0,2*difrad,-D(TEMP.Phi())-difrad+90);
	  
	  JETE->SetLineColor(kBlue);
	  JETE->SetLineWidth(2);
	  JETE->SetLineStyle(1);
	  JETE->SetFillStyle(0);
	  
	  JETE->DrawClone("SAME");
	}

	}




      /////////////////////////////////////////////////////
      ////          Next Radius Layer Drawn            ////
      /////////////////////////////////////////////////////
	
	  rad = 0.4*((double)1/(double)4)+0.00000001; // 2nd smallest radius
	  circ = new TEllipse(0.5,0.5,rad,rad,0,360,0);
	  circ->SetLineWidth(2);
	  circ->DrawClone("SAME");



      

      /////////////////////////////////////////////////////
      ////            Small R Jets 2nd Loop            ////
      /////////////////////////////////////////////////////

      

      for(int J = 0; J < Jets.size(); J++){ // Runs over all small R jets

	JETL = new TLine(0.5,0.5,0.5+0.4*sin(Jets.at(J).Phi()),0.5+0.4*cos(Jets.at(J).Phi()));
	JETL->SetLineWidth(3);
	if(J == 0 || J == 1){JETL->SetLineStyle(2);}
	JETL->DrawClone("SAME");

      }
      
                                                       /////////////////////////////////////////////////////
      ///////////////////////////////////////////////////////////                Draw MET Line                ///////////////////////////////////////////////////////////
                                                             /////////////////////////////////////////////////////
    
      
      TLine *METL = new TLine(0.5,0.5,0.5+0.4*sin(metFinalTrkPhi),0.5+0.4*cos(metFinalTrkPhi));
      METL->SetLineColor(kRed);
      METL->SetLineStyle(2);
      METL->SetLineWidth(3);
      METL->DrawClone();
    
      

      /////////////////////////////////////////////////////
      ////          Next Radius Layer Drawn            ////
      /////////////////////////////////////////////////////

      rad = 0.01;
      TArc *circa = new TArc(0.5,0.5,rad); // smallest radius
      circa->SetLineWidth(2);
      circa->SetLineStyle(1);
      circa->DrawClone("SAME");


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

      L[j][e] = new TLegend(0,0,1,1); 


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
      
      const int ntext = 3;
      TString TNAMES[ntext] = {"MET","Leading Jet Pt","Subleading Jet Pt"}; // The number's names
      float TVALS[ntext] = {metFinalTrk, (float)Jets.at(0).Pt(), (float)Jets.at(1).Pt()};                     // The actual numbers
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

//____________________________________________________________________________
/*!

\program Neutrino Energy reconstruction for CCQE events

\brief   selection of CCQE events : neutron momentum < 0.3 GeV

\author  R K pradhan

\created 2024

\cpright Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <string>

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TIterator.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TGraphSmooth.h>
#include <TCanvas.h>
#include <THStack.h>
#include <TVector3.h>
#include <TRandom3.h>

#include "Framework/EventGen/EventRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Ntuple/NtpMCFormat.h"
#include "Framework/Ntuple/NtpMCTreeHeader.h"
#include "Framework/Ntuple/NtpMCEventRecord.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/CmdLnArgParser.h"
#include "Framework/Interaction/ProcessInfo.h"
#include "Framework/Interaction/Interaction.h"
using std::string;
using namespace genie;
using namespace std;

void GetCommandLineArgs (int argc, char ** argv);

int    gOptNEvt;
string gOptInpFilename;

bool IsMeson(int pdgcode){    // function for meson
  int abscode = std::abs(pdgcode);
  std::string code = std::to_string(abscode);
  return(code.size()==3);
}

struct EveInfo{
  double nMom, pMom , TrueEnu, RecoEnu ;
};

EveInfo EventEnergy(const EventRecord& event)
{
  EveInfo info;
  
  GHepParticle* probe = event.Probe();
  GHepParticle * lep = event.FinalStatePrimaryLepton();
  GHepParticle * had = 0;

  double E_prob = probe->E();
  double El = lep->E();         //lepton energy
  double Ep = 0.0;              //proton energy  
  double px =0.0; double py=0.0; double pz=0.0;       //proton momentum 


  // calculating MX and  MY
  // for argon
  int Np = 18; int Nn = 22;       //number of proton and neutron 
  double mp = 0.93827208816; double mn = 0.9395654205;    //mass of proton and neutron
  double BE = 0.34384; double e = 0.02478;                // binding energy and separation energy
  double MX = (Np * mp) + (Nn * mn) - BE;                 // calculating M_A
  double MY = MX - mn + e;                                // M_(A-1)

  // // for carbon
  // int Np = 6; int Nn = 6;       //number of proton and neutron 
  // double mp = 0.93827208816; double mn = 0.9395654205;    //mass of proton and neutron
  // double BE = 0.09216; double e = 0.027139;                // binding energy and separation energy
  // double MX = (Np * mp) + (Nn * mn) - BE;                 // calculating M_A
  // double MY = MX - mn + e;                                // M_(A-1)


  TIter event_iter(&event);
  while((had =dynamic_cast<GHepParticle *>(event_iter.Next())))
  {
    if (had->Status() == kIStStableFinalState)
    {
      if (pdg::IsProton(had->Pdg()))
      {
        Ep = had->E();           // proton energy
        px = had->Px(); py = had->Py();  pz = had->Pz();   //proton momentum
      }
      }
      }

// momentum of muon
double lx = lep->Px(); double ly = lep->Py(); double lz = lep->Pz();
// TVector3 l(lx,ly,lz);

// mag of proton momentum

double pmom = sqrt(pow(px,2) + pow(py,2) + pow(pz,2));

// TVector3 p(px,py,pz);  //momentum of proton

// recostruction of neutron momentum
double Sq_nT = pow((lx+px),2) + pow((ly+py),2); // (perpendicular comp. of n (in paper pT))

double A = MX + lz + pz - El - Ep ;
double nL = (pow(A,2) - Sq_nT - pow(MY,2))/(2*A);


// total momentum

double Pn = std::sqrt(Sq_nT + pow(nL,2));

// neutrino energy reconstruction

double Enu = lz + pz - nL;

info.nMom = Pn;
info.pMom = pmom;
info.TrueEnu = E_prob;
info.RecoEnu =  Enu;

return info;
}

//___________________________________________________________________
int main(int argc, char ** argv)
{
  GetCommandLineArgs (argc, argv);

  //-- open the ROOT file and get the TTree & its header
  TTree *           tree = 0;
  NtpMCTreeHeader * thdr = 0;

  TFile file(gOptInpFilename.c_str(),"READ");

  tree = dynamic_cast <TTree *>           ( file.Get("gtree")  );
  thdr = dynamic_cast <NtpMCTreeHeader *> ( file.Get("header") );

  if(!tree) return 1;

// Historgrams

std::string filename = std::string(file.GetName());
std::string histofilename = filename.substr(0,filename.size()-5) +".out.root";
TFile *histofile = new TFile(histofilename.c_str(),"RECREATE");

TTree *tree1 = new TTree("neutron", "neutron");
TTree *tree2 = new TTree("selection","selection");

double pnFull ; double pnReco; int pID; int label; double Ereco; double Etrue; int sID; int EveID; double ppFull;

tree1->Branch("Process",&pID);
tree1->Branch("MomN",&pnFull);
tree1->Branch("MomP",&ppFull);

tree2->Branch("EventID",&EveID);
tree2->Branch("MomN",&pnReco);
tree2->Branch("label",&label);
tree2->Branch("ScatteringID",&sID);
tree2->Branch("RecoEnu",&Ereco);
tree2->Branch("TrueEnu",&Etrue);

TH1F* pnQE = new TH1F("QE","Reco Mom;Momentum (GeV/c);arb units",100,0.0,1.0);
TH1F* pnNonQE = new TH1F("Non-QE","Reco Mom;Momentum (GeV/c);arb units",100,0.0,1.0);




  NtpMCEventRecord * mcrec = 0;
  tree->SetBranchAddress("gmcrec", &mcrec);

  int nev = (gOptNEvt > 0) ?
        TMath::Min(gOptNEvt, (int)tree->GetEntries()) :
        (int) tree->GetEntries();

  int noSelection =0;
  int noSig =0;
  int noBkg=0;

  for(int i = 0; i < nev; i++)   //event loop
  {
    // get next tree entry
    tree->GetEntry(i);

    // get the GENIE event
    EventRecord &  event = *(mcrec->event);

    //Getting the process
    Interaction* summary  = event.Summary();
    const ProcessInfo & process = summary->ProcInfo();

    int Nproton = event.NEntries(2212,kIStStableFinalState,0);
    int Nmuon = event.NEntries(13,kIStStableFinalState,0);
    
    // counting number of pions in final state

    GHepParticle* pion = 0;
    TIter event_iterM(&event);
    int Npion = 0;
    while((pion =dynamic_cast<GHepParticle *>(event_iterM.Next())))
    {
      if (pion->Status() == kIStStableFinalState){
      if (pdg::IsPion(pion->Pdg()))
      // if (IsMeson(pion->Pdg()))
      {
        Npion++;
      }
    }}

    if (Nproton ==1 && Nmuon==1 && Npion==0)
    {
    // Kinematic calculation

    auto result = EventEnergy(event);
    
    // For the Full neutron momentum (1st tree)
    pID = process.ScatteringTypeId();
    pnFull = result.nMom;
    ppFull = result.pMom;

    if (process.IsQuasiElastic())
    {
      pnQE->Fill(result.nMom);
    }

    if (!process.IsQuasiElastic())
    {
      pnNonQE->Fill(result.nMom);
    }

    
    //Selection of Events
    // selection = Pn (initial) <= 0.3 GeV
    // final proton, p > 0.4 GeV

    if ((result.nMom < 0.300) && (result.pMom > 0.200))
    {
      noSelection++;

      // signal
      if (process.IsQuasiElastic())
      {
      noSig++;
      EveID=i;
      label = 1;
      pnReco = result.nMom;
      sID = process.ScatteringTypeId();
      Ereco = result.RecoEnu;
      Etrue = result.TrueEnu;
      // LOG("QE label",pNOTICE) <<"label: "<<label <<" event:" <<i;
      }

      // background
      else
      {
      noBkg++;
      EveID=i;
      label = 0;      
      pnReco = result.nMom;
      sID = process.ScatteringTypeId();
      Ereco = result.RecoEnu;
      Etrue = result.TrueEnu;
      // LOG("NonQE label",pNOTICE) <<"label: "<<label<<" event:" <<i <<" ID:" <<sID <<" pn:"<<pnReco<<" Ereco="<<Ereco;
      }

      tree2->Fill();

    }



tree1->Fill();
    }



    
// clear current mc event record
  mcrec->Clear();

  }//end loop over events

LOG("Events",pNOTICE)<<"Selected:"<<noSelection<<" sig:"<<noSig<<" bkg:"<<noBkg;

  tree1->Write();
  tree2->Write();

  // close input GHEP event file
  file.Close();

  histofile->Write();
  histofile->Close();

  LOG("myAnalysis", pNOTICE)  << "Done!";

  return 0;
}

//___________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("myAnalysis", pINFO) << "Parsing commad line arguments";

  CmdLnArgParser parser(argc,argv);

  // get GENIE event sample
  if( parser.OptionExists('f') ) {
    LOG("myAnalysis", pINFO) 
       << "Reading event sample filename";
    gOptInpFilename = parser.ArgAsString('f');
  } else {
    LOG("myAnalysis", pFATAL) 
        << "Unspecified input filename - Exiting";
    exit(1);
  }

  // number of events to analyse
  if( parser.OptionExists('n') ) {
    LOG("myAnalysis", pINFO) 
      << "Reading number of events to analyze";
    gOptNEvt = parser.ArgAsInt('n');
  } else {
    LOG("myAnalysis", pINFO)
      << "Unspecified number of events to analyze - Use all";
    gOptNEvt = -1;
  }
}
//_________________________________________________________________________________

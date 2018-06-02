// -*- C++ -*-
//
// Package:    SUSYAnalysis/TrigAnalyzerMiniAOD
// Class:      TrigAnalyzerMiniAOD
//
/**\class TrigAnalyzerMiniAOD TrigAnalyzerMiniAOD.cc SUSYAnalysis/TrigAnalyzerMiniAOD/plugins/TrigAnalyzerMiniAOD.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Oleksii Turkot
//         Created:  Tue, 22 May 2018 22:32:27 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
 #include "FWCore/Utilities/interface/InputTag.h"
 #include "DataFormats/TrackReco/interface/Track.h"
 #include "DataFormats/TrackReco/interface/TrackFwd.h"
 
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"

#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"


#include "SUSYAnalysis/TrigAnalyzerMiniAOD/interface/TrigAnalyzerMiniAOD.h"

// ROOT includes
#include "Math/VectorUtil.h"
#include "Math/Vector4D.h"

#include <cassert>

using namespace reco;
using namespace edm;

//
// constructors and destructor
//
//____________________________________________________________________________
TrigAnalyzerMiniAOD::TrigAnalyzerMiniAOD(const edm::ParameterSet& ps) 
{
  using namespace std;
  using namespace edm;

  processName_ = ps.getUntrackedParameter<std::string>("processName","HLT");
//  refTriggerName_ = ps.getUntrackedParameter<std::string>("refTriggerName","HLT_PFJet450_v20");
    refTriggerName_ = ps.getUntrackedParameter<std::string>("refTriggerName");
//  sigTriggerName_ = ps.getUntrackedParameter<std::string>("sigTriggerName","HLT_Ele15_IsoVVVL_PFHT450_v15");
      sigTriggerName_ = ps.getUntrackedParameter<std::string>("sigTriggerName");
//  triggerResultsToken_ = consumes<edm::TriggerResults> (ps.getUntrackedParameter<edm::InputTag>("triggerResultsTag", edm::InputTag("TriggerResults", "", "HLT")));
  pfMetToken_ = consumes<edm::View<pat::MET> >(ps.getUntrackedParameter<edm::InputTag>("pfMetInputTag_", edm::InputTag("slimmedMETs")));
  verbose_ = ps.getUntrackedParameter<bool>("verbose",false);

  
  electronCollection_  = consumes<edm::View<pat::Electron> > (ps.getParameter<edm::InputTag>("electronSrc"));
  muonCollection_      = consumes<edm::View<pat::Muon> >     (ps.getParameter<edm::InputTag>("muonSrc"));
  jetsCollection_      = consumes<edm::View<pat::Jet> >      (ps.getParameter<edm::InputTag>("jetSrc"));
  triggerResultsToken_ = consumes<edm::TriggerResults>       (ps.getParameter<edm::InputTag>("triggerResultsTag"));
  
  
  size_t idTmp = 0;
  // possible reference triggers:
  stTriggNames_.push_back("HLT_PFJet450");
        // new from Isabell
  stTriggNames_.push_back("HLT_PFJet500");
  stTriggNames_.push_back("HLT_PFJet550");
  stTriggNames_.push_back("HLT_PFJet60");
  
  // EleOR triggers:
  idTmp = stTriggNames_.size();     // index of the first EleOR trigger
  stTriggNames_.push_back("HLT_Ele115_CaloIdVT_GsfTrkIdT");
  stTriggNames_.push_back("HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165");
  stTriggNames_.push_back("HLT_Ele27_WPTight_Gsf");
  stTriggNames_.push_back("HLT_Ele15_IsoVVVL_PFHT450");
        // new from Isabell
  stTriggNames_.push_back("HLT_Ele15_IsoVVVL_PFHT450_CaloBTagCSV_4p5");
  stTriggNames_.push_back("HLT_Ele15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5"); // different name in later runs
  stTriggNames_.push_back("HLT_Ele15_IsoVVVL_PFHT450_PFMET50");
  stTriggNames_.push_back("HLT_Ele15_IsoVVVL_PFHT600");
  for(size_t i=idTmp; i < stTriggNames_.size(); ++i) {
    idTriggEleOR_.push_back(i);     // save indexes of EleOR triggers
  }
  
  // MuOR triggers:
  idTmp = stTriggNames_.size();     // index of the first MuOR trigger
  stTriggNames_.push_back("HLT_Mu50");
  stTriggNames_.push_back("HLT_IsoMu24");
  stTriggNames_.push_back("HLT_Mu15_IsoVVVL_PFHT450");
        // new from Isabell
  stTriggNames_.push_back("HLT_Mu15_IsoVVVL_PFHT450_CaloBTagCSV_4p5");
  stTriggNames_.push_back("HLT_Mu15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5"); // different name in later runs
  stTriggNames_.push_back("HLT_Mu15_IsoVVVL_PFHT450_PFMET50");
  stTriggNames_.push_back("HLT_Mu15_IsoVVVL_PFHT600");
  stTriggNames_.push_back("HLT_Mu50");
  stTriggNames_.push_back("HLT_IsoMu27_MET90");
  for(size_t i=idTmp; i < stTriggNames_.size(); ++i) {
    idTriggMuOR_.push_back(i);     // save indexes of MuOR triggers
  }
  
  // MetOR triggers:
  idTmp = stTriggNames_.size();     // index of the first MetOR trigger
  stTriggNames_.push_back("HLT_PFMET100_PFMHT100_IDTight_PFHT60");
  stTriggNames_.push_back("HLT_PFMET110_PFMHT110_IDTight");
  stTriggNames_.push_back("HLT_PFMET120_PFMHT120_IDTight");
        // new from Isabell
  stTriggNames_.push_back("HLT_PFHT500_PFMET100_PFMHT100_IDTight");
  stTriggNames_.push_back("HLT_PFHT500_PFMET110_PFMHT110_IDTight");
  stTriggNames_.push_back("HLT_PFMET200_HBHE_BeamHaloCleaned");
  stTriggNames_.push_back("HLT_PFMET250_HBHECleaned");
  stTriggNames_.push_back("HLT_PFMET300_HBHECleaned");
  for(size_t i=idTmp; i < stTriggNames_.size(); ++i) {
    idTriggMetOR_.push_back(i);     // save indexes of MetOR triggers
  }
  
  // EleOR, MuOR and MetOR themselves:
  idEleOR_ = stTriggNames_.size();
  stTriggNames_.push_back("HLT_EleOR");
  idMuOR_ = stTriggNames_.size();
  stTriggNames_.push_back("HLT_MuOR");
  idMetOR_ = stTriggNames_.size();
  stTriggNames_.push_back("HLT_MetOR");
  
  if (stTriggNames_.size() > nTriggMAX) {
    cout << "ERROR: " << stTriggNames_.size() << " triggers specified, while nTriggMAX = " << nTriggMAX << "." << endl;
    exit(1);
  }
  
  // Trigger bits:
  for(size_t i=0; i < stTriggNames_.size(); ++i) {
    blTriggBits_[i] = false;
  }
  
  if (verbose_) {
    cout << endl;
    cout << "EleOR triggers:" << endl;
    for(size_t i=0; i < idTriggEleOR_.size(); ++i) cout << "   " << stTriggNames_[idTriggEleOR_[i]] << endl;
    cout << "MuOR triggers:" << endl;
    for(size_t i=0; i < idTriggMuOR_.size(); ++i)  cout << "   " << stTriggNames_[idTriggMuOR_[i]] << endl;
    cout << "MetOR triggers:" << endl;
    for(size_t i=0; i < idTriggMetOR_.size(); ++i) cout << "   " << stTriggNames_[idTriggMetOR_[i]] << endl;
    cout << endl;
  }
  
  
  // check that there are no trigger ANDs in reference and signal triggers
  if ( (refTriggerName_.find("&") != string::npos) ||
       (sigTriggerName_.find("&") != string::npos) ) {
    cout << "ERROR: only trigger ORs are allowed in reference and signal triggers." << endl;
    exit(1);
  }
  // set ids for the reference and signal triggers
  for (size_t i = 0; i < stTriggNames_.size(); ++i) {
    if (refTriggerName_.find(stTriggNames_[i]) != string::npos) idRefTrigg_.push_back(i);
    if (sigTriggerName_.find(stTriggNames_[i]) != string::npos) idSigTrigg_.push_back(i);
  }
  // abort on invalid trigger name
  if (idRefTrigg_.empty()) {
    cout << "ERROR: " << refTriggerName_ << " - not found!" << endl;
    exit(1);
  }
  if (idSigTrigg_.empty()) {
    cout << "ERROR: " << sigTriggerName_ << " - not found!" << endl;
    exit(1);
  }
  
  // histogram setup
  edm::Service<TFileService> fs;
  tree_    = fs->make<TTree>("EventTree", "LW Trigger Event data 2018");
  hists_1d_["h_passreftrig"] = fs->make<TH1F>("h_passreftrig" , "; passed ref trigger" , 2 , 0. , 2. );
  hists_1d_["h_met_all"] = fs->make<TH1F>("h_met_all" , "; E_{T}^{miss} [GeV]" , 40, 100., 500. );
  hists_1d_["h_met_passtrig"] = fs->make<TH1F>("h_met_passtrig" , "; E_{T}^{miss} [GeV]" , 40, 100., 500. );

  branchesEvent(tree_);
  
}

//____________________________________________________________________________
TrigAnalyzerMiniAOD::~TrigAnalyzerMiniAOD()
{
}

//
// member functions
//
//____________________________________________________________________________
void
TrigAnalyzerMiniAOD::beginRun(edm::Run const & iRun, edm::EventSetup const& iSetup)
{
  using namespace std;
  using namespace edm;

  bool changed(true);
  if (hltConfig_.init(iRun,iSetup,processName_,changed)) {
    if (changed) {
/*
      const unsigned int n(hltConfig_.size());
      // check if trigger names in (new) config
      unsigned int refTriggerIndex(hltConfig_.triggerIndex(refTriggerName_));
      if (refTriggerIndex>=n) {
	cout << "TrigAnalyzerMiniAOD::analyze:"
	     << " TriggerName " << refTriggerName_ 
	     << " not available in config!" << endl;
      }
      unsigned int sigTriggerIndex(hltConfig_.triggerIndex(sigTriggerName_));
      if (sigTriggerIndex>=n) {
	cout << "TrigAnalyzerMiniAOD::analyze:"
	     << " TriggerName " << sigTriggerName_ 
	     << " not available in config!" << endl;
      }
*/
    } // if changed
  } else {
    cout << "TrigAnalyzerMiniAOD::analyze:"
	 << " config extraction failure with process name "
	 << processName_ << endl;
  }

}

//____________________________________________________________________________
// ------------ method called to produce the data  ------------
void
TrigAnalyzerMiniAOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace std;
  using namespace edm;
  using namespace reco;
  using namespace trigger;

  
  // fill event variables, check preselection result
  if (!fillEvent(iEvent, iSetup)) return;
  
  
  tree_->Fill();
  
  
  // modules on this trigger path
  bool refAccept = getTriggRes(idRefTrigg_);
  bool sigAccept = getTriggRes(idSigTrigg_);
  
  if (refAccept) hists_1d_["h_passreftrig"]->Fill(1);
  else {  
    // don't consider event if reference trigger didn't fire
    hists_1d_["h_passreftrig"]->Fill(0);
    return;
  }
  // fill MET distributions
  hists_1d_["h_met_all"]->Fill(MET_);
  if (sigAccept) hists_1d_["h_met_passtrig"]->Fill(MET_);
  
  
  return;
}

void TrigAnalyzerMiniAOD::branchesEvent(TTree* tree) {
  
  tree->Branch("run",         &run_);
  tree->Branch("event",       &event_);
  tree->Branch("lumis",       &lumis_);
  tree->Branch("isData",      &isData_);
  
  tree->Branch("nEle",        &nEle_);
  tree->Branch("elePt",       &elePt_);
  
  tree->Branch("nMu",         &nMu_);
  tree->Branch("muPt",        &muPt_);
  
  tree->Branch("nLep",        &nLep_);
  tree->Branch("lepPt",       &lepPt_);
  
  tree->Branch("nJet30",      &nJet_);
  tree->Branch("leadJet30Pt", &leadJetPt_);
  
  tree->Branch("MET",         &MET_);
  tree->Branch("LT",          &LT_);
  tree->Branch("HT",          &HT_);
  
  for(size_t i=0; i < stTriggNames_.size(); ++i) {
    tree->Branch(stTriggNames_[i].c_str(),      &(blTriggBits_[i]));
  }
  
}

bool TrigAnalyzerMiniAOD::fillEvent(const edm::Event& e, const edm::EventSetup& es) {
  
  using namespace std;
    
  run_       = e.id().run();
  event_     = e.id().event();
  lumis_     = e.luminosityBlock();
  isData_    = e.isRealData();
  
  
  nEle_      = 0;
  elePt_     = -10.0;
  ROOT::Math::PtEtaPhiEVector eleP4;   // used in jets selection
  edm::Handle<edm::View<pat::Electron> > electronHandle;
  e.getByToken(electronCollection_, electronHandle);
  if (!electronHandle.isValid()) {
    edm::LogWarning("TrigAnalyzerMiniAOD") << "no pat::Electrons in event";
    return false;
  }
  nEle_      = electronHandle->size();
  if (nEle_ > 1) return false;         // we are interested in nLep == 1.
  if (nEle_ > 0) {
    edm::View<pat::Electron>::const_iterator iEle = electronHandle->begin();
    elePt_   = iEle->pt();
    eleP4    = iEle->p4();
  }
  
  
  nMu_       = 0;
  muPt_      = -10.0;
  ROOT::Math::PtEtaPhiEVector muP4;   // used in jets selection
  edm::Handle<edm::View<pat::Muon> > muonHandle;
  e.getByToken(muonCollection_, muonHandle);
  if (!muonHandle.isValid()) {
    edm::LogWarning("TrigAnalyzerMiniAOD") << "no pat::Muons in event";
    return false;
  }
  nMu_       = muonHandle->size();
  if (nMu_ > 1) return false;          // we are interested in nLep == 1.
  if (nMu_ > 0) {
    edm::View<pat::Muon>::const_iterator iMu = muonHandle->begin();
    muPt_    = iMu->pt();
    muP4     = iMu->p4();
  }
  
  
  nLep_      = nEle_ + nMu_;
  ROOT::Math::PtEtaPhiEVector lepP4;   // used in jets selection
  if (nLep_ != 1) return false;        // preselection of nLep == 1.
  if (nEle_ > 0) {
    lepPt_ = elePt_;
    lepP4  = eleP4;
  } else {
    lepPt_ = muPt_;
    lepP4  = muP4;
  }
  if (lepPt_ < 15.0) return false;     // preselection of lepPt >= 15 GeV.
  
  
  nJet_      = 0;
  leadJetPt_ = -10.0;
  HT_        = 0.;
  edm::Handle<edm::View<pat::Jet> > jetHandle;
  e.getByToken(jetsCollection_, jetHandle);
  if (!jetHandle.isValid()) {
    edm::LogWarning("TrigAnalyzerMiniAOD") << "no pat::Jets in event";
    return false;
  }
  // central jets with pt > 30 GeV, |eta| < 2.4 and excluded lepton jet:
  std::vector<pat::Jet> centralJet30clean;
  // temporary variable to find jet nearest to te lepton:
  size_t idjNear = 0;
  float    dR    = 2.0;
  float    dRmin = 2.0;
  for(edm::View<pat::Jet>::const_iterator jet=jetHandle->begin(); jet!=jetHandle->end(); ++jet){
    if ( (     jet->pt()   > 30)  &&
         (fabs(jet->eta()) < 2.4) ) {
      centralJet30clean.push_back(*jet);
      dR = ROOT::Math::VectorUtil::DeltaR(jet->p4(), lepP4);
      if (dR < dRmin) {
        dRmin   = dR;
        idjNear = centralJet30clean.size() - 1;
      }
    }
  } // close loop over event jets
  if (dRmin < 0.4) centralJet30clean.erase(centralJet30clean.begin() + idjNear);  // remove one jet closest to the lepton with dR < 0.4
  
  if (centralJet30clean.empty()) return false; // preselection of at least 1 jet with Pt >= 30 GeV, |eta| < 2.4 and dR(jet, lep) > 0.4
  for (std::vector<pat::Jet>::iterator itj = centralJet30clean.begin(); itj != centralJet30clean.end(); ++itj) {
    nJet_++;
    if(itj->pt() > leadJetPt_) leadJetPt_ = itj->pt();
    // evaluation of HT
    HT_ += itj->pt();
  }
  
  
  // retrieve MET container
  Handle<edm::View<pat::MET> > pfMetHandle;
  e.getByToken(pfMetToken_, pfMetHandle);
  if (!pfMetHandle.isValid()) {
    edm::LogWarning("TrigAnalyzerMiniAOD") << "no pat::MET in event";
    return false;
  }
  MET_ = ( pfMetHandle->front() ).pt();
  
  
  // evaluation of LT
  LT_ = lepPt_ + MET_;
  
  
  // set trigger bits to false
  for(size_t i=0; i < stTriggNames_.size(); ++i) {
    blTriggBits_[i] = false;
  }
  // retrieve triggers container
  e.getByToken(triggerResultsToken_, triggerResultsHandle_);
  if (!triggerResultsHandle_.isValid()) {
    edm::LogWarning("TrigAnalyzerMiniAOD") << "no TriggerResults in event";
    return false;
  }
  // sanity check
  assert(triggerResultsHandle_->size()==hltConfig_.size());
  // loop over all event triggers
  const edm::TriggerNames &trgNames = e.triggerNames(*triggerResultsHandle_);
  for (size_t i = 0; i < trgNames.size(); ++i) {
    const string &name = trgNames.triggerName(i);
    // check if the event trigger correspond to any of our triggers
    for (size_t j = 0; j < stTriggNames_.size(); ++j) {
      if (name.find( stTriggNames_[j]+string("_v") ) != string::npos) blTriggBits_[j] = triggerResultsHandle_->accept(i);
    }
    // sanity check
    assert(i == hltConfig_.triggerIndex(name));
  }
  // fill bits of EleOR, MuOR and MetOR triggers.
  blTriggBits_[idEleOR_] = getTriggRes(idTriggEleOR_);
  blTriggBits_[idMuOR_]  = getTriggRes(idTriggMuOR_);
  blTriggBits_[idMetOR_] = getTriggRes(idTriggMetOR_);
  
  
  return true;
}

// additional function
bool TrigAnalyzerMiniAOD::getTriggRes(std::vector<size_t> v) {
  bool res = false;  
  for (std::vector<size_t>::iterator it = v.begin(); it != v.end(); ++it) res = (res || blTriggBits_[*it]);
  return res;
}

DEFINE_FWK_MODULE(TrigAnalyzerMiniAOD);



#ifndef HLTcore_TrigAnalyzerMiniAOD_h
#define HLTcore_TrigAnalyzerMiniAOD_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1.h"
#include "TH2.h"

#include "TTree.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

//
// class declaration
//
class TrigAnalyzerMiniAOD : public edm::EDAnalyzer {
  
  typedef math::XYZTLorentzVectorF LorentzVector;

 public:
  explicit TrigAnalyzerMiniAOD(const edm::ParameterSet&);
  ~TrigAnalyzerMiniAOD();

  virtual void beginRun(edm::Run const &, edm::EventSetup const&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  bool analyzeTrigger(const edm::Event&, const edm::EventSetup&, const std::string& triggerName);

 private:

  /// module config parameters
  std::string   processName_;
  std::string   refTriggerName_;
  std::string   sigTriggerName_;
  bool verbose_;
  
  std::vector<size_t>  idRefTrigg_;
  std::vector<size_t>  idSigTrigg_;
  
  bool getTriggRes(std::vector<size_t>);

  /// additional class data memebers
  edm::EDGetTokenT<edm::TriggerResults>       triggerResultsToken_;
  edm::Handle<edm::TriggerResults>            triggerResultsHandle_;
  HLTConfigProvider hltConfig_;

  TTree   *tree_;
  std::map<std::string,TH1F*> hists_1d_;
  
  Int_t       run_;
  Long64_t    event_;
  Int_t       lumis_;
  Bool_t      isData_;
  
  edm::EDGetTokenT<edm::View<pat::Electron> > electronCollection_;
  Int_t       nEle_;
  float       elePt_;
  
  edm::EDGetTokenT<edm::View<pat::Muon> >     muonCollection_;
  Int_t       nMu_;
  float       muPt_;
  
  Int_t       nLep_;
  float       lepPt_;
  
  edm::EDGetTokenT<edm::View<pat::Jet> >      jetsCollection_;
  Int_t       nJet_;
  float       leadJetPt_;
  
  edm::EDGetTokenT<edm::View<pat::MET> >      pfMetToken_;
  float       MET_;
  
  float       LT_;
  float       HT_;
  
  const static size_t nTriggMAX = 100;
  std::vector<std::string> stTriggNames_;
//  std::vector<bool>        blTriggBits_;
  bool                     blTriggBits_[nTriggMAX]; // blTriggBits_ elements should never reallocate so their references could be used for setting tree_ branches -> using not vector but array.
  
  size_t  idEleOR_, idMuOR_, idMetOR_;
  std::vector<size_t>      idTriggEleOR_;
  std::vector<size_t>      idTriggMuOR_;
  std::vector<size_t>      idTriggMetOR_;
  
  void branchesEvent(TTree*);
  bool fillEvent(const edm::Event&, const edm::EventSetup&);
  
};
#endif

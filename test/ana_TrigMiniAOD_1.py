import FWCore.ParameterSet.Config as cms

import FWCore.Utilities.FileUtils as FileUtils

process = cms.Process('HLTANALYZER')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")

process.maxEvents = cms.untracked.PSet(
#    input = cms.untracked.int32(-1)
    input = cms.untracked.int32(100000)
)
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

# Input from a file
listfiles = FileUtils.loadListFromFile ('JetHT_Run2018Av2_all.txt') 

# Input source
process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring(*listfiles),
)


process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.20 $'),
    annotation = cms.untracked.string('RelVal nevts:100'),
    name = cms.untracked.string('Applications')
)

# Output definition

process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string('ntup_TrigAnalyzer_2.root')
                                   )

### analyzer configuration

process.TrigAnalyzerMiniAOD = cms.EDAnalyzer("TrigAnalyzerMiniAOD")
process.TrigAnalyzerMiniAOD.refTriggerName = cms.untracked.string("HLT_PFJet500")
#process.TrigAnalyzerMiniAOD.sigTriggerName = cms.untracked.string("HLT_Ele27_WPTight_Gsf_v15")
#process.TrigAnalyzerMiniAOD.sigTriggerName = cms.untracked.string("HLT_Ele15_IsoVVVL_PFHT450_v15")
process.TrigAnalyzerMiniAOD.sigTriggerName = cms.untracked.string("HLT_EleOR||HLT_MetOR")
#process.TrigAnalyzerMiniAOD.verbose = cms.untracked.bool(False)
process.TrigAnalyzerMiniAOD.verbose = cms.untracked.bool(True)

process.TrigAnalyzerMiniAOD.electronSrc       = cms.InputTag("slimmedElectrons")
process.TrigAnalyzerMiniAOD.muonSrc           = cms.InputTag("slimmedMuons")
process.TrigAnalyzerMiniAOD.jetSrc            = cms.InputTag("slimmedJets")
process.TrigAnalyzerMiniAOD.triggerResultsTag = cms.InputTag("TriggerResults", "", "HLT")

#triggerResultsToken_ = consumes<edm::TriggerResults> (ps.getUntrackedParameter<edm::InputTag>("triggerResultsTag", edm::InputTag("TriggerResults", "", "HLT")));

#process.GlobalTag.globaltag = "80X_dataRun2_2016SeptRepro_v3"

# Path and EndPath definitions
process.HLTanalyzers = cms.Path(process.TrigAnalyzerMiniAOD)

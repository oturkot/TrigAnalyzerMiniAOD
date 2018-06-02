import FWCore.ParameterSet.Config as cms

import FWCore.Utilities.FileUtils as FileUtils

process = cms.Process('HLTANALYZER')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input from a file
listfiles = FileUtils.loadListFromFile ('JetHT_Run2018A_all.txt') 
#readFiles = cms.untracked.vstring( *listfiles)

# Input source
process.source = cms.Source("PoolSource",
#   fileNames = cms.untracked.vstring('/store/data/Run2016G/SingleElectron/MINIAOD/23Sep2016-v1/100000/004A7893-A990-E611-B29F-002590E7DE36.root'),
#   fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/data/Run2018A/JetHT/MINIAOD/PromptReco-v2/000/316/239/00000/02904165-4259-E811-AC97-FA163EB783E7.root'), # skimmed file on EOS at LPC
#    readFiles = cms.untracked.vstring( *listfiles),
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
                                       fileName = cms.string('histos_TrigAnalyzer_all.root')
                                   )

### analyzer configuration

process.TrigAnalyzerMiniAOD = cms.EDAnalyzer("TrigAnalyzerMiniAOD")
process.TrigAnalyzerMiniAOD.refTriggerName = cms.untracked.string("HLT_PFJet450_v20")
#process.TrigAnalyzerMiniAOD.sigTriggerName = cms.untracked.string("HLT_Ele27_WPTight_Gsf_v15")
process.TrigAnalyzerMiniAOD.sigTriggerName = cms.untracked.string("HLT_Ele15_IsoVVVL_PFHT450_v15")
process.TrigAnalyzerMiniAOD.verbose = cms.untracked.bool(False)

#process.GlobalTag.globaltag = "80X_dataRun2_2016SeptRepro_v3"

# Path and EndPath definitions
process.HLTanalyzers = cms.Path(process.TrigAnalyzerMiniAOD)

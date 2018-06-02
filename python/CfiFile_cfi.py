import FWCore.ParameterSet.Config as cms

demo = cms.EDAnalyzer('TrigAnalyzerMiniAOD'
     ,tracks = cms.untracked.InputTag('ctfWithMaterialTracks')
)

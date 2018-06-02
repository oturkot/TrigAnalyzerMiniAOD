
from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'test_LWTrigMiniAOD_01'
# from WorkBook:
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
#--------------#
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
# Name of the CMSSW configuration file
config.JobType.psetName = 'ana_TrigMiniAOD.py'

config.Data.inputDataset = '/JetHT/Run2018A-PromptReco-v2/MINIAOD'
# from WorkBook:
config.Data.inputDBS = 'global'
#--------------#
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 20
# for testing:
config.Data.totalUnits = 20
config.Data.publication = True
# This string is used to construct the output dataset name
config.Data.outputDatasetTag = 'CRAB3_Analysis_test_LWTrigMiniAOD_01'

# These values only make sense for processing data
#    Select input data based on a lumi mask
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions18/13TeV/PromptReco/Cert_314472-316723_13TeV_PromptReco_Collisions18_JSON.txt'
#    Select input data based on run-ranges
config.Data.runRange = '314472-316723'
# from WorkBook:
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
#--------------#

# Where the output files will be transmitted to
config.Site.storageSite = 'T2_DE_DESY'



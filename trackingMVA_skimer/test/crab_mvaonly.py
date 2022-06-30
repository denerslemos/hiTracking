if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException, UsernameException
    from httplib import HTTPException

    from CRABClient.UserUtilities import config, getUsername
    config = config()

    config.General.requestName ='trkMVA_ntuple_generalTk_mvaonly'
    config.General.workArea = 'jobStatus'
    config.General.transferOutputs = True
    config.General.transferLogs = False

    config.JobType.psetName = 'runTuple_mvaonly.py'
    config.JobType.outputFiles = ['output.root']
    config.JobType.inputFiles = ['extra_producer.py']
    config.JobType.maxMemoryMB = 4000
    config.JobType.numCores = 1

    config.Data.splitting = 'FileBased'
    config.Data.unitsPerJob = 1
    NJOBS = 100  # This is not a configuration parameter, but an auxiliary variable that we use in the next line.
    config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
    config.Data.publication = False
    config.Data.outLFNDirBase = '/store/user/%s/' % (getUsername())
    config.Data.inputDataset = '/MinBias_Hydjet_Drum5F_5p02TeV/Run3Winter22PbPbNoMixDIGI-122X_mcRun3_2021_realistic_HI_v10-v3/GEN-SIM-DIGI-RAW-HLTDEBUG'
    config.Data.inputDBS = 'global'
    config.Data.outputDatasetTag = config.General.requestName
    config.JobType.allowUndistributedCMSSW = True

    config.Site.storageSite = 'T2_CH_CERN'

    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException as hte:
            print "Failed submitting task: %s" % (hte.headers)
        except ClientException as cle:
            print "Failed submitting task: %s" % (cle)

    submit(config)

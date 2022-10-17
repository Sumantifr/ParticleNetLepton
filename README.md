# ParticleNetLepton
Lepton identification using ParticleNet. 

## Ntuple production 

- cmsrel CMSSW_10_6_29

- cd $CMSSW_BASE/src

- git clone git@github.com:Sumantifr/ParticleNetLepton.git . <br/> 
  OR <br/>
  git clone https://github.com/Sumantifr/ParticleNetLepton.git . <br/>
  *(Don't forget '.')*

- cmsenv

- scram b -j10

- cd DeepLeptonNtuples/NTuplizer/test/

- voms-proxy-init -rfc -voms cms -valid 96:00

- cmsRun PNet_MC_MINIAOD_cfg.py 


## For using XGBoost in CMSSW (needed to evaluate TopLeptonMVA) ##

 Taken from git clone https://github.com/simonepigazzini/XGBoostCMSSW.git

- cd $CMSSW_BASE/src

- cmsenv

- git clone https://github.com/simonepigazzini/XGBoostCMSSW.git

- cp XGBoostCMSSW/XGBoostInterface/toolbox/\*xml $CMSSW_BASE/config/toolbox/$SCRAM_ARCH/tools/selected/

- scram setup rabit

- scram setup xgboost

- scram b -j10

### Running crab jobs on MC samples 

- *If needed*, 'source /cvmfs/cms.cern.ch/crab3/crab.sh' OR 'source /cvmfs/cms.cern.ch/common/crab-setup.sh'

- check sample names (& nicknames) in generate_crab_miniaod_UL18.sh

- If okay, run './generate_crab_miniaod_UL18.sh' <br/>
  This will not submit crab jobs, but prepare all inputs needed for crab jobs & produce two scripts 'crab_submit.sh' and 'crab_monitor.sh', and cfg python files for each sample

- If happy, submit the jobs by running './crab_submit.sh'

- Monitor the job status by running './crab_monitor.sh'

- *Don't forget to remove crab_submit.sh and crab_monitor.sh before submitting new jobs*

- *Important* Running crab jobs with XGBoost will give problem due to shared libraries. A quick fix is the following:
  -- check the original libxgboost.so: 
  -- ls -ltr $CMSSW_BASE/external/$SCRAM_ARCH/lib/libxgboost.so
  -- it will show you a libxgboost.so file in /cvmfs/cms.cern.... path
  -- copy the original libxgboost.so to $CMSSW_BASE/lib/$SCRAM_ARCH

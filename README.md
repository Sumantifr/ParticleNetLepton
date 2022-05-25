# ParticleNetLepton
Lepton identification using ParticleNet. 

## Follow the steps below to set up the framework

- cmsrel CMSSW_10_6_29

- cd CMSSW_10_6_29/src

- git clone git@github.com:Sumantifr/ParticleNetLepton.git . <br/> 
  OR <br/>
  git clone https://github.com/Sumantifr/ParticleNetLepton.git .

- cmsenv

- scram b -j10

- cd DeepLeptonNtuples/NTuplizer/test/

*Create voms proxy*

- cmsRun PNet_MC_MINIAOD_cfg.py 

## Running crab jobs on MC samples 

- cd DeepLeptonNtuples/NTuplizer/test/

- cmsenv

- *If needed*, 'source /cvmfs/cms.cern.ch/crab3/crab.sh' OR 'source /cvmfs/cms.cern.ch/common/crab-setup.sh'

*Create voms proxy*

- check sample names (& nicknames) in generate_crab_miniaod_UL18.sh

- If okay, run './generate_crab_miniaod_UL18.sh' <br/>
  This will not submit crab jobs, but prepare all inputs needed for crab jobs & produce two scripts 'crab_submit.sh' and 'crab_monitor.sh', and cfg python files for each sample

- If happy, submit the jobs by running './crab_submit.sh'

- Monitor the job status by running './crab_monitor.sh'

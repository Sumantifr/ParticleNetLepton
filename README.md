# ParticleNetLepton
Lepton identification using ParticleNet. 

## Follow the steps below to set up the framework

- cmsrel CMSSW_10_6_29

- cd CMSSW_10_6_29/src

- git clone git@github.com:Sumantifr/ParticleNetLepton.git .

- cmsenv

- scram b -j10

- cd DeepLeptonNtuples/NTuplizer/test/

*Create voms proxy * 

- cmsRun PNet_MC_MINIAOD_cfg.py 

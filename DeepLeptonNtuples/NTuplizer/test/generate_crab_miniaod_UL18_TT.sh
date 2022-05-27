#!/bin/bash

config=PNet_MC_MINIAOD_cfg.py
publish=True
site=T2_AT_Vienna
DBS=global

sample_names=(
TT_Dilep TT_Semilep
)

sample_data=(
/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM
/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM
)
nsamples=${#sample_data[*]}
if [ $nsamples != ${#sample_names[*]} ]; 
then
	echo "No of names & samples are not same!! please check! (samples $nsamples names ${#sample_names[*]}"
	exit
fi

fil_list=crab_submit
mon_list=crab_monitor
truncate -s 0 ${fil_list}.sh
echo "#!/bin/bash" | cat >>${fil_list}.sh
truncate -s 0 ${mon_list}.sh

i=1
while [[ $i -le $nsamples ]]
do
	echo ${sample_data[i-1]} ${sample_names[i-1]}
	label=ParticleNetLepton_UL2018_${sample_names[i-1]}
	./crab_write.sh $label $config  ${sample_data[i-1]} $publish $site $DBS $username
	echo "crab submit -c crabfile_${label}.py" | cat >>${fil_list}.sh
	echo "crab status -d crab_${label}/crab_crab_${label}/" | cat >>${mon_list}.sh
	((i = i + 1))
done

chmod 744 ${fil_list}.sh
chmod 744 ${mon_list}.sh

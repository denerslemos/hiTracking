# Tracking code for cutbased studies 
(based on Xiao's code: https://github.com/Taburis/myProcesses/tree/hiTracking_12_2_X/hiTracking)

##CMSSW setup
```
cmsrel CMSSW_12_3_0
cd $CMSSW_BASE/src
cmsenv
git clone git@github.com:denerslemos/hiTracking.git
scram b -j 10
cd $CMSSW_BASE/src/hiTracking/trackingMVA_skimer/test/
```

Once compiled, you can run the codes for cutbased and mva only using:

```
cmsRun runTuple_cutbased.py
```
or
```
cmsRun runTuple_mvaonly.py
```

Crab configuration files are also available to run over full sample. To run, just use:
```
python crab_cutbased.py
```
or
```
python crab_mvaonly.py
```
please, remember to change the output storage location (```config.Site.storageSite```) to somewhere that your are allowed to write.

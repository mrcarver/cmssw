#!/bin/bash              
	#
	#PBS -r n 
	##Job settings                                                                                                                                                                       
	#PBS -N ROOT       
	#PBS -m a                                                                                                                                                                
	#PBS -o /cms/data/store/user/t2/users/mrcarver/Fall15AnalysisTuples/MC/WpWpJJ_EWK-QCD/HPCLogs_1.txt
	#PBS -e /cms/data/store/user/t2/users/mrcarver/Fall15AnalysisTuples/MC/WpWpJJ_EWK-QCD/HPCErrors_1.txt

	##Job Configuration                                                                                                                                                                  
	##Job Resources                                                                                                                                                                   
	#PBS -l walltime=00:07:00:00 
	#PBS -l nodes=1:ppn=1                                                                                                                                                                
	#PBS -l pmem=5200mb              
	# initialize environment for worker
	BASE=/lfs/scratch/mrcarver/Fall15AnalysisCode/CMSSW_7_4_12/src/SUSYAnalyzer/PatAnalyzer/test
	cd $BASE

	export SCRAM_ARCH=slc6_amd64_gcc491/
	export OSG_APP=/osg/app
	export VO_CMS_SW_DIR=${OSG_APP}/cmssoft/cms
	export CMS_PATH=${VO_CMS_SW_DIR}
	. /osg/app/cmssoft/cms/cmsset_default.sh;
	eval `scramv1 runtime -sh`


	cd $BASE
	export LD_LIBRARY_PATH="/osg/app/cmssoft/cms/osg/app/GLITE/GLITE_3_2_7-0/d-cache/dcap/lib:/lfs/scratch/mrcarver/Fall15AnalysisCode/CMSSW_7_4_12/biglib/slc6_amd64_gcc491:/lfs/scratch/mrcarver/Fall15AnalysisCode/CMSSW_7_4_12/lib/slc6_amd64_gcc491:/lfs/scratch/mrcarver/Fall15AnalysisCode/CMSSW_7_4_12/external/slc6_amd64_gcc491/lib:/cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/cmssw/CMSSW_7_4_12/biglib/slc6_amd64_gcc491:/cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/cmssw/CMSSW_7_4_12/lib/slc6_amd64_gcc491:/cvmfs/cms.cern.ch/slc6_amd64_gcc491/cms/cmssw/CMSSW_7_4_12/external/slc6_amd64_gcc491/lib:/cvmfs/cms.cern.ch/slc6_amd64_gcc491/external/llvm/3.6/lib:/cvmfs/cms.cern.ch/slc6_amd64_gcc491/external/gcc/4.9.1-cms/lib64:/cvmfs/cms.cern.ch/slc6_amd64_gcc491/external/gcc/4.9.1-cms/lib";

	# enter working area

cmsRun SSSync_Interactive.py input="file:/cms/data/store/mc/RunIISpring15DR74/WpWpJJ_EWK-QCD_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/50000/003F747A-D835-E511-BFC8-002590D9D976.root, file:/cms/data/store/mc/RunIISpring15DR74/WpWpJJ_EWK-QCD_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/50000/1C5348FA-142F-E511-A1B5-02163E011AE0.root, file:/cms/data/store/mc/RunIISpring15DR74/WpWpJJ_EWK-QCD_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/50000/54C90192-AC2E-E511-8E71-842B2B2AB616.root, file:/cms/data/store/mc/RunIISpring15DR74/WpWpJJ_EWK-QCD_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/50000/7459E177-AF2E-E511-8A53-F45214C748D2.root, file:/cms/data/store/mc/RunIISpring15DR74/WpWpJJ_EWK-QCD_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/50000/78253D89-FD34-E511-A83C-008CFA1113F4.root, file:/cms/data/store/mc/RunIISpring15DR74/WpWpJJ_EWK-QCD_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/50000/86BCDB16-A62E-E511-B644-0CC47A4DEE36.root, file:/cms/data/store/mc/RunIISpring15DR74/WpWpJJ_EWK-QCD_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/50000/9A54BF51-302F-E511-A41C-6C3BE5B5F0A0.root, file:/cms/data/store/mc/RunIISpring15DR74/WpWpJJ_EWK-QCD_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/50000/CE33AF75-AF2E-E511-AAF5-842B2B2B0D2E.root, file:/cms/data/store/mc/RunIISpring15DR74/WpWpJJ_EWK-QCD_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/50000/D0B751FF-B12E-E511-8B58-B8CA3A70A410.root, file:/cms/data/store/mc/RunIISpring15DR74/WpWpJJ_EWK-QCD_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/50000/DCE1B795-AC2E-E511-9CBF-000F530E47CC.root" output=/cms/data/store/user/t2/users/mrcarver/Fall15AnalysisTuples/MC/WpWpJJ_EWK-QCD/Job_1.root > /cms/data/store/user/t2/users/mrcarver/Fall15AnalysisTuples/MC/WpWpJJ_EWK-QCD/logs/Job_1.txt 2>  /cms/data/store/user/t2/users/mrcarver/Fall15AnalysisTuples/MC/WpWpJJ_EWK-QCD/errs/Job_1.txt


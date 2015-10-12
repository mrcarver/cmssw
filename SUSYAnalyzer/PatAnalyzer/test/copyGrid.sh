for i in `seq 2 1 21`;
do
xrdcp root://cms-xrd-global.cern.ch//store/user/cerati/TTWJets_Tune4C_13TeV-madgraph-tauola/Spring14dr-PU_S14_POSTLS170_V6AN1-miniAOD706p1/141013_122019/0000/miniAOD-prod_PAT_$i.root /scratch/lfs/lesya/TTW
done

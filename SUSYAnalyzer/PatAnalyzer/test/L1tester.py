import sys
import FWCore.ParameterSet.Config as cms
from RecoTauTag.RecoTau.PFRecoTauQualityCuts_cfi import PFTauQualityCuts


skim     ="NOSkim"
isMC     = bool(True)
isMSUGRA = bool(False)
isSMS    = bool(False)
useCrab  = bool(False)
doSusyTopProjection = bool(False)
inputFile=""
outputFile=""
#outputFile="SyncOutput/TriggerNamesTest.root"
def getVal(arg):
    i=0
    while i < len(arg) and arg[i] != "=": i+=1
    return arg[i+1:]

## loop over arguments
for i in range(1,len(sys.argv)):
    print "[arg "+str(i)+"] : ", sys.argv[i] 
    if "isMC" in sys.argv[i] :
        isMC=False if getVal(sys.argv[i]) == "False" else True
    elif "isMSUGRA" in sys.argv[i]:
        isMSUGRA=True
    elif "isSMS" in sys.argv[i]:
        isSMS=True
    elif "skim" in sys.argv[i]:
        skim=getVal(sys.argv[i])
    elif "output" in sys.argv[i]:
        outputFile=getVal(sys.argv[i])
    elif "input" in sys.argv[i] :
        inputFile=getVal(sys.argv[i])

        
if skim=="" : print "WARNING: No Skim Conditions have been provided \n"



#process = cms.Process("FakeElectrons")
process = cms.Process("pippo")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.MessageLogger.cerr.threshold = 'WARNING' # Options: INFO, WARNING, ERROR
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

if isMC:
    #process.GlobalTag.globaltag="GR_P_V40_AN1::All" # tag fo 53X 2012C Prompt Reco (FT_53_V10_AN1) in the future for some runs
    #process.GlobalTag.globaltag="FT_53_V6_AN1::All" # tag fo 53X 2012A/B 13Jul2012 rereco
	#process.GlobalTag.globaltag="PHYS14_25_V2::All" # tag fo PHYS14 25ns samples
	process.GlobalTag.globaltag="MCRUN2_74_V9::All" # tag for Run2 MC samples
else:
    #process.GlobalTag.globaltag="GR_P_V40_AN1::All" # tag fo 53X 2012C Prompt Reco (FT_53_V10_AN1) in the future for some runs
    process.GlobalTag.globaltag="FT_53_V6_AN1::All" # tag fo 53X 2012A/B 13Jul2012 rereco

process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.Reconstruction_cff')

process.source = cms.Source ("PoolSource",
                             # Disable duplicate event check mode because the run and event -numbers
                             # are incorrect in current Madgraph samples (Dec 16, 2008)
                             duplicateCheckMode = cms.untracked.string('noDuplicateCheck'), 
							 #eventsToProcess = cms.untracked.VEventRange('1:1557092937','1:1624289770','1:1137407433','1:1658883658'),                   
                             fileNames = cms.untracked.vstring(),      
                             )

#from  Data.ElectronHad_Run2012B_SSDiLepSkim_193752_194076 import *

if isMC:
    process.source.fileNames = cms.untracked.vstring(   
	
	
    #'/cms/data/store/mc/RunIISpring15DR74/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/00000/022B08C4-C702-E511-9995-D4856459AC30.root',
	 '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/029ED365-C209-E511-82BE-A0040420FE80.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/06FC4365-D209-E511-B648-0002C90A33FC.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/0AE236D5-E909-E511-AF51-0025901AEBD4.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/0C6566C1-C109-E511-9622-00259073E37A.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/0CA3B41E-E209-E511-B566-002590FD5A48.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/10A79AC4-C109-E511-9641-00259073E474.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/12369DB6-CD09-E511-9242-20CF3027A626.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/14A14638-C709-E511-9E91-0002C94CDA06.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/14ED3293-F209-E511-8B6E-002590FD5A48.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/162DE916-8609-E511-ADCC-00259073E4A2.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/16361F6A-CB09-E511-9DBC-00259073E3AE.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/16EE83C0-C109-E511-80B8-00259073E4FC.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/183E9A4D-DB09-E511-B078-002590FD5A78.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/22A60FA6-C309-E511-94AA-A0040420FE80.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/2A7D5B9B-CB09-E511-AE12-00259073E4A2.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/2C77B403-C309-E511-926A-A0040420FE80.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/2CBA6F3A-CA09-E511-8F5E-00259073E3AE.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/302C8971-D609-E511-AF1D-0CC47A0AD456.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/32DAACA6-C309-E511-BB5A-A0040420FE80.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/34738A1B-8609-E511-9FC7-00259074AE3C.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/362B3171-D409-E511-A0A8-00259029E66E.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/3678C713-8609-E511-9069-20CF305B050D.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/3697FDDF-8209-E511-8148-00259074AEAE.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/36E2F16A-D609-E511-B032-A0369F310120.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/3CEF6E12-E009-E511-BC0A-0025901F8740.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/469B8FDF-8209-E511-B760-00259073E3DA.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/48BD2B67-C209-E511-AEAE-A0040420FE80.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/4A851EA6-C609-E511-9E70-0002C90B73F0.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/508A6D30-DD09-E511-BCFB-0025901AEDA6.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/54301F1A-9809-E511-B3BC-00259073E4AC.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/54837A5F-6C0A-E511-964A-A0369F3102B6.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/566793CA-C109-E511-A3B1-0002C94D54D4.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/584F1920-6A0A-E511-91F4-00259073E4F6.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/5EFCAD49-DB09-E511-B874-002590FD5A78.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/6057D2C3-C109-E511-85DF-00074305CB73.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/62C103C7-C109-E511-B96B-A0040420FE80.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/66A2581C-DD09-E511-BA79-002590D9D956.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/66BC711E-DA09-E511-8E78-002590D9D990.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/6A40D3C0-C109-E511-94A8-002590574A44.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/6ACEE5CF-C109-E511-89EE-0002C94D5516.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/6CCBE516-8609-E511-9F51-00259073E506.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/6CF1B29D-E309-E511-840C-002590D9D8BC.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/6E04E317-0B0A-E511-B90B-0025907D2000.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/705EF2E8-8309-E511-9692-00259074AE40.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/70DF0D8B-C209-E511-AFBA-A0040420FE80.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/74A924A1-C209-E511-BCFE-0002C94CD150.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/76CCC6B1-C209-E511-BB37-0002C90A3698.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/78420AD4-C109-E511-83BB-0002C952E766.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/7A949CD8-9A09-E511-9BA3-00259073E4A2.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/7C5F3130-C309-E511-BBD5-A0040420FE80.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/7E7B1EC7-C109-E511-B980-A0040420FE80.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/824A87EB-D409-E511-BF43-0025901AEBD6.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/82C1D125-D409-E511-9A60-0CC47A0AD456.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/84F8A97B-9A09-E511-92AE-00259073E3E6.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/90EFED04-D909-E511-8DE6-0CC47A0AD63E.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/9231C502-6A0A-E511-AC25-002590D9D880.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/9462A307-C309-E511-A13B-00259073E37A.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/949BD125-D409-E511-96CC-0CC47A0AD456.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/989B9AC5-C109-E511-B9D5-0002C94DBB18.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/9C423BA6-C209-E511-9D63-0002C90A3452.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/9CCDBFEB-DB09-E511-BB67-002590D9D956.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/9CDE3535-CA09-E511-87D0-00259073E4A2.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/9E1EF167-080A-E511-B9FE-0025907D2000.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/9EF0D125-D409-E511-8DE4-0CC47A0AD456.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/A2FE8A19-9809-E511-BE52-002590747DDC.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/A4A315C8-C109-E511-829F-0002C94D5648.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/A644E8C5-C109-E511-9A47-00259073E4A2.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/AEF2D568-9909-E511-ADBE-00259073E3E6.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/B0ABB569-9909-E511-93DC-00259073E4CA.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/B2C2E618-C309-E511-BE1E-00259073E474.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/B4F7C647-D709-E511-9A98-A0369F310120.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/B85A7F2C-DD09-E511-ACE8-0025907D24E6.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/B8709B00-C309-E511-9F0C-00259073E53E.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/B899CDD3-C109-E511-9897-0002C952E766.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/B8CF3C39-DA09-E511-BD2F-0CC47A0AD6FC.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/BADED125-D409-E511-8C05-0CC47A0AD456.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/C03C8971-D609-E511-B73A-0CC47A0AD456.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/C0BB80CE-C109-E511-A07C-0002C94CD120.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/C0F63B50-E009-E511-83D3-0CC47A0AD456.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/C64695C4-C109-E511-9C01-00074305CB0D.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/C86F16C9-D409-E511-B117-0025901AEBDA.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/C8AA7E65-080A-E511-8AAC-00259019A418.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/C8D4F99C-C209-E511-BAFF-A0040420FE80.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/CE635BC2-DD09-E511-A27C-002590D9D966.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/D0CE2165-9909-E511-96C8-00259073E4AC.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/D0CFFC29-DA09-E511-B8BF-002590FD5A78.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/D866009A-9A09-E511-8A90-00259073E4CA.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/D86FD225-D409-E511-8532-0CC47A0AD456.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/E07D54DF-8209-E511-BE22-00259073E3DA.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/E24365AA-D909-E511-A56B-0002C94D552A.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/EA5BEC5F-D209-E511-BD78-A0040420FE80.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/EC382D26-E209-E511-91AB-002590D9D9F6.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/F027AFEE-DB09-E511-B34B-002590FD5122.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/F4400206-C309-E511-9DCB-A0040420FE80.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/F4BB3CEA-8309-E511-91F1-00259073E474.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/F4C68971-D609-E511-AC71-0CC47A0AD456.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/F69AFCB9-9A09-E511-BCA6-00259073E4AC.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/F6E794FF-D809-E511-A4EB-002590FD5A78.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/F6ECE6CF-E809-E511-A59A-002590D9D880.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/F8D49AC2-C109-E511-AE1A-A0040420FE80.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/F8F712FC-D309-E511-BD41-002590FD5A3A.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/FA607DBB-C109-E511-A798-00259073E53E.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/FC91CDD3-C109-E511-AAE9-0002C952E766.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/FCBD1EE1-8209-E511-86DC-00259074AE54.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/50000/FE4687EB-D409-E511-8A22-0025901AEBD6.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/026041C5-A609-E511-A0F5-0025907750A0.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/08A50251-9609-E511-B874-00259073E3C0.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/08FEF2FF-9809-E511-B75B-A0040420FE80.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/0E5C5B36-BE09-E511-9B95-0CC47A0AD668.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/12C1FF61-B709-E511-869F-0002C90F8030.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/12C5BAC2-7909-E511-9856-A0040420FE80.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/1658BAF0-7809-E511-BA96-0002C92DB4CC.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/16BEE98B-9309-E511-B76D-00259073E3C0.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/16FB765D-A709-E511-8844-A0040420FE80.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/1E58896C-9909-E511-A9EA-002590747DDC.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/26D0C154-7B09-E511-9B16-00074305CB2A.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/2AAA0188-9309-E511-938A-20CF305B053E.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/2AE26660-9609-E511-8A84-00259073E4AC.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/2C6F672A-9709-E511-9211-20CF305B053E.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/2CF67100-9909-E511-9AFC-00259073E506.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/2E47A16B-5B09-E511-ABB3-00259073E344.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/30CE84D2-7309-E511-B67F-A0040420FE80.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/36CD8801-AF09-E511-A329-002590D9D968.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/3AFBAFC9-7909-E511-A374-00074305CBAE.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/3C2A5F78-9A09-E511-A894-00259073E500.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/3C57D7D0-A609-E511-8E62-0025902BD8CE.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/3CD838C6-A609-E511-B532-00259073E4A2.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/3CFECD46-C009-E511-BA2D-00259029E7FC.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/3E72A3F5-9809-E511-80ED-A0040420FE80.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/421BDDA0-9809-E511-93FB-00259073E506.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/446EC59A-C009-E511-8B05-003048947BAA.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/4C69C9C6-8D09-E511-990D-0002C94CDAD2.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/4C6D6671-9709-E511-8929-00259073E506.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/4E0B1B88-8B09-E511-8089-00259073E3C0.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/4EBD528E-AA09-E511-9111-00259074AEC2.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/502D0D88-9309-E511-8570-20CF3027A5C9.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/540026F6-B709-E511-9DAC-0CC47A0AD742.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/56C6138F-AA09-E511-B47E-00259073E4AC.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/5A85985D-9909-E511-863F-00259073E500.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/5C2AE169-A709-E511-9F89-A0040420FE80.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/624EBBC5-A609-E511-A743-00074305D01A.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/6482CDC3-A609-E511-9668-00259074AE3E.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/6489658E-A609-E511-BF2E-002590D9D8A4.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/666B1A8B-9309-E511-9C6C-00259074AE8A.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/6C9BA4B2-7909-E511-81E5-A0369F3102B6.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/6E740F01-B809-E511-B266-002590D9D8BE.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/726B3A21-B809-E511-9E87-003048947CC6.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/74A3E9A9-6E09-E511-A2B9-20CF305B055B.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/74A935CF-A609-E511-AE8E-0025907DE266.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/764635FC-BC09-E511-B9F8-00259075D6B4.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/78506D5F-AC09-E511-9648-A0040420FE80.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/7A2ACBC1-A609-E511-AF24-A0040420FE80.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/7CD94CB2-7909-E511-B4AB-A0369F3016EC.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/7EDDCD89-9309-E511-9E17-00259073E3DA.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/80865F70-A709-E511-AB46-0025901AEBD8.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/8687B64F-8309-E511-BF20-B8CA3A709648.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/887698CC-A609-E511-9343-002590D9D990.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/8A7EE4C6-A609-E511-8F13-00074305CEEA.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/8C0BD265-9A09-E511-A49A-002590747DDC.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/8CE8BE89-9309-E511-A870-002590747DDC.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/9275A7C3-8D09-E511-A2AB-0002C90F8036.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/9624C097-C009-E511-840A-0025901AC3CE.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/9AACEE68-7B09-E511-9BB6-0002C92DB46C.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/9ED50E9B-A609-E511-9573-00259075D62E.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/A042B441-BE09-E511-BF8D-00259029E84C.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/A0E36CEB-8D09-E511-9A64-0002C94CD2C6.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/A82487D0-7C09-E511-8795-0002C92DB44E.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/AA1F6D85-7B09-E511-B14F-0002C92DB44E.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/B01ECC8C-A609-E511-8991-002590D9D9F6.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/B03A431A-B809-E511-A72C-002590D9D980.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/B8056A84-6C0A-E511-9AD9-0CC47A0AD6FC.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/B87EC4C2-A609-E511-BF0A-A0040420FE80.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/BA2B73EC-A609-E511-BFCE-0CC47A0AD6DA.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/BAB1E75D-6C0A-E511-A338-00259073E4D6.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/BC91441D-8309-E511-89A1-A0040420FE80.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/BCAC35B8-7C09-E511-BECB-A0369F310120.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/BE42C2F6-7809-E511-B16E-0002C90F80DA.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/BE86AA39-6909-E511-ADF4-20CF3027A629.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/BEEB7E68-5B09-E511-9202-00259074AE28.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/C07A7DB5-7909-E511-89F0-F45214C748C2.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/C2553D69-5B09-E511-8DA4-00259074AE44.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/C2AD3C08-5B09-E511-9EB9-00259073E4E6.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/C2B1FFC2-A609-E511-A1B4-00259074AEC2.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/C41DD28B-9309-E511-9323-20CF305B057E.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/C6D01749-6C0A-E511-885D-F45214C7B6BA.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/C8D7CC39-6909-E511-B285-20CF300E9EC3.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/CC7642E8-A609-E511-BD53-002590D9D8AE.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/CCDAAF4F-9609-E511-BAC6-20CF305B053E.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/CE9F856D-9509-E511-A795-00259073E4FC.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/D071B59D-7F09-E511-BCE4-A0369F3102B6.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/D2C2E2A9-7F09-E511-9BCA-F45214C748D0.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/D45D5813-9009-E511-BC68-0002C94D561C.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/D878A3F5-9809-E511-B699-A0040420FE80.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/DA1D508C-6D09-E511-BFDE-00259073E466.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/DA3A216A-5B09-E511-8090-00259074AEDE.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/DC23BD87-9309-E511-BE83-20CF3027A626.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/DEF39B36-C009-E511-B043-00304867FD3F.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/E26C922F-6C0A-E511-9046-0025907E33EA.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/E4143B5A-6C0A-E511-BCFD-002590FD5A48.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/E86AB31A-C109-E511-99D5-0025901ABB72.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/E8B40721-B809-E511-9199-0CC47A0AD498.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/EA6D5489-9309-E511-A0AD-20CF3027A5CE.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/EE49AD65-5B09-E511-97EC-20CF305B04F5.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/EEFB2923-6009-E511-92A1-00259073E3C0.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/F06CF22A-B809-E511-969D-0025901AEBD6.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/F0F21A13-9009-E511-8726-00074305CAFD.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/F4F23E04-B809-E511-B358-002590FD5A4C.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/F6C8C383-8B09-E511-AA3A-20CF305B057E.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/F8C57C9A-9909-E511-81D1-A0040420FE80.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/FA86E357-7B09-E511-AB65-B8CA3A709648.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/FC24D76C-AF09-E511-A324-00259074AE54.root',
       '/store/mc/RunIISpring15DR74/SingleNeutrino/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v2/60000/FE31D6A9-6E09-E511-88A8-20CF305616EA.root'
	 
	 
	)

else:
   import PhysicsTools.PythonAnalysis.LumiList as LumiList
#   process.source.lumisToProcess = LumiList.LumiList(filename = '/scratch/osg/lesya/CMSSW_5_3_11_patch6/src/SUSYAnalyzer/PatAnalyzer/test/fullJSON_SUSYFakes.json').getVLuminosityBlockRange()

   process.source.fileNames = cms.untracked.vstring( 
#    #'file:../../../CMSSW_5_2_3_PatTuple_test_PFiso_El_newMETCor.root'
#    'file:/raid/raid10/ddidar/CMSSW_5_2_3/src/CMSSW_5_2_3_PatTuple.root'
#    '/store/user/ddidar/DATA_52X/DoubleElectron_Run2012B-PromptReco-v1_AOD/TwoLepSkim/193754_196509/PatTuple_DoubleElectron_Run2012B-PromptReco-v1_AOD_TwoLepSkim_193754_196509_1813_1_Il0.root'
#     '/store/user/ddidar/DATA_52X/ElectronHad_Run2012B-PromptReco-v1_AOD/TwoLepSkim/193754_196509/PatTuple_ElectronHad_Run2012B-PromptReco-v1_AOD_TwoLepSkim_193754_196509_999_1_PPc.root'
#         '/cms/data/store/user/ddidar/DATA_53X_DCSONLY/MuHad_Run2012D-PromptReco-v1_AOD/OneLepOneBJetSkim/203768_207925/PatTuple_MuHad_Run2012D-PromptReco-v1_AOD_OneLepOneBJetSkim_203768_207925_214_1_YFH.root'
#         '/cms/data/store/user/ddidar/DATA_53X_DCSONLY/ElectronHad_Run2012A-13Jul2012-v1_AOD/OneLepOneBJetSkim/190456_193621/PatTuple_ElectronHad_Run2012A-13Jul2012-v1_AOD_OneLepOneBJetSkim_190456_193621_110_1_24w.root'
#         '/cms/data/store/user/ddidar/DATA_53X/ElectronHad_Run2012A-13Jul2012-v1_AOD/OneLepOneBJetSkim/190456_193754/PatTuple_ElectronHad_Run2012A-13Jul2012-v1_AOD_OneLepOneBJetSkim_190456_193754_101_1_aAA.root',
#        '/cms/data/store/user/ddidar/DATA_53X_DCSONLY/MuHad_Run2012B-13Jul2012-v1_AOD/OneLepOneBJetSkim/193833_196531/PatTuple_MuHad_Run2012B-13Jul2012-v1_AOD_OneLepOneBJetSkim_193833_196531_451_1_huk.root'
#         '/cms/data/store/user/ddidar/DATA_53X/MuHad_Run2012A-13Jul2012-v1_AOD/OneLepOneBJetSkim/190456_193754/PatTuple_MuHad_Run2012A-13Jul2012-v1_AOD_OneLepOneBJetSkim_190456_193754_332_1_G9q.root' 
#    '/scratch/lfs/lesya/FakeSync_jets_test_data_3rdShit.root'
#          '/cms/data/store/user/ddidar/DATA_53X_DCSONLY/DoubleMu_Run2012D-PromptReco-v1_AOD/NOSkim/203767_207926/PatTuple_DoubleMu_Run2012D-PromptReco-v1_AOD_NOSkim_203767_207926_85_1_CCE.root'
    #'/scratch/lfs/lesya/FakeSync_jets_test_data_Electrons.root'
    )
#    runOnData(process)
   

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(2000000))

#process.TFileService = cms.Service("TFileService", fileName = cms.string(outputFile) )
process.TFileService = cms.Service("TFileService", fileName = cms.string("L1ntuplizer.root") )

## input files
#for i in inputFile.split(","):
#   print "Adding: ", i
#   process.source.fileNames.append(i)

process.SyncExercise = cms.EDAnalyzer("L1RateNtuplizer",#"FakeMuonsYa",
                                       MuonLabel = cms.InputTag("slimmedMuons"),
                                       ElectronLabel = cms.InputTag("slimmedElectrons"),
                                       TauLabel = cms.InputTag("slimmedTaus"),
                                       #TauDiscriminatorLabel = cms.InputTag("recoPFTauDiscriminator"),
                                       JetLabel = cms.InputTag("slimmedJets"),
                                       BeamSpotLabel = cms.InputTag("offlineBeamSpot"),
                                       HLTResultsLabel = cms.InputTag("TriggerResults"),
									   MVAId = cms.InputTag("mvaTrigV050nsCSA14","","addMVAid"), #not used when running on PAT
                                       METLabel = cms.InputTag("slimmedMETs"),
                                       #METFilter = cms.InputTag("TriggerResults::PAT"),
                                       qualityCuts = PFTauQualityCuts,
                                       SampleLabel = cms.untracked.string("ElectronsMC") # defines a piece of code to run; helps to avoid code recompilation
                                       )

#FakeElectrons
#process.FakeElectrons = cms.EDAnalyzer("SSSync",#"FakeMuonsYa",
#                                       MuonLabel = cms.InputTag("slimmedMuons"),
#                                       ElectronLabel = cms.InputTag("slimmedElectrons"),
#                                       TauLabel = cms.InputTag("slimmedTaus"),
#                                       #TauDiscriminatorLabel = cms.InputTag("recoPFTauDiscriminator"),
#                                       JetLabel = cms.InputTag("slimmedJets"),
#                                      BeamSpotLabel = cms.InputTag("offlineBeamSpot"),
#                                      HLTResultsLabel = cms.InputTag("TriggerResults"),
#                                      METLabel = cms.InputTag("slimmedMETs"),
#                                      #METFilter = cms.InputTag("TriggerResults::PAT"),
#                                      qualityCuts = PFTauQualityCuts,
#                                       SampleLabel = cms.untracked.string("ElectronsMC") # defines a piece of code to run; helps to avoid code recompilation
#                                       )

#process.tauMCMatch = cms.EDProducer("MCTruthDeltaRMatcherNew",
#                                     src = cms.InputTag("slimmedTaus"),
#                                     matched = cms.InputTag("prunedGenParticles"),
#                                     distMin = cms.double(0.15),
#                                     matchPDGId = cms.vint32()
#                                     )

#process.tauMCMatch = cms.EDProducer("MCMatcher",
#    				src     = cms.InputTag("slimmedTaus"),
#    				matched = cms.InputTag("prunedGenParticles"),
#                   mcPdgId     = cms.vint32(),
#                   checkCharge = cms.bool(False),
#                   mcStatus = cms.vint32(),
#                   maxDeltaR = cms.double(0.15),
#                   maxDPtRel = cms.double(3.0),
#                   resolveAmbiguities = cms.bool(True),
#                   resolveByMatchQuality = cms.bool(False),
#)


#
# Example for a configuration of the MC match
# for taus (cuts are NOT tuned)
# (using old values from TQAF, january 2008)
#
process.tauMatch = cms.EDProducer("MCMatcher",
                          src         = cms.InputTag("slimmedTaus"),          # RECO objects to match
                          matched     = cms.InputTag("prunedGenParticles"),              # mc-truth particle collection
                          mcPdgId     = cms.vint32(15),                            # one or more PDG ID (15 = tau); absolute values (see below)
                          checkCharge = cms.bool(True),                            # True = require RECO and MC objects to have the same charge
                          mcStatus    = cms.vint32(2),                             # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
                          # NOTE that Taus can only be status 3 or 2, never 1!
                          maxDeltaR   = cms.double(999.9),                         # Minimum deltaR for the match.     By default any deltaR is allowed (why??)
                          maxDPtRel   = cms.double(999.9),                         # Minimum deltaPt/Pt for the match. By default anything is allowed   ( ""  )
                          resolveAmbiguities    = cms.bool(True),                  # Forbid two RECO objects to match to the same GEN object
                          resolveByMatchQuality = cms.bool(False),                 # False = just match input in order; True = pick lowest deltaR pair first
                          )

process.tauGenJetMatch = cms.EDProducer("GenJetMatcher",             # cut on deltaR, deltaPt/Pt; pick best by deltaR
                                src         = cms.InputTag("slimmedTaus"),          # RECO jets (any View<Jet> is ok)
                                matched     = cms.InputTag("tauGenJets"),  # GEN jets  (must be GenJetCollection) "tauGenJetsSelectorAllHadrons"
                                mcPdgId     = cms.vint32(),                              # n/a
                                mcStatus    = cms.vint32(),                              # n/a
                                checkCharge = cms.bool(False),                           # n/a
                                maxDeltaR   = cms.double(0.1),                           # Minimum deltaR for the match
                                maxDPtRel   = cms.double(3.0),                           # Minimum deltaPt/Pt for the match
                                resolveAmbiguities    = cms.bool(True),                  # Forbid two RECO objects to match to the same GEN object
                                resolveByMatchQuality = cms.bool(False),                 # False = just match input in order; True = pick lowest deltaR pair first
                                )
process.tauGenJets = cms.EDProducer(
                            "TauGenJetProducer",
                            prunedGenParticles =  cms.InputTag('prunedGenParticles'),
                            includeNeutrinos = cms.bool( False ),
                            verbose = cms.untracked.bool( False )
                            )

process.patMCTruth_Tau =  cms.Sequence ( process.tauMatch+
                                process.tauGenJets*
                                process.tauGenJetMatch )


process.electronMatch = cms.EDProducer("MCMatcher",
                                  src         = cms.InputTag("slimmedElectrons"),          # RECO objects to match
                                  matched     = cms.InputTag("prunedGenParticles"),              # mc-truth particle collection
                                  mcPdgId     = cms.vint32(),                            # one or more PDG ID (15 = tau); absolute values (see below)
                                  checkCharge = cms.bool(False),                            # True = require RECO and MC objects to have the same charge
                                  mcStatus    = cms.vint32(),                             # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
                                  # NOTE that Taus can only be status 3 or 2, never 1!
                                  maxDeltaR   = cms.double(0.25),                         # Minimum deltaR for the match.     By default any deltaR is allowed (why??)
                                  maxDPtRel   = cms.double(10.),                         # Minimum deltaPt/Pt for the match. By default anything is allowed   ( ""  )
                                  resolveAmbiguities    = cms.bool(False),                  # Forbid two RECO objects to match to the same GEN object
                                  resolveByMatchQuality = cms.bool(True),                 # False = just match input in order; True = pick lowest deltaR pair first
                                  )

process.muonMatch = cms.EDProducer("MCMatcher",
                                   src         = cms.InputTag("slimmedMuons"),              # RECO objects to match
                                   matched     = cms.InputTag("prunedGenParticles"),          # mc-truth particle collection
                                   mcPdgId     = cms.vint32(),                          # one or more PDG ID (15 = tau); absolute values (see below)
                                   checkCharge = cms.bool(False),                       # True = require RECO and MC objects to have the same charge
                                   mcStatus    = cms.vint32(),                          # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
                                   # NOTE that Taus can only be status 3 or 2, never 1!
                                   maxDeltaR   = cms.double(0.25),                         # Minimum deltaR for the match.     By default any deltaR is allowed
                                   maxDPtRel   = cms.double(10.),                         # Minimum deltaPt/Pt for the match. By default anything is allowed   ( ""  )
                                   resolveAmbiguities    = cms.bool(False),                  # Forbid two RECO objects to match to the same GEN object
                                   resolveByMatchQuality = cms.bool(True),                 # False = just match input in order; True = pick lowest deltaR pair first
                                   )

## reco-generator(parton) matching for jets
process.jetPartonMatch = cms.EDProducer("MCMatcher",      # cut on deltaR, deltaPt/Pt; pick best by deltaR
                                        src = cms.InputTag("slimmedJets"),      # RECO objects to match
                                        matched = cms.InputTag("prunedGenParticles"),     # mc-truth particle collection
                                        mcPdgId  = cms.vint32(1, 2, 3, 4, 5, -1, -2, -3, -4, -5, 21),# one or more PDG ID (quarks except top; gluons)
                                        mcStatus = cms.vint32(3),                   # PYTHIA status code (3 = hard scattering)
                                        checkCharge = cms.bool(False),              # False = any value of the charge of MC and RECO is ok
                                        maxDeltaR = cms.double(0.5),                # Minimum deltaR for the match
                                        maxDPtRel = cms.double(30.0),                # Minimum deltaPt/Pt for the match
                                        resolveAmbiguities = cms.bool(False),        # Forbid two RECO objects to match to the same GEN object
                                        resolveByMatchQuality = cms.bool(True),     # False = just match input in order; True = pick lowest deltaR pair first
                                        )

process.jetPartonMatch2 = cms.EDProducer("MCMatcher",      # cut on deltaR, deltaPt/Pt; pick best by deltaR
                                        src = cms.InputTag("slimmedJets"),      # RECO objects to match
                                        matched = cms.InputTag("prunedGenParticles"),     # mc-truth particle collection
                                        mcPdgId  = cms.vint32(1, 2, 3, 4, 5, -1, -2, -3, -4, -5, 21),# one or more PDG ID (quarks except top; gluons)
                                        mcStatus = cms.vint32(2),                   # PYTHIA status code (3 = hard scattering)
                                        checkCharge = cms.bool(False),              # False = any value of the charge of MC and RECO is ok
                                        maxDeltaR = cms.double(0.5),                # Minimum deltaR for the match
                                        maxDPtRel = cms.double(3.0),                # Minimum deltaPt/Pt for the match
                                        resolveAmbiguities = cms.bool(False),        # Forbid two RECO objects to match to the same GEN object
                                        resolveByMatchQuality = cms.bool(True),     # False = just match input in order; True = pick lowest deltaR pair first
                                        )

process.jetPartonMatchJets = cms.EDProducer("MCMatcher",      # cut on deltaR, deltaPt/Pt; pick best by deltaR
                                        src = cms.InputTag("slimmedJets"),      # RECO objects to match
                                        matched = cms.InputTag("ak5GenJets"),     # mc-truth particle collection
                                        mcPdgId  = cms.vint32(),# one or more PDG ID (quarks except top; gluons)
                                        mcStatus = cms.vint32(),                   # PYTHIA status code (3 = hard scattering)
                                        checkCharge = cms.bool(False),              # False = any value of the charge of MC and RECO is ok
                                        maxDeltaR = cms.double(0.5),                # Minimum deltaR for the match
                                        maxDPtRel = cms.double(3.0),                # Minimum deltaPt/Pt for the match
                                        resolveAmbiguities = cms.bool(True),        # Forbid two RECO objects to match to the same GEN object
                                        resolveByMatchQuality = cms.bool(True),     # False = just match input in order; True = pick lowest deltaR pair first
                                        )

process.electronPartonMatch = cms.EDProducer("MCMatcher",      # cut on deltaR, deltaPt/Pt; pick best by deltaR
                                             src = cms.InputTag("slimmedElectrons"),      # RECO objects to match
                                             matched = cms.InputTag("prunedGenParticles"),     # mc-truth particle collection
                                             mcPdgId  = cms.vint32(1, 2, 3, 4, 5, 21),# one or more PDG ID (quarks except top; gluons)
                                             mcStatus = cms.vint32(3),                   # PYTHIA status code (3 = hard scattering)
                                             checkCharge = cms.bool(False),              # False = any value of the charge of MC and RECO is ok
                                             maxDeltaR = cms.double(0.5),                # Minimum deltaR for the match
                                             maxDPtRel = cms.double(30.0),                # Minimum deltaPt/Pt for the match
                                             resolveAmbiguities = cms.bool(True),        # Forbid two RECO objects to match to the same GEN object
                                             resolveByMatchQuality = cms.bool(True),     # False = just match input in order; True = pick lowest deltaR pair first
                                             )

process.muonPartonMatch = cms.EDProducer("MCMatcher",      # cut on deltaR, deltaPt/Pt; pick best by deltaR
                                         src = cms.InputTag("slimmedMuons"),      # RECO objects to match
                                         matched = cms.InputTag("prunedGenParticles"),     # mc-truth particle collection
                                         mcPdgId  = cms.vint32(1, 2, 3, 4, 5, 21),# one or more PDG ID (quarks except top; gluons)
                                         mcStatus = cms.vint32(3),                   # PYTHIA status code (3 = hard scattering)
                                         checkCharge = cms.bool(False),              # False = any value of the charge of MC and RECO is ok
                                         maxDeltaR = cms.double(0.5),                # Minimum deltaR for the match
                                         maxDPtRel = cms.double(30.0),                # Minimum deltaPt/Pt for the match
                                         resolveAmbiguities = cms.bool(True),        # Forbid two RECO objects to match to the same GEN object
                                         resolveByMatchQuality = cms.bool(True),     # False = just match input in order; True = pick lowest deltaR pair first
                                         )

process.goodOfflinePrimaryVertices = cms.EDFilter( "PrimaryVertexObjectFilter" # checks for fake PVs automatically
                                                  , filterParams =cms.PSet(
                                                                           minNdof = cms.double( 4. )
                                                                           , maxZ    = cms.double( 24. )
                                                                           , maxRho  = cms.double( 2. ) )		
                                                  , filter       = cms.bool( False ) # use only as producer
                                                  , src          = cms.InputTag( 'offlineSlimmedPrimaryVertices' )
)

if isMC:
    process.p = cms.Path(
	process.goodOfflinePrimaryVertices
#   *process.tauMatch
#   *process.muonMatch
#   *process.electronMatch
#   *process.jetPartonMatch
#   *process.jetPartonMatch2
#   *process.electronPartonMatch
#   *process.muonPartonMatch
 #   *process.FakeElectrons
 	#*process.mvaNonTrigV025nsPHYS14
	*process.SyncExercise
        )
else:
   import PhysicsTools.PythonAnalysis.LumiList as LumiList
#   process.source.lumisToProcess = LumiList.LumiList(filename = '/scratch/hpc/lesya/CMSSW_5_3_6_patch1/src/SUSYAnalyzer/PatAnalyzer/test/JSON/Merged_190456-208686_8TeV_Collision12.json').getVLuminosityBlockRange()
   process.source.lumisToProcess = LumiList.LumiList(filename = '/scratch/osg/lesya/CMSSW_5_3_11_patch6/src/SUSYAnalyzer/PatAnalyzer/test/fullJSON_SUSYFakes.json').getVLuminosityBlockRange() 
   process.p = cms.Path(
	process.goodOfflinePrimaryVertices
       *process.FakeElectrons
   )
    


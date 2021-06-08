#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <sstream>
#include <map>
#include <algorithm>

#include "TFile.h" 
#include "TH1.h" 
#include "TH2.h"
#include "TGraph.h"
#include "TTree.h"
#include "TROOT.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TChain.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TError.h"
#include "TLorentzVector.h"
#include "TRandom.h"
#include "TSystem.h"

#include "TVector3.h"
#include "TMatrix.h"

#include "RooRealVar.h"
#include "RooWorkspace.h"

#include "DesyTauAnalyses/NTupleMaker/interface/Config.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AC1B.h"

#include "DesyTauAnalyses/NTupleMaker/interface/Synch17Tree.h"
#include "DesyTauAnalyses/NTupleMaker/interface/Synch17GenTree.h"
#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h"
#include "DesyTauAnalyses/NTupleMaker/interface/functions.h"
#include "DesyTauAnalyses/NTupleMaker/interface/leptau_jets_WIP.h"
#include "HTT-utilities/LepEffInterface/interface/ScaleFactor.h"
#include "DesyTauAnalyses/NTupleMaker/interface/PileUp.h"
#include "HTT-utilities/RecoilCorrections_KIT/interface/RecoilCorrector.h"
#include "DesyTauAnalyses/NTupleMaker/interface/functionsSynch2017.h"
#include "HiggsCPinTauDecays/IpCorrection/interface/IpCorrection.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>

#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "DesyTauAnalyses/NTupleMaker/interface/Systematics_WIP.h"
#include "DesyTauAnalyses/NTupleMaker/interface/LeptonScaleSys_WIP.h"
#include "DesyTauAnalyses/NTupleMaker/interface/ZPtWeightSys_WIP.h"
#include "DesyTauAnalyses/NTupleMaker/interface/TopPtWeightSys_WIP.h"
#include "DesyTauAnalyses/NTupleMaker/interface/JetEnergyScaleSys_WIP.h"
#include "DesyTauAnalyses/NTupleMaker/interface/PuppiMETSys_WIP.h"
#include "DesyTauAnalyses/NTupleMaker/interface/PFMETSys_WIP.h"
//#include "DesyTauAnalyses/NTupleMaker/interface/JESUncertainties.h"
#include "DesyTauAnalyses/NTupleMaker/interface/LepTauFakeRate.h"
#include "DesyTauAnalyses/NTupleMaker/interface/functionsCP.h"
#include "DesyTauAnalyses/NTupleMaker/interface/functionsCP_tt.h"
//#include "HTT-utilities/TauTriggerSFs2017/interface/TauTriggerSFs2017.h"
#include "TauPOG/TauIDSFs/interface/TauIDSFTool.h"
#include "TauAnalysisTools/TauTriggerSFs/interface/TauTriggerSFs2017.h"
#include "HTTutilities/Jet2TauFakes/interface/FakeFactor.h"
//#include "DesyTauAnalyses/NTupleMaker/interface/ImpactParameter.h"
#include "HiggsCPinTauDecays/ImpactParameter/interface/ImpactParameter.h"
#include "HTT-utilities/RecoilCorrections_KIT/interface/MEtSys.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include "RooFunctor.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LV;
typedef ROOT::Math::XYZPointD Point3D;
typedef ROOT::Math::XYZVectorD PV;

#define pi   3.14159265358979312
#define d2r  1.74532925199432955e-02
#define r2d  57.2957795130823229

#define electronMass 	 0.000511
#define muonMass 	 0.105658
#define tauMass 	 1.77682
#define pionMass 	 0.1396

#define expectedtauspinnerweights 5
void FillLeadingTau(const AC1B * analysisTree, Synch17Tree *otree, int leptonIndex, int tauIndex);
void FillSubLeadingTau(const AC1B * analysisTree, Synch17Tree *otree, int leptonIndex, int tauIndex);
void FillVertices(const AC1B * analysisTree,Synch17Tree *otree, const bool isData, int leptonIndex, int tauIndex, TString channel);
void initializeGenTree(Synch17GenTree *gentree);
void FillGenTree(const AC1B *analysisTree, Synch17GenTree *gentree);
float getEmbeddedWeight(const AC1B * analysisTree, RooWorkspace* WS);
int main(int argc, char * argv[]){

  string cmsswBase = (getenv("CMSSW_BASE"));
   if(argc < 4){
    std::cout << "RUN ERROR: wrong number of arguments"<< std::endl;
    std::cout << "Please run the code in the following way:"<< std::endl;
    std::cout << "SynchNTupleProducer_Run2 NameOfTheConfigurationFile FileList Channel" << std::endl;
    std::cout << "example: SynchNTupleProducer_tt_Run2 analysisMacroSynch_lept_mt_DATA.conf DATA_SingleMuon tt" << std::endl;
    exit(-1);
  }
  // **** configuration analysis  
  Config cfg(argv[1]);
  using namespace std;
  // configuration process
  const string sample = argv[2];
  const bool isData = cfg.get<bool>("isData");
  const string infiles = argv[2];
  TString ch = argv[3];
  std::string lep;
  if (ch != "tt") {
  	std::cout << " Channel " << ch << " is not a valid choice. Please try again with 'tt'. " << std::endl;
    exit(0);
  }
  lumi_json json;
  if (isData){ 
    const string json_name = cfg.get<string>("JSON");
    read_json(TString(TString(cmsswBase) + "/src/" + TString(json_name)).Data(), json);
  }
  const int era = cfg.get<int>("era");
  //  const bool Synch = cfg.get<bool>("Synch"); 
  const bool ApplyPUweight    = cfg.get<bool>("ApplyPUweight"); 
  const bool ApplyLepSF       = cfg.get<bool>("ApplyLepSF"); 
  const bool ApplyTrigger     = cfg.get<bool>("ApplyTrigger"); 
  const bool ApplyFastMTT       = cfg.get<bool>("ApplyFastMTT");
  const bool ApplySVFit       = cfg.get<bool>("ApplySVFit");
  const bool ApplyBTagScaling = cfg.get<bool>("ApplyBTagScaling");
  const bool ApplySystShift   = cfg.get<bool>("ApplySystShift");
  const bool ApplyMetFilters  = cfg.get<bool>("ApplyMetFilters");
  const bool usePuppiMET      = cfg.get<bool>("UsePuppiMET");
  const bool ApplyBTagCP5Correction = cfg.get<bool>("ApplyBTagCP5Correction");
  //const bool ApplyRecoilCorrections = cfg.get<bool>("ApplyRecoilCorrections");
  const bool ApplyIpCorrection = cfg.get<bool>("ApplyIpCorrection");
  const bool applyTauSpinnerWeights = cfg.get<bool>("applyTauSpinnerWeights");
  const float dRTrigMatch       = cfg.get<float>("dRTrigMatch");
  //pileup distrib
  const string pileUpInDataFile = cfg.get<string>("pileUpInDataFile");
  const string pileUpInMCFile = cfg.get<string>("pileUpInMCFile");
  const string pileUpforMC = cfg.get<string>("pileUpforMC");

  const string ipCorrFileName = cfg.get<string>("IpCorrFileName");
  TString IpCorrFileName(ipCorrFileName);

  //TauTriggerSFs2017 *tauTriggerSF = new TauTriggerSFs2017(cmsswBase + "/src/TauAnalysisTools/TauTriggerSFs/data/tauTriggerEfficiencies"+to_string(era) +".root","ditau",to_string(era),"tight",false);//MVAv2
  std::string year;
  if(era==2016)year = "2016Legacy";
  else if(era==2017)year = "2017ReReco";
  else if(era==2018)year = "2018ReReco";

  TauIDSFTool *tauIDSF_medium = new TauIDSFTool(year,"DeepTau2017v2p1VSjet","Medium",false);
  
  IpCorrection *ip = new IpCorrection(TString(cmsswBase) + "/src/" + IpCorrFileName);
  //svFit
  const string svFitPtResFile = TString(TString(cmsswBase)+"/src/"+TString(cfg.get<string>("svFitPtResFile"))).Data();

  //zptweight file 
  const string ZptweightFile = cfg.get<string>("ZptweightFile");
  //b-tag scale factors
  const string BTagAlgorithm = cfg.get<string>("BTagAlgorithm");
  const string BtagSfFile = cmsswBase + "/src/" + cfg.get<string>("BtagSfFile");
  if( ApplyBTagScaling && gSystem->AccessPathName( (TString) BtagSfFile) ){
    cout<<BtagSfFile<<" not found. Please check."<<endl;
    exit(-1);
  }
  // JER
  std::unique_ptr<JME::JetResolution> m_resolution_from_file;
  std::unique_ptr<JME::JetResolutionScaleFactor> m_scale_factor_from_file;
  if (era==2016) {
    m_resolution_from_file.reset(new JME::JetResolution(cmsswBase+"/src/DesyTauAnalyses/NTupleMaker/data/JER/Summer16_25nsV1_MC_PtResolution_AK4PFchs.txt"));
    m_scale_factor_from_file.reset(new JME::JetResolutionScaleFactor(cmsswBase+"/src/DesyTauAnalyses/NTupleMaker/data/JER/Summer16_25nsV1_MC_SF_AK4PFchs.txt"));
  }
  else if (era==2017) {
    m_resolution_from_file.reset(new JME::JetResolution(cmsswBase+"/src/DesyTauAnalyses/NTupleMaker/data/JER/Fall17_V3_MC_PtResolution_AK4PFchs.txt"));
    m_scale_factor_from_file.reset(new JME::JetResolutionScaleFactor(cmsswBase+"/src/DesyTauAnalyses/NTupleMaker/data/JER/Fall17_V3_MC_SF_AK4PFchs.txt"));
  }
  else {
    m_resolution_from_file.reset(new JME::JetResolution(cmsswBase+"/src/DesyTauAnalyses/NTupleMaker/data/JER/Autumn18_V7b_MC_PtResolution_AK4PFchs.txt"));
    m_scale_factor_from_file.reset(new JME::JetResolutionScaleFactor(cmsswBase+"/src/DesyTauAnalyses/NTupleMaker/data/JER/Autumn18_V7b_MC_SF_AK4PFchs.txt"));    
  }

  JME::JetResolution resolution = *m_resolution_from_file;
  JME::JetResolutionScaleFactor resolution_sf = *m_scale_factor_from_file;

  cout<<"using "<<BTagAlgorithm<<endl;
  BTagCalibration calib(BTagAlgorithm, BtagSfFile);
  BTagCalibrationReader reader_B(BTagEntry::OP_MEDIUM, "central",{"up","down"});
  BTagCalibrationReader reader_C(BTagEntry::OP_MEDIUM, "central");
  BTagCalibrationReader reader_Light(BTagEntry::OP_MEDIUM, "central");
  if(ApplyBTagScaling){
    reader_B.load(calib, BTagEntry::FLAV_B, "comb");
    reader_C.load(calib, BTagEntry::FLAV_C, "comb");
    reader_Light.load(calib, BTagEntry::FLAV_UDSG, "incl");
  }
    
  TString pathToTaggingEfficiencies = (TString) cmsswBase + "/src/" + cfg.get<string>("BtagMCeffFile");
  if (ApplyBTagScaling && gSystem->AccessPathName(pathToTaggingEfficiencies)){
    cout<<pathToTaggingEfficiencies<<" not found. Please check."<<endl;
    exit(-1);
  }
  TFile *fileTagging  = new TFile(pathToTaggingEfficiencies);
  TH2F  *tagEff_B     = 0;
  TH2F  *tagEff_C     = 0;
  TH2F  *tagEff_Light = 0;
  TH2F  *tagEff_B_nonCP5     = 0;
  TH2F  *tagEff_C_nonCP5     = 0;
  TH2F  *tagEff_Light_nonCP5 = 0;
  
  if(ApplyBTagScaling){
    tagEff_B     = (TH2F*)fileTagging->Get("btag_eff_b");
    tagEff_C     = (TH2F*)fileTagging->Get("btag_eff_c");
    tagEff_Light = (TH2F*)fileTagging->Get("btag_eff_oth");
  
    if (ApplyBTagCP5Correction) {
      TString pathToTaggingEfficiencies_nonCP5 = (TString) cmsswBase + "/src/" + cfg.get<string>("BtagMCeffFile_nonCP5");
      if (gSystem->AccessPathName(pathToTaggingEfficiencies_nonCP5)) {
	cout<<pathToTaggingEfficiencies_nonCP5<<" not found. Please check."<<endl;
	exit(-1);
      } 
      TFile *fileTagging_nonCP5  = new TFile(pathToTaggingEfficiencies_nonCP5);
      tagEff_B_nonCP5     = (TH2F*)fileTagging_nonCP5->Get("btag_eff_b");
      tagEff_C_nonCP5     = (TH2F*)fileTagging_nonCP5->Get("btag_eff_c");
      tagEff_Light_nonCP5 = (TH2F*)fileTagging_nonCP5->Get("btag_eff_oth");
    }
  }
  TRandom3 *rand = new TRandom3();
  const struct btag_scaling_inputs inputs_btag_scaling_medium = {reader_B, reader_C, reader_Light, tagEff_B, tagEff_C, tagEff_Light, tagEff_B_nonCP5, tagEff_C_nonCP5, tagEff_Light_nonCP5,rand};

  //FakeFactor Workspace
  TFile * ff_file = TFile::Open(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/fakefactors_ws_"+TString(year)+".root");
  if (ff_file->IsZombie()) {
    cout << "File " << TString(cmsswBase) << "/src/DesyTauAnalyses/NTupleMaker/data/fakefactors_ws_" << TString(year) << ".root not found" << endl;
    cout << "Quitting... " << endl;
    exit(-1);
  }
  // FakeFactor* ff = (FakeFactor*)ff_file->Get("ff_comb");
  // std::shared_ptr<RooWorkspace> ff_ws_;
  // std::map<std::string, std::shared_ptr<RooFunctor>> fns_;
  // ff_ws_ = std::shared_ptr<RooWorkspace>((RooWorkspace*)gDirectory->Get("w"));
  // fns_["ff_tt_medium_dmbins"] = std::shared_ptr<RooFunctor>(ff_ws_->function("ff_tt_medium_dmbins")->functor(ff_ws_->argSet("pt,dm,njets,pt_2,os,met_var_qcd")));
  // fns_["ff_tt_medium_mvadmbins"] = std::shared_ptr<RooFunctor>(ff_ws_->function("ff_tt_medium_mvadmbins")->functor(ff_ws_->argSet("pt,mvadm,njets,pt_2,os,met_var_qcd")));

  // MET Recoil Corrections
  const bool isDY = infiles.find("DY") != string::npos;
  const bool isWJets = (infiles.find("WJets") != string::npos) || (infiles.find("W1Jets") != string::npos) || (infiles.find("W2Jets") != string::npos) || (infiles.find("W3Jets") != string::npos) || (infiles.find("W4Jets") != string::npos) || (infiles.find("EWK") != string::npos);
  const bool isVBForGGHiggs = (infiles.find("VBFHTo")!= string::npos) || (infiles.find("GluGluHTo")!= string::npos);
  const bool isEWKZ =  infiles.find("EWKZ") != string::npos;
  const bool isMG = infiles.find("madgraph") != string::npos;
  const bool isMSSMsignal =  (infiles.find("SUSYGluGluToHToTauTau")!= string::npos) || (infiles.find("SUSYGluGluToBBHToTauTau")!= string::npos);
  const bool isTauSpinner = infiles.find("Uncorr") != string::npos;
  const bool isTTbar = infiles.find("TT") != string::npos;

  // bool applyTauSpinnerWeights = false;
  // if(isTauSpinner) applyTauSpinnerWeights = true;
  const bool isEmbedded = cfg.get<bool>("isEmbedded");//infiles.find("Embed") != string::npos;

  const bool ApplyRecoilCorrections = cfg.get<bool>("ApplyRecoilCorrections") && !isEmbedded && !isData && (isDY || isWJets || isVBForGGHiggs || isMSSMsignal);
  kit::RecoilCorrector recoilCorrector(cfg.get<string>("RecoilFilePath"));
  kit::MEtSys MetSys(cfg.get<string>("RecoilSysFilePath"));
  
  
  //tau Cuts 
  const float ptTauCut = cfg.get<float>("ptTauLowCut");
  const float etaTauCut = cfg.get<float>("etaTauCut");
  const float dzTauCut = cfg.get<float>("dzTauCut");

  //const bool Synch = cfg.get<bool>("Synch");
  
  // tau energy scale corrections
  const float shift_tes_1prong = cfg.get<float>("TauEnergyScaleShift_OneProng");
  const float shift_tes_1p1p0  = cfg.get<float>("TauEnergyScaleShift_OneProngOnePi0");
  const float shift_tes_3prong = cfg.get<float>("TauEnergyScaleShift_ThreeProng");
  const float shift_tes_3p1p0 = cfg.get<float>("TauEnergyScaleShift_ThreeProngOnePi0");
  
  const float shift_tes_1prong_e = cfg.get<float>("TauEnergyScaleShift_OneProng_Error");
  const float shift_tes_1p1p0_e  = cfg.get<float>("TauEnergyScaleShift_OneProngOnePi0_Error");
  const float shift_tes_3prong_e = cfg.get<float>("TauEnergyScaleShift_ThreeProng_Error");
  const float shift_tes_3prong1p0_e = cfg.get<float>("TauEnergyScaleShift_ThreeProngOnePi0_Error");

  const float shift_tes_1prong_efake = cfg.get<float>("TauEnergyScaleShift_OneProng_efake");
  const float shift_tes_1p1p0_efake  = cfg.get<float>("TauEnergyScaleShift_OneProngOnePi0_efake");

  const float shift_tes_1prong_mufake = cfg.get<float>("TauEnergyScaleShift_OneProng_mufake");
  const float shift_tes_1p1p0_mufake  = cfg.get<float>("TauEnergyScaleShift_OneProngOnePi0_mufake");

  const string CorrectionWorkspaceFileName = cfg.get<string>("CorrectionWorkspaceFileName");
  
  // *** end of configuration analysis

  //file list creation
  int ifile = 0;
  int jfile = -1;
  //create input files list 
  std::vector<string>fileList;
  int NumberofFile=0;
  if(infiles.find(".root")!=std::string::npos){
    ifile = 0;
    jfile = 1;
    fileList.push_back(infiles);
  }
  else{
    ifstream input;
    std::string infile;
    input.open(infiles);
    while(true){
      input>>infile;
      if(!input.eof()){
	if(infile.length()>0){
	    fileList.push_back(infile);
	    NumberofFile +=1;
	}
      }
      else
	break;
    }
    if(jfile<0)
      jfile = fileList.size();
  }
  if(NumberofFile < jfile)jfile = NumberofFile;
  for(int iF=ifile;iF<jfile;iF++)
    cout<<fileList[iF]<<endl;
  TString rootFileName(sample);
  std::string ntupleName("makeroottree/AC1B");
  std::string initNtupleName("initroottree/AC1B");
  std::string TauSpinnerWeightTreeName("icTauSpinnerProducer/TauSpinnerWeightTree");

  // PU reweighting - initialization
  PileUp *PUofficial = new PileUp();
  if(ApplyPUweight){
    TFile *filePUdistribution_data = new TFile(TString(cmsswBase) + "/src/" + TString(pileUpInDataFile), "read");
    TFile *filePUdistribution_MC = new TFile (TString(cmsswBase) + "/src/" + TString(pileUpInMCFile), "read");
    TH1D *PU_data = (TH1D *)filePUdistribution_data->Get("pileup");    
    TH1D *PU_mc = (TH1D *)filePUdistribution_MC->Get(TString(pileUpforMC));
    if (PU_mc == NULL) {
      std::cout << "Histogram " << pileUpforMC << " is not present in pileup file" << std::endl;
      exit(-1);
    }
    PUofficial->set_h_data(PU_data);
    PUofficial->set_h_MC(PU_mc);
  }  
  // Workspace with corrections
  TString workspace_filename = TString(cmsswBase) + "/src/" + CorrectionWorkspaceFileName;
  cout << "Taking correction workspace from " << workspace_filename << endl;
  TFile *f_workspace = new TFile(workspace_filename, "read");
  if (f_workspace->IsZombie()) {
    std::cout << " workspace file " << workspace_filename << " not found. Please check. " << std::endl;
    exit(-1);
  }
  RooWorkspace *w = (RooWorkspace*)f_workspace->Get("w");
  // Zpt reweighting for LO DY samples 
  TFile *f_zptweight = new TFile(TString(cmsswBase) + "/src/" + ZptweightFile, "read");
  TH2D *h_zptweight = (TH2D*)f_zptweight->Get("zptmass_histo");
  //output fileName with histograms
  rootFileName += "_";
  rootFileName += ifile;
  rootFileName += "_"+ch+"_Sync.root";

  cout<<rootFileName<<endl;
  
  TFile *file = new TFile(rootFileName,"recreate");
  file->cd("");
  
  TH1D *inputEventsH = new TH1D("inputEventsH", "", 1, -0.5, 0.5);
  TH1D *nWeightedEventsH = new TH1D("nWeightedEvents", "", 1, -0.5, 0.5);

  TTree *tree = new TTree("TauCheck", "TauCheck");
  TTree *gtree = new TTree("GenTauCheck", "GenTauCheck");
  Synch17Tree *otree = new Synch17Tree(tree);
  //initializeCPvar(otree);  
  Synch17GenTree *gentree = new Synch17GenTree(gtree);
  Synch17GenTree *gentreeForGoodRecoEvtsOnly = new Synch17GenTree(tree);
    
  int nTotalFiles = 0;
  int nEvents = 0;
  int selEvents = 0;
  int nFiles = 0;

  //svFit
  TH1::AddDirectory(false);  
  TFile *inputFile_visPtResolution = new TFile(svFitPtResFile.data());
  //Systematics init
  
  TauScaleSys *tauScaleSys = 0;
  TauOneProngScaleSys *tauOneProngScaleSys = 0;
  TauOneProngOnePi0ScaleSys *tauOneProngOnePi0ScaleSys = 0;
  TauThreeProngScaleSys *tauThreeProngScaleSys = 0;
  TauThreeProngOnePi0ScaleSys *tauThreeProngOnePi0ScaleSys = 0;

  //BtagSys * btagSys = 0;
  std::vector<JetEnergyScaleSys*> jetEnergyScaleSys;
  JESUncertainties * jecUncertainties = 0;
  std::vector<TString> metSysNames = {"CMS_scale_met_unclustered_13TeV"};
  std::vector<TString> recoilSysNames = {"CMS_htt_boson_reso_met_13TeV",
					 "CMS_htt_boson_scale_met_13TeV"};

  std::vector<PFMETSys*> metSys;
  std::vector<PuppiMETSys*> puppiMetSys;
  if((!isData||isEmbedded) && ApplySystShift){
    tauOneProngScaleSys = new TauOneProngScaleSys(otree);
    tauOneProngScaleSys->SetScale(shift_tes_1prong,shift_tes_1prong_e);
    tauOneProngScaleSys->SetSvFitVisPtResolution(inputFile_visPtResolution);
    tauOneProngScaleSys->SetUseSVFit(ApplySVFit);
    tauOneProngScaleSys->SetUseFastMTT(ApplyFastMTT);
    tauOneProngScaleSys->SetUsePuppiMET(usePuppiMET);
    
    tauOneProngOnePi0ScaleSys = new TauOneProngOnePi0ScaleSys(otree);
    tauOneProngOnePi0ScaleSys->SetScale(shift_tes_1p1p0,shift_tes_1p1p0_e);
    tauOneProngOnePi0ScaleSys->SetSvFitVisPtResolution(inputFile_visPtResolution);
    tauOneProngOnePi0ScaleSys->SetUseSVFit(ApplySVFit);
    tauOneProngOnePi0ScaleSys->SetUseFastMTT(ApplyFastMTT);
    tauOneProngOnePi0ScaleSys->SetUsePuppiMET(usePuppiMET);

    tauThreeProngScaleSys = new TauThreeProngScaleSys(otree);
    tauThreeProngScaleSys->SetScale(shift_tes_3prong,shift_tes_3prong_e);
    tauThreeProngScaleSys->SetSvFitVisPtResolution(inputFile_visPtResolution);
    tauThreeProngScaleSys->SetUseSVFit(ApplySVFit);
    tauThreeProngScaleSys->SetUseFastMTT(ApplyFastMTT);
    tauThreeProngScaleSys->SetUsePuppiMET(usePuppiMET);

    tauThreeProngOnePi0ScaleSys = new TauThreeProngOnePi0ScaleSys(otree);
    tauThreeProngOnePi0ScaleSys->SetScale(shift_tes_3prong,shift_tes_3prong_e);
    tauThreeProngOnePi0ScaleSys->SetSvFitVisPtResolution(inputFile_visPtResolution);
    tauThreeProngOnePi0ScaleSys->SetUseSVFit(ApplySVFit);
    tauThreeProngOnePi0ScaleSys->SetUseFastMTT(ApplyFastMTT);
    tauThreeProngOnePi0ScaleSys->SetUsePuppiMET(usePuppiMET);
    

    // if (!isEmbedded) {
    //   // if (!isDY && !isWJets && !isVBForGGHiggs) {
    //   // 	btagSys = new BtagSys(otree,TString("Btag"));
    //   // 	btagSys->SetConfig(&cfg);
    //   // 	btagSys->SetBtagScaling(&inputs_btag_scaling_medium);
    //   // }
    //   if (ApplyRecoilCorrections) {
    // 	if (usePuppiMET) {
    // 	  for (unsigned int i = 0; i < recoilSysNames.size(); ++i) {
    // 	    PuppiMETSys * puppiMetRecoilSys = new PuppiMETSys(otree,recoilSysNames[i]);
    // 	    puppiMetRecoilSys->SetMEtSys(&MetSys);
    // 	    puppiMetSys.push_back(puppiMetRecoilSys);
    // 	  }
    // 	}
    //   }
    //   else {
    // 	for (unsigned int i = 0; i<metSysNames.size(); ++i) {
    // 	  if (usePuppiMET)
    // 	    puppiMetSys.push_back(new PuppiMETSys(otree,metSysNames[i]));
    // 	  else
    // 	    metSys.push_back(new PFMETSys(otree,metSysNames[i]));
    // 	}
    //   }
    //   if (cfg.get<bool>("splitJES")){
    // 	JESUncertainties *jecUncertainties;
    // 	if (era==2016) 
    // 	  jecUncertainties = new JESUncertainties("DesyTauAnalyses/NTupleMaker/data/Regrouped_Summer16_07Aug2017_V11_MC_UncertaintySources_AK4PFchs.txt");
    // 	else if (era==2017)
    // 	  jecUncertainties = new JESUncertainties("DesyTauAnalyses/NTupleMaker/data/Regrouped_Fall17_17Nov2017_V32_MC_UncertaintySources_AK4PFchs.txt");
    // 	else 
    // 	  jecUncertainties = new JESUncertainties("DesyTauAnalyses/NTupleMaker/data/Regrouped_Autumn18_V19_MC_UncertaintySources_AK4PFchs.txt");
    // 	std::vector<std::string> JESnames = jecUncertainties->getUncertNames();
    // 	for (unsigned int i = 0; i < JESnames.size(); i++) std::cout << "i: "<< i << ", JESnames.at(i) : " << JESnames.at(i) << std::endl;
    // 	for (unsigned int i = 0; i < JESnames.size(); i++){
    // 	  JetEnergyScaleSys *aJESobject = new JetEnergyScaleSys(otree, TString(JESnames.at(i)));
    // 	  aJESobject->SetConfig(&cfg);
    // 	  aJESobject->SetBtagScaling(&inputs_btag_scaling_medium);
    // 	  aJESobject->SetJESUncertainties(jecUncertainties);
    // 	  jetEnergyScaleSys.push_back(aJESobject);
    // 	}  
    //   }
    //   else { // use JEC uncertainty from analysis tree
    // 	JetEnergyScaleSys *singleJES = new JetEnergyScaleSys(otree, TString("JES"));
    // 	singleJES->SetConfig(&cfg);
    // 	singleJES->SetBtagScaling(&inputs_btag_scaling_medium);
    // 	singleJES->SetJESUncertainties(jecUncertainties);
    // 	jetEnergyScaleSys.push_back(singleJES);
    //   }
    //   JetEnergyScaleSys * JERsys = new JetEnergyScaleSys(otree, TString("JER"));
    //   JERsys->SetConfig(&cfg);
    //   JERsys->SetBtagScaling(&inputs_btag_scaling_medium);
    //   JERsys->SetJESUncertainties(jecUncertainties);
    //   jetEnergyScaleSys.push_back(JERsys);
    // }

    if (!isEmbedded) {
      // if (!isDY && !isWJets && !isVBForGGHiggs) {
      // 	btagSys = new BtagSys(otree,TString("Btag"));
      // 	btagSys->SetConfig(&cfg);
      // 	btagSys->SetBtagScaling(&inputs_btag_scaling_medium);
      // }
      if (ApplyRecoilCorrections) {
	if (usePuppiMET) {
	  for (unsigned int i = 0; i < recoilSysNames.size(); ++i) {
	    PuppiMETSys * puppiMetRecoilSys = new PuppiMETSys(otree,recoilSysNames[i]);
	    puppiMetRecoilSys->SetMEtSys(&MetSys);
	    puppiMetSys.push_back(puppiMetRecoilSys);
	  }
	}
      }
      else {
	for (unsigned int i = 0; i<metSysNames.size(); ++i) {
	  if (usePuppiMET)
	    puppiMetSys.push_back(new PuppiMETSys(otree,metSysNames[i]));
	  else
	    metSys.push_back(new PFMETSys(otree,metSysNames[i]));
	}
      }
      if (cfg.get<bool>("splitJES")){
	JESUncertainties *jecUncertainties;
	if (era==2016) 
	  jecUncertainties = new JESUncertainties("DesyTauAnalyses/NTupleMaker/data/Regrouped_Summer16_07Aug2017_V11_MC_UncertaintySources_AK4PFchs.txt");
	else if (era==2017)
	  jecUncertainties = new JESUncertainties("DesyTauAnalyses/NTupleMaker/data/Regrouped_Fall17_17Nov2017_V32_MC_UncertaintySources_AK4PFchs.txt");
	else 
	  jecUncertainties = new JESUncertainties("DesyTauAnalyses/NTupleMaker/data/Regrouped_Autumn18_V19_MC_UncertaintySources_AK4PFchs.txt");
	std::vector<std::string> JESnames = jecUncertainties->getUncertNames();
	for (unsigned int i = 0; i < JESnames.size(); i++) std::cout << "i: "<< i << ", JESnames.at(i) : " << JESnames.at(i) << std::endl;
	for (unsigned int i = 0; i < JESnames.size(); i++){
	  JetEnergyScaleSys *aJESobject = new JetEnergyScaleSys(otree, TString(JESnames.at(i)));
	  aJESobject->SetConfig(&cfg);
	  aJESobject->SetBtagScaling(&inputs_btag_scaling_medium);
	  aJESobject->SetJESUncertainties(jecUncertainties);
	  jetEnergyScaleSys.push_back(aJESobject);
	}  
      }
      else { // use JEC uncertainty from analysis tree
	JetEnergyScaleSys *singleJES = new JetEnergyScaleSys(otree, TString("JES"));
	singleJES->SetConfig(&cfg);
	singleJES->SetBtagScaling(&inputs_btag_scaling_medium);
	singleJES->SetJESUncertainties(jecUncertainties);
	jetEnergyScaleSys.push_back(singleJES);
      }
      JetEnergyScaleSys * JERsys = new JetEnergyScaleSys(otree, TString("JER"));
      JERsys->SetConfig(&cfg);
      JERsys->SetBtagScaling(&inputs_btag_scaling_medium);
      JERsys->SetJESUncertainties(jecUncertainties);
      jetEnergyScaleSys.push_back(JERsys);
    }

  }//shape Systematics
  //list of met filters from config
  std::vector<TString>met_filters_list;
  for(unsigned int i =1;i < (unsigned int)cfg.get<int>("num_met_filters")+1;i++){
    met_filters_list.push_back(cfg.get<string>("met_filter_"+std::to_string(i)));
  }
  int file_count=0;
  ofstream check;
  TString checkfile = "check.txt";
  // if(file_count%50==0)
  //   checkfile = "check_"+std::to_string(file_count)+".txt";
  check.open(checkfile, std::ofstream::out | std::ofstream::app);
  check<<"checking...";
  int counter[20] = {0};
  ///////////////FILE LOOP///////////////

  for (int iF = ifile; iF < jfile; ++iF) {
     std::cout << "file " << iF + 1 << " out of " << fileList.size() << " filename : " << fileList[iF] << std::endl;
     file_count += iF;
    TFile *file_ = TFile::Open(fileList[iF].data());
    TTree *_tree = NULL;
    _tree = (TTree*)file_->Get(TString(ntupleName));  
    if (_tree == NULL) continue;
    
    TH1D *histoInputEvents = NULL;
    histoInputEvents = (TH1D*)file_->Get("makeroottree/nEvents");
    if (histoInputEvents == NULL) continue;
    int NE = int(histoInputEvents->GetEntries());
    std::cout << "      number of input events    = " << NE << std::endl;
    for (int iE = 0; iE < NE; ++iE)
      inputEventsH->Fill(0.);

    AC1B analysisTree(_tree, isData);
    
    double * TSweight = new double[expectedtauspinnerweights];
    TTree  * _treeTauSpinnerWeights = NULL;
    
    TTree * _inittree = NULL;
    _inittree = (TTree*)file_->Get(TString(initNtupleName));
    if (_inittree!=NULL) {
      Float_t genweight;
      if (!isData)
	_inittree->SetBranchAddress("genweight",&genweight);
      Long64_t numberOfEntriesInitTree = _inittree->GetEntries();
      std::cout << "      number of entries in Init Tree = " << numberOfEntriesInitTree << std::endl;
      for (Long64_t iEntry=0; iEntry<numberOfEntriesInitTree; iEntry++) {
	_inittree->GetEntry(iEntry);
	if (isData)
	  nWeightedEventsH->Fill(0.,1.);
	else
	  nWeightedEventsH->Fill(0.,genweight);
      }
    }
    delete _inittree;
    
  

    vector<string> filterDiTau;// = cfg.get<vector<string>>("filterDiTaus");                                                                                                    
    if(era==2018){
      if(isData && !isEmbedded){
        if(analysisTree.event_run >= 317509)
          filterDiTau = cfg.get<vector<string>>("filterDiTauswithHPS");
        else
          filterDiTau = cfg.get<vector<string>>("filterDiTausNoHPS");
      }
      else{
	filterDiTau = cfg.get<vector<string>>("filterDiTaus");
      }
    }
    else{
      filterDiTau = cfg.get<vector<string>>("filterDiTaus");
    }
    cout<<"Number of Ditau trigger legs = "<<filterDiTau.size()<<endl;


    Long64_t numberOfEntries = analysisTree.GetEntries();
    
    for(Long64_t iEntry = 0; iEntry < numberOfEntries; iEntry++){
      counter[0]++;
      analysisTree.GetEntry(iEntry);
      nEvents++;
      bool passed_all_metFilters = passedAllMetFilters(&analysisTree,met_filters_list);
      //if(ApplyMetFilters && !passed_all_metFilters)continue;
      otree->passedAllMetFilters = passed_all_metFilters;
      counter[1]++;
     
      if(!applyTauSpinnerWeights){
	for(int tsindex = 0; tsindex < expectedtauspinnerweights; tsindex++) 
          TSweight[tsindex] = 1;

	otree->TauSpinnerWeightsEven = TSweight[0];

	gentreeForGoodRecoEvtsOnly->sm_htt125 = TSweight[0];
	gentree->sm_htt125 = TSweight[0];
	otree->TauSpinnerWeightsMaxMix = TSweight[1];
	gentree->mm_htt125 = TSweight[1];
	gentreeForGoodRecoEvtsOnly->mm_htt125 = TSweight[1];
	gentree->ps_htt125 = TSweight[2];
	otree->TauSpinnerWeightsOdd = TSweight[2];
     
	gentreeForGoodRecoEvtsOnly->ps_htt125 = TSweight[2];

	otree->TauSpinnerWeightsMinusMaxMix = TSweight[3];
	gentree->minusmm_htt125 = TSweight[3];
	gentreeForGoodRecoEvtsOnly->minusmm_htt125 = TSweight[3];

	otree->TauSpinnerWeightsMix0p375 = TSweight[4];
	gentree->mix0p375_htt125 = TSweight[4];
	gentreeForGoodRecoEvtsOnly->mix0p375_htt125 = TSweight[4];
      
      	}
      else{

      	otree->TauSpinnerWeightsEven = analysisTree.TauSpinnerWeight[0];

      	gentreeForGoodRecoEvtsOnly->sm_htt125 = analysisTree.TauSpinnerWeight[0];

      	otree->TauSpinnerWeightsMaxMix = analysisTree.TauSpinnerWeight[1];
       
      	gentreeForGoodRecoEvtsOnly->mm_htt125 = analysisTree.TauSpinnerWeight[1];

      	otree->TauSpinnerWeightsOdd = analysisTree.TauSpinnerWeight[2];

      	gentreeForGoodRecoEvtsOnly->ps_htt125 = analysisTree.TauSpinnerWeight[2];

      	if(analysisTree.TauSpinAngles_count>=5){
      	  otree->TauSpinnerWeightsMinusMaxMix = analysisTree.TauSpinnerWeight[3];
      	  gentree->minusmm_htt125 = analysisTree.TauSpinnerWeight[3];
      	  gentreeForGoodRecoEvtsOnly->minusmm_htt125 = analysisTree.TauSpinnerWeight[3];
	    
      	  otree->TauSpinnerWeightsMix0p375 = analysisTree.TauSpinnerWeight[4];
      	  gentree->mix0p375_htt125 = analysisTree.TauSpinnerWeight[4];
      	  gentreeForGoodRecoEvtsOnly->mix0p375_htt125 = analysisTree.TauSpinnerWeight[4];
       	}
      }

      
       if (!isData){
       	 FillGenTree(&analysisTree,gentree);
       	 //cout<<"Gen Tauspinor weight  "<<gentree->sm_htt125<<endl;
       	 gentree->Fill();
       }   
       

     
      //  if (!isData){
      // 	FillGenTree(&analysisTree,gentree);
      // 	//cout<<"Gen Tauspinor weight  "<<gentree->sm_htt125<<endl;
      // 	gentree->Fill();
      // }   


      
      vector<int> nDiTauTrig(filterDiTau.size(),-1);      
      if(ApplyTrigger){
	vector<bool> checkFilterDiTauTrig(filterDiTau.size(), false);
      	unsigned int nfilters = analysisTree.run_hltfilters->size();
	for (unsigned int i=0; i<nfilters; ++i) {
	  TString HLTFilter(analysisTree.run_hltfilters->at(i));
	  for(unsigned int i_trig=0; i_trig<filterDiTau.size(); i_trig++){
	    if (HLTFilter==filterDiTau.at(i_trig)){ nDiTauTrig.at(i_trig) = i; checkFilterDiTauTrig.at(i_trig) = true;}
	  } 
	}
      }

      cout<<"OK1"<<endl;
      if (nEvents % 10000 == 0) 
      	cout << "      processed " << nEvents << " events" << endl; 
      
      otree->run  = analysisTree.event_run;
      otree->lumi = analysisTree.event_luminosityblock;
      otree->evt  = analysisTree.event_nr;
      
      otree->npv = analysisTree.primvertex_count;
      otree->npu = analysisTree.numtruepileupinteractions;// numpileupinteractions;
      otree->rho = analysisTree.rho;

      // embedded weight
      otree->embweight = 1;
      if (isEmbedded) {
	otree->embweight = getEmbeddedWeight(&analysisTree, w);
	if (otree->embweight>10.0)
	  cout << "warning : embedding weight = " << otree->embweight << endl;
      }
      initializeGenTree(gentree);
      //GenTree Filling
   
      //LOOP OVER TAUS
      vector<int> taus; taus.clear();
      for (unsigned int it = 0; it < analysisTree.tau_count; ++it) { 

	cout<<"OK2"<<endl;        
	float ptTauCut_corr = ptTauCut;

	if (analysisTree.tau_decayMode[it]==0 && analysisTree.tau_genmatch[it]==5) ptTauCut_corr = ptTauCut/(1+shift_tes_1prong);
	if (analysisTree.tau_decayMode[it]>=1 && analysisTree.tau_decayMode[it]<=3 && analysisTree.tau_genmatch[it]==5) ptTauCut_corr = ptTauCut/(1+shift_tes_1p1p0);
	if (analysisTree.tau_decayMode[it]==10 && analysisTree.tau_genmatch[it]==5) ptTauCut_corr = ptTauCut/(1+shift_tes_3prong);
	if (analysisTree.tau_decayMode[it]==11 && analysisTree.tau_genmatch[it]==5) ptTauCut_corr = ptTauCut/(1+shift_tes_3p1p0);
	// if (analysisTree.tau_decayMode[it]==0 && (analysisTree.tau_genmatch[it]==1 || analysisTree.tau_genmatch[it]==3)) ptTauCut_corr = ptTauCut/shift_tes_1prong_efake;
	// if (analysisTree.tau_decayMode[it]==1 && (analysisTree.tau_genmatch[it]==1 || analysisTree.tau_genmatch[it]==3)) ptTauCut_corr = ptTauCut/shift_tes_1p1p0_efake;
	// if (analysisTree.tau_decayMode[it]==0 && (analysisTree.tau_genmatch[it]==2 || analysisTree.tau_genmatch[it]==4)) ptTauCut_corr = ptTauCut/shift_tes_1prong_mufake;
	// if (analysisTree.tau_decayMode[it]==1 && (analysisTree.tau_genmatch[it]==2 || analysisTree.tau_genmatch[it]==4)) ptTauCut_corr = ptTauCut/shift_tes_1p1p0_mufake;

	if (analysisTree.tau_pt[it] <= ptTauCut_corr) continue;
        if (fabs(analysisTree.tau_eta[it]) >= etaTauCut) continue;
        if (fabs(analysisTree.tau_leadchargedhadrcand_dz[it]) >= dzTauCut) continue;
        if (fabs(fabs(analysisTree.tau_charge[it]) - 1) > 0.001) continue;
      	if (analysisTree.tau_byVVVLooseDeepTau2017v2p1VSjet[it] < 0.5) continue;
      	if (analysisTree.tau_byVVVLooseDeepTau2017v2p1VSe[it] < 0.5) continue;
      	if (analysisTree.tau_byVLooseDeepTau2017v2p1VSmu[it] < 0.5) continue;
        if (analysisTree.tau_decayModeFindingNewDMs[it] < 0.5) continue; //always true, cut applied in NTupleMaker
        if (analysisTree.tau_decayMode[it] == 5 || analysisTree.tau_decayMode[it] == 6) continue;
	//if (analysisTree.tau_MVADM2017v1[it] < -1) continue;
        taus.push_back(it);
      }
      counter[2]++;
      //loop over pair of taus
      int tauIndex_1 = -1;
      int tauIndex_2 = -1;
      float isoTauMin = 1e+10;
      float isoTauMax = -1;
      float lep_pt_max = -1;
      float tau_pt_max = -1;

      cout<<"OK3"<<endl;
      for (unsigned int it=0; it<taus.size(); ++it) {
        for (unsigned int itn=it+1; itn<taus.size(); ++itn) {
	  unsigned int tIndex1 = taus.at(it);
          unsigned int tIndex2 = taus.at(itn);
	  float dR = deltaR(analysisTree.tau_eta[tIndex1],analysisTree.tau_phi[tIndex1],analysisTree.tau_eta[tIndex2],analysisTree.tau_phi[tIndex2]);
	  if(dR < 0.5) continue;
	  float isoTau1 = analysisTree.tau_byDeepTau2017v2p1VSjetraw[tIndex1];  //by AN
          float isoTau2 = analysisTree.tau_byDeepTau2017v2p1VSjetraw[tIndex2]; //by AN
	  float tau_pt1 = analysisTree.tau_pt[tIndex1];
	  float tau_pt2 = analysisTree.tau_pt[tIndex2];
	  counter[3]++;
	  if (analysisTree.tau_decayMode[tIndex1]==0 && analysisTree.tau_genmatch[tIndex1]==5) tau_pt1 = tau_pt1*(1+shift_tes_1prong);
	  if (analysisTree.tau_decayMode[tIndex1]>=1 && analysisTree.tau_decayMode[tIndex1]<=3 && analysisTree.tau_genmatch[tIndex1]==5) tau_pt1 = tau_pt1*(1+shift_tes_1p1p0);
	  if (analysisTree.tau_decayMode[tIndex1]==10 && analysisTree.tau_genmatch[tIndex1]==5) tau_pt1 = tau_pt1*(1+shift_tes_3prong);
	  if (analysisTree.tau_decayMode[tIndex1]==11 && analysisTree.tau_genmatch[tIndex1]==5) tau_pt1 = tau_pt1*(1+shift_tes_3p1p0);

	  if (analysisTree.tau_decayMode[tIndex2]==0 && analysisTree.tau_genmatch[tIndex2]==5) tau_pt2 = tau_pt2*(1+shift_tes_1prong);
	  if (analysisTree.tau_decayMode[tIndex2]>=1 && analysisTree.tau_decayMode[tIndex2]<=3 && analysisTree.tau_genmatch[tIndex2]==5) tau_pt2 = tau_pt2*(1+shift_tes_1p1p0);
	  if (analysisTree.tau_decayMode[tIndex2]==10 && analysisTree.tau_genmatch[tIndex2]==5) tau_pt2 = tau_pt2*(1+shift_tes_3prong);
	  if (analysisTree.tau_decayMode[tIndex2]==11 && analysisTree.tau_genmatch[tIndex2]==5) tau_pt2 = tau_pt2*(1+shift_tes_3p1p0);
	  //change pair
	  bool changepair=false;
	  if(isoTau1 < isoTauMin)
	    changepair = true;
	  else if(fabs(isoTau1-isoTauMin)<1.e-5){
	    if(tau_pt1>lep_pt_max)
	      changepair=true;
	    else if(fabs(tau_pt1 - lep_pt_max)<1.e-5){
	      if(isoTau2 > isoTauMax)
	  	changepair=true;
	      else if(fabs(isoTau2-isoTauMax)<1.e-5){
	  	if(tau_pt2 > tau_pt_max)
	  	  changepair = true;
	      }
	    }
	  }
	  counter[4]++;
	  if(changepair=true){
	    isoTauMin=isoTau1;
	    isoTauMax=isoTau2;
	    lep_pt_max = tau_pt1;
	    tau_pt_max = tau_pt2;
	    tauIndex_1 = tIndex1;
	    tauIndex_2 = tIndex2;
	  }
	    
	}//tau1
      }//tau2
      if(tauIndex_1<0)continue;
      if(tauIndex_2<0)continue;
      counter[5]++;

      cout<<"OK4"<<endl;



      //////////////////////////////////////
      //Trigger Matching
      /////////////////////////////////////

      bool isTau1Trig = false;
      bool isTau2Trig = false;
      bool isDiTauTrig = false;
      otree->trg_doubletau = false;
      for(unsigned int iT=0;iT < analysisTree.trigobject_count;++iT){
	float dRtrigTau1 = deltaR(analysisTree.tau_eta[tauIndex_1],analysisTree.tau_phi[tauIndex_1],
				  analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	float dRtrigTau2 = deltaR(analysisTree.tau_eta[tauIndex_2],analysisTree.tau_phi[tauIndex_2],
				  analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	if(dRtrigTau1 < dRTrigMatch){
	  for(unsigned int i_trig=0;i_trig<filterDiTau.size();i_trig++){
	    if (nDiTauTrig.at(i_trig)==-1)continue;
	    if(analysisTree.trigobject_filters[iT][nDiTauTrig.at(i_trig)])isTau1Trig = true;
	  }
	}
	if(dRtrigTau2 < dRTrigMatch){
	  for(unsigned int i_trig=0;i_trig<filterDiTau.size();i_trig++){
	    if (nDiTauTrig.at(i_trig)==-1)continue;
	    if(analysisTree.trigobject_filters[iT][nDiTauTrig.at(i_trig)])isTau2Trig = true;
	  }
	}
      }
      //L1 matching for tt-channel
      bool is_l1tau_match1=false,is_l1tau_match2=false;
      if(era == 2017){

	for(unsigned int iT=0;iT < analysisTree.trigobject_count;++iT){
	  TLorentzVector l1Tau_LV(analysisTree.l1tau_px[iT],analysisTree.l1tau_py[iT],analysisTree.l1tau_pz[iT],
				  sqrt(pow(analysisTree.l1tau_px[iT],2)+pow(analysisTree.l1tau_py[iT],2)+pow(analysisTree.l1tau_pz[iT],2)));
	  float dRL1trigTau1 = deltaR(analysisTree.tau_eta[tauIndex_1],analysisTree.tau_phi[tauIndex_1],
				      l1Tau_LV.Eta(),l1Tau_LV.Phi());
	  float dRL1trigTau2 = deltaR(analysisTree.tau_eta[tauIndex_2],analysisTree.tau_phi[tauIndex_2],
				      l1Tau_LV.Eta(),l1Tau_LV.Phi());
	  if(dRL1trigTau1 < dRTrigMatch){
	    if(analysisTree.l1tau_pt[iT]>32.0)is_l1tau_match1 =true;
	  }
	  if(dRL1trigTau2 < dRTrigMatch){
	    if(analysisTree.l1tau_pt[iT]>32.0)is_l1tau_match2 =true;
	  }
	  if( is_l1tau_match1 && is_l1tau_match2)break;
	}
      }
      //all criterua passed, we fill vertices here;
      FillVertices(&analysisTree, otree, isData, tauIndex_1, tauIndex_2, ch);
      (era == 2017)?otree->trg_doubletau = (isTau1Trig && isTau2Trig && is_l1tau_match1 && is_l1tau_match2):otree->trg_doubletau = (isTau1Trig && isTau2Trig);
      otree->trg_doubletau = (isTau1Trig && isTau2Trig);
      //if(otree->trg_doubletau)continue;
      FillLeadingTau(&analysisTree,otree,tauIndex_1,tauIndex_2);
      FillSubLeadingTau(&analysisTree,otree,tauIndex_1,tauIndex_2);

      cout<<"OK5"<<endl;
      counter[6]++;
      
      //////////////////////////////////////////////////////////
      //JETS
      /////////////////////////////////////////////////////////
      jets::counting_jets(&analysisTree, otree, &cfg, &inputs_btag_scaling_medium);

      counter[7]++;
      ////////////////////////////////////////////////////////////
      // ID/Iso and Trigger Corrections
      ////////////////////////////////////////////////////////////

      double eff_data_trig_lt_tau1 = 1;
      double eff_mc_trig_lt_tau1   = 1;
      double eff_data_trig_lt_tau2 = 1;
      double eff_mc_trig_lt_tau2 = 1;
      // setting weights to 1
      otree->trkeffweight = 1;
      otree->trigweight_1 = 1;
      otree->trigweight_2 = 1;
      otree->idisoweight_1 = 1;
      otree->idisoweight_antiiso_1 = 1;
      otree->idisoweight_2 = 1;
      otree->idisoweight_antiiso_2 = 1;
      otree->trigweight = 1;
      otree->effweight = 1;
      otree->puweight = 1; 
      otree->mcweight = 1;

      if (ApplyPUweight) 
        otree->puweight = float(PUofficial->get_PUweight(double(analysisTree.numtruepileupinteractions)));
    

      cout<<"OK5.1"<<endl;  
      if(!isData){
      
	otree->mcweight = analysisTree.genweight;
        otree->gen_noutgoing = analysisTree.genparticles_noutgoing;
	if (isEmbedded&&otree->mcweight>1.0)
	  otree->mcweight = 0.0;
	
	TString suffix = "";
	TString suffixRatio = "ratio";
	if (isEmbedded) {suffix = "_embed"; suffixRatio = "embed_ratio";}

	if(analysisTree.tau_genmatch[tauIndex_1] == 5){
	  //w->var("t_dm")->setVal(analysisTree.tau_decayMode[tauIndex_1]);
	  if (analysisTree.tau_MVADM2017v1[tauIndex_1]<0.)
	    w->var("t_mvadm")->setVal(analysisTree.tau_decayMode[tauIndex_1]);
	  else
	    w->var("t_mvadm")->setVal(analysisTree.tau_MVADM2017v1[tauIndex_1]);
	  

	  otree->idisoweight_1 = w->function("t_deeptauid_mvadm"+suffix+"_medium_highpt")->getVal();

	  otree->idisoweight_1 = w->function("t_deeptauid_mvadm"+suffix+"_medium")->getVal();
	  
	}
	if(analysisTree.tau_genmatch[tauIndex_2] == 5){
	  //w->var("t_dm")->setVal(analysisTree.tau_decayMode[tauIndex_2]);
	  if (analysisTree.tau_MVADM2017v1[tauIndex_2]<0.)
	    w->var("t_mvadm")->setVal(analysisTree.tau_decayMode[tauIndex_2]);
	  else
	    w->var("t_mvadm")->setVal(analysisTree.tau_MVADM2017v1[tauIndex_2]);

	  otree->idisoweight_2 = w->function("t_deeptauid_mvadm"+suffix+"_medium_lowpt")->getVal();
	
	}
	cout<<"OK5.2"<<endl;
	TString suffix_trg = "_mc";
	if (isEmbedded) {suffix_trg = "_embed_data";}
      	///Tau triger SF
      	TString tau1mvadm = TString::Itoa(analysisTree.tau_MVADM2017v1[tauIndex_1],10);
      	TString tau2mvadm = TString::Itoa(analysisTree.tau_MVADM2017v1[tauIndex_2],10);
      	if (analysisTree.tau_MVADM2017v1[tauIndex_1]<0.0)
      	  tau1mvadm = TString::Itoa(analysisTree.tau_decayMode[tauIndex_1],10);
      	if (analysisTree.tau_MVADM2017v1[tauIndex_2]<0.0)
      	  tau2mvadm = TString::Itoa(analysisTree.tau_decayMode[tauIndex_2],10);
	
      	w->var("t_pt")->setVal(analysisTree.tau_pt[tauIndex_1]);
      	w->var("t_eta")->setVal(analysisTree.tau_eta[tauIndex_1]);
      	w->var("t_phi")->setVal(analysisTree.tau_phi[tauIndex_1]);
      	//w->var("t_dm")->setVal(analysisTree.tau_decayMode[tauIndex_1]);
      	if (analysisTree.tau_MVADM2017v1[tauIndex_1]<0.)
      	  w->var("t_mvadm")->setVal(analysisTree.tau_decayMode[tauIndex_1]);
      	else
      	  w->var("t_mvadm")->setVal(analysisTree.tau_MVADM2017v1[tauIndex_1]);
	//	t_trg_ic_deeptau_medium_ditau_mvadm1_mc(
      	eff_data_trig_lt_tau1 = w->function("t_trg_ic_deeptau_medium_ditau_mvadm"+tau1mvadm+"_data")->getVal();
      	eff_mc_trig_lt_tau1   = w->function("t_trg_ic_deeptau_medium_ditau_mvadm"+tau1mvadm+suffix_trg)->getVal();
	
	//otree->trigweight_1 = w->function("t_trg_ic_deeptau_medium_ditau_mvadm"+tau1mvadm+"_"+suffixRatio)->getVal();
	
      	w->var("t_pt")->setVal(analysisTree.tau_pt[tauIndex_2]);
      	w->var("t_eta")->setVal(analysisTree.tau_eta[tauIndex_2]);
      	w->var("t_phi")->setVal(analysisTree.tau_phi[tauIndex_2]);
      	//w->var("t_dm")->setVal(analysisTree.tau_decayMode[tauIndex_2]);
      	
	if (analysisTree.tau_MVADM2017v1[tauIndex_2]<0.)
      	  w->var("t_mvadm")->setVal(analysisTree.tau_decayMode[tauIndex_2]);
      	else
      	  w->var("t_mvadm")->setVal(analysisTree.tau_MVADM2017v1[tauIndex_2]);
	
      	eff_data_trig_lt_tau2 = w->function("t_trg_ic_deeptau_medium_ditau_mvadm"+tau2mvadm+"_data")->getVal();
      	eff_mc_trig_lt_tau2 = w->function("t_trg_ic_deeptau_medium_ditau_mvadm"+tau2mvadm + suffix_trg)->getVal();

      	//otree->trigweight_2 = w->function("t_trg_ic_deeptau_medium_ditau_mvadm"+tau2mvadm+"_"+suffixRatio)->getVal();
	
      	otree->trigweight_1 = eff_data_trig_lt_tau1/eff_mc_trig_lt_tau1;
      	otree->trigweight_2 = eff_data_trig_lt_tau2/eff_mc_trig_lt_tau2;;//(eff_data_trig_lt_tau1*eff_data_trig_lt_tau2)/(eff_mc_trig_lt_tau1*eff_mc_trig_lt_tau2);
      
	
	otree->trigweight = (eff_data_trig_lt_tau1*eff_data_trig_lt_tau2)/(eff_mc_trig_lt_tau1*eff_mc_trig_lt_tau2);
	otree->effweight = otree->idisoweight_1 * otree->idisoweight_2 * otree->trigweight_1 * otree->trigweight_2;
	otree->weight = otree->effweight * otree->puweight * otree->mcweight; 
	cout<<"OK6"<<endl;

	otree->idisoweight_2 = w->function("t_deeptauid_mvadm"+suffix+"_medium")->getVal();
	
      }

      
     
	
      otree->effweight = otree->idisoweight_1 * otree->idisoweight_2 * otree->trigweight_1 * otree->trigweight_2;
      otree->weight = otree->effweight * otree->puweight * otree->mcweight; 
      ////////////////////////////////////////////////////////////
      // Z pt weight
      ////////////////////////////////////////////////////////////
      
      TLorentzVector genV( 0., 0., 0., 0.);
      TLorentzVector genL( 0., 0., 0., 0.);

      otree->zptweight = 1.;
      if (!isData && isDY){
        genV = genTools::genV(analysisTree); // gen Z boson ?
	float bosonMass = genV.M();
	float bosonPt = genV.Pt();

        //Merijn determine here some min and max values:
        double massxmin = h_zptweight->GetXaxis()->GetXmin();
        double massxmax = h_zptweight->GetXaxis()->GetXmax();

        double ptxmin = h_zptweight->GetYaxis()->GetXmin();
        double ptxmax = h_zptweight->GetYaxis()->GetXmax();

	//Merijn 2019 6 13: adjust to T/M functions, to get boundaries right. Otherwise, for 2017 data we get few outliers that screw up the weight histogram dramatically.
	Float_t zptmassweight = 1;
	if (bosonMass > 50.0) {
          float bosonMassX = bosonMass;
          float bosonPtX = bosonPt;
          if (bosonMassX > massxmax) bosonMassX = massxmax - h_zptweight->GetXaxis()->GetBinWidth(h_zptweight->GetYaxis()->GetNbins())*0.5;//Merijn: if doesn't work, lower by 1/2 binwidth..
          if (bosonPtX < ptxmin)     bosonPtX = ptxmin + h_zptweight->GetYaxis()->GetBinWidth(1)*0.5;
          if (bosonPtX > ptxmax)     bosonPtX = ptxmax - h_zptweight->GetYaxis()->GetBinWidth(h_zptweight->GetYaxis()->GetNbins())*0.5;
          zptmassweight = h_zptweight->GetBinContent(h_zptweight->GetXaxis()->FindBin(bosonMassX), h_zptweight->GetYaxis()->FindBin(bosonPtX));
	}
	otree->zptweight = zptmassweight;
      }
      
     
     
      //Theory uncertainties for CP analysis
      
      // otree->weight_CMS_scale_gg_13TeVUp   = analysisTree.weightScale4;
      // otree->weight_CMS_scale_gg_13TeVDown = analysisTree.weightScale8;

      // otree->weight_CMS_PS_ISR_ggH_13TeVUp   = 1.;
      // otree->weight_CMS_PS_ISR_ggH_13TeVDown = 1.;
      // otree->weight_CMS_PS_FSR_ggH_13TeVUp   = 1.;
      // otree->weight_CMS_PS_FSR_ggH_13TeVDown = 1.;

      // if(isVBForGGHiggs){
      // 	otree->weight_CMS_PS_ISR_ggH_13TeVUp   = analysisTree.gen_pythiaweights[6];
      // 	otree->weight_CMS_PS_ISR_ggH_13TeVDown = analysisTree.gen_pythiaweights[8];
      // 	otree->weight_CMS_PS_FSR_ggH_13TeVUp   = analysisTree.gen_pythiaweights[7];
      // 	otree->weight_CMS_PS_FSR_ggH_13TeVDown = analysisTree.gen_pythiaweights[9];
      // }

      // //Prefiring weights for CP analysis
      // otree->prefiringweight     = analysisTree.prefiringweight;
      // otree->prefiringweightUp   = analysisTree.prefiringweightup;
      // otree->prefiringweightDown = analysisTree.prefiringweightdown;


      cout<<"OK7"<<endl;
      counter[8]++;
      ////////////////////////////////////////////////////////////
      // MET and Recoil Corrections
      ////////////////////////////////////////////////////////////
      
      otree->met = TMath::Sqrt(analysisTree.pfmetcorr_ex*analysisTree.pfmetcorr_ex + analysisTree.pfmetcorr_ey*analysisTree.pfmetcorr_ey);
      otree->metphi = TMath::ATan2(analysisTree.pfmetcorr_ey,analysisTree.pfmetcorr_ex);
      otree->metcov00 = analysisTree.pfmetcorr_sigxx;
      otree->metcov01 = analysisTree.pfmetcorr_sigxy;
      otree->metcov10 = analysisTree.pfmetcorr_sigyx;
      otree->metcov11 = analysisTree.pfmetcorr_sigyy;
      otree->puppimet = TMath::Sqrt(analysisTree.puppimet_ex*analysisTree.puppimet_ex + analysisTree.puppimet_ey*analysisTree.puppimet_ey);
      otree->puppimetphi = TMath::ATan2(analysisTree.puppimet_ey,analysisTree.puppimet_ex);
      otree->puppimetcov00 = analysisTree.puppimet_sigxx;
      otree->puppimetcov01 = analysisTree.puppimet_sigxy;
      otree->puppimetcov10 = analysisTree.puppimet_sigyx;
      otree->puppimetcov11 = analysisTree.puppimet_sigyy;

      otree->met_uncorr = otree->puppimet;
      otree->metphi_uncorr = otree->puppimetphi;
      otree->njetshad = otree->njets;
      if(ApplyRecoilCorrections){        
	TLorentzVector genV( 0., 0., 0., 0.);
	TLorentzVector genL( 0., 0., 0., 0.);
	genV = genTools::genV(analysisTree);
	genL = genTools::genL(analysisTree);

	genTools::KITRecoilCorrections( recoilCorrector, ApplyRecoilCorrections, // pass the value != 0 to apply corrections
					otree->puppimet, otree->puppimetphi,
					genV.Px(), genV.Py(),
					genL.Px(), genL.Py(),
					otree->njetshad,
					otree->met_rcmr, otree->metphi_rcmr
					);
        
        // overwriting with recoil-corrected values 
        otree->puppimet = otree->met_rcmr;
        otree->puppimetphi = otree->metphi_rcmr;   
      }
      // using PF MET
      TLorentzVector metLV, puppimetLV; 
      metLV.SetXYZT(otree->met*TMath::Cos(otree->metphi), otree->met*TMath::Sin(otree->metphi), 0,
                    TMath::Sqrt( otree->met*TMath::Sin(otree->metphi)*otree->met*TMath::Sin(otree->metphi) +
				 otree->met*TMath::Cos(otree->metphi)*otree->met*TMath::Cos(otree->metphi)));
      puppimetLV.SetXYZT(otree->puppimet*TMath::Cos(otree->puppimetphi), otree->puppimet*TMath::Sin(otree->puppimetphi), 0,
			 TMath::Sqrt( otree->puppimet*TMath::Sin(otree->puppimetphi)*otree->puppimet*TMath::Sin(otree->puppimetphi) +
	
				      otree->puppimet*TMath::Cos(otree->puppimetphi)*otree->puppimet*TMath::Cos(otree->puppimetphi)));
      cout<<"OK8"<<endl;
      counter[9]++;


      //extra lepton veto
      otree->extraelec_veto = extra_electron_veto(tauIndex_1,"tt", &cfg, &analysisTree,era, isEmbedded);
      otree->extramuon_veto = extra_muon_veto(tauIndex_1, "tt", &cfg, &analysisTree, isData);
    
      /////////////////////////////////////////
      //Di-Tau System
      ////////////////////////////////////////
      otree->os = (otree->q_1 * otree->q_2) < 0.0;
      TLorentzVector tauLV1; tauLV1.SetXYZM(analysisTree.tau_px[tauIndex_1],
					  analysisTree.tau_py[tauIndex_1],
					  analysisTree.tau_pz[tauIndex_1],
					  analysisTree.tau_mass[tauIndex_1]);
      TLorentzVector tauLV2; tauLV2.SetXYZM(analysisTree.tau_px[tauIndex_2],
					  analysisTree.tau_py[tauIndex_2],
					  analysisTree.tau_pz[tauIndex_2],
					  analysisTree.tau_mass[tauIndex_2]);


      //////////////////////////////////////////////
      //TauES + MET propagation ////////////////////
      /////////////////////////////////////////////
      if (!isData||isEmbedded) {
      	bool isOneProng1 = false;
      	bool isOneProng2 = false;
      	float shift_tes_1 = 0.0, shift_tes_2 = 0.0;
      	if (otree->gen_match_2 >= 5){
      	  if (otree->tau_decay_mode_2 == 0){
            shift_tes_2 = shift_tes_1prong; 
            isOneProng2 = true;
          }
      	  else if (otree->tau_decay_mode_2 == 1) shift_tes_2 = shift_tes_1p1p0; 
      	  else if (otree->tau_decay_mode_2 == 10) shift_tes_2 = shift_tes_3prong;
      	}

      	if (otree->gen_match_1 >= 5){
      	  if (otree->tau_decay_mode_1 == 0){
            shift_tes_1 = shift_tes_1prong; 
            isOneProng1 = true;
          }
      	  else if (otree->tau_decay_mode_1 == 1) shift_tes_1 = shift_tes_1p1p0; 
      	  else if (otree->tau_decay_mode_1 == 10) shift_tes_1 = shift_tes_3prong;
      	}
      correctTauES(tauLV1, puppimetLV, shift_tes_1, isOneProng1);
      correctTauES(tauLV2, puppimetLV, shift_tes_2, isOneProng2);
      }//isData
      otree->pt_1 = tauLV1.Pt();
      otree->eta_1 = tauLV1.Eta();
      otree->phi_1 = tauLV1.Phi();
      otree->pt_2 = tauLV2.Pt();
      otree->eta_2 = tauLV2.Eta();
      otree->phi_2 = tauLV2.Phi();

      TLorentzVector dileptonLV = tauLV1 + tauLV2;
      otree->m_vis = dileptonLV.M();
      otree->Prompt_pT = dileptonLV.Pt();
      otree->pt_tt = (dileptonLV+metLV).Pt();   
      
      
      TVector3 tau1_vec = tauLV1.Vect();
      TVector3 met_vec = puppimetLV.Vect();
      double deltaPhi = tau1_vec.Dot(met_vec)/(tau1_vec.Mag()*met_vec.Mag());
      otree->met_var_qcd = otree->puppimet*TMath::Cos(deltaPhi)/tauLV1.Pt();
    
      if (usePuppiMET)
	otree->pt_tt = (dileptonLV+puppimetLV).Pt();
      

      otree->m_sv   = -10;
      otree->pt_sv  = -10;
      otree->eta_sv = -10;
      otree->phi_sv = -10;
      otree->met_sv = -10;
      otree->mt_sv = -10;
      otree->m_fast = -10;
      otree->mt_fast = -10;
      otree->pt_fast = -10;
      otree->phi_fast = -10;
      otree->eta_fast = -10;
      
      bool Synch = otree->trg_doubletau>0.5&&otree->byMediumDeepTau2017v2p1VSjet_1>0.5;
      if((ApplySVFit||ApplyFastMTT)&& otree->trg_doubletau>0.5&&otree->byMediumDeepTau2017v2p1VSjet_2>0.5){
	svfit_variables(ch, &analysisTree, otree, &cfg, inputFile_visPtResolution);
      }
      cout<<"OK9"<<endl;
      counter[10]++;
      // ***********************************
      // ** IPSignificance calibration ->
      // ***********************************


      //acott_Impr_tt(&analysisTree, otree, tauIndex_1, tauIndex_2);
      counter[11]++;
      if ((!isData||isEmbedded) && ApplySystShift) {
	tauOneProngScaleSys->Eval(utils::TAUTAU);
	tauOneProngOnePi0ScaleSys->Eval(utils::TAUTAU);
	tauThreeProngScaleSys->Eval(utils::TAUTAU);
	tauThreeProngOnePi0ScaleSys->Eval(utils::TAUTAU);
      }

      cout<<"OK10"<<endl;
      //acott_Impr_tt(&analysisTree, otree, tauIndex_1, tauIndex_2);
      
      counter[11]++;

      otree->Fill();
     
      ++selEvents;  
    }//event loop
    delete _tree;
    file_->Close();
    delete file_;
    check.close();
  }//file closing
  std::cout<<"COUNTERS "<<std::endl;

  for(int iC = 0; iC < 20; iC++) std::cout << "Counter " << iC << ":    " << counter[iC] << std::endl;

    for(int iC = 0; iC < 20; iC++) std::cout << "Counter " << iC << ":    " << counter[iC] << std::endl;

  std::cout << std::endl;
  std::cout << "Total number of input events    = " << int(inputEventsH->GetEntries()) << std::endl;
  std::cout << "Total number of events in Tree  = " << nEvents << std::endl;
  std::cout << "Total number of selected events = " << selEvents << std::endl;
  std::cout << std::endl;
  
  file->cd("");
  file->Write();
  // delete systematics objects

  if(tauScaleSys != 0){
    tauScaleSys->Write();
    delete tauScaleSys;
  }
  if(tauOneProngScaleSys != 0){
    tauOneProngScaleSys->Write();
    delete tauOneProngScaleSys;
  }
  if(tauOneProngOnePi0ScaleSys != 0){
    tauOneProngOnePi0ScaleSys->Write();
    delete tauOneProngOnePi0ScaleSys;
  }

  if(tauThreeProngScaleSys != 0){
    tauThreeProngScaleSys->Write();
    delete tauThreeProngScaleSys;
  }

  if(tauThreeProngOnePi0ScaleSys != 0){
    tauThreeProngOnePi0ScaleSys->Write();
    delete tauThreeProngOnePi0ScaleSys;
  }
  if(jetEnergyScaleSys.size() > 0){
    for (unsigned int i = 0; i < jetEnergyScaleSys.size(); i++){
      (jetEnergyScaleSys.at(i))->Write();
      delete jetEnergyScaleSys.at(i);
    }
  }

  if (metSys.size() > 0){
    for (unsigned int i = 0; i < metSys.size(); i++ ) {
      (metSys.at(i))->Write();
      delete metSys.at(i);
    }
  }

  if (puppiMetSys.size() > 0){
    for (unsigned int i = 0; i < puppiMetSys.size(); i++ ) {
      (puppiMetSys.at(i))->Write();
      delete puppiMetSys.at(i);
    }
  }
  file->Close();
  delete file;

}//closing int main()
//fill the otree with the tau variables 
void FillSubLeadingTau(const AC1B *analysisTree, Synch17Tree *otree, int leptonIndex, int tauIndex){
  otree->pt_2 = analysisTree->tau_pt[tauIndex];
  otree->eta_2 = analysisTree->tau_eta[tauIndex];
  otree->phi_2 = analysisTree->tau_phi[tauIndex];
  otree->q_2 = analysisTree->tau_charge[tauIndex];
  otree->gen_match_2 = analysisTree->tau_genmatch[tauIndex];
  otree->mva_2 = analysisTree->tau_byTightIsolationMVArun2017v2DBoldDMwLT2017[tauIndex];
  otree->mva17_2= analysisTree->tau_byTightIsolationMVArun2017v2DBoldDMwLT2017[tauIndex];
  otree->d0_2 = analysisTree->tau_leadchargedhadrcand_dxy[tauIndex];
  otree->dZ_2 = analysisTree->tau_leadchargedhadrcand_dz[tauIndex];      
  otree->iso_2 = analysisTree->tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[tauIndex];
  otree->m_2 = analysisTree->tau_mass[tauIndex];
  otree->tau_decay_mode_2 = analysisTree->tau_decayMode[tauIndex];
  otree->dm_2 = analysisTree->tau_decayMode[tauIndex];
  otree->dmMVA_2 = analysisTree->tau_MVADM2017v1[tauIndex];

  std::vector<float> PV_with_BS_cov_components = {};
  for(auto i:  analysisTree->primvertexwithbs_cov) PV_with_BS_cov_components.push_back(i);
  TVector3 PV_with_BS (analysisTree->primvertexwithbs_x, analysisTree->primvertexwithbs_y, analysisTree->primvertexwithbs_z );
  //otree->IP_signif_PV_with_BS_2 = IP_significance_helix_tauh(analysisTree, tauIndex, PV_with_BS, PV_with_BS_cov_components);

  TLorentzVector constituents_P4 = charged_constituents_P4(analysisTree, tauIndex);
  TLorentzVector tau_P4;
  tau_P4.SetXYZM(analysisTree->tau_px[tauIndex],
                 analysisTree->tau_py[tauIndex],
                 analysisTree->tau_pz[tauIndex],
                 analysisTree->tau_mass[tauIndex]);
  otree->chpt_2 = constituents_P4.Pt();
  otree->cheta_2 = constituents_P4.Eta();
  otree->chphi_2 = constituents_P4.Phi();
  otree->chm_2 = constituents_P4.M();
  
  otree->npt_2 = (tau_P4 - constituents_P4).Pt();
  otree->neta_2 = (tau_P4 - constituents_P4).Eta();
  otree->nphi_2 = (tau_P4 - constituents_P4).Phi();
  otree->nm_2 = (tau_P4 - constituents_P4).M();
  
  otree->tau_pca2D_x_2 = analysisTree->tau_pca2D_x[tauIndex];
  otree->tau_pca2D_y_2 = analysisTree->tau_pca2D_y[tauIndex];
  otree->tau_pca2D_z_2 = analysisTree->tau_pca2D_z[tauIndex];
  otree->tau_pca3D_x_2 = analysisTree->tau_pca3D_x[tauIndex];
  otree->tau_pca3D_y_2 = analysisTree->tau_pca3D_y[tauIndex];
  otree->tau_pca3D_z_2 = analysisTree->tau_pca3D_z[tauIndex];
  otree->tau_SV_x_2 = analysisTree->tau_SV_x[tauIndex];
  otree->tau_SV_y_2 = analysisTree->tau_SV_y[tauIndex];
  otree->tau_SV_z_2 = analysisTree->tau_SV_z[tauIndex];
  otree->tau_SV_covxx_2 = analysisTree->tau_SV_cov[tauIndex][0];
  otree->tau_SV_covyx_2 = analysisTree->tau_SV_cov[tauIndex][1];
  otree->tau_SV_covzx_2 = analysisTree->tau_SV_cov[tauIndex][2];
  otree->tau_SV_covyy_2 = analysisTree->tau_SV_cov[tauIndex][3];
  otree->tau_SV_covzy_2 = analysisTree->tau_SV_cov[tauIndex][4];
  otree->tau_SV_covzz_2 = analysisTree->tau_SV_cov[tauIndex][5];

  otree->deepTauVsEleRaw_2                = analysisTree->tau_byDeepTau2017v2p1VSeraw[tauIndex];
  otree->deepTauVsJetRaw_2                = analysisTree->tau_byDeepTau2017v2p1VSjetraw[tauIndex];
  otree->deepTauVsMuRaw_2                 = analysisTree->tau_byDeepTau2017v2p1VSmuraw[tauIndex];
  otree->byLooseDeepTau2017v2p1VSe_2      = analysisTree->tau_byLooseDeepTau2017v2p1VSe[tauIndex];
  otree->byLooseDeepTau2017v2p1VSjet_2    = analysisTree->tau_byLooseDeepTau2017v2p1VSjet[tauIndex];
  otree->byLooseDeepTau2017v2p1VSmu_2     = analysisTree->tau_byLooseDeepTau2017v2p1VSmu[tauIndex];
  otree->byMediumDeepTau2017v2p1VSe_2     = analysisTree->tau_byMediumDeepTau2017v2p1VSe[tauIndex];
  otree->byMediumDeepTau2017v2p1VSjet_2   = analysisTree->tau_byMediumDeepTau2017v2p1VSjet[tauIndex];
  otree->byMediumDeepTau2017v2p1VSmu_2    = analysisTree->tau_byMediumDeepTau2017v2p1VSmu[tauIndex];
  otree->byTightDeepTau2017v2p1VSe_2      = analysisTree->tau_byTightDeepTau2017v2p1VSe[tauIndex];
  otree->byTightDeepTau2017v2p1VSjet_2    = analysisTree->tau_byTightDeepTau2017v2p1VSjet[tauIndex];
  otree->byTightDeepTau2017v2p1VSmu_2     = analysisTree->tau_byTightDeepTau2017v2p1VSmu[tauIndex];
  otree->byVLooseDeepTau2017v2p1VSe_2     = analysisTree->tau_byVLooseDeepTau2017v2p1VSe[tauIndex];
  otree->byVLooseDeepTau2017v2p1VSjet_2   = analysisTree->tau_byVLooseDeepTau2017v2p1VSjet[tauIndex];
  otree->byVLooseDeepTau2017v2p1VSmu_2    = analysisTree->tau_byVLooseDeepTau2017v2p1VSmu[tauIndex];
  otree->byVTightDeepTau2017v2p1VSe_2     = analysisTree->tau_byVTightDeepTau2017v2p1VSe[tauIndex];
  otree->byVTightDeepTau2017v2p1VSjet_2   = analysisTree->tau_byVTightDeepTau2017v2p1VSjet[tauIndex];
  otree->byVVLooseDeepTau2017v2p1VSe_2    = analysisTree->tau_byVVLooseDeepTau2017v2p1VSe[tauIndex];
  otree->byVVLooseDeepTau2017v2p1VSjet_2  = analysisTree->tau_byVVLooseDeepTau2017v2p1VSjet[tauIndex];
  otree->byVVTightDeepTau2017v2p1VSe_2    = analysisTree->tau_byVVTightDeepTau2017v2p1VSe[tauIndex];
  otree->byVVTightDeepTau2017v2p1VSjet_2  = analysisTree->tau_byVVTightDeepTau2017v2p1VSjet[tauIndex];
  otree->byVVVLooseDeepTau2017v2p1VSe_2   = analysisTree->tau_byVVVLooseDeepTau2017v2p1VSe[tauIndex];
  otree->byVVVLooseDeepTau2017v2p1VSjet_2 = analysisTree->tau_byVVVLooseDeepTau2017v2p1VSjet[tauIndex];

  otree->MVADM2017v1DM0raw_2 = analysisTree->tau_MVADM2017v1DM0raw[tauIndex];
  otree->MVADM2017v1DM10raw_2 = analysisTree->tau_MVADM2017v1DM10raw[tauIndex];
  otree->MVADM2017v1DM11raw_2 = analysisTree->tau_MVADM2017v1DM11raw[tauIndex];
  otree->MVADM2017v1DM1raw_2 = analysisTree->tau_MVADM2017v1DM1raw[tauIndex];
  otree->MVADM2017v1DM2raw_2 = analysisTree->tau_MVADM2017v1DM2raw[tauIndex];
  otree->MVADM2017v1DMotherraw_2 = analysisTree->tau_MVADM2017v1DMotherraw[tauIndex];

  otree->byCombinedIsolationDeltaBetaCorrRaw3Hits_2 = analysisTree->tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[tauIndex];
  otree->byLooseCombinedIsolationDeltaBetaCorr3Hits_2 = analysisTree->tau_byLooseCombinedIsolationDeltaBetaCorr3Hits[tauIndex];
  otree->byMediumCombinedIsolationDeltaBetaCorr3Hits_2 = analysisTree->tau_byMediumCombinedIsolationDeltaBetaCorr3Hits[tauIndex];
  otree->byTightCombinedIsolationDeltaBetaCorr3Hits_2 = analysisTree->tau_byTightCombinedIsolationDeltaBetaCorr3Hits[tauIndex];
  otree->againstMuonLoose3_2 = analysisTree->tau_againstMuonLoose3[tauIndex];
  otree->againstMuonTight3_2 = analysisTree->tau_againstMuonTight3[tauIndex];
  otree->againstElectronLooseMVA6_2 = analysisTree->tau_againstElectronLooseMVA6[tauIndex];
  otree->againstElectronMediumMVA6_2 = analysisTree->tau_againstElectronMediumMVA6[tauIndex];
  otree->againstElectronTightMVA6_2 = analysisTree->tau_againstElectronTightMVA6[tauIndex];
  otree->againstElectronVLooseMVA6_2 = analysisTree->tau_againstElectronVLooseMVA6[tauIndex];
  otree->againstElectronVTightMVA6_2 = analysisTree->tau_againstElectronVTightMVA6[tauIndex];

  otree->byVLooseIsolationMVArun2017v2DBoldDMwLT2017_2 = analysisTree->tau_byVLooseIsolationMVArun2017v2DBoldDMwLT2017[tauIndex];
  otree->byLooseIsolationMVArun2017v2DBoldDMwLT2017_2 = analysisTree->tau_byLooseIsolationMVArun2017v2DBoldDMwLT2017[tauIndex];
  otree->byMediumIsolationMVArun2017v2DBoldDMwLT2017_2 = analysisTree->tau_byMediumIsolationMVArun2017v2DBoldDMwLT2017[tauIndex];
  otree->byTightIsolationMVArun2017v2DBoldDMwLT2017_2 = analysisTree->tau_byTightIsolationMVArun2017v2DBoldDMwLT2017[tauIndex];
  otree->byVTightIsolationMVArun2017v2DBoldDMwLT2017_2 = analysisTree->tau_byVTightIsolationMVArun2017v2DBoldDMwLT2017[tauIndex];
  otree->byVVTightIsolationMVArun2017v2DBoldDMwLT2017_2 = analysisTree->tau_byVVTightIsolationMVArun2017v2DBoldDMwLT2017[tauIndex];
  otree->byIsolationMVArun2017v2DBoldDMwLTraw2017_2 = analysisTree->tau_byIsolationMVArun2017v2DBoldDMwLTraw2017[tauIndex];
  otree->byIsolationMVA3oldDMwLTraw_2 = analysisTree->tau_byIsolationMVArun2017v2DBoldDMwLTraw2017[tauIndex];
  otree->chargedIsoPtSum_2 = analysisTree->tau_chargedIsoPtSum[tauIndex];
  otree->neutralIsoPtSum_2 = analysisTree->tau_neutralIsoPtSum[tauIndex];
  otree->puCorrPtSum_2 = analysisTree->tau_puCorrPtSum[tauIndex];
  //otree->isolationGammaCands_size_2 = analysisTree->tau_isolationGammaCands_size[tauIndex];
  //otree->signalGammaCands_size_2 = analysisTree->tau_signalGammaCands_size[tauIndex];

  otree->correction_againstElectronVLooseMVA6_2 = 1;
  //otree->correction_againstElectronLooseMVA6_2 = 1;
  //otree->correction_againstElectronMediumMVA6_2 = 1;
  otree->correction_againstElectronTightMVA6_2 = 1;
  //otree->correction_againstElectronVTightMVA6_2 = 1;
  otree->correction_againstMuonLoose3_2 = 1;
  otree->correction_againstMuonTight3_2 = 1;

  if (analysisTree->tau_genmatch[tauIndex]==1){
    if (abs(otree->eta_2)<1.460){
      otree->correction_againstElectronVLooseMVA6_2 = 1.09;
      //		otree->correction_againstElectronLooseMVA6_2 = 1.17;
      //		otree->correction_againstElectronMediumMVA6_2 = 1.40;
      otree->correction_againstElectronTightMVA6_2 = 1.80;
      //		otree->correction_againstElectronVTightMVA6_2 = 1.96;
    }else if (abs(otree->eta_2)>1.558){
      otree->correction_againstElectronVLooseMVA6_2 = 1.19;
      //		otree->correction_againstElectronLooseMVA6_2 = 1.25;
      //		otree->correction_againstElectronMediumMVA6_2 = 1.21;
      otree->correction_againstElectronTightMVA6_2 = 1.53;
      //		otree->correction_againstElectronVTightMVA6_2 = 1.66;
    }
  }
  
  if (analysisTree->tau_genmatch[tauIndex]==2){
    if ((0 < abs(otree->eta_2))&&(abs(otree->eta_2)<= 0.4)){
      otree->correction_againstMuonLoose3_2 = 1.06;
      otree->correction_againstMuonTight3_2 = 1.17;
    }else if ((0.4 < abs(otree->eta_2))&&(abs(otree->eta_2) <= 0.8)){
      otree->correction_againstMuonLoose3_2 = 1.02;
      otree->correction_againstMuonTight3_2 = 1.14;
    }else if ((0.8 < abs(otree->eta_2))&&(abs(otree->eta_2) <= 1.2)){
      otree->correction_againstMuonLoose3_2 = 1.10;
      otree->correction_againstMuonTight3_2 = 1.14;
    }else if ((1.2 < abs(otree->eta_2))&&(abs(otree->eta_2) <= 1.7)){
      otree->correction_againstMuonLoose3_2 = 1.03;
      otree->correction_againstMuonTight3_2 = 0.93;
    }else if ((1.7 < abs(otree->eta_2))&&(abs(otree->eta_2) <= 2.3)){
      otree->correction_againstMuonLoose3_2 = 1.94;
      otree->correction_againstMuonTight3_2 = 1.61;
    }
  }
  ///////////////////////////////////////////////////////NEW
  otree->efficiency_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017_2 = 1;
  //  otree->efficiency_byLooseIsolationMVArun2017v2DBoldDMwLT2017_2 = 1;
  //  otree->efficiency_byMediumIsolationMVArun2017v2DBoldDMwLT2017_2 = 1;
  otree->efficiency_byTightIsolationMVArun2017v2DBoldDMwLT2017_2 = 1;
  //  otree->efficiency_byVTightIsolationMVArun2017v2DBoldDMwLT2017_2 = 1;
  //  otree->efficiency_byVVTightIsolationMVArun2017v2DBoldDMwLT2017_2 = 1;
  if (analysisTree->tau_genmatch[tauIndex]==5){
    otree->efficiency_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017_2 = 0.88;
    //  otree->efficiency_byLooseIsolationMVArun2017v2DBoldDMwLT2017_2 = 0.89;
    //  otree->efficiency_byMediumIsolationMVArun2017v2DBoldDMwLT2017_2 = 0.89;
    otree->efficiency_byTightIsolationMVArun2017v2DBoldDMwLT2017_2 = 0.89;
    //  otree->efficiency_byVTightIsolationMVArun2017v2DBoldDMwLT2017_2 = 0.86;
    //  otree->efficiency_byVVTightIsolationMVArun2017v2DBoldDMwLT2017_2 = 0.84;
  }

}
void FillLeadingTau(const AC1B *analysisTree, Synch17Tree *otree, int leptonIndex, int tauIndex){
  otree->pt_1 = analysisTree->tau_pt[leptonIndex];
  otree->eta_1 = analysisTree->tau_eta[leptonIndex];
  otree->phi_1 = analysisTree->tau_phi[leptonIndex];
  otree->q_1 = analysisTree->tau_charge[leptonIndex];
  otree->gen_match_1 = analysisTree->tau_genmatch[leptonIndex];
  otree->mva_1 = analysisTree->tau_byTightIsolationMVArun2017v2DBoldDMwLT2017[leptonIndex];
  otree->mva17_1= analysisTree->tau_byTightIsolationMVArun2017v2DBoldDMwLT2017[leptonIndex];
  otree->d0_1 = analysisTree->tau_leadchargedhadrcand_dxy[leptonIndex];
  otree->dZ_1 = analysisTree->tau_leadchargedhadrcand_dz[leptonIndex];      
  otree->iso_1 = analysisTree->tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[leptonIndex];
  otree->m_1 = analysisTree->tau_mass[leptonIndex];
  otree->tau_decay_mode_1 = analysisTree->tau_decayMode[leptonIndex];
  otree->dm_1 = analysisTree->tau_decayMode[leptonIndex];
  otree->dmMVA_1 = analysisTree->tau_MVADM2017v1[leptonIndex];

  std::vector<float> PV_with_BS_cov_components = {};
  for(auto i:  analysisTree->primvertexwithbs_cov) PV_with_BS_cov_components.push_back(i);
  TVector3 PV_with_BS (analysisTree->primvertexwithbs_x, analysisTree->primvertexwithbs_y, analysisTree->primvertexwithbs_z );
  //otree->IP_signif_PV_with_BS_1 = IP_significance_helix_mu(analysisTree, leptonIndex, PV_with_BS, PV_with_BS_cov_components);
  TLorentzVector constituents_P4 = charged_constituents_P4(analysisTree, tauIndex);
  
  TLorentzVector tau_P4;
  tau_P4.SetXYZM(analysisTree->tau_px[leptonIndex],
                 analysisTree->tau_py[leptonIndex],
                 analysisTree->tau_pz[leptonIndex],
                 analysisTree->tau_mass[leptonIndex]);
  otree->chpt_1 = constituents_P4.Pt();
  otree->cheta_1 = constituents_P4.Eta();
  otree->chphi_1 = constituents_P4.Phi();
  otree->chm_1 = constituents_P4.M();
  
  otree->npt_1 = (tau_P4 - constituents_P4).Pt();
  otree->neta_1 = (tau_P4 - constituents_P4).Eta();
  otree->nphi_1 = (tau_P4 - constituents_P4).Phi();
  otree->nm_1 = (tau_P4 - constituents_P4).M();
  
  otree->tau_pca2D_x_1 = analysisTree->tau_pca2D_x[leptonIndex];
  otree->tau_pca2D_y_1 = analysisTree->tau_pca2D_y[leptonIndex];
  otree->tau_pca2D_z_1 = analysisTree->tau_pca2D_z[leptonIndex];
  otree->tau_pca3D_x_1 = analysisTree->tau_pca3D_x[leptonIndex];
  otree->tau_pca3D_y_1 = analysisTree->tau_pca3D_y[leptonIndex];
  otree->tau_pca3D_z_1 = analysisTree->tau_pca3D_z[leptonIndex];
  otree->tau_SV_x_1 = analysisTree->tau_SV_x[leptonIndex];
  otree->tau_SV_y_1 = analysisTree->tau_SV_y[leptonIndex];
  otree->tau_SV_z_1 = analysisTree->tau_SV_z[leptonIndex];
  otree->tau_SV_covxx_1 = analysisTree->tau_SV_cov[leptonIndex][0];
  otree->tau_SV_covyx_1 = analysisTree->tau_SV_cov[leptonIndex][1];
  otree->tau_SV_covzx_1 = analysisTree->tau_SV_cov[leptonIndex][2];
  otree->tau_SV_covyy_1 = analysisTree->tau_SV_cov[leptonIndex][3];
  otree->tau_SV_covzy_1 = analysisTree->tau_SV_cov[leptonIndex][4];
  otree->tau_SV_covzz_1 = analysisTree->tau_SV_cov[leptonIndex][5];

  otree->deepTauVsEleRaw_1                = analysisTree->tau_byDeepTau2017v2p1VSeraw[leptonIndex];
  otree->deepTauVsJetRaw_1                = analysisTree->tau_byDeepTau2017v2p1VSjetraw[leptonIndex];
  otree->deepTauVsMuRaw_1                 = analysisTree->tau_byDeepTau2017v2p1VSmuraw[leptonIndex];
  otree->byLooseDeepTau2017v2p1VSe_1      = analysisTree->tau_byLooseDeepTau2017v2p1VSe[leptonIndex];
  otree->byLooseDeepTau2017v2p1VSjet_1    = analysisTree->tau_byLooseDeepTau2017v2p1VSjet[leptonIndex];
  otree->byLooseDeepTau2017v2p1VSmu_1     = analysisTree->tau_byLooseDeepTau2017v2p1VSmu[leptonIndex];
  otree->byMediumDeepTau2017v2p1VSe_1     = analysisTree->tau_byMediumDeepTau2017v2p1VSe[leptonIndex];
  otree->byMediumDeepTau2017v2p1VSjet_1   = analysisTree->tau_byMediumDeepTau2017v2p1VSjet[leptonIndex];
  otree->byMediumDeepTau2017v2p1VSmu_1    = analysisTree->tau_byMediumDeepTau2017v2p1VSmu[leptonIndex];
  otree->byTightDeepTau2017v2p1VSe_1      = analysisTree->tau_byTightDeepTau2017v2p1VSe[leptonIndex];
  otree->byTightDeepTau2017v2p1VSjet_1    = analysisTree->tau_byTightDeepTau2017v2p1VSjet[leptonIndex];
  otree->byTightDeepTau2017v2p1VSmu_1     = analysisTree->tau_byTightDeepTau2017v2p1VSmu[leptonIndex];
  otree->byVLooseDeepTau2017v2p1VSe_1     = analysisTree->tau_byVLooseDeepTau2017v2p1VSe[leptonIndex];
  otree->byVLooseDeepTau2017v2p1VSjet_1   = analysisTree->tau_byVLooseDeepTau2017v2p1VSjet[leptonIndex];
  otree->byVLooseDeepTau2017v2p1VSmu_1    = analysisTree->tau_byVLooseDeepTau2017v2p1VSmu[leptonIndex];
  otree->byVTightDeepTau2017v2p1VSe_1     = analysisTree->tau_byVTightDeepTau2017v2p1VSe[leptonIndex];
  otree->byVTightDeepTau2017v2p1VSjet_1   = analysisTree->tau_byVTightDeepTau2017v2p1VSjet[leptonIndex];
  otree->byVVLooseDeepTau2017v2p1VSe_1    = analysisTree->tau_byVVLooseDeepTau2017v2p1VSe[leptonIndex];
  otree->byVVLooseDeepTau2017v2p1VSjet_1  = analysisTree->tau_byVVLooseDeepTau2017v2p1VSjet[leptonIndex];
  otree->byVVTightDeepTau2017v2p1VSe_1    = analysisTree->tau_byVVTightDeepTau2017v2p1VSe[leptonIndex];
  otree->byVVTightDeepTau2017v2p1VSjet_1  = analysisTree->tau_byVVTightDeepTau2017v2p1VSjet[leptonIndex];
  otree->byVVVLooseDeepTau2017v2p1VSe_1   = analysisTree->tau_byVVVLooseDeepTau2017v2p1VSe[leptonIndex];
  otree->byVVVLooseDeepTau2017v2p1VSjet_1 = analysisTree->tau_byVVVLooseDeepTau2017v2p1VSjet[leptonIndex];

  otree->MVADM2017v1DM0raw_1 = analysisTree->tau_MVADM2017v1DM0raw[leptonIndex];
  otree->MVADM2017v1DM10raw_1 = analysisTree->tau_MVADM2017v1DM10raw[leptonIndex];
  otree->MVADM2017v1DM11raw_1 = analysisTree->tau_MVADM2017v1DM11raw[leptonIndex];
  otree->MVADM2017v1DM1raw_1 = analysisTree->tau_MVADM2017v1DM1raw[leptonIndex];
  otree->MVADM2017v1DM2raw_1 = analysisTree->tau_MVADM2017v1DM2raw[leptonIndex];
  otree->MVADM2017v1DMotherraw_1 = analysisTree->tau_MVADM2017v1DMotherraw[leptonIndex];

  otree->byCombinedIsolationDeltaBetaCorrRaw3Hits_1 = analysisTree->tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[leptonIndex];
  otree->byLooseCombinedIsolationDeltaBetaCorr3Hits_1 = analysisTree->tau_byLooseCombinedIsolationDeltaBetaCorr3Hits[leptonIndex];
  otree->byMediumCombinedIsolationDeltaBetaCorr3Hits_1 = analysisTree->tau_byMediumCombinedIsolationDeltaBetaCorr3Hits[leptonIndex];
  otree->byTightCombinedIsolationDeltaBetaCorr3Hits_1 = analysisTree->tau_byTightCombinedIsolationDeltaBetaCorr3Hits[leptonIndex];
  otree->againstMuonLoose3_1 = analysisTree->tau_againstMuonLoose3[leptonIndex];
  otree->againstMuonTight3_1 = analysisTree->tau_againstMuonTight3[leptonIndex];
  otree->againstElectronLooseMVA6_1 = analysisTree->tau_againstElectronLooseMVA6[leptonIndex];
  otree->againstElectronMediumMVA6_1 = analysisTree->tau_againstElectronMediumMVA6[leptonIndex];
  otree->againstElectronTightMVA6_1 = analysisTree->tau_againstElectronTightMVA6[leptonIndex];
  otree->againstElectronVLooseMVA6_1 = analysisTree->tau_againstElectronVLooseMVA6[leptonIndex];
  otree->againstElectronVTightMVA6_1 = analysisTree->tau_againstElectronVTightMVA6[leptonIndex];

  otree->byVLooseIsolationMVArun2017v2DBoldDMwLT2017_1 = analysisTree->tau_byVLooseIsolationMVArun2017v2DBoldDMwLT2017[leptonIndex];
  otree->byLooseIsolationMVArun2017v2DBoldDMwLT2017_1 = analysisTree->tau_byLooseIsolationMVArun2017v2DBoldDMwLT2017[leptonIndex];
  otree->byMediumIsolationMVArun2017v2DBoldDMwLT2017_1 = analysisTree->tau_byMediumIsolationMVArun2017v2DBoldDMwLT2017[leptonIndex];
  otree->byTightIsolationMVArun2017v2DBoldDMwLT2017_1 = analysisTree->tau_byTightIsolationMVArun2017v2DBoldDMwLT2017[leptonIndex];
  otree->byVTightIsolationMVArun2017v2DBoldDMwLT2017_1 = analysisTree->tau_byVTightIsolationMVArun2017v2DBoldDMwLT2017[leptonIndex];
  otree->byVVTightIsolationMVArun2017v2DBoldDMwLT2017_1 = analysisTree->tau_byVVTightIsolationMVArun2017v2DBoldDMwLT2017[leptonIndex];
  otree->byIsolationMVArun2017v2DBoldDMwLTraw2017_1 = analysisTree->tau_byIsolationMVArun2017v2DBoldDMwLTraw2017[leptonIndex];
  otree->byIsolationMVA3oldDMwLTraw_1 = analysisTree->tau_byIsolationMVArun2017v2DBoldDMwLTraw2017[leptonIndex];
  otree->chargedIsoPtSum_1 = analysisTree->tau_chargedIsoPtSum[leptonIndex];
  otree->neutralIsoPtSum_1 = analysisTree->tau_neutralIsoPtSum[leptonIndex];
  otree->puCorrPtSum_1 = analysisTree->tau_puCorrPtSum[leptonIndex];
  //otree->isolationGammaCands_size_1 = analysisTree->tau_isolationGammaCands_size[leptonIndex];
  //otree->signalGammaCands_size_1 = analysisTree->tau_signalGammaCands_size[leptonIndex];

  otree->correction_againstElectronVLooseMVA6_1 = 1;
  //otree->correction_againstElectronLooseMVA6_1 = 1;
  //otree->correction_againstElectronMediumMVA6_1 = 1;
  otree->correction_againstElectronTightMVA6_1 = 1;
  //otree->correction_againstElectronVTightMVA6_1 = 1;
  otree->correction_againstMuonLoose3_1 = 1;
  otree->correction_againstMuonTight3_1 = 1;

  if (analysisTree->tau_genmatch[leptonIndex]==1){
    if (abs(otree->eta_1)<1.460){
      otree->correction_againstElectronVLooseMVA6_1 = 1.09;
      //		otree->correction_againstElectronLooseMVA6_1 = 1.17;
      //		otree->correction_againstElectronMediumMVA6_1 = 1.40;
      otree->correction_againstElectronTightMVA6_1 = 1.80;
      //		otree->correction_againstElectronVTightMVA6_1 = 1.96;
    }else if (abs(otree->eta_1)>1.558){
      otree->correction_againstElectronVLooseMVA6_1 = 1.19;
      //		otree->correction_againstElectronLooseMVA6_1 = 1.25;
      //		otree->correction_againstElectronMediumMVA6_1 = 1.21;
      otree->correction_againstElectronTightMVA6_1 = 1.53;
      //		otree->correction_againstElectronVTightMVA6_1 = 1.66;
    }
  }
  
  if (analysisTree->tau_genmatch[leptonIndex]==2){
    if ((0 < abs(otree->eta_1))&&(abs(otree->eta_1)<= 0.4)){
      otree->correction_againstMuonLoose3_1 = 1.06;
      otree->correction_againstMuonTight3_1 = 1.17;
    }else if ((0.4 < abs(otree->eta_1))&&(abs(otree->eta_1) <= 0.8)){
      otree->correction_againstMuonLoose3_1 = 1.02;
      otree->correction_againstMuonTight3_1 = 1.14;
    }else if ((0.8 < abs(otree->eta_1))&&(abs(otree->eta_1) <= 1.2)){
      otree->correction_againstMuonLoose3_1 = 1.10;
      otree->correction_againstMuonTight3_1 = 1.14;
    }else if ((1.2 < abs(otree->eta_1))&&(abs(otree->eta_1) <= 1.7)){
      otree->correction_againstMuonLoose3_1 = 1.03;
      otree->correction_againstMuonTight3_1 = 0.93;
    }else if ((1.7 < abs(otree->eta_1))&&(abs(otree->eta_1) <= 2.3)){
      otree->correction_againstMuonLoose3_1 = 1.94;
      otree->correction_againstMuonTight3_1 = 1.61;
    }
  }
  ///////////////////////////////////////////////////////NEW
  otree->efficiency_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017_1 = 1;
  //  otree->efficiency_byLooseIsolationMVArun2017v2DBoldDMwLT2017_1 = 1;
  //  otree->efficiency_byMediumIsolationMVArun2017v2DBoldDMwLT2017_1 = 1;
  otree->efficiency_byTightIsolationMVArun2017v2DBoldDMwLT2017_1 = 1;
  //  otree->efficiency_byVTightIsolationMVArun2017v2DBoldDMwLT2017_1 = 1;
  //  otree->efficiency_byVVTightIsolationMVArun2017v2DBoldDMwLT2017_1 = 1;
  if (analysisTree->tau_genmatch[leptonIndex]==5){
    otree->efficiency_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017_1 = 0.88;
    //  otree->efficiency_byLooseIsolationMVArun2017v2DBoldDMwLT2017_1 = 0.89;
    //  otree->efficiency_byMediumIsolationMVArun2017v2DBoldDMwLT2017_1 = 0.89;
    otree->efficiency_byTightIsolationMVArun2017v2DBoldDMwLT2017_1 = 0.89;
    //  otree->efficiency_byVTightIsolationMVArun2017v2DBoldDMwLT2017_1 = 0.86;
    //  otree->efficiency_byVVTightIsolationMVArun2017v2DBoldDMwLT2017_1 = 0.84;
  }
  
}
////FILLING FUNCTIONS//////

void FillVertices(const AC1B *analysisTree, Synch17Tree *otree, const bool isData, int leptonIndex, int tauIndex, TString channel){

  otree->RecoVertexX = analysisTree->primvertex_x;
  otree->RecoVertexY = analysisTree->primvertex_y;
  otree->RecoVertexZ = analysisTree->primvertex_z;
  
  bool is_refitted_PV_with_BS = true;
  TVector3 vertex_refitted_BS = get_refitted_PV_with_BS(analysisTree, leptonIndex, tauIndex,"tt",is_refitted_PV_with_BS);
  otree->pvx = vertex_refitted_BS.X();
  otree->pvy = vertex_refitted_BS.Y();
  otree->pvz = vertex_refitted_BS.Z();
  otree->is_refitted_PV_with_BS = is_refitted_PV_with_BS;

  otree->pvx_bs = analysisTree->primvertexwithbs_x;
  otree->pvy_bs = analysisTree->primvertexwithbs_y;
  otree->pvz_bs = analysisTree->primvertexwithbs_z;

  if(!isData){
    for (unsigned int igen = 0; igen < analysisTree->genparticles_count; ++igen) {

      //here fill the generator vertices to have the gen information present in tree PER GOOD RECO EVENT
      //Note: we may want to add constraint that the W and Z are prompt. If we remove these, may get in trouble with a DY or W MC sample..

      if ((analysisTree->genparticles_pdgid[igen] == 23 || analysisTree->genparticles_pdgid[igen] == 24 ||
	   analysisTree->genparticles_pdgid[igen] == 25 || analysisTree->genparticles_pdgid[igen] == 35 || analysisTree->genparticles_pdgid[igen] == 36) && 
	  analysisTree->genparticles_isLastCopy[igen] == 1 && analysisTree->genparticles_isPrompt[igen] == 1) {
        otree->GenVertexX = analysisTree->genparticles_vx[igen];
        otree->GenVertexY = analysisTree->genparticles_vy[igen];
        otree->GenVertexZ = analysisTree->genparticles_vz[igen];
        break;
      }
    }
  }
  else {//if it is data, fill with something recognisable nonsensible
    otree->GenVertexX = -9999;
    otree->GenVertexY = -9999;
    otree->GenVertexZ = -9999;
  }
}
void FillGenTree(const AC1B *analysisTree, Synch17GenTree *gentree){
  int ntaus=analysisTree->gentau_count;
  int npart=analysisTree->genparticles_count;
  int leptonid=15;
  TLorentzVector Tau1,Tau2,Tau;
  TLorentzVector Lepton;
  TLorentzVector lvector;
  int tauHIndex=-1;
  int tauLIndex=-1;
  int LeadingtauIndex=-1;
  int TrailingtauIndex=-1;
  double taumaxpt=-1;
  
  for(int itau=0;itau<ntaus;itau++){
    if(analysisTree->gentau_isLastCopy[itau]==1&&analysisTree->gentau_isPrompt[itau]==1){
      if(analysisTree->gentau_visible_pt[itau]>=taumaxpt) {
	LeadingtauIndex=itau; 
	taumaxpt=analysisTree->gentau_visible_pt[itau];
      }
    }
  }

  taumaxpt=-1; 
  for(int itau=0;itau<ntaus;itau++){
    if(analysisTree->gentau_isLastCopy[itau]==1&&analysisTree->gentau_isPrompt[itau]==1&&itau!=LeadingtauIndex){
      if(analysisTree->gentau_visible_pt[itau]>=taumaxpt) {
	TrailingtauIndex=itau; 
	taumaxpt=analysisTree->gentau_visible_pt[itau];
      }
    }
  }

  TLorentzVector genTauVis1; genTauVis1.SetXYZT(0,0,0,0);
  TLorentzVector genTauVis2; genTauVis2.SetXYZT(0,0,0,0);
  gentree->decaymode_1 = -1;
  gentree->decaymode_2 = -1;
  if (LeadingtauIndex>-1) {
    genTauVis1.SetXYZT(analysisTree->gentau_visible_px[LeadingtauIndex],
		       analysisTree->gentau_visible_py[LeadingtauIndex],
		       analysisTree->gentau_visible_pz[LeadingtauIndex],
		       analysisTree->gentau_visible_e[LeadingtauIndex]);
    gentree->decaymode_1 = analysisTree->gentau_decayMode[LeadingtauIndex];
  }
  if (TrailingtauIndex>-1) {
    genTauVis2.SetXYZT(analysisTree->gentau_visible_px[TrailingtauIndex],
		       analysisTree->gentau_visible_py[TrailingtauIndex],
		       analysisTree->gentau_visible_pz[TrailingtauIndex],
		       analysisTree->gentau_visible_e[TrailingtauIndex]);
    gentree->decaymode_2 = analysisTree->gentau_decayMode[TrailingtauIndex];
  }
  gentree->pt_1 = genTauVis1.Pt();
  gentree->eta_1 = genTauVis1.Eta();
  gentree->phi_1 = genTauVis1.Phi();

  gentree->pt_2 = genTauVis2.Pt();
  gentree->eta_2 = genTauVis2.Eta();
  gentree->phi_2 = genTauVis2.Phi();

  double dR;
  // double dR1 = deltaR(otree->eta_1,otree->phi_1,gentree->eta_1,gentree->phi_1);
  // double dR2 = deltaR(otree->eta_2,otree->phi_2,gentree->eta_2,gentree->phi_2);
  const double dRcut=0.3;
  //TLorentzVector genHP4;genHP4.clear();
  for(int ipart=0;ipart<npart;ipart++){
    if((abs(analysisTree->genparticles_pdgid[ipart])==25||
	abs(analysisTree->genparticles_pdgid[ipart])==35||
	abs(analysisTree->genparticles_pdgid[ipart])==36)&&
       analysisTree->genparticles_isLastCopy[ipart]==1){
      TLorentzVector Higgs;      
      Higgs.SetPxPyPzE(analysisTree->genparticles_px[ipart],
		       analysisTree->genparticles_py[ipart],
		       analysisTree->genparticles_pz[ipart],
		       analysisTree->genparticles_e[ipart]);
      gentree->Higgs_pt=Higgs.Pt();
      gentree->Higgs_eta=Higgs.Eta();
      gentree->Higgs_phi=Higgs.Phi();
      gentree->Higgs_mass=Higgs.M();
    }
  }

  //ROOT::Math::Boost boost_to_rf_h(genHP4_.BoostToCM());
  if (LeadingtauIndex>-1&&TrailingtauIndex>-1){
    //if(dR1<0.4 && dR2<0.4 ){
      gen_acott(analysisTree,gentree,LeadingtauIndex,TrailingtauIndex);
      gentree->sm_htt125 = analysisTree->TauSpinnerWeight[0];
      gentree->mm_htt125 = analysisTree->TauSpinnerWeight[1];
      gentree->ps_htt125 = analysisTree->TauSpinnerWeight[2];
      //}
  }
  
  gentree->a1polarization_1=gen_A1Polarization(analysisTree,LeadingtauIndex);
  gentree->a1polarization_2=gen_A1Polarization(analysisTree,TrailingtauIndex);

  //here fill the generator vertices to have the information present in tree
  //Note: we may want to add constraint that the W and Z are prompt. If we remove these, may get in trouble with a DY or W MC sample..


  for (unsigned int igen=0; igen<analysisTree->genparticles_count; ++igen) {
    if ((analysisTree->genparticles_pdgid[igen]==23||analysisTree->genparticles_pdgid[igen]==24||
	 analysisTree->genparticles_pdgid[igen]==25||analysisTree->genparticles_pdgid[igen]==35||analysisTree->genparticles_pdgid[igen]==36)&&analysisTree->genparticles_isLastCopy[igen]==1&&analysisTree->genparticles_isPrompt[igen]==1) {
      gentree->VertexX=analysisTree->genparticles_vx[igen];
      gentree->VertexY=analysisTree->genparticles_vy[igen];
      gentree->VertexZ=analysisTree->genparticles_vz[igen];
      break;
    }
  }
};
void initializeGenTree(Synch17GenTree *gentree){
  gentree->Higgs_pt=-9999;
  gentree->Higgs_eta=-9999;
  gentree->Higgs_phi=-9999;
  gentree->Higgs_mass=-9999;
  gentree->pt_1=-9999;
  gentree->eta_1=-9999;
  gentree->phi_1=-9999;
  gentree->pt_2=-9999;
  gentree->eta_2=-9999;
  gentree->phi_2=-9999;
  gentree->acotautau_00 = -9999;
  gentree->acotautau_10 = -9999;
  gentree->acotautau_01 = -9999;
  gentree->acotautau_11 = -9999;
  gentree->acotautau_02=-9999;
  gentree->acotautau_20=-9999;
  gentree->acotautau_12=-9999;
  gentree->acotautau_21=-9999;
  gentree->acotautau_22=-9999;

  //Merijn added the angle psi, currently for debugging purpose. Later may extend to 3-prong..
  gentree->acotautauPsi_00=-9999;
  gentree->acotautauPsi_01=-9999;
  gentree->acotautauPsi_10=-9999;
  gentree->acotautauPsi_11=-9999;

  //init the new vertex variables to something nonsensible also
  gentree->VertexX=-9999;
  gentree->VertexY=-9999;
  gentree->VertexZ=-99999;

  gentree->VxConstitTau1=-9999;
  gentree->VyConstitTau1=-9999;
  gentree->VzConstitTau1=-9999;
  
  gentree->VxConstitTau2=-9999;
  gentree->VyConstitTau2=-9999;
  gentree->VzConstitTau2=-9999;

  //Merijn add initialiser for
  gentree->chconst_1_pt=-9999;
  gentree->chconst_1_eta=-9999;
  gentree->chconst_1_phi=-9999;
  
  gentree->chconst_2_pt=-9999;
  gentree->chconst_2_eta=-9999;
  gentree->chconst_2_phi=-9999;
  gentree->alphaminus=-9999;

};
float getEmbeddedWeight(const AC1B *analysisTree, RooWorkspace * wEm) {

  std::vector<TLorentzVector> taus; taus.clear();
  float emWeight = 1;
  for (unsigned int igentau = 0; igentau < analysisTree->gentau_count; ++igentau) {
    TLorentzVector tauLV; tauLV.SetXYZT(analysisTree->gentau_px[igentau], 
					analysisTree->gentau_py[igentau],
					analysisTree->gentau_pz[igentau],
					analysisTree->gentau_e[igentau]);
    if (analysisTree->gentau_isPrompt[igentau]&&analysisTree->gentau_isFirstCopy[igentau]) {
      taus.push_back(tauLV);
    }
  }

  //  std::cout << "n taus = " << taus.size() << "  :  wEm = " << wEm << std::endl;

  if (taus.size() == 2) {
    double gt1_pt  = taus[0].Pt();
    double gt1_eta = taus[0].Eta();
    double gt2_pt  = taus[1].Pt();
    double gt2_eta = taus[1].Eta();
    wEm->var("gt_pt")->setVal(gt1_pt);
    wEm->var("gt_eta")->setVal(gt1_eta);
    double id1_embed = wEm->function("m_sel_id_ic_ratio")->getVal();
    wEm->var("gt_pt")->setVal(gt2_pt);
    wEm->var("gt_eta")->setVal(gt2_eta);
    double id2_embed = wEm->function("m_sel_id_ic_ratio")->getVal();
    wEm->var("gt1_pt")->setVal(gt1_pt);
    wEm->var("gt2_pt")->setVal(gt2_pt);
    wEm->var("gt1_eta")->setVal(gt1_eta);
    wEm->var("gt2_eta")->setVal(gt2_eta);
    double trg_emb = wEm->function("m_sel_trg_ic_ratio")->getVal();
    emWeight = id1_embed * id2_embed * trg_emb;
  }

  return emWeight;

}

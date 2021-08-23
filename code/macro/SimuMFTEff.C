#if !defined(__CLING__) || defined(__ROOTCLING__)

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THStack.h"
#include "SimulationDataFormat/MCTrack.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "SimulationDataFormat/MCEventHeader.h"
#include "ITSMFTSimulation/Hit.h"
#include "ITSMFTReconstruction/ChipMappingMFT.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsMFT/TrackMFT.h"

#endif

void SetStyle(Bool_t graypalette=kFALSE);

bool DEBUG_VERBOSE = false;

void SimuMFTEff(Double_t pMax = 11.0,
		       Double_t pMin = 0.0,
		       Double_t etaMin = -5,
		       Double_t etaMax = 0.,
		       const Char_t *kineFileName = "o2sim_Kine.root",
		       const Char_t *hitsFileName = "o2sim_HitsMFT.root",
		       const Char_t *clsFileName = "mftclusters.root",
		       const Char_t *trkFileName = "mfttracks.root")
{
  //SetStyle();
  using o2::itsmft::Hit;
  using o2::itsmft::CompClusterExt;
  using o2::MCTrack;
  
  using eventFoundTracks = std::vector<bool>;
  vector<eventFoundTracks> allFoundTracksLTF, allFoundTracksCA; // True for reconstructed tracks - one vector of bool per event
  using trackHasHitsInMFTDisks = std::array<bool,5>; // Disks with hits from a MFT track

  o2::itsmft::ChipMappingMFT mftChipMapper;
 
  std::unique_ptr<TH1F> MCTrackspT = std::make_unique<TH1F> ("MCTracks_Pt", "MCTracks_Pt", 200, pMin, 0.2 * pMax);
  MCTrackspT->GetXaxis()->SetTitle("#it{p}_{T}");
  MCTrackspT->GetYaxis()->SetTitle("N_{entries}");
  std::unique_ptr<TH1F> MCTracksp = std::make_unique<TH1F> ("MCTracks_P", "MCTracks_P", 200, pMin, pMax);
  MCTracksp->GetXaxis()->SetTitle("Total #it{p}");
  MCTracksp->GetYaxis()->SetTitle("N_{entries}");
  std::unique_ptr<TH1F> MCTrackEta = std::make_unique<TH1F> ("MCTracks_eta", "MCTracks_eta", 200, etaMin, etaMax);
  MCTrackEta->GetXaxis()->SetTitle("#eta");
  std::unique_ptr<TH1F> MCTrackVz = std::make_unique<TH1F> ("MCTracks_vz", "MCTracks_vz", 200, -15, 16);
  MCTrackVz->GetXaxis()->SetTitle("V_{z}");
  MCTrackVz->GetYaxis()->SetTitle("N_{entries}");

  std::unique_ptr<TH1F> MFTTrackspT = std::make_unique<TH1F> ("MFTTracks_Pt", "MFTTracks_Pt", 200, pMin, 0.2 * pMax);
  MFTTrackspT->GetXaxis()->SetTitle("#it{p}_{T}");
  MFTTrackspT->GetYaxis()->SetTitle("N_{entries}");
  std::unique_ptr<TH1F> MFTTracksp = std::make_unique<TH1F> ("MFTTracks_P", "MFTTracks_P", 200, pMin, pMax);
  MFTTracksp->GetXaxis()->SetTitle("Total #it{p}");
  MFTTracksp->GetYaxis()->SetTitle("N_{entries}");
  std::unique_ptr<TH1F> MFTTrackEta = std::make_unique<TH1F> ("MFTTracks_eta", "MFTTracks_eta", 200, etaMin, etaMax);
  MFTTrackEta->GetXaxis()->SetTitle("#eta");
  MFTTrackEta->GetYaxis()->SetTitle("N_{entries}");
  std::unique_ptr<TH1F> MFTTrackVz = std::make_unique<TH1F> ("MFTTracks_vz", "MFTTracks_vz", 200, -15, 16);
  MFTTrackVz->GetXaxis()->SetTitle("V_{z}");
  MFTTrackVz->GetYaxis()->SetTitle("N_{entries}");

  std::unique_ptr<TH1F> LTFTrackspT = std::make_unique<TH1F> ("LTFTracks_Pt", "LTFTracks_Pt", 200, pMin, 0.2 * pMax);
  LTFTrackspT->GetXaxis()->SetTitle("#it{p}_{T}");
  LTFTrackspT->GetYaxis()->SetTitle("N_{entries}");
  std::unique_ptr<TH1F> LTFTracksp = std::make_unique<TH1F> ("LTFTracks_P", "LTFTracks_P", 200, pMin, pMax);
  LTFTracksp->GetXaxis()->SetTitle("Total #it{p}");
  LTFTracksp->GetYaxis()->SetTitle("N_{entries}");
  std::unique_ptr<TH1F> LTFTracksEta = std::make_unique<TH1F> ("LTFTracks_eta", "LTFTracks_eta", 200, etaMin, etaMax);
  LTFTracksEta->GetXaxis()->SetTitle("#eta");
  LTFTracksEta->GetYaxis()->SetTitle("N_{entries}");
  std::unique_ptr<TH1F> LTFTracksVz = std::make_unique<TH1F> ("LTFTracks_Vz", "LTFTracks_Vz", 200, -15, 16);
  LTFTracksVz->GetXaxis()->SetTitle("V_{z}");
  LTFTracksVz->GetYaxis()->SetTitle("N_{entries}");

  std::unique_ptr<TH1F> CATrackspT = std::make_unique<TH1F> ("CATracks_Pt", "CATracks_Pt", 200, pMin, 0.2 * pMax);
  CATrackspT->GetXaxis()->SetTitle("#it{p}_{T}");
  CATrackspT->GetYaxis()->SetTitle("N_{entries}");
  std::unique_ptr<TH1F> CATracksp = std::make_unique<TH1F> ("CATracks_P", "CATracks_P", 200, pMin, pMax);
  CATracksp->GetXaxis()->SetTitle("Total #it{p}");
  CATracksp->GetYaxis()->SetTitle("N_{entries}");
  std::unique_ptr<TH1F> CATracksEta = std::make_unique<TH1F> ("CATracks_eta", "CATracks_eta", 200, etaMin, etaMax);
  CATracksEta->GetXaxis()->SetTitle("#eta");
  CATracksEta->GetYaxis()->SetTitle("N_{entries}");
  std::unique_ptr<TH1F> CATracksVz = std::make_unique<TH1F> ("CATracks_Vz", "CATracks_Vz", 200, -15, 16);
  CATracksVz->GetXaxis()->SetTitle("V_{z}");
  CATracksVz->GetYaxis()->SetTitle("N_{entries}");

  std::unique_ptr<TH1F> MCTracksEta5 = std::make_unique<TH1F> ("MCTracks_Vz5_eta", "-5 cm < V_{z} < 5 cm", 200, etaMin, etaMax);
  MCTracksEta5->GetXaxis()->SetTitle("#eta");
  MCTracksEta5->GetYaxis()->SetTitle("N_{entries}");
  std::unique_ptr<TH1F> MCTracksEtaOut5 = std::make_unique<TH1F> ("MCTracks_VzOut5_eta", "|V_{z}|>5 cm", 200, etaMin, etaMax);
  MCTracksEtaOut5->GetXaxis()->SetTitle("#eta");
  MCTracksEtaOut5->GetYaxis()->SetTitle("N_{entries}");

  std::unique_ptr<TH1F> MFTTracksEta5 = std::make_unique<TH1F> ("MFTTracks_Vz5_eta", "-5 cm < V_{z} < 5 cm", 200, etaMin, etaMax);
  MFTTracksEta5->GetXaxis()->SetTitle("#eta");
  MFTTracksEta5->GetYaxis()->SetTitle("N_{entries}");
  std::unique_ptr<TH1F> MFTTracksEtaOut5 = std::make_unique<TH1F> ("MFTTracks_VzOut5_eta", "|V_{z}| >5 cm", 200, etaMin, etaMax);
  MFTTracksEtaOut5->GetXaxis()->SetTitle("#eta");
  MFTTracksEtaOut5->GetYaxis()->SetTitle("N_{entries}");


  std::unique_ptr<TH1F> MCTracksp5 = std::make_unique<TH1F> ("MCTracks_Vz5_p", "-5 cm < V_{z} < 5 cm", 200, pMin, pMax);
  MCTracksp5->GetXaxis()->SetTitle("Total #it{p}");
  MCTracksp5->GetYaxis()->SetTitle("N_{entries}");

  std::unique_ptr<TH1F> MFTTracksp5 = std::make_unique<TH1F> ("MFTTracks_Vz5_p", "-5 cm < V_{z} < 5 cm", 200, pMin, pMax);
  MFTTracksp5->GetXaxis()->SetTitle("Total #it{p}");
  MFTTracksp5->GetYaxis()->SetTitle("N_{entries}");

  std::unique_ptr<TH1I> Trackablility = std::make_unique<TH1I> ("Trackablility", "In how many disks the tracks has hits", 6, 0, 6);
  Trackablility->GetXaxis()->SetTitle("Number of disks");
  Trackablility->GetYaxis()->SetTitle("N_{entries}");

  //Histos for Missed (missed tracks that could be tracked)
  std::unique_ptr<TH1F> MissedlepT = std::make_unique<TH1F> ("MissedTracks_Pt", "MissedTracks_Pt", 200, pMin, 0.2 * pMax);
  std::unique_ptr<TH1F> Missedp = std::make_unique<TH1F> ("MissedTracks_P", "MissedTracks_P", 200, pMin, pMax);
  std::unique_ptr<TH1F> MissedEta = std::make_unique<TH1F> ("MissedTracks_eta", "MissedTracks_eta", 200, etaMin, etaMax);
  std::unique_ptr<TH1F> MissedVz = std::make_unique<TH1F> ("MissedTracks_Vz", "MissedTracks_Vz", 200, -15, 16);

  //Histos for Trackables
  std::unique_ptr<TH1F> TrackablepT = std::make_unique<TH1F> ("TrackablesTracks_Pt", "TrackablesTracks_Pt", 200, pMin, 0.2 * pMax);
  TrackablepT->GetXaxis()->SetTitle("#it{p}_{T}");
  TrackablepT->GetYaxis()->SetTitle("N_{entries}");
  std::unique_ptr<TH1F> Trackablep = std::make_unique<TH1F> ("TrackablesTracks_P", "TrackablesTracks_P", 200, pMin, pMax);
  Trackablep->GetXaxis()->SetTitle("Total #it{p}");
  Trackablep->GetYaxis()->SetTitle("N_{entries}");
  std::unique_ptr<TH1F> TrackableEta = std::make_unique<TH1F> ("TrackablesTracks_eta", "TrackablesTracks_eta", 200, etaMin, etaMax);
  TrackableEta->GetXaxis()->SetTitle("#eta");
  TrackableEta->GetYaxis()->SetTitle("N_{entries}");
  std::unique_ptr<TH1F> TrackableVz = std::make_unique<TH1F> ("TrackablesTracks_Vz", "TrackablesTracks_Vz", 200, -15, 16);
  TrackableVz->GetXaxis()->SetTitle("V_{z}");
  TrackableVz->GetYaxis()->SetTitle("N_{entries}");

  //2D Histos
  std::unique_ptr<TH2F> MFTTrackedEtaZ = std::make_unique<TH2F> ("MFT_Tracked_eta_z", "ReconstructedTracks", 31, -15, 16, 25, etaMin, etaMax);
  MFTTrackedEtaZ->GetXaxis()->SetTitle("V_{z} [cm]");
  MFTTrackedEtaZ->GetYaxis()->SetTitle("N_{Reconstructed}");
  std::unique_ptr<TH2F> MFTTrackablesEtaZ = std::make_unique<TH2F> ("MFT_Trackables_eta_z", "MFTTrackables:", 31, -15, 16, 25, etaMin, etaMax);
  MFTTrackablesEtaZ->GetXaxis()->SetTitle("V_{z} [cm]");
  MFTTrackablesEtaZ->GetYaxis()->SetTitle("N_{Trackables}");
  std::unique_ptr<TH2F> MCTracksEtaZ = std::make_unique<TH2F> ("MCTracks_eta_z", "MC Tracks: Pseudorapidity vs V_{z}", 31, -15, 16, 25, etaMin, etaMax);
  MCTracksEtaZ->GetXaxis()->SetTitle("V_{z} [cm]");
  MCTracksEtaZ->GetYaxis()->SetTitle("N_{MCTracks}");
  
  // MC tracks
  TFile kineFile(kineFileName);
  TTree *kineTree = (TTree*)kineFile.Get("o2sim");
  std::vector<o2::MCTrack> mcTrkVec, *mcTrkVecP = &mcTrkVec;
  kineTree->SetBranchAddress("MCTrack",&mcTrkVecP);
  o2::dataformats::MCEventHeader* eventHeader = nullptr;
  kineTree->SetBranchAddress("MCEventHeader.", &eventHeader);

  // hits
  TFile hitsFile(hitsFileName);
  TTree* hitTree = (TTree*)hitsFile.Get("o2sim");
  std::vector<Hit> hitVec, *hitVecP = &hitVec;
  hitTree->SetBranchAddress("MFTHit", &hitVecP);

  Int_t nEvents = hitTree->GetEntries();
  printf("Number of events %d \n", nEvents);

  // clusters
  TFile clsFile(clsFileName);
  TTree *clsTree = (TTree*)clsFile.Get("o2sim");
  std::vector<CompClusterExt> clsVec, *clsVecP = &clsVec;
  clsTree->SetBranchAddress("MFTClusterComp", &clsVecP);
  o2::dataformats::MCTruthContainer<o2::MCCompLabel>* clsLabels = nullptr;
  if (clsTree->GetBranch("MFTClusterMCTruth")) {
    clsTree->SetBranchAddress("MFTClusterMCTruth", &clsLabels);
  } else {
    printf("No Monte-Carlo information in this file\n");
    return;
  }
  
  clsTree->GetEntry(0);
  Int_t nClusters = clsVec.size();
  printf("Number of clusters %d \n", nClusters);
  
  // tracks
  TFile trkFile(trkFileName);
  TTree *trkTree = (TTree*)trkFile.Get("o2sim");
  std::vector<o2::mft::TrackMFT> trackVec, *trackVecP = &trackVec;
  trkTree->SetBranchAddress("MFTTrack", &trackVecP);
  o2::dataformats::MCTruthContainer<o2::MCCompLabel>* trkLabels = nullptr;
  if (trkTree->GetBranch("MFTTrackMCTruth")) {
    trkTree->SetBranchAddress("MFTTrackMCTruth", &trkLabels);
  } else {
    printf("No Monte-Carlo information in this file\n");
    return;
  }
  std::vector<int> trackExtClsVec, *trackExtClsVecP = &trackExtClsVec;
  trkTree->SetBranchAddress("MFTTrackClusIdx", &trackExtClsVecP);

  trkTree->GetEntry(0);

  // MFTTrackerChecker.C
  // Resize vector to accomodate found status of all tracks in all events
  for (auto event = 0 ; event < nEvents ; event++) { 
    kineTree->GetEntry(event);
    hitTree->GetEntry(event);
    auto nTracks = eventHeader->getMCEventStats().getNKeptTracks();
    if(DEBUG_VERBOSE) std::cout << "Resizing allFoundTracks for event " << event <<  " with ntracks = " << nTracks << std::endl;
    eventFoundTracks tempFoundTracks(nTracks, false);
    allFoundTracksLTF.push_back(tempFoundTracks); // reserve size and initialize
    allFoundTracksCA.push_back(tempFoundTracks); // reserve size and initialize
  }

  // MFTTrackerChecker.C
  // Part 1: Quality of reconstructed MFT tracks
  //   - Loop over reconstructed tracks to identify clean and mixed/noise tracks
  //   - Clean tracks have at least 80% of its clusters from the same track
  //   - If track is not clear it is a mixed/noise track
  // Part 2: MC hits and tracks
  //   - Identify trackable tracks (clusters in at least 4 disks)
  //   - Identify successfully reconstructed tracks
  // Part 3: Calculate Efficiencies

  // Part 1: Quality of reconstructed MFT tracks

  Int_t nCleanTracksLTF = 0, nCleanTracksCA = 0, nInvalidTracksLTF = 0, nInvalidTracksCA = 0, nMFTTrackable = 0;

  int srcID, trkID, evnID;
  bool fake;
  Int_t iTrack = 0;
  for (auto &track : trackVec) {
    
    auto ncls = track.getNumberOfPoints();
    auto offset = track.getExternalClusterIndexOffset();
    std::map<Int_t, Int_t> trkIDs;
    for (int icls = 0; icls < ncls; ++icls) {
      auto clsEntry = trackExtClsVec[offset + icls];
      auto cluster = clsVec[clsEntry];
      auto& clsLabel = (clsLabels->getLabels(clsEntry))[0];
      clsLabel.get(trkID, evnID, srcID, fake);
      if (!clsLabel.isNoise()) {
        trkIDs[trkID] = trkIDs[trkID] + 1;
      }
    }
    
    Int_t thisEvnID = -1, thisSrcID = -1, thisTrkID = -1, thisEventIDLabel = -1, nPoints = track.getNumberOfPoints();
    for (int icls = 0; icls < ncls; ++icls) {
      auto clsEntry = trackExtClsVec[offset + icls];
      auto cluster = clsVec[clsEntry];
      auto& clsLabel = (clsLabels->getLabels(clsEntry))[0];
      clsLabel.get(trkID, evnID, srcID, fake);
      if (!clsLabel.isNoise()) {
        if (((Float_t)(trkIDs[trkID]) / (Float_t)(nPoints)) >= 0.8) { // Must have at least 80% of its clusters from the same MC Track
          thisTrkID = trkID;
          thisSrcID = srcID;
          thisEvnID = evnID;
          thisEventIDLabel = icls;
        }
      }
    }
    
    auto eventID = thisEvnID;
    
    if ((thisTrkID >= 0) & (thisTrkID != 0x7FFFFFF) & (eventID < nEvents)) { // If is good match and not noise ...
      if (!track.isCA()) {
        allFoundTracksLTF[eventID][thisTrkID] = true;
        nCleanTracksLTF++;
      } else {
        allFoundTracksCA[eventID][thisTrkID] = true;
        nCleanTracksCA++;
      }
    } else {
      if (!track.isCA()) {
        nInvalidTracksLTF++;
      } else {
        nInvalidTracksCA++;
      }
    }
    if(DEBUG_VERBOSE) std::cout << "This TrackLTF ID = " << thisTrkID << " from sourceID = " << thisSrcID << " in eventID = " << eventID << " nPoints = " << nPoints << " nCleanTracksLTF = " << nCleanTracksLTF << std::endl;
   
    ++iTrack;
  }

  // Part 2: MC hits and tracks

  for (Int_t event = 0; event < nEvents ; event++) {

    hitTree->GetEntry(event);
    kineTree->GetEntry(event);

    Int_t nMFTHits = hitVec.size();
    std::vector<trackHasHitsInMFTDisks> mcTrackHasHitsInMFTDisks(eventHeader->getMCEventStats().getNKeptTracks(), {0, 0, 0, 0, 0});

    if(DEBUG_VERBOSE) std::cout << "Loop over " << nMFTHits << " mfthits to identify trackable MFT tracks in event " <<  event << std::endl;

    for (Int_t n_hit = 0 ; n_hit < nMFTHits; n_hit++) { // Loop over mfthits to identify trackable tracks
      
      Hit* hitp = &(hitVec).at(n_hit);
      Int_t trkID = hitp->GetTrackID(); // ID of the tracks having given the hit
      Float_t z = hitp->GetZ(); // Z position of the hit => Identify MFT disk
      mcTrackHasHitsInMFTDisks.at(trkID)[mftChipMapper.chip2Layer(hitp->GetDetectorID()) / 2] = true;
    }
    
    for (Int_t trkID = 0 ; trkID < eventHeader->getMCEventStats().getNKeptTracks(); trkID++) { // Loop on MC tracks

      //fill MC histograms
      MCTrack* thisTrack =  &(mcTrkVec)[trkID];
      if (!thisTrack->isPrimary()) continue;
      
      auto z = thisTrack->GetStartVertexCoordinatesZ();
      auto pt = thisTrack->GetPt();
      auto p = thisTrack->GetP();
      auto eta = atanh (thisTrack->GetStartVertexMomentumZ()/p); // eta;
      MCTrackspT->Fill(pt);
      MCTracksp->Fill(p);
      MCTrackEta->Fill(eta);
      MCTrackVz->Fill(z);
      MCTracksEtaZ->Fill(z,eta);
      if( (z > -5) & (z < 5) ) {
        MCTracksEta5->Fill(eta);
        MCTracksp5->Fill(p);
      }
      else{
        MCTracksEtaOut5->Fill(eta);
      }
      // Count disks "touched" by the track
      int nMFTDisksHasHits = 0;
      for(auto disk: {0, 1, 2, 3, 4}) {
        nMFTDisksHasHits+= int(mcTrackHasHitsInMFTDisks[trkID][disk]);
      }
      Trackablility->Fill(nMFTDisksHasHits);
      //std::cout << "nMFTDisksHasHits = " << nMFTDisksHasHits; // << std::endl;
      if (nMFTDisksHasHits >=4) {   //Track is trackable if has left hits on at least 4 disks
        nMFTTrackable++;
        MFTTrackablesEtaZ->Fill(z,eta);
        TrackablepT->Fill(pt);
        Trackablep->Fill(p);
        TrackableEta->Fill(eta);
        TrackableVz->Fill(z);
        bool wasFound = allFoundTracksLTF[event][trkID] | allFoundTracksCA[event][trkID];
        if(wasFound) {	  
          MFTTrackspT->Fill(pt);
          MFTTracksp->Fill(p);
          MFTTrackEta->Fill(eta);
          MFTTrackVz->Fill(z);
          if( (z > -5) & (z < 5) ) {
            MFTTracksEta5->Fill(eta);
            MFTTracksp5->Fill(p);
          }
          else{
            MFTTracksEtaOut5->Fill(eta);
          }
          if(allFoundTracksLTF[event][trkID]) {
            LTFTrackspT->Fill(thisTrack->GetPt());
            LTFTracksp->Fill(thisTrack->GetP());
            LTFTracksEta->Fill(eta);
            LTFTracksVz->Fill(thisTrack->GetStartVertexCoordinatesZ());
            MFTTrackedEtaZ->Fill(thisTrack->GetStartVertexCoordinatesZ(),eta);
          }
          if(allFoundTracksCA[event][trkID]) {
            CATrackspT->Fill(thisTrack->GetPt());
            CATracksp->Fill(thisTrack->GetP());
            CATracksEta->Fill(eta);
            CATracksVz->Fill(thisTrack->GetStartVertexCoordinatesZ());
          }
        }
      } else {  // Fill histograms for Missed Tracks
        MissedlepT->Fill(thisTrack->GetPt());
        Missedp->Fill(thisTrack->GetP());
        MissedEta->Fill(eta);
        MissedVz->Fill(z);
      }
    } // end loop on tracks   
  } // end loop on events

  // Part 3: Calculate Efficiencies
  //std::cout << "Building efficiencies histos..." << std::endl;
  TH1F MFTEfficiencypT = (*MFTTrackspT)/ (*TrackablepT);
  TH1F MFTEfficiencyp = (*MFTTracksp) / (*Trackablep);
  TH1F MFTEfficiencyEta = (*MFTTrackEta) / (*TrackableEta);
  TH1F MFTEfficiencyVz = (*MFTTrackVz) / (*TrackableVz);
  TH2F MFTTrackerEfficiency = (*MFTTrackedEtaZ) / (*MFTTrackablesEtaZ);
  TH2F MFTEfficiency2D = (*MFTTrackedEtaZ) / (*MCTracksEtaZ);
  TH2F MFTAcceptance = (*MFTTrackablesEtaZ) / (*MCTracksEtaZ);
  
  TH1F MFTEffsEta5 = (*MFTTracksEta5)/(*MCTracksEta5); 
  TH1F MFTEffsEtaOut5 = (*MFTTracksEtaOut5)/(*MCTracksEtaOut5); 
  TH1F MFTEffsp5 = (*MFTTracksp5)/(*MCTracksp5);

  TH1F MFTAlgoRatioPt = (*LTFTrackspT) /(*CATrackspT);
  TH1F MFTAlgoRatioEta = (*LTFTracksEta) /(*CATracksEta);
  TH1F MFTAlgoRatioP = (*LTFTracksp) /(*CATracksp);
  TH1F MFTAlgoRatioVz = (*LTFTracksVz) /(*CATracksVz);

  TH1F MFTAlgoPercLTFPt = (*LTFTrackspT) /(*CATrackspT + *LTFTrackspT);
  TH1F MFTAlgoPercLTFEta = (*LTFTracksEta) /(*CATracksEta + *LTFTracksEta);
  TH1F MFTAlgoPercLTFP = (*LTFTracksp) /(*CATracksp + *LTFTracksp);
  TH1F MFTAlgoPercLTFVz = (*LTFTracksVz) /(*CATracksVz + *LTFTracksVz);


  TH1F MFTLTFEfficiencypT = (*LTFTrackspT)/ (*TrackablepT);
  TH1F MFTLTFEfficiencyp = (*LTFTracksp) / (*Trackablep);
  TH1F MFTLTFEfficiencyEta = (*LTFTracksEta) / (*TrackableEta);
  TH1F MFTLTFEfficiencyVz = (*LTFTracksVz) / (*TrackableVz);

  // TH1F MFTAlgoRatioPt_5_5 = (*LTFTrackspT) /(*CATrackspT);
  // TH1F MFTAlgoRatioEta_5_5 = (*LTFTracksEta) /(*CATracksEta);
  // TH1F MFTAlgoRatioP_5_5 = (*LTFTracksp) /(*CATracksp);
  
  MFTEfficiencypT.SetNameTitle("MFT_Efficiency_pT", "MFT_Efficiency_pT");
  MFTEfficiencypT.GetXaxis()->SetTitle("#it{p}_{T}");
  MFTEfficiencypT.GetYaxis()->SetTitle("Efficiency");
  MFTEfficiencyp.SetNameTitle("MFT_Efficiency_p", "MFT_Efficiency_p");
  MFTEfficiencyp.GetXaxis()->SetTitle("#it{p}");
  MFTEfficiencyp.GetYaxis()->SetTitle("Efficiency");
  MFTEfficiencyEta.SetNameTitle("MFT_Efficiency_eta", "MFT_Efficiency_Pseudorapidity");
  MFTEfficiencyEta.GetXaxis()->SetTitle("#eta");
  MFTEfficiencyEta.GetYaxis()->SetTitle("Efficiency");
  MFTEfficiencyVz.SetNameTitle("MFT_Efficiency_Vz", "MFT_Efficiency_Vz");
  MFTEfficiencyVz.GetXaxis()->SetTitle("V_{z} [cm]");
  MFTEfficiencyVz.GetYaxis()->SetTitle("Efficiency");
  MFTTrackerEfficiency.SetNameTitle("MFT Tracker Efficiency", "MFT Tracker Efficiency");
  MFTTrackerEfficiency.GetXaxis()->SetTitle("V_{z} [cm]");
  MFTTrackerEfficiency.GetYaxis()->SetTitle("Efficiency");
  
  MFTEfficiency2D.SetNameTitle("MFT_Efficiency", "MFT_Efficiency");
  MFTEfficiency2D.GetXaxis()->SetTitle("Vertex PosZ [cm]");
  
  MFTAcceptance.SetNameTitle("MFT_Acceptance", "MFT_Acceptance");
  MFTAcceptance.GetXaxis()->SetTitle("V_{z} [cm]");
  
  MFTEffsEta5.SetNameTitle("MFT_Eta_Efficiency5_5", "-5 cm < z < 5 cm");
  MFTEffsEta5.GetXaxis()->SetTitle("#eta");
  
  MFTEffsp5.SetNameTitle("MFT_P_Efficiency5_5", "-5 cm < z < 5 cm");
  MFTEffsp5.GetXaxis()->SetTitle("P (GeV)");
  
  MFTEffsEtaOut5.SetNameTitle("MFT_Eta_EfficiencyOut5_5", "|z| > 5 cm");
  MFTEffsEtaOut5.GetXaxis()->SetTitle("#eta");

  MFTAlgoRatioPt.SetNameTitle("RatioLTFtoCA_pt","Ratio of tracks found by the LTF algo to the CA algo vs #it{p}_{T}");
  MFTAlgoRatioPt.GetXaxis()->SetTitle("#it{p}_{T}");
  MFTAlgoRatioPt.GetYaxis()->SetTitle("Ratio LTF/CA");
  MFTAlgoRatioP.SetNameTitle("RatioLTFtoCA_p","Ratio of tracks found by the LTF algo to the CA algo vs #it{p}");
  MFTAlgoRatioP.GetXaxis()->SetTitle("#it{p}");
  MFTAlgoRatioP.GetYaxis()->SetTitle("Ratio LTF/CA");
  MFTAlgoRatioEta.SetNameTitle("RatioLTFtoCA_eta","Ratio of tracks found by the LTF algo to the CA algo vs #eta");
  MFTAlgoRatioEta.GetXaxis()->SetTitle("#eta");
  MFTAlgoRatioEta.GetYaxis()->SetTitle("Ratio LTF/CA");
  MFTAlgoRatioVz.SetNameTitle("RatioLTFtoCA_Vz","Ratio of tracks found by the LTF algo to the CA algo vs V_{z}");
  MFTAlgoRatioVz.GetXaxis()->SetTitle("V_{z}");
  MFTAlgoRatioVz.GetYaxis()->SetTitle("Ratio LTF/CA");


  MFTLTFEfficiencypT.SetNameTitle("MFT_LTF_Efficiency_pT", "Efficiency of the LTF algo vs #it{p}_{T}");
  MFTLTFEfficiencypT.GetXaxis()->SetTitle("#it{p}_{T}");
  MFTLTFEfficiencypT.GetYaxis()->SetTitle("Efficiency");
  MFTLTFEfficiencyp.SetNameTitle("MFT_LTF_Efficiency_p", "Efficiency of the LTF algo vs #it{p}");
  MFTLTFEfficiencyp.GetXaxis()->SetTitle("#it{p}");
  MFTLTFEfficiencyp.GetYaxis()->SetTitle("Efficiency");
  MFTLTFEfficiencyEta.SetNameTitle("MFT_LTF_Efficiency_eta", "Efficiency of the LTF algo vs #eta");
  MFTLTFEfficiencyEta.GetXaxis()->SetTitle("#eta");
  MFTLTFEfficiencyEta.GetYaxis()->SetTitle("Efficiency");
  MFTLTFEfficiencyVz.SetNameTitle("MFT_LTF_Efficiency_Vz", "Efficiency of the LTF algo vs V_{z}");
  MFTLTFEfficiencyVz.GetXaxis()->SetTitle("V_{z} [cm]");
  MFTLTFEfficiencyVz.GetYaxis()->SetTitle("Efficiency");


  MFTAlgoPercLTFPt.SetNameTitle("PercLTF_pt","Percentage of tracks found by the LTF algo to the (LTF+CA) total vs #it{p}_{T}");
  MFTAlgoPercLTFPt.GetXaxis()->SetTitle("#it{p}_{T}");
  MFTAlgoPercLTFPt.GetYaxis()->SetTitle("LTF percentage");
  MFTAlgoPercLTFP.SetNameTitle("PercLTF_p","Percentage of tracks found by the LTF algo to the (LTF+CA) total vs #it{p}");
  MFTAlgoPercLTFP.GetXaxis()->SetTitle("#it{p}");
  MFTAlgoPercLTFP.GetYaxis()->SetTitle("LTF percentage");
  MFTAlgoPercLTFEta.SetNameTitle("PercLTF_eta","Percentage of tracks found by the LTF algo to the (LTF+CA) total vs #eta");
  MFTAlgoPercLTFEta.GetXaxis()->SetTitle("#eta");
  MFTAlgoPercLTFEta.GetYaxis()->SetTitle("LTF percentage");
  MFTAlgoPercLTFVz.SetNameTitle("PercLTF_Vz","Percentage of tracks found by the LTF algo to the (LTF+CA) total vs V_{z}");
  MFTAlgoPercLTFVz.GetXaxis()->SetTitle("V_{z}");
  MFTAlgoPercLTFVz.GetYaxis()->SetTitle("LTF percentage");

  // Write histograms to file
  //std::cout << "Writting histograms to file..." << std::endl;
  TFile outFile("MFTEffCheck.root","RECREATE");
  
  MCTrackspT->Write();
  MCTracksp->Write();
  MCTrackEta->Write();
  MCTrackVz->Write();
  
  MFTTrackspT->Write();
  MFTTracksp->Write();
  MFTTrackEta->Write();
  MFTTrackVz->Write();
  
  LTFTrackspT->Write();
  LTFTracksp->Write();
  LTFTracksEta->Write();
  LTFTracksVz->Write();
  
  CATrackspT->Write();
  CATracksp->Write();
  CATracksEta->Write();
  CATracksVz->Write();

  MFTAlgoRatioPt.Write();
  MFTAlgoRatioP.Write();
  MFTAlgoRatioEta.Write();
  MFTAlgoRatioVz.Write();

  MFTAlgoPercLTFPt.Write();
  MFTAlgoPercLTFP.Write();
  MFTAlgoPercLTFEta.Write();
  MFTAlgoPercLTFVz.Write();
  
  MFTEfficiencypT.Write();
  MFTEfficiencyp.Write();
  MFTEfficiencyEta.Write();
  MFTEfficiencyVz.Write();
  
  MFTLTFEfficiencypT.Write();
  MFTLTFEfficiencyp.Write();
  MFTLTFEfficiencyEta.Write();
  MFTLTFEfficiencyVz.Write();
  
  MCTracksEta5->Write();
  MCTracksEtaOut5->Write();
  MFTTracksEta5->Write();
  MFTTracksEtaOut5->Write();
  MCTracksp5->Write();
  MFTTracksp5->Write();
  MFTEffsEta5.Write();
  MFTEffsEtaOut5.Write();
  MFTEffsp5.Write();
  
  MissedlepT->Write();
  Missedp->Write();
  MissedEta->Write();
  MissedVz->Write();
  
  Trackablility->Write();
  
  TrackablepT->Write();
  Trackablep->Write();
  TrackableEta->Write();
  TrackableVz->Write();
  
  MCTracksEtaZ->SetOption("CONT4");
  MCTracksEtaZ->Write();
  
  MFTTrackablesEtaZ->SetOption("CONT4");
  MFTTrackablesEtaZ->Write();
  
  MFTTrackedEtaZ->SetOption("CONT4");
  MFTTrackedEtaZ->Write();
  
  MFTTrackerEfficiency.SetOption("CONT4");
  MFTTrackerEfficiency.Write();
  
  MFTEfficiency2D.SetOption("CONT4");
  MFTEfficiency2D.Write();
  
  MFTAcceptance.SetOption("CONT4");
  MFTAcceptance.Write();
  
  outFile.Close();
  
  Int_t totalRecoMFTTracks = nCleanTracksLTF + nCleanTracksCA + nInvalidTracksLTF + nInvalidTracksCA;
  std::cout << "---------------------------------------------------" << std::endl;
  
  std::cout << "Number of MFT trackables = " << nMFTTrackable << std::endl;
  std::cout << "Number of reconstructed MFT Tracks = " << totalRecoMFTTracks << std::endl;
  std::cout << "Number of clean MFT Tracks = " << nCleanTracksLTF + nCleanTracksCA << std::endl;
  std::cout << "Number of mixed MFT Tracks = " << nInvalidTracksLTF + nInvalidTracksCA << std::endl;
  std::cout << "---------------------------------------------------" << std::endl;
  std::cout << "nCleanTracksLTF = " << nCleanTracksLTF << std::endl;
  std::cout << "nCleanTracksCA = " << nCleanTracksCA << std::endl;
  std::cout << "nInvalidTracksLTF = " << nInvalidTracksLTF  << " (" << 100.f*nInvalidTracksLTF/(nCleanTracksLTF+nInvalidTracksLTF) << " %)" << std::endl;
  std::cout << "nInvalidTracksCA = " << nInvalidTracksCA << " (" << 100.f*nInvalidTracksCA/(nCleanTracksCA+nInvalidTracksCA) << " %)" << std::endl;
  std::cout << "---------------------------------------------------" << std::endl;

}
//__________________________________________________________________________
void SetStyle(Bool_t graypalette) {
  cout << "Setting style!" << endl;

  gStyle->Reset("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  if(graypalette) gStyle->SetPalette(8,0);
  else gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetPadColor(10);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetHistLineWidth(1);
  gStyle->SetHistLineColor(kRed);
  gStyle->SetFuncWidth(2);
  gStyle->SetFuncColor(kGreen);
  gStyle->SetLineWidth(2.);
  gStyle->SetLabelSize(0.045,"xyz");
  gStyle->SetLabelOffset(0.01,"y");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetTitleOffset(1.,"y");
  gStyle->SetTitleOffset(0.9,"x");
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTextSizePixels(26);
  gStyle->SetTextFont(42);
  gStyle->SetTextAlign(11);
  //  gStyle->SetTickLength(0.04,"X");  gStyle->SetTickLength(0.04,"Y");

  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  //  gStyle->SetFillColor(kWhite);
  gStyle->SetLegendFont(42);

}
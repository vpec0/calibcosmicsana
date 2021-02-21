////////////////////////////////////////////////////////////////////////
// Class:       CalibCosmicsAna
// Plugin Type: analyzer (art v3_02_06)
// File:        CalibCosmicsAna_module.cc
//
// Generated at Tue Jul  2 08:25:23 2019 by Viktor Pec using cetskelgen
// from cetlib version v3_07_02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

#include "art_root_io/TFileService.h"

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataalg/DetectorInfo/DetectorClocksData.h"


#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"

#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"

#include "CalibCosmicsEvent.h"

class CalibCosmicsAna;

using namespace std;

class CalibCosmicsAna : public art::EDAnalyzer {
public:
  explicit CalibCosmicsAna(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CalibCosmicsAna(CalibCosmicsAna const&) = delete;
  CalibCosmicsAna(CalibCosmicsAna&&) = delete;
  CalibCosmicsAna& operator=(CalibCosmicsAna const&) = delete;
  CalibCosmicsAna& operator=(CalibCosmicsAna&&) = delete;

  // Required functions.
    virtual void analyze(art::Event const& e) override;

    virtual void beginJob() override;

public:
    static const size_t fMaxTrajPoints = 3000;

private:
    bool insideTPC(const TVector3& pos);
    void GetTPClimits();
    double length(const simb::MCParticle& p,
		  TLorentzVector& start, TLorentzVector& end,
		  unsigned int &starti, unsigned int &endi,
		  float *dedx, float* pos);

    int analyzeMC(detinfo::DetectorClocksData const& clockData, art::Event const& e);
    int analyzeReco(detinfo::DetectorClocksData const& clockData, art::Event const& e);

    const std::vector<const recob::Hit*>
    GetRecoTrackHits(const recob::Track &track,
		     art::Event const &evt) const;

    std::vector<std::pair<const simb::MCParticle*, double> >
    GetMCParticleListFromRecoTrack(detinfo::DetectorClocksData const& clockData,
				   const recob::Track &track,
    				   art::Event const & evt) const;

private:
    // Declare member data here.
    // labels
    std::string fSimModuleLabel;
    std::string fHitModuleLabel;
    std::string fRecoTrackModuleLabel;
    std::string fCaloModuleLabel;


    int fEvent;
    int fRun;
    int fSubRun;

    /// Flat tree to store info of primary muon which enters TPC
    TTree* fTree;

    CalibCosmicsEvent fCosmicEvent;

    geo::GeometryCore const* fGeometryService;

    double ActiveBounds[6]; // Cryostat boundaries ( neg x, pos x, neg y, pos y, neg z, pos z )

    double fXmin;
    double fXmax;
    double fYmin;
    double fYmax;
    double fZmin;
    double fZmax;
};


CalibCosmicsAna::CalibCosmicsAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  // More initializers here.
    fSimModuleLabel   (p.get< std::string >("SimModuleLabel")),
    fHitModuleLabel  (p.get< std::string >("HitModuleLabel")),
    fRecoTrackModuleLabel  (p.get< std::string >("RecoTrackModuleLabel")),
    fCaloModuleLabel   (p.get< std::string >("CaloModuleLabel"))
{

    // Get a pointer to the geometry service provider.
    fGeometryService = lar::providerFrom<geo::Geometry>();

    GetTPClimits();

    // Call appropriate consumes<>() for any products to be retrieved by this module.
    // Since art 2.8, you can and should tell beforehand, here in the constructor,
    // all the data the module is going to read ("consumes") or might read
    // ("may_consume"). Diligence here will in the future help the framework
    // execute modules in parallel, making sure the order is correct.
    consumes<std::vector<simb::MCParticle>>(fSimModuleLabel);
    // consumes<std::vector<sim::SimChannel>>(fSimulationProducerLabel);
    consumes<art::Assns<simb::MCTruth, simb::MCParticle>>(fSimModuleLabel);
    consumes<std::vector<recob::Hit>>(fHitModuleLabel);
    consumes<std::vector<recob::Track>>(fRecoTrackModuleLabel);
    consumes<std::vector<anab::Calorimetry>>(fCaloModuleLabel);
    consumes<art::Assns<recob::Track, anab::Calorimetry>>(fCaloModuleLabel);
    consumes<art::Assns<recob::Track, recob::Hit>>(fRecoTrackModuleLabel);
    // consumes<std::vector<recob::Hit>>(fHitProducerLabel);
    // consumes<std::vector<recob::Cluster>>(fClusterProducerLabel);
    // consumes<art::Assns<recob::Cluster, recob::Hit>>(fHitProducerLabel);
}

void CalibCosmicsAna::beginJob()
{
    art::ServiceHandle<art::TFileService const> tfs;

    fTree = tfs->make<TTree>("muon_tree", "Truth data for all cosmic muons entering the TPC");

    createBranches(fTree, &fCosmicEvent);
}


void CalibCosmicsAna::analyze(art::Event const& e)
{
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(e);
    auto const detProp =  art::ServiceHandle<detinfo::DetectorPropertiesService>()->DataFor(e, clockData);

  // Implementation of required member function here.
    fEvent = e.id().event();
    fRun = e.run();
    fSubRun = e.subRun();

    fCosmicEvent.event = fEvent;
    fCosmicEvent.subRun = fSubRun;
    fCosmicEvent.run = fRun;

    if (! analyzeMC(clockData, e) ) // primary muon did not make it into TPC
	return;
    if (! analyzeReco(clockData, e) ) // no reconstruction?
	return;

    // Fill MC and Reco info of this event
    fTree->Fill();
}

int CalibCosmicsAna::analyzeMC(detinfo::DetectorClocksData const& clockData, art::Event const& e)
{
    art::Handle< std::vector<simb::MCParticle> > particleHandle;
    // Then tell the event to fill the vector with all the objects of
    // that type produced by a particular producer.
    //
    // Note that if I don't test whether this is successful, and there
    // aren't any simb::MCParticle objects, then the first time we
    // access particleHandle, art will display a "ProductNotFound"
    // exception message and, depending on art settings, it may skip
    // all processing for the rest of this event (including any
    // subsequent analysis steps) or stop the execution.
    if (!e.getByLabel(fSimModuleLabel, particleHandle))
	{
	    // If we have no MCParticles at all in an event, then we're in
	    // big trouble. Throw an exception to force this module to
	    // stop. Try to include enough information for the user to
	    // figure out what's going on. Note that when we throw a
	    // cet::exception, the run and event number will automatically
	    // be displayed.
	    //
	    // __LINE__ and __FILE__ are values computed by the compiler.
	    throw cet::exception("AnalysisExample")
		<< " No simb::MCParticle objects in this event - "
		<< " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
	}

    int total_muons = 0;
    int tpc_crossers = 0;
    for ( auto const& particle: (*particleHandle) ) {
	// check if this is the primary muon
	if ( particle.Process() != "primary" ) continue;
	if ( TMath::Abs(particle.PdgCode()) != 13 ) continue;
	//if ( insideTPC(particle.Position().Vect()) ) continue;

	total_muons++;

	auto traj = particle.Trajectory();

	TLorentzVector mcstart, mcend;
	unsigned int pstarti, pendi; //mcparticle indices for starts and ends in tpc or drifted volumes
	// check for trajectory's length within TPC, and get start and end position, and index of 1st and last point in TPC
	double plen = length(particle, mcstart, mcend, pstarti, pendi,
			     fCosmicEvent.mc_dedx_tpcAV, (float*)fCosmicEvent.mc_pos_tpcAV);
	fCosmicEvent.mc_ntrajpoints_tpcAV = pendi - pstarti + 1;

	// find whether it has entered the TPC and where
	auto it = traj.begin(); // iterator is a pointer to a pair of 2 TLorentzVectors, pos, mom
	// store generated position and momentum
	it->first.GetXYZT(fCosmicEvent.gen_pos);
	it->second.GetXYZT(fCosmicEvent.gen_mom);

	if ( plen == 0. ) // muon did not enter the TPC
	    return 0;

	tpc_crossers++;

	// save position of muon's entry and exit of TPC
	mcstart.GetXYZT(fCosmicEvent.mc_startPoint_tpcAV);
	mcstart.GetXYZT(fCosmicEvent.mc_startMom_tpcAV);
	mcend.GetXYZT(fCosmicEvent.mc_endPoint_tpcAV);
	mcend.GetXYZT(fCosmicEvent.mc_endMom_tpcAV);

	// test last position in TPC
	auto endit = traj.begin() + pendi + 1;
	if (endit != traj.end()) {
	    fCosmicEvent.mc_endProcess_tpcAV = "exit";
	    fCosmicEvent.mc_exited_tpcAV = 1;
	} else {
	    fCosmicEvent.mc_endProcess_tpcAV = particle.EndProcess();
	    fCosmicEvent.mc_exited_tpcAV = 0;
	}

	// fill trajectory length
	fCosmicEvent.mc_length_tpcAV = plen;

	// found primary muon and it entered the TPC
	return 1;
    }
    return 0;
}


int CalibCosmicsAna::analyzeReco(detinfo::DetectorClocksData const& clockData, art::Event const& evt)
{
    /// Retrieves reconstructed tracks
    art::Handle< std::vector<recob::Track> >  trackListHandle;
    std::vector<art::Ptr<recob::Track> > tracklist;

    if ( evt.getByLabel(fRecoTrackModuleLabel, trackListHandle) )
	art::fill_ptr_vector(tracklist, trackListHandle);

    size_t ntracks = tracklist.size();

    // test number of track in order not to overflow allocated memory
    if (ntracks > MAX_RECO_TRACKS) {
	cout<<"CalibCosmicsAna: WARNING: Event "
	    <<" has "<<ntracks<<" reconstructed tracks. Will only store "<<MAX_RECO_TRACKS<<endl;
	ntracks = MAX_RECO_TRACKS;
    }

    fCosmicEvent.reco_nTracks = ntracks;

    // Fill in basic info about each track
    for (size_t i = 0; i < ntracks; i++) {
	auto track = tracklist.at(i);
	auto start = track->Start();
	auto end = track->End();

	start.GetCoordinates(fCosmicEvent.reco_startPoint[i]);
	end.GetCoordinates(fCosmicEvent.reco_endPoint[i]);

	fCosmicEvent.reco_length[i] = track->Length();


	std::vector<std::pair<const simb::MCParticle*, double> > mc_parts =
	    GetMCParticleListFromRecoTrack(clockData, *track, evt);
	if (!mc_parts.size())
	    fCosmicEvent.reco_mcpdg[i] = -999;
	else
	    fCosmicEvent.reco_mcpdg[i] = mc_parts[0].first->PdgCode();
    }



    /// Get hits list of all hits
    art::Handle< std::vector<recob::Hit> >  hitListHandle;
    std::vector<art::Ptr<recob::Hit> > hitList;
    if ( evt.getByLabel(fHitModuleLabel, hitListHandle) )
    	art::fill_ptr_vector(hitList, hitListHandle);

    // Get calorimetry for each track
    art::FindManyP<anab::Calorimetry> fmcal(trackListHandle, evt, fCaloModuleLabel);
    if (fmcal.isValid()) {
	// loop over calorimetry by track
	for (size_t i = 0; i < fmcal.size(); i++) {
	    auto calos = fmcal.at(i);
	    // loop over calorimetry by plane
	    for (auto calo : calos) {
	    	size_t end = calo->TpIndices().size();
		int plane = calo->PlaneID().Plane;

		// test size of calorimetry points, not to overflow allocated memory
		if (end > MAX_TRACK_POINTS) {
		    cout<<"CalibCosmicsAna: WARNING: Track "<<i<<", in plane"<<plane
			<<" has "<<end<<" calorimetry points. Will only store "<<MAX_TRACK_POINTS<<endl;
		    end = MAX_TRACK_POINTS;
		}

		fCosmicEvent.reco_nPoints[i][plane] = end;
		// loop over calo points
	    	for (size_t ipt = 0; ipt < end; ++ipt) {
		    // hit related info
	    	    size_t tpidx = calo->TpIndices().at(ipt);
	    	    auto hit = hitList[tpidx]; // assuming trajectory point index is equal to index in hit list...
		    int tpc = hit->WireID().TPC;
		    float thit = hit->PeakTime();

		    // calorimetry info
		    float dqdx = calo->dQdx().at(ipt);
		    float dedx = calo->dEdx().at(ipt);
		    float resrange = calo->ResidualRange().at(ipt);
		    float pitch = calo->TrkPitchVec().at(ipt);

		    // save calo info into the event struct
		    fCosmicEvent.reco_dqdx[i][plane][ipt] = dqdx;
		    fCosmicEvent.reco_dedx[i][plane][ipt] = dedx;
		    fCosmicEvent.reco_resrange[i][plane][ipt] = resrange;
		    fCosmicEvent.reco_pitch[i][plane][ipt] = pitch;
		    // set position of the calo point and time of the hit
		    calo->XYZ().at(ipt).GetCoordinates(fCosmicEvent.reco_point[i][plane][ipt]);
		    fCosmicEvent.reco_point[i][plane][ipt][3] = thit;
		    fCosmicEvent.reco_pointtpc[i][plane][ipt] = tpc;
	    	}
	    }
	}
    }

    fCosmicEvent.reco_nTracks = ntracks;


    return 1;
}

bool CalibCosmicsAna::insideTPC(const TVector3& pos)
{
    // check if this point is inside active volume of the TPC
    if (pos.X() < fXmin || pos.X() > fXmax)
	return false;
    if (pos.Y() < fYmin || pos.Y() > fYmax)
	return false;
    if (pos.Z() < fZmin || pos.Z() > fZmax)
	return false;

    return true;
}

void CalibCosmicsAna::GetTPClimits()
{
    // Don't remember where I got this definition of the active region.
    std::cout << "----> HERE!!!!: N TPCs: " << fGeometryService->NTPC() << std::endl;
    std::cout << fGeometryService->DetectorName() << std::endl;

    if (fGeometryService->DetectorName() == "protodune") {
	auto box1 = fGeometryService->TPC(1).ActiveBoundingBox(); // beam-side upstream tpc
	auto box2 = fGeometryService->TPC(2).ActiveBoundingBox(); // rack-side upstream tpc
	auto box3 = fGeometryService->TPC(9).ActiveBoundingBox(); // beam-side downstream tpc

	fXmin = box1.MinX(); fXmax = box2.MaxX();
	fYmin = box1.MinY(); fYmax = box2.MaxY();
	fZmin = box1.MinZ(); fZmax = box3.MaxZ();
    } else { // assume 10kt FD
	auto box1 = fGeometryService->TPC(1).ActiveBoundingBox(); // beam-right bottom upstream tpc
	auto box2 = fGeometryService->TPC(10).ActiveBoundingBox(); // mean-left top upstream tpc
	auto box3 = fGeometryService->TPC(289).ActiveBoundingBox(); // beam-right bottom downstream tpc

	fXmin = box1.MinX(); fXmax = box2.MaxX();
	fYmin = box1.MinY(); fYmax = box2.MaxY();
	fZmin = box1.MinZ(); fZmax = box3.MaxZ();
	// fXmin = -300; fXmax = 300;
	// fYmin = 0; fYmax = 600;
	// fZmin = 0; fZmax = 690;
    }

    std::cout << "min X: " << fXmin << "max X: " << fXmax << std::endl
	      << "min Y: " << fYmin << "max Y: " << fYmax << std::endl
	      << "min Z: " << fZmin << "max Z: " << fZmax << std::endl;


    // Taken from dune::AnalysisTree, this definition of active region
    // extends about 8 cm further along drift coordinate on both sides
    //
    // Build my Cryostat boundaries array...Taken from Tyler Alion in
    // Geometry Core. Should still return the same values for uBoone.
    ActiveBounds[0] = ActiveBounds[2] = ActiveBounds[4] = DBL_MAX;
    ActiveBounds[1] = ActiveBounds[3] = ActiveBounds[5] = -DBL_MAX;
    // assume single cryostats
    auto const* geom = lar::providerFrom<geo::Geometry>();
    for (geo::TPCGeo const& TPC: geom->IterateTPCs()) {
	// get center in world coordinates
	double origin[3] = {0.};
	double center[3] = {0.};
	TPC.LocalToWorld(origin, center);
	double tpcDim[3] = {TPC.HalfWidth(), TPC.HalfHeight(), 0.5*TPC.Length() };

	if( center[0] - tpcDim[0] < ActiveBounds[0] ) ActiveBounds[0] = center[0] - tpcDim[0];
	if( center[0] + tpcDim[0] > ActiveBounds[1] ) ActiveBounds[1] = center[0] + tpcDim[0];
	if( center[1] - tpcDim[1] < ActiveBounds[2] ) ActiveBounds[2] = center[1] - tpcDim[1];
	if( center[1] + tpcDim[1] > ActiveBounds[3] ) ActiveBounds[3] = center[1] + tpcDim[1];
	if( center[2] - tpcDim[2] < ActiveBounds[4] ) ActiveBounds[4] = center[2] - tpcDim[2];
	if( center[2] + tpcDim[2] > ActiveBounds[5] ) ActiveBounds[5] = center[2] + tpcDim[2];
    } // for all TPC
    std::cout << "Active Boundaries: "
	      << "\n\tx: " << ActiveBounds[0] << " to " << ActiveBounds[1]
	      << "\n\ty: " << ActiveBounds[2] << " to " << ActiveBounds[3]
	      << "\n\tz: " << ActiveBounds[4] << " to " << ActiveBounds[5]
	      << std::endl;
}

// Length of MC particle, trajectory by trajectory (with out the manual shifting for x correction)
// and fill in dE/dx and positions at each point in TPC
// Taken from dune::AnalysisTree
double CalibCosmicsAna::length(const simb::MCParticle& p,
			       TLorentzVector& start, TLorentzVector& end,
			       unsigned int &starti, unsigned int &endi,
			       float *dedx_arr, float* pos_arr)
{
    double result = 0.;
    TVector3 disp;
    bool first = true;

    TLorentzVector lastPos; // energy in GeV
    double lastEnergy = 0.; // energy in GeV
    double dedx = 0.;

    int tpindex = 0;

    unsigned int size = p.NumberTrajectoryPoints();

    for(unsigned int i = 0; i < size; ++i) {
	// test number of points to be stored.
	if (tpindex >= MAX_TRAJECTORY_POINTS) {
	    cout<<"CalibCosmicsAna: WARNING: MCParticle "<<p.TrackId()<<"/"<<p.PdgCode()
		<<" has "<<size<<" trajectory points. And more than "<<MAX_TRAJECTORY_POINTS
		<<" are in TPC active. Will only store "<<MAX_TRAJECTORY_POINTS<<endl;
	    break;
	}

	// check if the particle is inside a TPC
	if (p.Vx(i) >= ActiveBounds[0] && p.Vx(i) <= ActiveBounds[1] && p.Vy(i) >= ActiveBounds[2] && p.Vy(i) <= ActiveBounds[3] && p.Vz(i) >= ActiveBounds[4] && p.Vz(i) <= ActiveBounds[5]){
	    if(first){
		start = p.Position(i);
		first = false;
		starti = i;
		lastPos = start;
		lastEnergy = p.E(i);
	    }else{
		disp -= p.Position(i).Vect();
		result += disp.Mag();
	    }
	    // calculate truth based dE/dx
	    double de = lastEnergy - p.E(i); // GeV
	    TLorentzVector pos (p.Position(i));
	    double dx = (lastPos.Vect() - pos.Vect()).Mag(); // cm
	    if (dx == 0.)
		dedx = 0.;
	    else
		dedx = de*1e3/dx; // MeV/cm

	    lastPos = pos;
	    lastEnergy = p.E(i);

	    dedx_arr[tpindex] = dedx;
	    pos.GetXYZT(pos_arr+4*tpindex);
	    tpindex++;

	    disp = p.Position(i).Vect();
	    end = p.Position(i);
	    endi = i;
	}
    }
    return result;
}

// code taken from protoana::ProtoDUNETrackUtils
const std::vector<const recob::Hit*>
CalibCosmicsAna::GetRecoTrackHits(const recob::Track &track,
				  art::Event const &evt) const
{

    auto recoTracks = evt.getValidHandle<std::vector<recob::Track> >(fRecoTrackModuleLabel);
    art::FindManyP<recob::Hit> findHits(recoTracks,evt,fRecoTrackModuleLabel);
    std::vector<art::Ptr<recob::Hit>> inputHits = findHits.at(track.ID());

    std::vector<const recob::Hit*> trackHits;

    for(const art::Ptr<recob::Hit> hit : inputHits){

	trackHits.push_back(hit.get());

    }

    return trackHits;

}

// code taken from protoana::ProtoDUNETruthUtils: GetMCParticleListFromRecoTrack and GetMCParticleListFromTrackHits
std::vector<std::pair<const simb::MCParticle*, double> >
CalibCosmicsAna::GetMCParticleListFromRecoTrack(detinfo::DetectorClocksData const& clockData,
						const recob::Track &track,
						art::Event const & evt) const
{
    // We must have MC for this module to make sense
    if(evt.isRealData()) return {};

    // Get the reconstructed track hits
    std::vector<const recob::Hit*> trackHits = GetRecoTrackHits(track, evt);


   using weightedMCPair = std::pair<const simb::MCParticle*, double>;
   std::vector<weightedMCPair> outVec;

   // Loop over all hits in the input vector and record the contributing MCParticles.
   art::ServiceHandle<cheat::BackTrackerService> bt_serv;
   art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
   std::unordered_map<const simb::MCParticle*, double> mcEMap;
   double hitTotalE = 0;
   for(const recob::Hit* hitp : trackHits) {
       for(const sim::TrackIDE& ide : bt_serv->HitToTrackIDEs(clockData, *hitp)) {
	   const simb::MCParticle* curr_part = pi_serv->TrackIdToParticle_P(ide.trackID);
	   mcEMap[curr_part] += ide.energy;
	   hitTotalE += ide.energy;
       }
   }

   // Fill and sort the output vector
   for (weightedMCPair const& p : mcEMap) {
       outVec.push_back(p);
   }
   std::sort(outVec.begin(), outVec.end(),
	     [](weightedMCPair a, weightedMCPair b){ return a.second > b.second;});

   // Normalise the weights by the total track energy.
   if (hitTotalE < 1e-5) { hitTotalE = 1; } // Protect against zero division
   for (weightedMCPair& p : outVec) {
       p.second /= hitTotalE;
   }

   return outVec;

}

DEFINE_ART_MODULE(CalibCosmicsAna)

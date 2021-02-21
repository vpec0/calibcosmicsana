#ifndef CALIBCOSMICSEVENT_H
#define CALIBCOSMICSEVENT_H


#define EVENT_LIST				\
    INT(run, "Run number");			\
    INT(subRun, "Sub-run number");		\
    INT(event, "Event number")			\

#define MC_TRUTH_LIST							\
    /* Primary muon */							\
    /* Generator level */						\
    FLOAT_4VEC(gen_mom, "Four momentum");				\
    FLOAT_4VEC(gen_pos, "Position");					\
									\
    /* Active volume level */						\
    FLOAT_4VEC(mc_startPoint_tpcAV, "Start point in TPC Active Volume"); \
    FLOAT_4VEC(mc_endPoint_tpcAV, "End point in TPC Active Volume");	\
    STRING(mc_endProcess_tpcAV, "Process at end point in TPC Active Volume"); \
    FLOAT_4VEC(mc_startMom_tpcAV, "Start-point 4-momentum in TPC Active Volume"); \
    FLOAT_4VEC(mc_endMom_tpcAV, "End-point 4-momentum in TPC Active Volume"); \
    FLOAT(mc_length_tpcAV, "Length of muon trajectory in TPC Active Volume [cm]"); \
    INT(mc_exited_tpcAV, "Whether muon exited TPC.");			\
    INT(mc_ntrajpoints_tpcAV, "Number of trajectory points in TPC.");	\
    FLOAT_ARR(mc_dedx_tpcAV, "dE/dx for segment between trajectory points in TPC", mc_ntrajpoints_tpcAV); \
    FLOAT_ARR2D(mc_pos_tpcAV, "Position of trajectory point in TPC", mc_ntrajpoints_tpcAV, 4) \
									\

#define RECO_TRACK_LIST							\
    /* Info on reconstructed tracks */					\
    INT(reco_nTracks, "Number of Reco Tracks");				\
    INT_ARR(reco_mcpdg, "MC Truth PDG code", reco_nTracks);		\
    INT_ARR(reco_g4id, "MC G4 ID", reco_nTracks);			\
    FLOAT_ARR2D(reco_startPoint, "Start point of reco track", reco_nTracks, 3); \
    FLOAT_ARR2D(reco_endPoint, "End point of reco track", reco_nTracks, 3); \
    FLOAT_ARR(reco_length, "Length of reco track [cm]", reco_nTracks);	\
    /* Info on individual track points */				\
    INT_ARR2D(reco_nPoints, "Number of track points in each plane", reco_nTracks, 3); \
    FLOAT_ARR3D(reco_dqdx, "Track point dQ/dx", reco_nTracks, 3, MAX_TRACK_POINTS); \
    FLOAT_ARR3D(reco_dedx, "Track point dE/dx", reco_nTracks, 3, MAX_TRACK_POINTS); \
    FLOAT_ARR3D(reco_resrange, "Track point residual range", reco_nTracks, 3, MAX_TRACK_POINTS); \
    FLOAT_ARR3D(reco_pitch, "Track pitch at respective point", reco_nTracks, 3, MAX_TRACK_POINTS); \
    FLOAT_ARR4D(reco_point, "Track point location", reco_nTracks, 3, MAX_TRACK_POINTS, 4); \
    SHORT_ARR3D(reco_pointtpc, "Track point TPC ID", reco_nTracks, 3, MAX_TRACK_POINTS) \



const int MAX_RECO_TRACKS = 100;
const int MAX_TRACK_POINTS = 5000; // 25 TPC along Z times 1000 wires each
const int MAX_TRAJECTORY_POINTS = 2000; // looks like trajectory point every ~3 cm


struct CalibCosmicsEvent {
#define INT(var, title) int var
#define FLOAT(var, title) float var

    EVENT_LIST;

#define STRING(var, title) std::string var
#define FLOAT_ARR(var, title, size) float var[MAX_TRAJECTORY_POINTS]
#define FLOAT_4VEC(var, title) float var[4]
#define FLOAT_ARR2D(var, title, size1, size2) float var[MAX_TRAJECTORY_POINTS][size2]

    MC_TRUTH_LIST;

#define INT_ARR(var, title, size) int var[MAX_RECO_TRACKS]
#define INT_ARR2D(var, title, size1, size2) int var[MAX_RECO_TRACKS][size2]
#undef FLOAT_ARR
#define FLOAT_ARR(var, title, size) float var[MAX_RECO_TRACKS]
#undef FLOAT_ARR2D
#define FLOAT_ARR2D(var, title, size1, size2) float var[MAX_RECO_TRACKS][size2]
#define FLOAT_ARR3D(var, title, size1, size2, size3) float var[MAX_RECO_TRACKS][size2][size3]
#define FLOAT_ARR4D(var, title, size1, size2, size3, size4) float var[MAX_RECO_TRACKS][size2][size3][size4]
#define SHORT_ARR3D(var, title, size1, size2, size3) short var[MAX_RECO_TRACKS][size2][size3]

    RECO_TRACK_LIST;

/// delete the macro definitions
#undef INT
#undef INT_ARR
#undef INT_ARR2D
#undef FLOAT
#undef FLOAT_4VEC
#undef FLOAT_ARR
#undef FLOAT_ARR2D
#undef FLOAT_ARR3D
#undef FLOAT_ARR4D
#undef STRING
#undef SHORT_ARR3D
}; // struct CalibCosmicsEvent


/// Create branches
TTree* createBranches(TTree* tree, CalibCosmicsEvent* evt)
{
#define INT(var, title) tree->Branch(#var, &evt->var, #var"/I")
#define FLOAT(var, title) tree->Branch(#var, &evt->var, #var"/F")
#define STRING(var, title) tree->Branch(#var, &evt->var)
#define FLOAT_4VEC(var, title) tree->Branch(#var, evt->var, #var"[4]/F")
#define FLOAT_ARR(var, title, size) tree->Branch(#var, evt->var, #var"["#size"]/F")
#define FLOAT_ARR2D(var, title, size1, size2) tree->Branch(#var, evt->var, Form(#var"["#size1"][%d]/F", size2))
#define INT_ARR(var, title, size) tree->Branch(#var, evt->var, #var"["#size"]/I")
#define INT_ARR2D(var, title, size1, size2) tree->Branch(#var, evt->var, Form(#var"["#size1"][%d]/I", size2))
#define FLOAT_ARR3D(var, title, size1, size2, size3) tree->Branch(#var, evt->var, Form(#var"["#size1"][%d][%d]/F", size2, size3))
#define FLOAT_ARR4D(var, title, size1, size2, size3, size4) tree->Branch(#var, evt->var, Form(#var"["#size1"][%d][%d][%d]/F", size2, size3, size4))
#define SHORT_ARR3D(var, title, size1, size2, size3) tree->Branch(#var, evt->var, Form(#var"["#size1"][%d][%d]/S", size2, size3))

    EVENT_LIST;
    MC_TRUTH_LIST;
    RECO_TRACK_LIST;

    return tree;
}


/// delete the macro definitions
#undef INT
#undef INT_ARR
#undef INT_ARR2D
#undef FLOAT
#undef FLOAT_4VEC
#undef FLOAT_ARR
#undef FLOAT_ARR2D
#undef FLOAT_ARR3D
#undef FLOAT_ARR4D
#undef STRING
#undef SHORT_ARR3D


/// set addresses for tree branches
TTree* setAddresses(TTree* tree, CalibCosmicsEvent* evt)
{
#define INT(var, title) tree->SetBranchAddress(#var, &evt->var)
#define FLOAT(var, title) tree->SetBranchAddress(#var, &evt->var)
#define STRING(var, title) tree->SetBranchAddress(#var, &evt->var)
#define FLOAT_4VEC(var, title) tree->SetBranchAddress(#var, evt->var)
#define FLOAT_ARR(var, title, size) tree->SetBranchAddress(#var, evt->var)
#define FLOAT_ARR2D(var, title, size1, size2) tree->SetBranchAddress(#var, evt->var)
#define INT_ARR(var, title, size) tree->SetBranchAddress(#var, evt->var)
#define INT_ARR2D(var, title, size1, size2) tree->SetBranchAddress(#var, evt->var)
#define FLOAT_ARR3D(var, title, size1, size2, size3) tree->SetBranchAddress(#var, evt->var)
#define FLOAT_ARR4D(var, title, size1, size2, size3, size4) tree->SetBranchAddress(#var, evt->var)
#define SHORT_ARR3D(var, title, size1, size2, size3) tree->SetBranchAddress(#var, evt->var)

    EVENT_LIST;
    MC_TRUTH_LIST;
    RECO_TRACK_LIST;

    return tree;
}


/// delete the macro definitions
#undef INT
#undef INT_ARR
#undef INT_ARR2D
#undef FLOAT
#undef FLOAT_4VEC
#undef FLOAT_ARR
#undef FLOAT_ARR2D
#undef FLOAT_ARR3D
#undef FLOAT_ARR4D
#undef STRING
#undef SHORT_ARR3D



#endif

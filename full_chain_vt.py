#!/usr/bin/env python3
from pathlib import Path
import os

import shutil
import os
import argparse
import pathlib
import acts
import acts.examples
import time


from acts import UnitConstants as u

from acts.examples.simulation import (
    addParticleReader,
    addFatras,
    addGeant4,
    ParticleSelectorConfig,
    addDigitization,
)
from acts.examples.reconstruction import (
    addNewSeeding,
    TruthSeedRanges,
    SeedFinderConfigArgNA60,
    SeedFinderOptionsArgNA60,
    SeedFilterConfigArgNA60,

    SpacePointGridConfigArg,
    SeedingAlgorithmConfigArg,
    addCKFTracks,
    addUsedMeasurementsFilter,
    TrackSelectorConfig,
    addAmbiguityResolution,
    addContainerMerger,
    AmbiguityResolutionConfig,
    addVertexFitting,
    addSecondaryVertexFitting,
    VertexFinder,
    CkfConfig,
    addMatching
)
from acts.examples import TGeoDetector




def getArgumentParser():
    """Get arguments from command line"""
    parser = argparse.ArgumentParser(
        description="Command line arguments for CKF")

    parser.add_argument(
        "-o",
        "--output",
        dest="outdir",
        help="Output directory for new ntuples",
        default=None,
    )
    parser.add_argument(
        "-n", "--nEvents", dest="nEvts", help="Number of events to run over", default=1, type=int
    )

    parser.add_argument(
        '-g',
        '--geant4',
        help='Use Geant4',
        action='store_true'
    )
    
    parser.add_argument(
        '-ve',
        '--verbose',
        help='Use Geant4',
        action='store_true'
    )
    
    parser.add_argument(
        '-d',
        '--dead',
        help='Apply dead zones',
        action='store_true'
    )
    parser.add_argument(
        '-f',
        '--fast',
        help='Apply fast sim selections',
        action='store_true'
    )

    parser.add_argument(
        '-st',
        '--primary',
        help='No primary particles',
        action='store_true'
    )

    parser.add_argument(
        '-nd',
        '--secondary',
        help='No secondary particles',
        action='store_true'
    )

    parser.add_argument(
        '-gu',
        '--noguessing',
        help='No secondary particles',
        action='store_true'
    )


    parser.add_argument(
        '-sec',
        '--sec',
        help='No secondary particles',
        action='store_true'
    )

    parser.add_argument(
        '-osec',
        '--osec',
        help='No secondary particles',
        action='store_true'
    )


    parser.add_argument(
        '-zz',
        '--zz',
        help='No secondary particles',
        action='store_true'
    )


    parser.add_argument(
        '-t0',
        '--t0',
        help='No secondary particles',
        action='store_true'
    )


    parser.add_argument(
        '--se',
        '--seed',
        dest="seed",
        help='Use particle gun',
        type=int,
        default=4,
    )

    ##########################################
    ##
    ## eint-> collision energy
    ##
    ## pdg code
    ##
    ##########################################

    parser.add_argument(
        '-e',
        '--eint',
        dest="eint",
        help='Collision energy',
        type=int,
        default=40,
    )    
    
    parser.add_argument(
        '-fullev',
        '--fullev',
        dest="fullev",
        help='Events with also K0S and Lambda0s',
        action='store_true'
    )

    parser.add_argument(
        '-bkg',
        '--bkg',
        dest="bkg",
        help='Remove bkg',
        action='store_true'
    )

    parser.add_argument(
        '-trk',
        '--trkVtxOnly',
        dest="trkVtxOnly",
        help='Remove bkg',
        action='store_true'
    )

    parser.add_argument(
        '-del',
        '--addDeltas',
        dest="addDeltas",
        help='Remove bkg',
        action='store_true'
    )
    
    parser.add_argument(
        '-per',
        '--periferal',
        dest="periferal",
        type=float,
        default=1,
    )

    parser.add_argument(
        '-pro',
        '--projective',
        dest="projective",
        help='Remove bkg',
        action='store_true'
    )
    ##########################################
    ##
    ##
    ##
    #################
    """
    python3 full_chain_vt.py -n1000 -d -per 0.01 > test_per0.01.out &&

    python3 full_chain_vt.py -n1000 -d -sec -per 0.01 > test_sec_per0.01.out &&

    python3 full_chain_vt.py -n1000 -d -per 0.02 > test_per0.02_rej0.114.out &&
    python3 full_chain_vt.py -n1000 -d -per 0.05 > test_per0.05_rej0.114.out &&
    python3 full_chain_vt.py -n1000 -d -per 0.1 > test_per0.1_rej0.114.out &&
    python3 full_chain_vt.py -n1000 -d -per 0.2 > test_per0.2_rej0.114.out &&
    python3 full_chain_vt.py -n1000 -d -per 0.4 > test_per0.4_rej0.114.out &&
    python3 full_chain_vt.py -n1000 -d -per 0.6 > test_per0.6_rej0.114.out &&
    python3 full_chain_vt.py -n1000 -d -per 0.8 > test_per0.8_rej0.114.out &&
    python3 full_chain_vt.py -n1000 -d         > test_central_rej0.114.out

    python3 full_chain_vt.py -n1000 -d -sec -per 0.02 > test_sec_per0.02_rej0.114.out &&
    python3 full_chain_vt.py -n1000 -d -sec -per 0.05 > test_sec_per0.05_rej0.114.out &&
    python3 full_chain_vt.py -n1000 -d -sec -per 0.1 > test_sec_per0.1_rej0.114.out &&
    python3 full_chain_vt.py -n1000 -d -sec -per 0.2 > test_sec_per0.2_rej0.114.out &&
    python3 full_chain_vt.py -n1000 -d -sec -per 0.4 > test_sec_per0.4_rej0.114.out &&
    python3 full_chain_vt.py -n1000 -d -sec -per 0.6 > test_sec_per0.6_rej0.114.out &&
    python3 full_chain_vt.py -n1000 -d -sec -per 0.8 > test_sec_per0.8_rej0.114.out &&
    python3 full_chain_vt.py -n1000 -d -sec         > test_sec_central_rej0.114.out

    python3 full_chain_vt.py -n1000 -d --sf_rejectedFraction 0 -per 0.02 > test_per0.02_rej0.out &&
    python3 full_chain_vt.py -n1000 -d --sf_rejectedFraction 0 -per 0.05 > test_per0.05_rej0.out &&
    python3 full_chain_vt.py -n1000 -d --sf_rejectedFraction 0 -per 0.1 > test_per0.1_rej0.out &&
    python3 full_chain_vt.py -n1000 -d --sf_rejectedFraction 0 -per 0.2 > test_per0.2_rej0.out &&
    python3 full_chain_vt.py -n1000 -d --sf_rejectedFraction 0 -per 0.4 > test_per0.4_rej0.out &&
    python3 full_chain_vt.py -n1000 -d --sf_rejectedFraction 0 -per 0.6 > test_per0.6_rej0.out &&
    python3 full_chain_vt.py -n1000 -d --sf_rejectedFraction 0 -per 0.8 > test_per0.8_rej0.out &&
    python3 full_chain_vt.py -n1000 -d --sf_rejectedFraction 0         > test_central_rej0.out

    python3 full_chain_vt.py -n1000 -d --sf_rejectedFraction 0 -sec -per 0.02 > test_sec_per0.02_rej0.out &&
    python3 full_chain_vt.py -n1000 -d --sf_rejectedFraction 0 -sec -per 0.05 > test_sec_per0.05_rej0.out &&
    python3 full_chain_vt.py -n1000 -d --sf_rejectedFraction 0 -sec -per 0.1 > test_sec_per0.1_rej0.out &&
    python3 full_chain_vt.py -n1000 -d --sf_rejectedFraction 0 -sec -per 0.2 > test_sec_per0.2_rej0.out &&
    python3 full_chain_vt.py -n1000 -d --sf_rejectedFraction 0 -sec -per 0.4 > test_sec_per0.4_rej0.out &&
    python3 full_chain_vt.py -n1000 -d --sf_rejectedFraction 0 -sec -per 0.6 > test_sec_per0.6_rej0.out &&
    python3 full_chain_vt.py -n1000 -d --sf_rejectedFraction 0 -sec -per 0.8 > test_sec_per0.8_rej0.out &&
    python3 full_chain_vt.py -n1000 -d --sf_rejectedFraction 0 -sec         > test_sec_central_rej0.out

    python3 full_chain_vt.py -n1000 -d -sec -per 0.02 > test.out &&


    """
    #########################
    parser.add_argument(
        '-ef',
        '--eff',
        dest="eff",
        help='Collision energy',
        type=float,
        default=None,
    )

    # arguments for the seeding
    parser.add_argument(
        "--sf_maxSeedsPerSpM_primary",
        dest="sf_maxSeedsPerSpM_primary",
        help="Number of compatible seeds considered for middle seed in the secondary ste",
        type=int,
        default=1,
    )
    
    # arguments for the seeding
    parser.add_argument(
        "--sf_maxSeedsPerSpM_secondary",
        dest="sf_maxSeedsPerSpM_secondary",
        help="Number of compatible seeds considered for middle seed in the secondary step",
        type=int,
        default=20
    )

    parser.add_argument(
        "--sf_impactMax",
        dest="sf_impactMax",
        help="max impact parameter in mm",
        type=float,
        default=1,
    )

    parser.add_argument(
        "--sf_seedConfirmation_primary",
        dest="sf_seedConfirmation_primary",
        help="Use seed confirmation",
        type=bool,
        default=False
    )

    parser.add_argument(
        "--sf_seedConfirmation_secondary",
        dest="sf_seedConfirmation_secondary",
        help="Use seed confirmation",
        type=bool,
        default=False
    )

    parser.add_argument(
        "--sf_deltaZMax",
        dest="sf_deltaZMax",
        help="Use seed confirmation",
        type=float,
        default=50,
    )

    parser.add_argument(
        "--sf_numMeasurementsCutOff",
        dest="sf_numMeasurementsCutOff",
        help="maximum value for deltaR separation in mm",
        type=int,
        default=1,
    )

    
    parser.add_argument(
        "-z",
        dest="zperigee",
        help="Use seed confirmation",
        type=float,
        default=0,
    )

    
    parser.add_argument(
        "-t",
        dest="target",
        help="Use seed confirmation",
        type=int,
        default=-1,
    )

    #for vtx optimization
    #  significanceCutSeeding
    #  maximumChi2cutForSeeding
    #  maxVertices
    #  SplitVertices
    #  splitVerticesTrkInvFraction
    #  reassignTracksAfterFirstFit
    #  doMaxTracksCut
    #  maxTracks
    #  cutOffTrackWeight
    #  cutOffTrackWeightReassign

    parser.add_argument(
        "--sf_cutOffTrackWeight",
        dest="sf_cutOffTrackWeight",
        help="",
        type=float,
        default=0.01,
    )

    parser.add_argument(
        "--sf_cutOffTrackWeightReassign",
        dest="sf_cutOffTrackWeightReassign",
        help="",
        type=float,
        default=1,
    )

    parser.add_argument(
        "--sf_significanceCutSeeding",
        dest="sf_significanceCutSeeding",
        help="",
        type=float,
        default=10, #0.06758286918225764,
    )

    parser.add_argument(
        "--sf_maximumChi2cutForSeeding",
        dest="sf_maximumChi2cutForSeeding",
        help="",
        type=float,
        default=36,
    )

    parser.add_argument(
        "--sf_maxVertices",
        dest="sf_maxVertices",
        help="",
        type=int,
        default=1,
    )

    parser.add_argument(
        "--sf_createSplitVertices",
        dest="sf_createSplitVertices",
        help="",
        type=bool,
        default=False,
    )

    parser.add_argument(
        "--sf_reassignTracksAfterFirstFit",
        dest="sf_reassignTracksAfterFirstFit",
        help="",
        type=bool,
        default=False,
    )

    parser.add_argument(
        "--sf_doMaxTracksCut",
        dest="sf_doMaxTracksCut",
        help="",
        type=bool,
        default=False,
    )

    parser.add_argument(
        "--sf_splitVerticesTrkInvFraction",
        dest="sf_splitVerticesTrkInvFraction",
        help="",
        type=int,
        default=2,
    )

    parser.add_argument(
        "--sf_maxTracks",
        dest="sf_maxTracks",
        help="",
        type=int,
        default=5000,
    )

    parser.add_argument(
        "--sf_rejectedFraction",
        dest="sf_rejectedFraction",
        help="",
        type=float,
        default=0.114,
    )



    return parser


def runVTchain(
    detector,
    trackingGeometry,
    inputDir: Path,
    outputDir: Path,
    jsonDigi="geometry/geomVTNA60p/digismear.json",
    jsonSeed="geometry/geomVTNA60p/seed_config.json",
    useGeant4=False,
    truthSeeding=False,
    applyDeadZones=False,
    applyFastSimSelections=False,
    applyReadout=False,
    applyBackbone=False,
    applyHole=True,
    applyEndcapShort=False,
    applyEndcapLong=False,
    NumEvents=1,
    CotThetaMax=5.809222379141632,
    SigmaScattering=7.293572849910582,
    RadLengthPerSeed=0.00015,
    MaxPtScattering=10.0,
    MaxSeedsPerSpMPrimary=1,
    MaxSeedsPerSpMSecondary=1,
    SeedConfirmationPrimary=False,
    SeedConfirmationSecondary=False,
    ImpactMax=2,
    DeltaZMax=5,
    NumMeasurementsCutOff=5,
    verbose=False,
    Efficiency=None,
    KeepPrimary=True,
    KeepSecondary=True,
    noGuessing=False,
    trkVtxOnly=False,
    addDeltas=True,
    projective=False,
    zPerigee=0,
    significanceCutSeeding=0,
    maximumChi2cutForSeeding=0,
    maxVertices=0,
    createSplitVertices=0,
    splitVerticesTrkInvFraction=0,
    reassignTracksAfterFirstFit=0,
    doMaxTracksCut=0,
    maxTracks=0,
    cutOffTrackWeight=0,
    cutOffTrackWeightReassign=0,
    rejectedFraction=0.9,

):

    field = acts.ConstantBField(acts.Vector3(
        0.0, 1.5 * u.T, 0.0))  # ~dipole field

    rnd = acts.examples.RandomNumbers(seed=44)

    s = acts.examples.Sequencer(
        events=NumEvents,
        # events=1,
        numThreads=1,
        outputDir=str(outputDir),
        # to remove mechanism to catch FPE during algorithm execution (see Mattermost 18/7/23 Igor)
        trackFpes=False,
    )

    addParticleReader(
        s,
        inputDir=inputDir,
        outputDirRoot=outputDir,
    )

    addFatras(
        s,
        trackingGeometry,
        field,
        rnd,
        preSelectParticles=ParticleSelectorConfig(
            eta=(None, None),
            pt=(100 * u.MeV, None),
            removeNeutral=True,
        ),
        outputDirRoot=outputDir,
    )

    
    addDigitization(
        s,
        trackingGeometry,
        field,
        digiConfigFile=jsonDigi,
        outputDirRoot=outputDir,
        efficiency=1 if Efficiency is not None else 0.99 if applyDeadZones or applyFastSimSelections else 1,
        applyDeadAreas=applyDeadZones,
        applyFastSimSelections=applyFastSimSelections,
         applyReadout=applyReadout,
        applyBackbone=applyBackbone,
        applyHole=applyHole,
        applyEndcapShort=applyEndcapShort,
        applyEndcapLong=applyEndcapLong,
        rnd=rnd,
    )

    # rotation of B field, so that B = Bz, needed for seeding
    field2 = field if truthSeeding else acts.ConstantBField(acts.Vector3(0.0, 0.0, -1.5 * u.T))

    TRUTHSEEDRANGE = TruthSeedRanges(pt=(100 * u.MeV, None),
                     eta=(None, None),
                     nHits=(4, None),
                     keep=(KeepPrimary,KeepSecondary))
    
    #SEEDFINDEROPTIONSARG = 
    
    SPACEPOINTGRIDCONFIGARG = SpacePointGridConfigArg(#NO NEED TO CHANGE THIS
                                    rMax=400 * u.mm,
                                    zBinEdges=[-150.,-75.,75.,150.]

                                )
    
    CKFCONFIG = CkfConfig(
            chi2CutOff=15,
            numMeasurementsCutOff=int(NumMeasurementsCutOff),
            maxSteps=None,
            seedDeduplication=False,
            stayOnSeed=True,
            )
    TRACKSELECTORCONFIG = TrackSelectorConfig(
            pt=(100 * u.MeV, None),
            absEta=(None, None),
            nMeasurementsMin=4,
        )
    
    AMBIGUITYRESOLUTIONCONFIG = AmbiguityResolutionConfig(
            maximumSharedHits=1,#MaxSharedHits,
            nMeasurementsMin=4,
            maximumIterations=1000000,
        )
    
    SEEDINGALGORITHMCONFIGARG = SeedingAlgorithmConfigArg(#NO NEED TO CHANGE THIS
            zBinNeighborsTop=[[0,1],[-1,1],[-1,0]],
            zBinNeighborsBottom=[[0,1],[-1,1],[-1,0]]
        )
    
    SEEDFINDEROPTIONSARG = SeedFinderOptionsArgNA60(#NO NEED TO CHANGE THIS
                                                beamPos=(0 * u.mm, 0 * u.mm),
                                                bFieldInZ=1.5 * u.T,
                                            )

    print(jsonSeed)
    #########
    # STEP 1
    #########
    addNewSeeding(
        s = s,
        trackingGeometry = trackingGeometry,
        field = field2,
        geoSelectionConfigFile = jsonSeed,
        truthSeedRanges = TRUTHSEEDRANGE,
        initialSigmas = None,
        initialVarInflation = None,
        seedFinderConfigArg = SeedFinderConfigArgNA60(
            maxSeedsPerSpM=MaxSeedsPerSpMPrimary,
            cotThetaMax=CotThetaMax,
            sigmaScattering=SigmaScattering,
            radLengthPerSeed=RadLengthPerSeed,
            maxPtScattering=MaxPtScattering,
            # min and max R between Middle and Top SP
            deltaYTopSP=(10 * u.mm, 1000 * u.mm),
            # min and max R between Middle and Bottom SP
            impactMax=ImpactMax * u.mm,
            deltaZMax=5 * u.mm,  # was 5
            minPt=100 * u.MeV,
            #interactionPointCut=True,
            interactionPointCut=True,
            verbose=verbose,
            collisionRegion=(-ImpactMax * u.mm, ImpactMax * u.mm),  # 0.5
            # NOT USED??? seems to be used in Orthogonal seeding
            y=(0 * u.mm, 400 * u.mm),
            yMiddle=(139 * u.mm, 170 * u.mm), #use only layer 2 and 3
            seedConfirmation=False, #WE ONLY HAVE 1 REGION
            seedConfirmationRange=acts.SeedConfirmationRangeConfig(  # it should not be used....
                zMinSeedConf=-150 * u.mm,#no need to change this
                zMaxSeedConf=150 * u.mm,#no need to change this
                rMaxSeedConf=140 * u.mm,#r first layer
                nTopForLargeR=2,
                nTopForSmallR=1,
                seedConfMinBottomRadius = 60. * u.mm,
                seedConfMaxZOrigin = 250. * u.mm,
                minImpactSeedConf = 10. * u.mm,
            )
        ),
        seedFinderOptionsArg = SEEDFINDEROPTIONSARG,
        seedFilterConfigArg = SeedFilterConfigArgNA60(
            seedConfirmation=False,
            maxSeedsPerSpMConf=1,
            maxQualitySeedsPerSpMConf=1,
            compatSeedLimit=2,  # added 4/10/23 (value not tuned)
            impactWeightFactor = 200.,
            zOriginWeightFactor = 200.,
            compatSeedWeight = 200.,
            seedWeightIncrement = 0,
            deltaYMin= 115* u.mm,

            verbose=verbose
        ),

        spacePointGridConfigArg = SPACEPOINTGRIDCONFIGARG,
        seedingAlgorithmConfigArg = SEEDINGALGORITHMCONFIGARG,
        outputDirRoot=outputDir,
        verbose=True,
        noGuessing=noGuessing,
        trkVtxOnly = trkVtxOnly,
        addDeltas = addDeltas,
        projective=projective,
        zPerigee = zPerigee
    )

    addCKFTracks(
        s,
        trackingGeometry,
        field,
        TRACKSELECTORCONFIG,
        CKFCONFIG,
        outputDirRoot=outputDir,
        suffixIn = "",
        suffixOut = ""
    )

    addAmbiguityResolution(
        s,
        AMBIGUITYRESOLUTIONCONFIG,
        outputDirRoot=outputDir,
        suffixIn = "",
        suffixOut = "ambi"
    )

    
    #########
    # STEP 2
    #########

    addUsedMeasurementsFilter(
        s,
        inputSourceLinks="sourcelinks",
        outputSourceLinks="llsourcelinks",     
        inputTracks="ambitracks"
    )

    addNewSeeding(
        s = s,
        trackingGeometry = trackingGeometry,
        field = field2,
        geoSelectionConfigFile = jsonSeed,
        #truthSeedRanges = TRUTHSEEDRANGE,
        initialSigmas = None,
        initialVarInflation = None,
        seedFinderConfigArg = SeedFinderConfigArgNA60(
            maxSeedsPerSpM=MaxSeedsPerSpMPrimary,
            cotThetaMax=CotThetaMax,
            sigmaScattering=SigmaScattering,
            radLengthPerSeed=RadLengthPerSeed,
            maxPtScattering=MaxPtScattering,
            # min and max R between Middle and Top SP
            deltaYTopSP=(10 * u.mm, 1000 * u.mm),
            # min and max R between Middle and Bottom SP
            deltaYBottomSP=(10 * u.mm, 100 * u.mm),
            impactMax=ImpactMax * u.mm,
            deltaZMax=5 * u.mm,  # was 5
            minPt=100 * u.MeV,
            #interactionPointCut=True,
            interactionPointCut=True,
            verbose=verbose,
            collisionRegion=(-ImpactMax * u.mm, ImpactMax * u.mm),  # 0.5
            # NOT USED??? seems to be used in Orthogonal seeding
            y=(0 * u.mm, 400 * u.mm),
            yMiddle=(160 * u.mm, 220 * u.mm), #use only layer 2 and 3
            seedConfirmation=False, #WE ONLY HAVE 1 REGION
            seedConfirmationRange=acts.SeedConfirmationRangeConfig(  # it should not be used....
                zMinSeedConf=-150 * u.mm,#no need to change this
                zMaxSeedConf=150 * u.mm,#no need to change this
                rMaxSeedConf=140 * u.mm,#r first layer
                nTopForLargeR=2,
                nTopForSmallR=1,
                seedConfMinBottomRadius = 60. * u.mm,
                seedConfMaxZOrigin = 250. * u.mm,
                minImpactSeedConf = 10. * u.mm,)
        ),
        seedFinderOptionsArg = SEEDFINDEROPTIONSARG,
                                            
        seedFilterConfigArg = SeedFilterConfigArgNA60(
            seedConfirmation=False,
            maxSeedsPerSpMConf=MaxSeedsPerSpMPrimary,
            maxQualitySeedsPerSpMConf=MaxSeedsPerSpMPrimary,
            compatSeedLimit=2,  # added 4/10/23 (value not tuned)
            impactWeightFactor = 200.,
            zOriginWeightFactor = 200.,
            compatSeedWeight = 200.,
            seedWeightIncrement = 0,
            deltaYMin= 115* u.mm,

            verbose=verbose
        ),

        spacePointGridConfigArg = SPACEPOINTGRIDCONFIGARG,
        seedingAlgorithmConfigArg = SEEDINGALGORITHMCONFIGARG,
        inputMeasurements="measurements",
        inputSourceLinks="llsourcelinks",
        outputDirRoot=outputDir,
        suffix = "ll",
        doTrkVtx = False #switch of the vertexing with the tracklets
    )

    addCKFTracks(
        s,
        trackingGeometry,
        field,
        CKFCONFIG,
        TRACKSELECTORCONFIG,
        outputDirRoot=outputDir,
        suffixIn = "ll",
        suffixOut = "ll"
    )

    addAmbiguityResolution(
        s,
        AMBIGUITYRESOLUTIONCONFIG,
        outputDirRoot=outputDir,
        suffixIn = "ll",
        suffixOut = "ambill"
    )
    
    addContainerMerger(
        s,
        inputTrackParameters=["ambitrackpars","ambilltrackpars"],
        outputTrackParameters="mergedtrackpars"   
    )
    
    addVertexFitting(
        s,
        field,
        vertexFinder=VertexFinder.Iterative,
        outputDirRoot=outputDir,
        suffixIn= "ambi",
        trackParameters="mergedtrackpars",
        outputProtoVertices= "protoverticesmerged",
        outputVertices= "fittedVerticesmerged",
        suffixOut = "_merged",
        logLevel=acts.logging.Level.DEBUG,
        significanceCutSeeding=significanceCutSeeding,
        maximumChi2cutForSeeding=maximumChi2cutForSeeding,
        maxVertices=maxVertices,
        createSplitVertices=createSplitVertices,
        splitVerticesTrkInvFraction=splitVerticesTrkInvFraction,
        reassignTracksAfterFirstFit=reassignTracksAfterFirstFit,
        doMaxTracksCut=doMaxTracksCut,
        maxTracks=maxTracks,
        cutOffTrackWeight=cutOffTrackWeight,
        cutOffTrackWeightReassign=cutOffTrackWeightReassign,
        rejectedFraction=rejectedFraction
        
    )
    

    
    
    #########
    # STEP 3
    #########

    addUsedMeasurementsFilter(
        s,
        inputSourceLinks="llsourcelinks",
        outputSourceLinks="llllsourcelinks",     
        inputTracks="ambilltracks"
    )

    addNewSeeding(
        s = s,
        trackingGeometry = trackingGeometry,
        field = field2,
        geoSelectionConfigFile = jsonSeed,
        #truthSeedRanges = TRUTHSEEDRANGE,
        initialSigmas = None,
        initialVarInflation = None,
        seedFinderConfigArg = SeedFinderConfigArgNA60(
            maxSeedsPerSpM=MaxSeedsPerSpMSecondary,
            cotThetaMax=CotThetaMax,
            sigmaScattering=SigmaScattering,
            radLengthPerSeed=RadLengthPerSeed,
            maxPtScattering=MaxPtScattering,
            # min and max R between Middle and Top SP
            deltaYTopSP=(10 * u.mm, 1000 * u.mm),
            # min and max R between Middle and Bottom SP
            deltaYBottomSP=(10 * u.mm, 300 * u.mm),
            impactMax=10000 * u.mm,
            deltaZMax=10000 * u.mm,  # was 5
            minPt=100 * u.MeV,
            interactionPointCut=False,
            verbose=verbose,
            collisionRegion=(-10000 * u.mm, 10000 * u.mm),  # 0.5
            # NOT USED??? seems to be used in Orthogonal seeding
            y=(0 * u.mm, 400 * u.mm),
            yMiddle=(139 * u.mm, 170 * u.mm), #use only layer 2 and 3
            seedConfirmation=SeedConfirmationSecondary, #WE ONLY HAVE 1 REGION
            seedConfirmationRange=acts.SeedConfirmationRangeConfig(  # it should not be used....
                zMinSeedConf=-150 * u.mm,#no need to change this
                zMaxSeedConf=150 * u.mm,#no need to change this
                rMaxSeedConf=140 * u.mm,#r first layer
                nTopForLargeR=2,
                nTopForSmallR=1,
                seedConfMinBottomRadius = 60. * u.mm,
                seedConfMaxZOrigin = 250. * u.mm,
                minImpactSeedConf = 10. * u.mm,
            )
        ),
        seedFinderOptionsArg = SEEDFINDEROPTIONSARG,
        seedFilterConfigArg = SeedFilterConfigArgNA60(
            seedConfirmation=SeedConfirmationSecondary,
            maxSeedsPerSpMConf=MaxSeedsPerSpMSecondary,
            maxQualitySeedsPerSpMConf=MaxSeedsPerSpMSecondary,
            compatSeedLimit=2,  # added 4/10/23 (value not tuned)
            impactWeightFactor = 0.,
            zOriginWeightFactor = 0.,
            compatSeedWeight = 200.,
            seedWeightIncrement = 0,
            deltaYMin= 115* u.mm,

            verbose=verbose
        ),

        spacePointGridConfigArg = SPACEPOINTGRIDCONFIGARG,
        seedingAlgorithmConfigArg = SEEDINGALGORITHMCONFIGARG,
        inputMeasurements="measurements",
        inputSourceLinks="llllsourcelinks",
        outputDirRoot=outputDir,
        suffix = "llll",
        doTrkVtx = False #switch of the vertexing with the tracklets
    )

    addCKFTracks(
        s,
        trackingGeometry,
        field,
        CKFCONFIG,
        TRACKSELECTORCONFIG,
        outputDirRoot=outputDir,
        suffixIn = "llll",
        suffixOut = "llll"
    )

    addAmbiguityResolution(
        s,
        AMBIGUITYRESOLUTIONCONFIG,
        outputDirRoot=outputDir,
        suffixIn = "llll",
        suffixOut = "ambillll"
    )

    #########
    # STEP 4
    #########

    addUsedMeasurementsFilter(
        s,
        inputSourceLinks="llllsourcelinks",
        outputSourceLinks="llllllsourcelinks",     
        inputTracks="ambilllltracks"
    )

    addNewSeeding(
        s = s,
        trackingGeometry = trackingGeometry,
        field = field2,
        geoSelectionConfigFile = jsonSeed,
        #truthSeedRanges = TRUTHSEEDRANGE,
        initialSigmas = None,
        initialVarInflation = None,
        seedFinderConfigArg = SeedFinderConfigArgNA60(
            maxSeedsPerSpM=MaxSeedsPerSpMSecondary,
            cotThetaMax=CotThetaMax,
            sigmaScattering=SigmaScattering,
            radLengthPerSeed=RadLengthPerSeed,
            maxPtScattering=MaxPtScattering,
            # min and max R between Middle and Top SP
            deltaYTopSP=(10 * u.mm, 1000 * u.mm),
            # min and max R between Middle and Bottom SP
            deltaYBottomSP=(10 * u.mm, 300 * u.mm),
            impactMax=10000 * u.mm,
            deltaZMax=10000 * u.mm,  # was 5
            minPt=100 * u.MeV,
            interactionPointCut=False,
            verbose=verbose,
            collisionRegion=(-10000 * u.mm, 10000 * u.mm),  # 0.5
            # NOT USED??? seems to be used in Orthogonal seeding
            y=(0 * u.mm, 400 * u.mm),
            yMiddle=(170 * u.mm, 220 * u.mm), #use only layer 2 and 3
            seedConfirmation=SeedConfirmationSecondary, #WE ONLY HAVE 1 REGION
            seedConfirmationRange=acts.SeedConfirmationRangeConfig(  # it should not be used....
                zMinSeedConf=-150 * u.mm,#no need to change this
                zMaxSeedConf=150 * u.mm,#no need to change this
                rMaxSeedConf=140 * u.mm,#r first layer
                nTopForLargeR=2,
                nTopForSmallR=1,
                seedConfMinBottomRadius = 60. * u.mm,
                seedConfMaxZOrigin = 250. * u.mm,
                minImpactSeedConf = 10. * u.mm,
            )
        ),
        seedFinderOptionsArg = SEEDFINDEROPTIONSARG,
        seedFilterConfigArg = SeedFilterConfigArgNA60(
            seedConfirmation=SeedConfirmationSecondary,
            maxSeedsPerSpMConf=MaxSeedsPerSpMSecondary,
            maxQualitySeedsPerSpMConf=MaxSeedsPerSpMSecondary,
            compatSeedLimit=2,  # added 4/10/23 (value not tuned)
            impactWeightFactor = 0.,
            zOriginWeightFactor = 0.,
            compatSeedWeight = 200.,
            seedWeightIncrement = 0,
            deltaYMin= 115* u.mm,

            verbose=verbose
        ),

        spacePointGridConfigArg = SPACEPOINTGRIDCONFIGARG,
        seedingAlgorithmConfigArg = SEEDINGALGORITHMCONFIGARG,
        inputMeasurements="measurements",
        inputSourceLinks="llllllsourcelinks",
        outputDirRoot=outputDir,
        suffix = "llllll",
        doTrkVtx = False #switch of the vertexing with the tracklets
    )

    addCKFTracks(
        s,
        trackingGeometry,
        field,
        CKFCONFIG,
        TRACKSELECTORCONFIG,
        outputDirRoot=outputDir,
        suffixIn = "llllll",
        suffixOut = "llllll"
    )

    addAmbiguityResolution(
        s,
        AMBIGUITYRESOLUTIONCONFIG,
        outputDirRoot=outputDir,
        suffixIn = "llllll",
        suffixOut = "ambillllll"
    )
    

    addVertexFitting(
        s,
        field,
        vertexFinder=VertexFinder.Iterative,
        outputDirRoot=outputDir,
        suffixIn= "ambi",
        trackParameters="ambitrackpars",
        outputProtoVertices= "protoverticesstep1",
        outputVertices= "fittedVerticesstep1",
        suffixOut = "_step1",
        logLevel=acts.logging.Level.DEBUG,
        significanceCutSeeding=significanceCutSeeding,
        maximumChi2cutForSeeding=maximumChi2cutForSeeding,
        maxVertices=maxVertices,
        createSplitVertices=createSplitVertices,
        splitVerticesTrkInvFraction=splitVerticesTrkInvFraction,
        reassignTracksAfterFirstFit=reassignTracksAfterFirstFit,
        doMaxTracksCut=doMaxTracksCut,
        maxTracks=maxTracks,
        cutOffTrackWeight=cutOffTrackWeight,
        cutOffTrackWeightReassign=cutOffTrackWeightReassign,
        rejectedFraction=rejectedFraction
        
    )
    """


    addContainerMerger(
        s,
        inputTrackParameters=["ambitrackpars","ambilltrackpars","ambilllltrackpars","ambilllllltrackpars"],
        outputTrackParameters="mergedlltrackpars"   
    )
    addSecondaryVertexFitting(
        s,
        field,
        trackParameters="mergedlltrackpars",
        outputProtoVertices = "protoverticestest",
        outputMasses = "masses",
        outputDirRoot=outputDir,
        logLevel=acts.logging.Level.DEBUG
    )
    """
    return s

def copy_py_file(source_file, destination_dir):
    # Check if the source file is a Python file
    if not source_file.endswith('.py'):
        print("Source file is not a Python file (.py)")
        return
    
    # Check if the destination directory exists
    if not os.path.exists(destination_dir):
        print("Destination directory does not exist.")
        return
    
    # Get the base filename of the source file
    file_name = os.path.basename(source_file)
    
    # Construct the destination path
    destination_path = os.path.join(destination_dir, file_name)
    
    # Copy the file
    shutil.copy(source_file, destination_path)
    print(f"File '{file_name}' copied to '{destination_dir}' successfully.")

if "__main__" == __name__:
    options = getArgumentParser().parse_args()

    # "events_full_sim"#"fullchain_input"
    Eint = options.eint

    event_type = "" if not options.bkg else "_noBkg"

    if options.fullev:
        event_type += "_FullEv"

    dir = "event_generation/events_"+str(Eint)+"GeV"+event_type


    target_option = "_onlyTarget_"+str(options.target) if options.target!=-1 else "_realTarget"
    has_sec = ""
    if options.sec:
        has_sec = "_Sec"
    if options.osec:
        has_sec = "_noBkg_Sec"
    if options.zz:
        beamsize = "0.000000"
    else:
        beamsize = "0.500000"
    if options.periferal != 1:
        dir = "event_generation/events_40GeV"+has_sec+target_option+"_beamSigma_"+beamsize+"_periferal_factor_"+str(options.periferal)+("00000" if options.periferal >= 0.1 else "0000")
    else:
        dir = "event_generation/events_40GeV"+has_sec+target_option+"_beamSigma_"+beamsize+""

    #dir ="event_generation/events_40GeV_Sec_periferal_factor_0.100000"

    inputDir = pathlib.Path.cwd() / dir
    suffix = "_newSeeding"
    suffix += "_standardSeeding"
        
    if options.primary:
        suffix += "_noPrimary"
    if options.secondary:
        suffix += "_noSecondary"


    if options.geant4:
        suffix += "_geant4"
    if options.dead:
        suffix += "_deadZones"
    if options.fast:
        suffix += "_fastSim"

    suffix += "_maxSeedSpMPrim"+str(options.sf_maxSeedsPerSpM_primary)
    suffix += "_maxSeedSpMSec"+str(options.sf_maxSeedsPerSpM_secondary)

    if options.sf_seedConfirmation_primary:
        suffix += "_confirmation_primary"
    if options.sf_seedConfirmation_secondary:
        suffix += "_confirmation_secondary"

    suffix += "_ImpMax"+str(options.sf_impactMax)
    suffix += "_dZMax"+str(options.sf_deltaZMax)


    if options.sf_numMeasurementsCutOff:
        suffix += "_branch"+str(options.sf_numMeasurementsCutOff)
    if options.eff is not None:
        suffix += "_eff"+str(options.eff)
    if options.noguessing:
        suffix += "_noGuessing"
        suffix+= "_z"+str(options.zperigee)
    if options.trkVtxOnly:
        suffix += "_trkVtxOnly"
    if options.addDeltas:
        suffix += "_addDeltas"
    if options.projective:
        suffix += "_projective"
    if options.periferal != 1:
        suffix += "_periferal_factor_"+str(options.periferal)
    suffix += "_onlyTarget_"+str(options.target) if options.target!=-1 else "_realTarget"
    if options.sec:
        suffix += "_Sec"
    if options.osec:
        suffix += "_onlySec"
    if options.zz:
        suffix += "_beam0.0"
    else:
        suffix += "_beam0.5"
    if options.periferal != 1:
        dir = "event_generation/events_40GeV"+has_sec+target_option+"_beamSigma_"+beamsize+"_periferal_factor_"+str(options.periferal)+("00000" if options.periferal >= 0.1 else "0000")
    else:
        dir = "event_generation/events_40GeV"+has_sec+target_option+"_beamSigma_"+beamsize+""

    if options.sf_maxVertices != 1:
        suffix += "_nVtxMax"+str(options.sf_maxVertices)

    suffix += "_twosteps_rej"+str(options.sf_rejectedFraction)
    suffix += "_suffix"

    current_dir = pathlib.Path.cwd()
    if options.outdir:
        outputDir = options.outdir
    else:
        outputDir = str(current_dir / ("output/output_"+str(Eint)+"GeV" + event_type + suffix))

    matDeco = acts.IMaterialDecorator.fromFile(
        "geometry/geomVTNA60p/material-map_VTNA60p.json")
    jsonFile = "geometry/geomVTNA60p/tgeo-config_VTNA60p.json"
    tgeo_fileName = "geometry/geomVTNA60p/geom_VTNA60p.root"
    jsonDigi = "geometry/geomVTNA60p/digismear0.005.json"

    jsonSeed = "geometry/geomVTNA60p/seed_config_2_4_6_8_10.json"
    print(jsonSeed)
    logLevel = acts.logging.INFO
    customLogLevel = acts.examples.defaultLogging(logLevel=logLevel)

    detector, trackingGeometry, decorators = TGeoDetector.create(jsonFile=str(jsonFile),
                                                                 fileName=str(tgeo_fileName),
                                                                 surfaceLogLevel=customLogLevel(),
                                                                 layerLogLevel=customLogLevel(),
                                                                 volumeLogLevel=customLogLevel(),
                                                                 mdecorator=matDeco,
                                                                )
    print(outputDir if options.outdir == None else options.outdir)

    start_time = time.time()

    runVTchain(
        detector,
        trackingGeometry,
        inputDir=inputDir,
        outputDir=outputDir if options.outdir == None else options.outdir,
        jsonDigi=jsonDigi,
        jsonSeed=jsonSeed,
        useGeant4=options.geant4,
        applyDeadZones=options.dead,
        applyFastSimSelections=options.fast,
        truthSeeding=False,
        NumEvents=options.nEvts,

        MaxSeedsPerSpMPrimary=options.sf_maxSeedsPerSpM_primary,
        MaxSeedsPerSpMSecondary=options.sf_maxSeedsPerSpM_secondary,
        SeedConfirmationPrimary=options.sf_seedConfirmation_primary,
        SeedConfirmationSecondary=options.sf_seedConfirmation_secondary,
        ImpactMax=options.sf_impactMax,
        DeltaZMax=options.sf_deltaZMax,

        NumMeasurementsCutOff=options.sf_numMeasurementsCutOff,
        verbose=options.verbose,
        Efficiency=options.eff,
        KeepPrimary=not options.primary,
        KeepSecondary=not options.secondary,
        noGuessing=options.noguessing,
        trkVtxOnly=options.trkVtxOnly,
        addDeltas=options.addDeltas,
        projective=options.projective,
        zPerigee=options.zperigee,
        significanceCutSeeding=options.sf_significanceCutSeeding,
        maximumChi2cutForSeeding=options.sf_maximumChi2cutForSeeding,
        maxVertices=options.sf_maxVertices,
        createSplitVertices=options.sf_createSplitVertices,
        splitVerticesTrkInvFraction=options.sf_splitVerticesTrkInvFraction,
        reassignTracksAfterFirstFit=options.sf_reassignTracksAfterFirstFit,
        doMaxTracksCut=options.sf_doMaxTracksCut,
        maxTracks=options.sf_maxTracks,
        cutOffTrackWeight=options.sf_cutOffTrackWeight,
        cutOffTrackWeightReassign=options.sf_cutOffTrackWeightReassign,
        rejectedFraction=options.sf_rejectedFraction

    ).run()
    end_time = time.time()

    execution_time = end_time - start_time
    print("NA60+_Summary_Run time", execution_time)

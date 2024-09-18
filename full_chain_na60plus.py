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
    addSeeding,
    TruthSeedRanges,
    SeedFinderConfigArgNA60,
    SeedFinderOptionsArgNA60,
    SeedFilterConfigArgNA60,
    SeedFinderConfigArg,
    SeedFinderOptionsArg,
    SeedFilterConfigArg,
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
    TruthEstimatedSeedingAlgorithmConfigArg,
    SeedingAlgorithm,
    addMatching,
)
from acts.examples import TGeoDetector


def getArgumentParser():
    """Get arguments from command line"""
    parser = argparse.ArgumentParser(description="Command line arguments for CKF")

    parser.add_argument(
        "-o",
        "--output",
        dest="outdir",
        help="Output directory for new ntuples",
        default=None,
    )
    parser.add_argument(
        "-n",
        "--nEvents",
        dest="nEvts",
        help="Number of events to run over",
        default=1,
        type=int,
    )

    parser.add_argument("-g", "--geant4", help="Use Geant4", action="store_true")

    parser.add_argument("-ve", "--verbose", help="Use Geant4", action="store_true")

    parser.add_argument("-d", "--dead", help="Apply dead zones", action="store_true")
    parser.add_argument(
        "-f", "--fast", help="Apply fast sim selections", action="store_true"
    )

    parser.add_argument(
        "-st", "--primary", help="No primary particles", action="store_true"
    )

    parser.add_argument(
        "-nd", "--secondary", help="No secondary particles", action="store_true"
    )

    parser.add_argument(
        "-gu", "--noguessing", help="No secondary particles", action="store_true"
    )

    parser.add_argument(
        "-sec", "--sec", help="No secondary particles", action="store_true"
    )

    parser.add_argument(
        "-osec", "--osec", help="No secondary particles", action="store_true"
    )

    parser.add_argument(
        "-zz", "--zz", help="No secondary particles", action="store_true"
    )

    parser.add_argument(
        "-t0", "--t0", help="No secondary particles", action="store_true"
    )

    parser.add_argument(
        "--se",
        "--seed",
        dest="seed",
        help="Use particle gun",
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
        "-e",
        "--eint",
        dest="eint",
        help="Collision energy",
        type=int,
        default=40,
    )

    parser.add_argument(
        "-vz",
        "--vz",
        dest="vz",
        help="perigee z for matching",
        type=float,
        default=400,
    )

    parser.add_argument(
        "-fullev",
        "--fullev",
        dest="fullev",
        help="Events with also K0S and Lambda0s",
        action="store_true",
    )

    parser.add_argument(
        "-bkg", "--bkg", dest="bkg", help="Remove bkg", action="store_true"
    )

    parser.add_argument(
        "-trk",
        "--trkVtxOnly",
        dest="trkVtxOnly",
        help="Remove bkg",
        action="store_true",
    )

    parser.add_argument(
        "-del", "--addDeltas", dest="addDeltas", help="Remove bkg", action="store_true"
    )

    parser.add_argument(
        "-per",
        "--periferal",
        dest="periferal",
        type=float,
        default=1,
    )

    parser.add_argument(
        "-pro",
        "--projective",
        dest="projective",
        help="Remove bkg",
        action="store_true",
    )
    ##########################################
    ##
    ##
    ##
    #################
    #########################
    parser.add_argument(
        "-ef",
        "--eff",
        dest="eff",
        help="Collision energy",
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
        default=20,
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
        default=False,
    )

    parser.add_argument(
        "--sf_seedConfirmation_secondary",
        dest="sf_seedConfirmation_secondary",
        help="Use seed confirmation",
        type=bool,
        default=False,
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

    # for vtx optimization
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
        default=10,  # 0.06758286918225764,
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


def runFullChain(
    detectorVT,
    trackingGeometryVT,
    detectorMS,
    trackingGeometryMS,
    inputDirVT: Path,
    inputDirMS: Path,
    outputDir: Path,
    jsonDigiVT="jsonDigiVT.json",
    jsonSeedVT="jsonSeedVT.json",
    jsonDigiMS="jsonDigiMS.json",
    jsonSeedMS="jsonSeedMS.json",
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
    CotThetaMax=5.8,
    SigmaScattering=7.3,
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
    perigeeZ = 400
):

    field = acts.examples.MagneticFieldMapXyz("bfield/BFieldNA60plus_longsetup.txt")
    field2 = acts.examples.MagneticFieldMapXyz(
        "bfield/BFieldNA60plus_ZRotated_longsetup.txt"
    )

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
        inputDir=inputDirVT,
        outputDirRoot=outputDir,
        #det_suffix = ""
    )
    addParticleReader(
        s,
        inputDir=inputDirMS,
        outputDirRoot=outputDir,
        det_suffix = "_ms"
    )

    addFatras(
        s,
        trackingGeometryVT,
        field,
        rnd,
        preSelectParticles=ParticleSelectorConfig(
            eta=(None, None),
            pt=(100 * u.MeV, None),
            removeNeutral=True,
        ),
        outputDirRoot=outputDir,
        outputParticlesInitial="particles_initial_vt",
        outputParticlesFinal="particles_final_vt",
        outputSimHits="simhits_vt",
        det_suffix="_vt",
    )

    addFatras(
        s,
        trackingGeometryMS,
        field,
        rnd,
        preSelectParticles=ParticleSelectorConfig(
            eta=(None, None),
            pt=(100 * u.MeV, None),
            removeNeutral=True,
        ),
        outputDirRoot=outputDir,
        inputParticles = "particles_input_ms",
        outputParticlesInitial="particles_initial_ms",
        outputParticlesFinal="particles_final_ms",
        outputSimHits="simhits_ms",
        det_suffix="_ms",
    )

    addDigitization(
        s,
        trackingGeometryVT,
        field,
        digiConfigFile=jsonDigiVT,
        outputDirRoot=outputDir,
        efficiency=(
            1
            if Efficiency is not None
            else 0.99 if applyDeadZones or applyFastSimSelections else 1
        ),
        applyDeadAreas=applyDeadZones,
        applyFastSimSelections=applyFastSimSelections,
        applyReadout=applyReadout,
        applyBackbone=applyBackbone,
        applyHole=applyHole,
        applyEndcapShort=applyEndcapShort,
        applyEndcapLong=applyEndcapLong,
        rnd=rnd,
        suffix="_vt",
    )

    addDigitization(
        s,
        trackingGeometryMS,
        field,
        digiConfigFile=jsonDigiMS,
        outputDirRoot=outputDir,
        efficiency=(
            1
            if Efficiency is not None
            else 0.99 if applyDeadZones or applyFastSimSelections else 1
        ),
        applyDeadAreas=applyDeadZones,
        applyFastSimSelections=applyFastSimSelections,
        applyReadout=applyReadout,
        applyBackbone=applyBackbone,
        applyHole=applyHole,
        applyEndcapShort=applyEndcapShort,
        applyEndcapLong=applyEndcapLong,
        rnd=rnd,
        suffix="_ms",
    )

    TRUTHSEEDRANGE = TruthSeedRanges(
        pt=(100 * u.MeV, None),
        eta=(None, None),
        nHits=(4, None),
        keep=(KeepPrimary, KeepSecondary),
    )

    SPACEPOINTGRIDCONFIGARG = SpacePointGridConfigArg(  # NO NEED TO CHANGE THIS
        rMax=400 * u.mm, zBinEdges=[-150.0, -75.0, 75.0, 150.0]
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
        maximumSharedHits=1,
        nMeasurementsMin=4,
        maximumIterations=1000000,
    )

    SEEDINGALGORITHMCONFIGARG = SeedingAlgorithmConfigArg(  # NO NEED TO CHANGE THIS
        zBinNeighborsTop=[[0, 1], [-1, 1], [-1, 0]],
        zBinNeighborsBottom=[[0, 1], [-1, 1], [-1, 0]],
    )

    SEEDFINDEROPTIONSARG = SeedFinderOptionsArgNA60(  # NO NEED TO CHANGE THIS
        beamPos=(0 * u.mm, 0 * u.mm),
        bFieldInZ=1.5 * u.T,
    )

    #########
    # STEP 1
    #########
    """
    
            deltaYTopSP=(10 * u.mm, 1000 * u.mm),
            # min and max R between Middle and Bottom SP
            deltaYBottomSP=(10 * u.mm, 100 * u.mm),
            impactMax=ImpactMax * u.mm,
            deltaZMax=5 * u.mm,  # was 5
            minPt=100 * u.MeV,
            # interactionPointCut=True,
            interactionPointCut=True,
            verbose=verbose,
            collisionRegion=(-ImpactMax * u.mm, ImpactMax * u.mm),  # 0.5
            # NOT USED??? seems to be used in Orthogonal seeding
            y=(0 * u.mm, 400 * u.mm),
            yMiddle=(160 * u.mm, 220 * u.mm),  # use only layer 2 and 3

    """
    addNewSeeding(
        s=s,
        trackingGeometry=trackingGeometryVT,
        field=field2,
        geoSelectionConfigFile=jsonSeedVT,
        truthSeedRanges=TRUTHSEEDRANGE,
        initialSigmas=None,
        initialVarInflation=None,
        seedFinderConfigArg=SeedFinderConfigArgNA60(
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
            # interactionPointCut=True,
            interactionPointCut=True,
            verbose=verbose,
            collisionRegion=(-ImpactMax * u.mm, ImpactMax * u.mm),  # 0.5
            # NOT USED??? seems to be used in Orthogonal seeding
            y=(0 * u.mm, 400 * u.mm),
            yMiddle=(139 * u.mm, 170 * u.mm),  # use only layer 2 and 3
            #yMiddle=(139 * u.mm, 220 * u.mm),  # use only layer 2 and 3
            seedConfirmation=False,  # WE ONLY HAVE 1 REGION
            seedConfirmationRange=acts.SeedConfirmationRangeConfig(  # it should not be used....
                zMinSeedConf=-150 * u.mm,  # no need to change this
                zMaxSeedConf=150 * u.mm,  # no need to change this
                rMaxSeedConf=140 * u.mm,  # r first layer
                nTopForLargeR=2,
                nTopForSmallR=1,
                seedConfMinBottomRadius=60.0 * u.mm,
                seedConfMaxZOrigin=250.0 * u.mm,
                minImpactSeedConf=10.0 * u.mm,
            ),
        ),
        seedFinderOptionsArg=SEEDFINDEROPTIONSARG,
        seedFilterConfigArg=SeedFilterConfigArgNA60(
            seedConfirmation=False,
            maxSeedsPerSpMConf=1,
            maxQualitySeedsPerSpMConf=1,
            compatSeedLimit=2,  # added 4/10/23 (value not tuned)
            impactWeightFactor=200.0,
            zOriginWeightFactor=200.0,
            compatSeedWeight=200.0,
            seedWeightIncrement=0,
            deltaYMin=115 * u.mm,
            verbose=verbose,
        ),
        spacePointGridConfigArg=SPACEPOINTGRIDCONFIGARG,
        seedingAlgorithmConfigArg=SEEDINGALGORITHMCONFIGARG,
        outputDirRoot=outputDir,
        verbose=False,
        noGuessing=False,  # noGuessing,
        trkVtxOnly=trkVtxOnly,
        addDeltas=addDeltas,
        projective=projective,
        zPerigee=zPerigee,
        inputSourceLinks="sourcelinks_vt",
        inputMeasurements="measurements_vt",
        det_suffix="_vt",
        inputParticles="particles_vt",
    )

    addCKFTracks(
        s,
        trackingGeometryVT,
        field,
        TRACKSELECTORCONFIG,
        CKFCONFIG,
        outputDirRoot=outputDir,
        inputSourceLinks="sourcelinks_vt",
        inputMeasurements="measurements_vt",
        suffixIn="",
        suffixOut="",
        det_suffix="_vt",
    )

    addAmbiguityResolution(
        s,
        AMBIGUITYRESOLUTIONCONFIG,
        outputDirRoot=outputDir,
        suffixIn="",
        suffixOut="ambi",
        det_suffix="_vt",
    )

    #########
    # STEP 2
    #########

    addUsedMeasurementsFilter(
        s,
        inputSourceLinks="sourcelinks_vt",
        outputSourceLinks="llsourcelinks_vt",
        inputTracks="ambitracks",
    )

    addNewSeeding(
        s=s,
        trackingGeometry=trackingGeometryVT,
        field=field2,
        geoSelectionConfigFile=jsonSeedVT,
        # truthSeedRanges = TRUTHSEEDRANGE,
        initialSigmas=None,
        initialVarInflation=None,
        seedFinderConfigArg=SeedFinderConfigArgNA60(
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
            # interactionPointCut=True,
            interactionPointCut=True,
            verbose=verbose,
            collisionRegion=(-ImpactMax * u.mm, ImpactMax * u.mm),  # 0.5
            # NOT USED??? seems to be used in Orthogonal seeding
            y=(0 * u.mm, 400 * u.mm),
            yMiddle=(160 * u.mm, 220 * u.mm),  # use only layer 2 and 3
            seedConfirmation=False,  # WE ONLY HAVE 1 REGION
            seedConfirmationRange=acts.SeedConfirmationRangeConfig(  # it should not be used....
                zMinSeedConf=-150 * u.mm,  # no need to change this
                zMaxSeedConf=150 * u.mm,  # no need to change this
                rMaxSeedConf=140 * u.mm,  # r first layer
                nTopForLargeR=2,
                nTopForSmallR=1,
                seedConfMinBottomRadius=60.0 * u.mm,
                seedConfMaxZOrigin=250.0 * u.mm,
                minImpactSeedConf=10.0 * u.mm,
            ),
        ),
        seedFinderOptionsArg=SEEDFINDEROPTIONSARG,
        seedFilterConfigArg=SeedFilterConfigArgNA60(
            seedConfirmation=False,
            maxSeedsPerSpMConf=MaxSeedsPerSpMPrimary,
            maxQualitySeedsPerSpMConf=MaxSeedsPerSpMPrimary,
            compatSeedLimit=2,  # added 4/10/23 (value not tuned)
            impactWeightFactor=200.0,
            zOriginWeightFactor=200.0,
            compatSeedWeight=200.0,
            seedWeightIncrement=0,
            deltaYMin=115 * u.mm,
            verbose=verbose,
        ),
        spacePointGridConfigArg=SPACEPOINTGRIDCONFIGARG,
        seedingAlgorithmConfigArg=SEEDINGALGORITHMCONFIGARG,
        inputSourceLinks="llsourcelinks_vt",
        outputDirRoot=outputDir,
        suffix="ll",
        inputMeasurements="measurements_vt",
        det_suffix="_vt",
        doTrkVtx=False,  # switch of the vertexing with the tracklets
        inputParticles="particles_vt",
    )

    addCKFTracks(
        s,
        trackingGeometryVT,
        field,
        TRACKSELECTORCONFIG,
        CKFCONFIG,
        outputDirRoot=outputDir,
        inputSourceLinks="llsourcelinks_vt",
        inputMeasurements="measurements_vt",
        suffixIn="ll",
        suffixOut="ll",
        det_suffix="_vt",
    )

    addAmbiguityResolution(
        s,
        AMBIGUITYRESOLUTIONCONFIG,
        outputDirRoot=outputDir,
        suffixIn="ll",
        suffixOut="ambill",
        det_suffix="_vt",
    )

    #########
    # STEP 3
    #########


    addUsedMeasurementsFilter(
        s,
        inputSourceLinks="llsourcelinks_vt",
        outputSourceLinks="llllsourcelinks_vt",
        inputTracks="ambilltracks",
    )

    addNewSeeding(
        s=s,
        trackingGeometry=trackingGeometryVT,
        field=field2,
        geoSelectionConfigFile=jsonSeedVT,
        # truthSeedRanges = TRUTHSEEDRANGE,
        initialSigmas=None,
        initialVarInflation=None,
        seedFinderConfigArg=SeedFinderConfigArgNA60(
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
            yMiddle=(139 * u.mm, 170 * u.mm),  # use only layer 2 and 3
            seedConfirmation=SeedConfirmationSecondary,  # WE ONLY HAVE 1 REGION
            seedConfirmationRange=acts.SeedConfirmationRangeConfig(  # it should not be used....
                zMinSeedConf=-150 * u.mm,  # no need to change this
                zMaxSeedConf=150 * u.mm,  # no need to change this
                rMaxSeedConf=140 * u.mm,  # r first layer
                nTopForLargeR=2,
                nTopForSmallR=1,
                seedConfMinBottomRadius=60.0 * u.mm,
                seedConfMaxZOrigin=250.0 * u.mm,
                minImpactSeedConf=10.0 * u.mm,
            ),
        ),
        seedFinderOptionsArg=SEEDFINDEROPTIONSARG,
        seedFilterConfigArg=SeedFilterConfigArgNA60(
            seedConfirmation=SeedConfirmationSecondary,
            maxSeedsPerSpMConf=MaxSeedsPerSpMSecondary,
            maxQualitySeedsPerSpMConf=MaxSeedsPerSpMSecondary,
            compatSeedLimit=2,  # added 4/10/23 (value not tuned)
            impactWeightFactor=0.0,
            zOriginWeightFactor=0.0,
            compatSeedWeight=200.0,
            seedWeightIncrement=0,
            deltaYMin=115 * u.mm,
            verbose=verbose,
        ),
        spacePointGridConfigArg=SPACEPOINTGRIDCONFIGARG,
        seedingAlgorithmConfigArg=SEEDINGALGORITHMCONFIGARG,
        inputSourceLinks="llllsourcelinks_vt",
        outputDirRoot=outputDir,
        suffix="llll",
        inputMeasurements="measurements_vt",
        det_suffix="_vt",
        doTrkVtx=False,  # switch of the vertexing with the tracklets
        inputParticles="particles_vt",
    )

    addCKFTracks(
        s,
        trackingGeometryVT,
        field,
        CKFCONFIG,
        TRACKSELECTORCONFIG,
        outputDirRoot=outputDir,
        inputSourceLinks="llllsourcelinks_vt",
        inputMeasurements="measurements_vt",
        suffixIn="llll",
        suffixOut="llll",
        det_suffix="_vt",
    )

    addAmbiguityResolution(
        s,
        AMBIGUITYRESOLUTIONCONFIG,
        outputDirRoot=outputDir,
        suffixIn="llll",
        suffixOut="ambillll",
        det_suffix="_vt",
    )

    #########
    # STEP 4
    #########

    addUsedMeasurementsFilter(
        s,
        inputSourceLinks="llllsourcelinks_vt",
        outputSourceLinks="llllllsourcelinks_vt",
        inputTracks="ambilllltracks",
    )

    addNewSeeding(
        s=s,
        trackingGeometry=trackingGeometryVT,
        field=field2,
        geoSelectionConfigFile=jsonSeedVT,
        # truthSeedRanges = TRUTHSEEDRANGE,
        initialSigmas=None,
        initialVarInflation=None,
        seedFinderConfigArg=SeedFinderConfigArgNA60(
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
            yMiddle=(170 * u.mm, 220 * u.mm),  # use only layer 2 and 3
            seedConfirmation=SeedConfirmationSecondary,  # WE ONLY HAVE 1 REGION
            seedConfirmationRange=acts.SeedConfirmationRangeConfig(  # it should not be used....
                zMinSeedConf=-150 * u.mm,  # no need to change this
                zMaxSeedConf=150 * u.mm,  # no need to change this
                rMaxSeedConf=140 * u.mm,  # r first layer
                nTopForLargeR=2,
                nTopForSmallR=1,
                seedConfMinBottomRadius=60.0 * u.mm,
                seedConfMaxZOrigin=250.0 * u.mm,
                minImpactSeedConf=10.0 * u.mm,
            ),
        ),
        seedFinderOptionsArg=SEEDFINDEROPTIONSARG,
        seedFilterConfigArg=SeedFilterConfigArgNA60(
            seedConfirmation=SeedConfirmationSecondary,
            maxSeedsPerSpMConf=MaxSeedsPerSpMSecondary,
            maxQualitySeedsPerSpMConf=MaxSeedsPerSpMSecondary,
            compatSeedLimit=2,  # added 4/10/23 (value not tuned)
            impactWeightFactor=0.0,
            zOriginWeightFactor=0.0,
            compatSeedWeight=200.0,
            seedWeightIncrement=0,
            deltaYMin=115 * u.mm,
            verbose=verbose,
        ),
        spacePointGridConfigArg=SPACEPOINTGRIDCONFIGARG,
        seedingAlgorithmConfigArg=SEEDINGALGORITHMCONFIGARG,
        inputSourceLinks="llllllsourcelinks_vt",
        outputDirRoot=outputDir,
        suffix="llllll",
        inputMeasurements="measurements_vt",
        det_suffix="_vt",
        doTrkVtx=False,  # switch of the vertexing with the tracklets
        inputParticles="particles_vt",
    )

    addCKFTracks(
        s,
        trackingGeometryVT,
        field,
        CKFCONFIG,
        TRACKSELECTORCONFIG,
        outputDirRoot=outputDir,
        suffixIn="llllll",
        suffixOut="llllll",
        inputSourceLinks="llllllsourcelinks_vt",
        inputMeasurements="measurements_vt",
        det_suffix="_vt",
    )

    addAmbiguityResolution(
        s,
        AMBIGUITYRESOLUTIONCONFIG,
        outputDirRoot=outputDir,
        suffixIn="llllll",
        suffixOut="ambillllll",
        det_suffix="_vt",
    )

    #########
    # MUON RECONSTRUCTION
    #########

    addSeeding(
        s,
        trackingGeometryMS,
        field2,
        TruthSeedRanges(
            pt=(100 * u.MeV, None),
            eta=(None, None),
            nHits=(6, None),
            keep=(KeepPrimary, KeepSecondary),
        ),
        SeedFinderConfigArg(
            maxSeedsPerSpM=5,
            cotThetaMax=1000000,  # 0.25,
            deltaZMax=100000,
            sigmaScattering=5.0,
            radLengthPerSeed=0.01,
            minPt=100 * u.MeV,
            impactMax=1000000 * u.mm,
            interactionPointCut=False,
            useVariableMiddleSPRange=False,
            deltaRMiddleSPRange=(
                0 * u.mm,
                10000000 * u.mm,
            ),  # not useful if useVariableMiddleSPRange=False
            deltaRTopSP=(
                350 * u.mm,
                40000 * u.mm,
            ),  # min and max R between Middle and Top SP
            deltaRBottomSP=(
                350 * u.mm,
                40000 * u.mm,
            ),  # min and max R between Middle and Bottom SP
            collisionRegion=(-1000 * u.mm, 1000 * u.mm),
            r=(20 * u.mm, 9000 * u.mm),
            rMiddle=(20 * u.mm, 9000 * u.mm),
            seedConfirmation=False,
            forwardSeedConfirmationRange=acts.SeedConfirmationRangeConfig(
                zMinSeedConf=150 * u.mm,
                zMaxSeedConf=0 * u.mm,
                rMaxSeedConf=9000 * u.mm,
                nTopForLargeR=1,
                nTopForSmallR=1,
            ),
            centralSeedConfirmationRange=acts.SeedConfirmationRangeConfig(
                zMinSeedConf=7 * u.mm,
                zMaxSeedConf=-7 * u.mm,
                rMaxSeedConf=9000 * u.mm,
                nTopForLargeR=1,
                nTopForSmallR=1,
            ),
        ),
        SeedFinderOptionsArg(
            beamPos=(0 * u.mm, 0 * u.mm),
            bFieldInZ=0 * u.T,
        ),
        SeedFilterConfigArg(  # not used, why?
            seedConfirmation=False,
            maxSeedsPerSpMConf=2,
            maxQualitySeedsPerSpMConf=1,
            compatSeedLimit=0,  # added 4/10/23 (value not tuned)
        ),
        SpacePointGridConfigArg(
            rMax=9000 * u.mm,
            zBinEdges=[-3200.0, -1600.0, 1600.0, 3200.0],
            impactMax=0.1 * u.mm,
            phiBinDeflectionCoverage=1,
        ),
        SeedingAlgorithmConfigArg(
            zBinNeighborsTop=[[0, 1], [-1, 1], [-1, 0]],
            zBinNeighborsBottom=[[0, 1], [-1, 1], [-1, 0]],
            numPhiNeighbors=1,
        ),
        #      seedingAlgorithm=SeedingAlgorithm.TruthEstimated,
        #    seedingAlgorithm=SeedingAlgorithm.Orthogonal,
        seedingAlgorithm=SeedingAlgorithm.Default,
        geoSelectionConfigFile=jsonSeedMS,
        outputDirRoot=outputDir,
        logLevel=acts.logging.DEBUG,
        inputSourceLinks="sourcelinks_ms",
        inputMeasurements="measurements_ms",
        suffix="ms",
        det_suffix="_ms",
        suffixSpacepoint="ms",
        doTrkVtx=False,  # switch of the vertexing with the tracklets
        inputParticles="particles_ms",
        verbose=True,
        initialVarInflation = [1,1,1,1,1000,1]
    )

    addCKFTracks(
        s,
        trackingGeometryMS,
        field,
        TrackSelectorConfig(
            pt=(100 * u.MeV, None),
            absEta=(None, None),
            nMeasurementsMin=5,
        ),
        CkfConfig(
            chi2CutOff=15,
            numMeasurementsCutOff=1,
            maxSteps=None,
            seedDeduplication=False,
            stayOnSeed=True,
        ),
        outputDirRoot=outputDir,
        inputSourceLinks="sourcelinks_ms",
        inputMeasurements="measurements_ms",
        suffixIn="ms",
        suffix="ms",
        suffixOut="ms",
        det_suffix="_ms",
    )

    addAmbiguityResolution(
        s,
        AMBIGUITYRESOLUTIONCONFIG,
        outputDirRoot=outputDir,
        suffixIn="ms",
        suffixOut="ambims",
        det_suffix="_ms",
    )

    #########
    # MATCHING
    #########
    """
    addContainerMerger(
        s,
        inputTrackParameters=["ambitrackpars","ambilltrackpars","ambilllltrackpars","ambilllllltrackpars"],
        outputTrackParameters="mergedtrackpars",   
        inputTracks=["ambitracks","ambilltracks","ambilllltracks","ambilllllltracks"],
        outputTracks="mergedtracks" 
    )
    """
    addMatching(
        s,
        trackingGeometry=trackingGeometryVT,
        magneticField=field,
        inputTrackParametersMS=["ambimstrackpars"],
        inputTrackParametersVT=["ambitrackpars","ambilltrackpars","ambilllltrackpars","ambilllllltrackpars"],
        inputTrackContainerMS=["ambimstracks"],
        inputTrackContainerVT=["ambitracks","ambilltracks","ambilllltracks","ambilllllltracks"], #"ambitracks",
        outputTrackParameters="outputTrackParameters",
        outputTracks="outputTracks",
        useRecVtx=False,
        px=0,
        py=0,
        pz=perigeeZ,
        chi2max=10000000,
    )

    #########
    # VERTEXING AND PROPAGATION TO THE PV
    #########

    addVertexFitting(
        s,
        field,
        vertexFinder=VertexFinder.Iterative,
        outputDirRoot=outputDir,
        suffixIn="ambi",
        trackParameters="ambitrackpars",
        outputProtoVertices="protoverticesstep1",
        outputVertices="fittedVerticesstep1",
        inputParticles = "particles",
        selectedParticles = "particles_selected_vt",
        inputVertices = "vertices_input",
        suffixOut="_step1",
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
        rejectedFraction=rejectedFraction,
    )

    return s
    """

    addVertexFitting(
        s,
        field,
        vertexFinder=VertexFinder.Truth,
        outputDirRoot=outputDir,
        suffixIn= "ambi",
        trackParameters="ambitrackpars",
        outputProtoVertices= "protoverticestruth",
        outputVertices= "fittedVerticestruth",
        suffixOut = "_truth",
        logLevel=acts.logging.Level.DEBUG,
        
    )
    
    """
    """
    
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
    

    """
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


if "__main__" == __name__:
    options = getArgumentParser().parse_args()

    # "events_full_sim"#"fullchain_input"
    Eint = options.eint

    event_type = "" if not options.bkg else "_noBkg"

    if options.fullev:
        event_type += "_FullEv"

    dir = "event_generation/events_" + str(Eint) + "GeV" + event_type

    target_option = (
        "_onlyTarget_" + str(options.target) if options.target != -1 else "_realTarget"
    )
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
        dir = (
            "event_generation/events_40GeV"
            + has_sec
            + target_option
            + "_beamSigma_"
            + beamsize
            + "_periferal_factor_"
            + str(options.periferal)
            + ("00000" if options.periferal >= 0.1 else "0000")
        )
    else:
        dir = (
            "event_generation/events_40GeV"
            + has_sec
            + target_option
            + "_beamSigma_"
            + beamsize
            + ""
        )

    # dir ="event_generation/events_40GeV_Sec_periferal_factor_0.100000"

    dirMS = "event_generation/events_40GeV_noBkg_jpsi_muons"
    dirVT = "event_generation/events_40GeV_noBkg_jpsi_muons"

    dirMS = "event_generation/events_40GeV_Sec_realTarget_beamSigma_0.500000_jpsi_muons"
    dirVT = "event_generation/events_40GeV_Sec_realTarget_beamSigma_0.500000_jpsi_muons"


    dirMS = "event_generation/events_40GeV_Sec_realTarget_beamSigma_0.500000_jpsi_muons"
    dirVT = "event_generation/events_40GeV_Sec_realTarget_beamSigma_0.500000_jpsi"
    

    inputDirVT = pathlib.Path.cwd() / dirVT
    inputDirMS = pathlib.Path.cwd() / dirMS
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

    suffix += "_maxSeedSpMPrim" + str(options.sf_maxSeedsPerSpM_primary)
    suffix += "_maxSeedSpMSec" + str(options.sf_maxSeedsPerSpM_secondary)

    if options.sf_seedConfirmation_primary:
        suffix += "_confirmation_primary"
    if options.sf_seedConfirmation_secondary:
        suffix += "_confirmation_secondary"

    suffix += "_ImpMax" + str(options.sf_impactMax)
    suffix += "_dZMax" + str(options.sf_deltaZMax)

    if options.sf_numMeasurementsCutOff:
        suffix += "_branch" + str(options.sf_numMeasurementsCutOff)
    if options.eff is not None:
        suffix += "_eff" + str(options.eff)
    if options.noguessing:
        suffix += "_noGuessing"
        suffix += "_z" + str(options.zperigee)
    if options.trkVtxOnly:
        suffix += "_trkVtxOnly"
    if options.addDeltas:
        suffix += "_addDeltas"
    if options.projective:
        suffix += "_projective"
    if options.periferal != 1:
        suffix += "_periferal_factor_" + str(options.periferal)
    suffix += (
        "_onlyTarget_" + str(options.target) if options.target != -1 else "_realTarget"
    )
    if options.sec:
        suffix += "_Sec"
    if options.osec:
        suffix += "_onlySec"
    if options.zz:
        suffix += "_beam0.0"
    else:
        suffix += "_beam0.5"
    if options.periferal != 1:
        dir = (
            "event_generation/events_40GeV"
            + has_sec
            + target_option
            + "_beamSigma_"
            + beamsize
            + "_periferal_factor_"
            + str(options.periferal)
            + ("00000" if options.periferal >= 0.1 else "0000")
        )
    else:
        dir = (
            "event_generation/events_40GeV"
            + has_sec
            + target_option
            + "_beamSigma_"
            + beamsize
            + ""
        )

    if options.sf_maxVertices != 1:
        suffix += "_nVtxMax" + str(options.sf_maxVertices)

    suffix += "_twosteps_rej" + str(options.sf_rejectedFraction)
    suffix += "_perigeeZ"+str(options.vz)
    suffix += "_suffix"

    current_dir = pathlib.Path.cwd()
    if options.outdir:
        outputDir = options.outdir
    else:
        outputDir = str(
            current_dir / ("output_" + str(Eint) + "GeV" + event_type + suffix)
        )

    matDeco = acts.IMaterialDecorator.fromFile("geomVTNA60p/material-map_VTNA60p.json")
    jsonFile = "geomVTNA60p/tgeo-config_VTNA60p.json"
    tgeo_fileName = "geomVTNA60p/geom_VTNA60p.root"
    jsonDigiVT = "geomVTNA60p/digismear0.005.json"
    jsonSeedVT = "geomVTNA60p/seed_config_2_4_6_8_10.json"
    logLevel = acts.logging.INFO
    customLogLevel = acts.examples.defaultLogging(logLevel=logLevel)

    detectorVT, trackingGeometryVT, decoratorsVT = TGeoDetector.create(
        jsonFile=str(jsonFile),
        fileName=str(tgeo_fileName),
        surfaceLogLevel=customLogLevel(),
        layerLogLevel=customLogLevel(),
        volumeLogLevel=customLogLevel(),
        mdecorator=matDeco,
    )

    matDeco = acts.IMaterialDecorator.fromFile(
        "geomMuons/material-map_muons_longsetup.json"
    )
    jsonFile = "geomMuons/tgeo-config_muons_longsetup.json"
    tgeo_fileName = "geomMuons/geom_muons_longsetup.root"

    jsonDigiMS = "geomMuons/digismearMS.json"
    jsonSeedMS = "geomMuons/seed_configMS.json"

    detectorMS, trackingGeometryMS, decoratorsMS = TGeoDetector.create(
        jsonFile=str(jsonFile),
        fileName=str(tgeo_fileName),
        surfaceLogLevel=customLogLevel(),
        layerLogLevel=customLogLevel(),
        volumeLogLevel=customLogLevel(),
        mdecorator=matDeco,
    )

    start_time = time.time()

    runFullChain(
        detectorVT,
        trackingGeometryVT,
        detectorMS,
        trackingGeometryMS,
        inputDirMS=inputDirMS,
        inputDirVT=inputDirVT,
        outputDir=outputDir if options.outdir == None else options.outdir,
        jsonDigiVT=jsonDigiVT,
        jsonSeedVT=jsonSeedVT,
        jsonDigiMS=jsonDigiMS,
        jsonSeedMS=jsonSeedMS,
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
        rejectedFraction=options.sf_rejectedFraction,
        perigeeZ = options.vz
    ).run()
    end_time = time.time()

    execution_time = end_time - start_time
    print("NA60+_Summary_Run time", execution_time)

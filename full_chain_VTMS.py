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


from acts.examples import (
    GaussianVertexGenerator,
    ParametricParticleGenerator,
    FixedMultiplicityGenerator,
    EventGenerator,
)

from acts import UnitConstants as u

from acts.examples.simulation import (
    addParticleReader,
    addParticleGun,
    addFatras,
    addSimHitsReader,
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
    addPropagation,
    VertexFinder,
    CkfConfig,
    TruthEstimatedSeedingAlgorithmConfigArg,
    SeedingAlgorithm,
    addMatching,
    addGx2fTracks,
    addKalmanTracks,
    addTrajectoriesToPrototracks,
    addTracksToTrajectories
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

    parser.add_argument("-ve", "--verbose", help="Verobese output", action="store_true")

    parser.add_argument("-b", "--bfield", help="", action="store_true")
    parser.add_argument("-m", "--muons", help="", action="store_true")
    parser.add_argument("-def", "--seeddef", help="", action="store_true")

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
        "-trk",
        "--trkVtxOnly",
        dest="trkVtxOnly",
        help="Remove bkg",
        action="store_true",
    )

    #########################
    parser.add_argument(
        "-ef",
        "--eff",
        dest="eff",
        help="Detector efficiency",
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
        default=10,
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
        "-fa",
        dest="fatras",
        help="Use fatras",
        action="store_true",
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
        "-dil",
        dest="dilepton",
        help="Use seed confirmation",
        type=int,
        default=0,
    )

    parser.add_argument(
        "-geo",
        dest="geometry",
        help="Use seed confirmation",
        type=int,
        default=0,
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
    trackingGeometry,
    inputDir: Path,
    outputDir: Path,
    jsonDigi="jsonDigiVT.json",
    jsonSeedVT1="jsonSeedVT.json",
    jsonSeedVT2="jsonSeedVT.json",
    jsonSeedMS="jsonSeedMS.json",
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
    useFatras = False
):

    field = acts.examples.MagneticFieldMapXyz("bfield/NewBFieldNA60plus_longsetup.txt")
    field_rotated = acts.examples.MagneticFieldMapXyz("bfield/NewBFieldNA60plus_longsetupRotated.txt")

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
        outputParticles="particles_input" if useFatras else "particles_selected",
    )

    if useFatras:
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
            enableInteractions=True,
            outputDirRoot=outputDir,
            # outputDirCsv=outputDir,
            outputParticlesInitial="particles_initial",
            outputParticlesFinal="particles_final",
            outputSimHits="simhits",
            inputSimHits= None,
            det_suffix="",
            logLevel=acts.logging.Level.MAX
        )
    else:
        addSimHitsReader(
            s,
            inputDir=inputDir,
            outputSimHits="simhits",
            outputDirRoot=outputDir,

        )

    addDigitization(
        s,
        trackingGeometry,
        field,
        digiConfigFile=jsonDigi,
        outputDirRoot=outputDir,
        efficiency=1,
        applyHole=False,
        rnd=rnd,
        suffix="",
    )


    TRUTHSEEDRANGE = TruthSeedRanges(
        pt=(100* u.MeV, None),
        eta=(None, None),
        nHits=(4, None),
        nHitsVT=(4, None),
        nHitsMS=(None, None),
        isMuon=False,
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
        stayOnSeed=False,
    )

    TRACKSELECTORCONFIG = TrackSelectorConfig(
        pt=(100 * u.MeV, None),
        absEta=(None, None),
        nMeasurementsMin=4,
    )

    AMBIGUITYRESOLUTIONCONFIG = AmbiguityResolutionConfig(
        maximumSharedHits=1, # if > 0 it will break durin the matching
        nMeasurementsMin=4,
        maximumIterations=100000000,
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

    addNewSeeding(
        s=s,
        trackingGeometry=trackingGeometry,
        field=field_rotated,
        geoSelectionConfigFile=jsonSeedVT1,
        truthSeedRanges=TRUTHSEEDRANGE,
        initialSigmas=None,
        initialVarInflation=None,
        seedFinderConfigArg=SeedFinderConfigArgNA60(
            maxSeedsPerSpM=1, #MaxSeedsPerSpMPrimary,
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
            yMiddle=(0 * u.mm, 400 * u.mm),
            #yMiddle=(139 * u.mm, 170 * u.mm),  # use only layer 2 and 3
            # yMiddle=(139 * u.mm, 220 * u.mm),  # use only layer 2 and 3
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
        noGuessing=noGuessing,
        trkVtxOnly=trkVtxOnly,
        projective=False,
        zPerigee=0,
        inputSourceLinks="sourcelinksVT",
        inputMeasurements="measurements",
        det_suffix="_vt",
        inputParticles="particles",
    )
    
    addCKFTracks(
        s,
        trackingGeometry,
        field,
        TRACKSELECTORCONFIG,
        CKFCONFIG,
        outputDirRoot=outputDir,
        inputSourceLinks="sourcelinksVT",
        inputMeasurements="measurements",
        suffixIn="",
        suffixOut="",
        det_suffix="_vt",
        logLevel=acts.logging.Level.DEBUG,
        twoWay = True,
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
        inputSourceLinks="sourcelinksVT",
        outputSourceLinks="llsourcelinks_vt",
        inputTracks="ambitracks",
    )

    addNewSeeding(
        s=s,
        trackingGeometry=trackingGeometry,
        field=field_rotated,
        geoSelectionConfigFile=jsonSeedVT2,
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
        inputMeasurements="measurements",
        det_suffix="_vt",
        doTrkVtx=False,  # switch of the vertexing with the tracklets
        inputParticles="particles",
    )

    addCKFTracks(
        s,
        trackingGeometry,
        field,
        TRACKSELECTORCONFIG,
        CKFCONFIG,
        outputDirRoot=outputDir,
        inputSourceLinks="llsourcelinks_vt",
        inputMeasurements="measurements",
        suffixIn="ll",
        suffixOut="ll",
        det_suffix="_vt",
        logLevel=acts.logging.Level.MAX,
        twoWay = True,
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
        trackingGeometry=trackingGeometry,
        field=field_rotated,
        geoSelectionConfigFile=jsonSeedVT1,
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
            deltaYBottomSP=(10 * u.mm, 1000 * u.mm),
            impactMax=10 * u.cm,
            deltaZMax=10000 * u.cm,  # was 5
            minPt=100 * u.MeV,
            interactionPointCut=False,
            verbose=verbose,
            collisionRegion=(-10 * u.cm, 10 * u.cm),  # 0.5
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
        inputMeasurements="measurements",
        det_suffix="_vt",
        doTrkVtx=False,  # switch of the vertexing with the tracklets
        inputParticles="particles",
        
    )

    addCKFTracks(
        s,
        trackingGeometry,
        field,
        CKFCONFIG,
        TRACKSELECTORCONFIG,
        outputDirRoot=outputDir,
        inputSourceLinks="llllsourcelinks_vt",
        inputMeasurements="measurements",
        suffixIn="llll",
        suffixOut="llll",
        det_suffix="_vt",
        logLevel=acts.logging.Level.MAX,
        twoWay = True,
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
        trackingGeometry=trackingGeometry,
        field=field_rotated,
        geoSelectionConfigFile=jsonSeedVT2,
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
            impactMax=10 * u.cm,
            deltaZMax=10000 * u.cm,  # was 5
            minPt=100 * u.MeV,
            interactionPointCut=False,
            verbose=verbose,
            collisionRegion=(-10 * u.cm, 10 * u.cm),  # 0.5
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
        inputMeasurements="measurements",
        det_suffix="_vt",
        doTrkVtx=False,  # switch of the vertexing with the tracklets
        inputParticles="particles",
    )

    addCKFTracks(
        s,
        trackingGeometry,
        field,
        CKFCONFIG,
        TRACKSELECTORCONFIG,
        outputDirRoot=outputDir,
        suffixIn="llllll",
        suffixOut="llllll",
        inputSourceLinks="llllllsourcelinks_vt",
        inputMeasurements="measurements",
        det_suffix="_vt",
        logLevel=acts.logging.Level.MAX,
        twoWay = True,
    )

    addAmbiguityResolution(
        s,
        AMBIGUITYRESOLUTIONCONFIG,
        outputDirRoot=outputDir,
        suffixIn="llllll",
        suffixOut="ambillllll",
        det_suffix="_vt",
    )
    """
    addVertexFitting(
        s,
        field,
        vertexFinder=VertexFinder.Iterative,
        outputDirRoot=outputDir,
        suffixIn="ambi",
        trackParameters="ambitrackpars",
        outputProtoVertices="protoverticesstep1",
        outputVertices="fittedVerticesstep1",
        inputParticles="particles",
        selectedParticles="particles_selected",
        inputVertices="vertices_input",
        suffixOut="_step1",
        det_suffix="_vt",
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
    """
    #########
    # MUON RECONSTRUCTION
    #########

    addSeeding(
        s,
        trackingGeometry,
        field_rotated,
        TruthSeedRanges(pt=(0, None), nHits=(5, None), nHitsVT=(None, None), nHitsMS=(5, None), isMuon=True, keep=(True, True)),
        SeedFinderConfigArg(
            maxSeedsPerSpM=5,
            cotThetaMax=0.25,
            deltaZMax=500,  # was 5
            sigmaScattering=50.0,
            radLengthPerSeed=0.01,
            minPt=100 * u.MeV,
            impactMax=1000 * u.mm,  # was 0.1 0.15 //very impactful cut
            interactionPointCut=False,
            useVariableMiddleSPRange=False,  # MODIFICATO 22/5/23
            deltaRMiddleSPRange=(
                0 * u.mm,
                0 * u.mm,
            ),  # not useful if useVariableMiddleSPRange=False
            deltaRTopSP=(
                350 * u.mm,
                4000 * u.mm,
            ),  # min and max R between Middle and Top SP
            deltaRBottomSP=(
                350 * u.mm,
                4000 * u.mm,
            ),  # min and max R between Middle and Bottom SP
            collisionRegion=(-1000 * u.mm, 1000 * u.mm),
            r=(20 * u.mm, 9000 * u.mm),
            seedConfirmation=False,
            forwardSeedConfirmationRange=acts.SeedConfirmationRangeConfig(  # it should not be used....
                zMinSeedConf=150 * u.mm,
                zMaxSeedConf=0 * u.mm,
                rMaxSeedConf=9000 * u.mm,
                nTopForLargeR=1,
                nTopForSmallR=1,
            ),
            centralSeedConfirmationRange=acts.SeedConfirmationRangeConfig(  # it should not be used....
                zMinSeedConf=7 * u.mm,
                zMaxSeedConf=-7 * u.mm,
                rMaxSeedConf=9000 * u.mm,
                nTopForLargeR=1,
                nTopForSmallR=1,
            ),
            verbose=False,
        ),
        SeedFinderOptionsArg(
            beamPos=(0 * u.mm, 0 * u.mm),
            bFieldInZ=1.5 * u.T,  # WHY????????37?
        ),  # why should I give the b field? to compute phi bins in SpacePointGrid.ipp
        SeedFilterConfigArg(  # not used, why?
            seedConfirmation=False,
            maxSeedsPerSpMConf=2,
            maxQualitySeedsPerSpMConf=1,
            compatSeedLimit=0,  # added 4/10/23 (value not tuned)
            verbose=False,
        ),
        SpacePointGridConfigArg(
            rMax=9000 * u.mm,
            zBinEdges=[
                -3200.0,
                -1600.0,
                1600.0,
                3200.0,
            ],  # valori sensati???? forse devo stare in =- 3790/2
            impactMax=0.1
            * u.mm,  # not used if bfieldZ is 0, otherwise it's used to compute number of phi bins
            phiBinDeflectionCoverage=1,  # not used if bfieldZ is 0, otherwise it's used to compute number of phi bins
        ),
        SeedingAlgorithmConfigArg(
            zBinNeighborsTop=[[0, 1], [-1, 1], [-1, 0]],
            zBinNeighborsBottom=[[0, 1], [-1, 1], [-1, 0]],
            numPhiNeighbors=1,
        ),
        TruthEstimatedSeedingAlgorithmConfigArg(
            deltaR=(0, 10000),
        ),
        # seedingAlgorithm=SeedingAlgorithm.Default,
        seedingAlgorithm=(
            SeedingAlgorithm.Default #if seeddef else SeedingAlgorithm.TruthSmeared
        ),
        geoSelectionConfigFile=jsonSeedMS,
        particleHypothesis=acts.ParticleHypothesis.muon,
        outputDirRoot=outputDir,
        logLevel = acts.logging.INFO,

        #logLevel=acts.logging.MAX,
        inputSourceLinks="sourcelinksMS",
        inputMeasurements="measurements",
        suffix="ms",
        det_suffix="_ms",
        suffixSpacepoint="ms",
        doTrkVtx=False,  # switch of the vertexing with the tracklets
        inputParticles="particles",
        verbose=False,
        initialVarInflation=(
            [10, 10, 10, 10, 10, 10] #[50.0, 50.0, 10.0, 10.0, 10.0, 10.0] #if seeddef else 
        ),
    )
    
    addCKFTracks(
        s,
        trackingGeometry,
        field,
        TrackSelectorConfig(
            pt=(0 * u.MeV, None),
            absEta=(None, None),
            nMeasurementsMin=5,
        ),
        CkfConfig(chi2CutOff=15, numMeasurementsCutOff=2),
        outputDirRoot=outputDir,
        inputSourceLinks="sourcelinksMS",
        inputMeasurements="measurements",
        suffixIn="ms",
        suffixOut="ms",
        det_suffix="_ms",
        logLevel=acts.logging.Level.MAX,
        twoWay = True,
    )


    addAmbiguityResolution(
        s,
        AmbiguityResolutionConfig(
            maximumSharedHits=1,  
            maximumIterations=1000000,
            nMeasurementsMin=5,
        ),
        outputDirRoot=outputDir,
        suffixIn="ms",
        suffixOut="ambims",
        det_suffix="_ms",
    )

    return s

    addVertexFitting(
        s,
        field,
        vertexFinder=VertexFinder.Truth,
        outputDirRoot=outputDir,
        suffixIn="ambims",
        trackParameters="ambitrackpars",
        outputProtoVertices="protoverticestruth",
        outputVertices="fittedVerticestruth",
        inputParticles="particles",
        inputVertices="vertices_input",
        selectedParticles="particles_selected",
        suffixOut="_truth",
        #logLevel=acts.logging.Level.VERBOSE,
        doConstrainedFit = True
    )

    addContainerMerger(
        s,
        inputTrackParameters=[
            "ambitrackpars",
            #"ambilltrackpars",
            #"ambilllltrackpars",
            #"ambilllltrackpars"#,
            #"ambilllllltrackpars"
        ],
        outputTrackParameters="mergedtrackpars",
        inputTracks=[
            "ambitracks",
            #"ambilltracks",
            #"ambilllltracks",
            #"ambilllllltracks"
        ],
        outputTracks="mergedtracks",
        outputDirRoot=outputDir,
    )
    addMatching(
        s,
        trackingGeometry=trackingGeometry,
        magneticField=field,
        inputVertices="fittedVerticesstep1",
        inputTrackParametersMS=["ambimstrackpars"],
        inputTrackParametersVT=["mergedtrackpars"],
        inputTrackContainerMS=["ambimstracks"],
        inputTrackContainerVT=["mergedtracks"], 
        outputTrackParameters="outputTrackParameters",
        inputMeasurementParticlesMapVT = "measurement_particles_mapVT",
        inputMeasurementParticlesMapMS = "measurement_particles_map",
    
        #outputTracks="outputTracks",
        useRecVtx=False,
        px=0,
        py=0,
        pz=0,
        chi2max=10000000,
        outputDirRoot=outputDir,
    ) 

    return s
    addMatching(
        s,
        trackingGeometry=trackingGeometry,
        magneticField=field,
        inputVertices="fittedVerticesstep1",
        inputTrackParametersMS=["ambimstrackpars"],
        inputTrackParametersVT=[
            "ambitrackpars",
            "ambilltrackpars",
            #"ambilllltrackpars",
            "ambilllllltrackpars",
        ],
        inputTrackContainerMS=["ambimstracks"],
        inputTrackContainerVT=[
            "ambitracks",
            "ambilltracks",
            #"ambilllltracks",
            "ambilllllltracks",
        ],  # "ambitracks",
        outputTrackParameters="outputTrackParameters",
        inputMeasurementParticlesMapVT = "measurement_particles_mapVT",
        inputMeasurementParticlesMapMS = "measurement_particles_map",
    
        #outputTracks="outputTracks",
        useRecVtx=False,
        px=0,
        py=0,
        pz=0,
        chi2max=10000000,
        outputDirRoot=outputDir,
    ) 

    return s

    #########
    # MATCHING
    #########
    addPropagation(
        s,
        trackingGeometry,
        field,
        inputTrackParameters=[
            "ambitrackpars",
            "ambilltrackpars",
            "ambilllltrackpars",
            "ambilllllltrackpars",
        ],
        # inputTracks=["ambitracks","ambilltracks","ambilllltracks","ambilllllltracks"],
        inputVertices="fittedVerticesstep1",
        outputTracks="mergetracks",
        outputDirRoot=outputDir,
        useRecVtx=False,
        suffixOut="merged",
        det_suffix="_vt",
    )
    """
    #########
    # VERTEXING AND PROPAGATION TO THE PV
    #########

    """


    addVertexFitting(
        s,
        field,
        vertexFinder=VertexFinder.Iterative,
        outputDirRoot=outputDir,
        suffixIn="ambi",
        trackParameters="mergedtrackpars",
        outputProtoVertices="protoverticesmerged",
        outputVertices="fittedVerticesmerged",
        suffixOut="_merged",
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

    """
    addContainerMerger(
        s,
        inputTrackParameters=["ambitrackpars","ambilltrackpars","ambilllltrackpars","ambilllllltrackpars"],
        outputTrackParameters="mergedlltrackpars"   
    )
    """

if "__main__" == __name__:
    options = getArgumentParser().parse_args()

    dir = "event_generation/simulatedEvents/events_40GeV_noBkg_omega2Body_muons"
    dir = "event_generation/simulatedEvents/events_40GeV_noBkg_jpsi"
    if options.dilepton == 0:
        dir = "event_generation/simulatedEvents/events_40GeV_Sec_omega2Body"
        suffix="omega_full"
    if options.dilepton == 1:
        dir = "event_generation/simulatedEvents/events_40GeV_Sec_phi"
        suffix="phi_full"
    if options.dilepton == 2:
        dir = "event_generation/simulatedEvents/events_40GeV_Sec_jpsi"
        suffix="jpsi_full"
    
    if options.geometry == 0:
        # RUBEN
        matDeco = acts.IMaterialDecorator.fromFile("geometry/geoRubenAll/material-map.json")
        jsonFile = "geometry/geoRubenAll/tgeoRubenVol.json"
        tgeo_fileName = "geometry/geoRubenAll/geometry_Ruben.root"
        suffix += "_ruben"

    if options.geometry == 1:
        # ROBERTA
        matDeco = acts.IMaterialDecorator.fromFile("geometry/geomVTMSLong/material-map_VT_MSlongsetup.json")
        jsonFile = "geometry/geomVTMSLong/tgeo-config_VT_MSlongsetup.json"
        tgeo_fileName = "geometry/geomVTMSLong/geom_VT_MSlongsetup.root"
        suffix += "_roberta_hole"

    if options.geometry == 2:
        # ROBERTA NO ABS
        matDeco = acts.IMaterialDecorator.fromFile("geometry/geomVTMSLong/material-map_VTMSlongsetup_noabs.json")
        jsonFile = "geometry/geomVTMSLong/tgeo-config_VTMSlongsetup_noabs.json"
        tgeo_fileName = "geometry/geomVTMSLong/geom_VT_MSlongsetup_noabs.root"
        suffix += "_roberta_noabs"

    suffix += "_lay12"
    dir ="event_generation/simulatedEvents/rubenxprino_40GeV_jpsi_filtered"
    suffix = "rubenxprino_40GeV_jpsi"
    inputDir = pathlib.Path.cwd() / dir

    current_dir = pathlib.Path.cwd()
    if options.outdir:
        outputDir = options.outdir
    else:
        outputDir = str(
            current_dir / ("output/output_" + suffix)
        )


    jsonDigi = "geometry/geomVTMSLong/digismearVTMS3.json"
    jsonSeedVT1 = "geometry/geomVTMSLong/seed_configVTMS.json"
    jsonSeedVT2 = "geometry/geomVTMSLong/seed_configVTMS.json"
    jsonSeedMS = "geometry/geomVTMSLong/seed_configMS.json"

    logLevel = acts.logging.INFO
    customLogLevel = acts.examples.defaultLogging(logLevel=logLevel) 
    detector, trackingGeometry, decorators = TGeoDetector.create(
        jsonFile=str(jsonFile),
        fileName=str(tgeo_fileName),
        surfaceLogLevel=customLogLevel(),
        layerLogLevel=customLogLevel(),
        volumeLogLevel=customLogLevel(),
        mdecorator=matDeco,
    )
    runFullChain(
        trackingGeometry,
        inputDir=inputDir,
        outputDir=outputDir if options.outdir == None else options.outdir,
        jsonDigi=jsonDigi,
        jsonSeedVT1=jsonSeedVT1,
        jsonSeedVT2=jsonSeedVT2,
        jsonSeedMS=jsonSeedMS,
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
        useFatras = False, #options.fatras
    ).run()
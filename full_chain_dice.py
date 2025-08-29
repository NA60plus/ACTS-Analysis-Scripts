#!/usr/bin/env python3

import acts
import acts.examples
import logging
# Also enable Python logging
logging.basicConfig(level=logging.DEBUG)
import sys
sys.path.append('/home/giacomo/acts_for_NA60+/acts/Examples/Scripts/Python')
import argparse
import pathlib

from dice_geometry import buildDICEgeometry


from acts import UnitConstants as u

from acts.examples.simulation import (
    addFatras,
    MomentumConfig,
    EtaConfig,
    PhiConfig,
    ParticleConfig,
    ParticleSelectorConfig,
    addParticleReader,
    addDigitization,
    addDigiParticleSelection
)


from acts.examples.reconstruction import (
    addSeeding,
    SeedingAlgorithm,
    addCKFTracks,
    CkfConfig,
    addAmbiguityResolution,
    AmbiguityResolutionConfig,
    addVertexFitting,
    VertexFinder,
    addTrackletVertexing,
    SeedFinderOptionsArg,
    SeedingAlgorithmConfigArg,
    SeedFinderConfigArg,
    SeedFilterConfigArg,
    SpacePointGridConfigArg,
    TrackSelectorConfig,
    AmbiguityResolutionConfig,
    addUsedMeasurementsFilter,
)

def getArgumentParser():
    """Get arguments from command line"""
    parser = argparse.ArgumentParser(description="Command line arguments for reconstruction chain")

    parser.add_argument(
        "-n",
        "--nEvents",
        dest="nEvts",
        help="Number of events to run over",
        default=1,
        type=int,
    )


    return parser


def runFullChain(
    NumEvents=10
):
    

    TRACKSELECTORCONFIG = TrackSelectorConfig(
        pt=(0 * u.MeV, None),
        absEta=(None, None),
        nMeasurementsMin=4,
    )

    AMBIGUITYRESOLUTIONCONFIG = AmbiguityResolutionConfig(
        maximumSharedHits=1,
        nMeasurementsMin=4,
        maximumIterations=1000000,
    )


    s = acts.examples.Sequencer(
        events=NumEvents,
        numThreads=1,
        outputDir=str(outputDir),
        trackFpes=False,
    )

    dirVT = "event_generation/simulatedEvents/events_40GeV"

    #field = acts.examples.MagneticFieldMapXyz("bfield/field_map.txt")
    #field_rotated = acts.examples.MagneticFieldMapXyz("bfield/rotated_field_map.txt")


    field = acts.ConstantBField(acts.Vector3(
        0.0, 1.5 * u.T,  0.0))  # ~dipole field
    field_rotated = acts.ConstantBField(acts.Vector3(0.0, 0.0, -1.5 * u.T))    
    
    inputDirVT = pathlib.Path.cwd() / dirVT
    addParticleReader(
        s,
        inputDir=inputDirVT,
        outputDirRoot=outputDir,
        det_suffix ="",
        outputParticles = "particles_generated_selected"

    )

    addFatras(
        s,
        trackingGeometry,
        field,
        rnd,
        outputDirRoot=outputDir
    )

    addDigitization(
        s,
        trackingGeometry,
        field,
        digiConfigFile="geometry/digismear.json",
        outputDirRoot=outputDir,
        rnd=rnd,
    )

    addDigiParticleSelection(
        s,
        ParticleSelectorConfig(
            pt=(0.1 * u.GeV, None),
            eta=(0, None),
            measurements=(5, None),
            removeNeutral=True,
        ),
    )

    addSeeding(
        s,
        trackingGeometry,
        field_rotated,
        seedingAlgorithm=SeedingAlgorithm.Default,
        geoSelectionConfigFile="geometry/seed_config_VS.json",
        outputDirRoot=outputDir,
        logLevel=acts.logging.ERROR,
        seedFinderOptionsArg=SeedFinderOptionsArg(  # NO NEED TO CHANGE THIS
            beamPos=(0 * u.mm, 0 * u.mm, 0 * u.mm),
            bFieldInZ=1.5 * u.T,
        ),
        seedFinderConfigArg=SeedFinderConfigArg(
            maxSeedsPerSpM=1,
            cotThetaMax=0.25,  # 0.25,
            deltaZMax=10,
            sigmaScattering=5.0,
            radLengthPerSeed=10,
            minPt=100 * u.MeV,
            impactMax=2 * u.mm,
            interactionPointCut=True,
            deltaRTopSP=(
                50 * u.mm,
                150 * u.mm,
            ),  # min and max R between Middle and Top SP
            deltaRBottomSP=(
                50 * u.mm,
                150 * u.mm,
            ),  # min and max R between Middle and Bottom SP
            collisionRegion=(-2 * u.mm, 2 * u.mm),
            r=(0 * u.mm, 500 * u.mm),
            rMiddle=(139 * u.mm, 220 * u.mm)
        ),
        seedFilterConfigArg=SeedFilterConfigArg(
            seedConfirmation=False,
            maxSeedsPerSpMConf=2,
            maxQualitySeedsPerSpMConf=1,
            compatSeedLimit=0,
        ),
        spacePointGridConfigArg=SpacePointGridConfigArg(
            rMax=500 * u.mm,
            zBinEdges=[-150.,-75.,75.,150.],
            impactMax=1 * u.mm,
            phiBinDeflectionCoverage=1,
        ),
        seedingAlgorithmConfigArg=SeedingAlgorithmConfigArg(
            zBinNeighborsTop=[[0,1],[-1,1],[-1,0]],
            zBinNeighborsBottom=[[0,1],[-1,1],[-1,0]],
            numPhiNeighbors=1,
            allowSeparateRMax=True,
        ),

        particleHypothesis=acts.ParticleHypothesis.pion,
    )

    addTrackletVertexing(
        s,
        outputDirRoot=outputDir,
    )

    addCKFTracks(
        s,
        trackingGeometry,
        field,
        ckfConfig=CkfConfig(
            seedDeduplication=True,
            stayOnSeed=True,
            constrainToVolumes=[3,5]
        ),
        trackSelectorConfig = TrackSelectorConfig(
            pt=(0 * u.MeV, None),
            absEta=(None, None),
            nMeasurementsMin=4,
        ),
        twoWay = True,
        outputDirRoot=outputDir,
        writeTrackSummary=False
    )


    addAmbiguityResolution(
        s,
        AmbiguityResolutionConfig(
            maximumSharedHits=1,
            maximumIterations=1000000,
            nMeasurementsMin=5,
        ),
        outputDirRoot=outputDir,
    )

    addUsedMeasurementsFilter(
        s,
        logLevel = acts.logging.VERBOSE,
        outputDirRoot=outputDir,
        trackingGeometry=trackingGeometry
    )
    return s

if "__main__" == __name__:

    options = getArgumentParser().parse_args()

    logLevel = acts.logging.INFO
    customLogLevel = acts.examples.defaultLogging(logLevel=logLevel)

    matDeco = (
        acts.IMaterialDecorator.fromFile("geometry/fullgeo/material-map.json")
    )
    detector = buildDICEgeometry(matDeco=matDeco)
    trackingGeometry = detector.trackingGeometry()
    decorators = detector.contextDecorators()

    outputDir = "test_geometry_output"


    rnd = acts.examples.RandomNumbers(seed=44)

    runFullChain(
        NumEvents=options.nEvts
    ).run()

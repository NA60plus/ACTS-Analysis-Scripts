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
    addSimHitsReader,
    addGenParticleSelection,
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
    addMatching,
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
"""
What should be added:
- FATRAS with genereated tracks and hits from Geant4


"""

def getArgumentParser():
    """Get arguments from command line"""
    parser = argparse.ArgumentParser(description="Command line arguments for reconstruction chain")

    parser.add_argument(
        "-n",
        "--nev",
        dest="nev",
        help="Number of events to run over",
        default=1,
        type=int,
    )

    parser.add_argument("-f","--fatras", action="store_true", help="Run FATRAS simulation")
    parser.add_argument("-rvs","--remove-vs", action="store_true", help="Remove Vertex Spectrometer from the geometry")
    parser.add_argument("-rms","--remove-ms", action="store_true", help="Remove Muon Spectrometer from the geometry")
    parser.add_argument("--outputDir", type=str, default="output", help="Output directory")
    parser.add_argument("--inputDir", type=str, default="event_generation/simulatedEvents/events_40GeV_noBkg_jpsi", help="Input directory")

    return parser


def runFullChain(
    nEvents: int = 10,
    useFATRAS : bool = False,
    inputDir=pathlib.Path("."),
    outputDir=pathlib.Path("output"),
    seedConfig="geometry/seed_config_VS_layers234.json",
    constrainToVolumesVS = None,
    endOfWorldVolumesVS = None
):
    

    TRACKSELECTORCONFIG = TrackSelectorConfig(
        pt=(0 * u.MeV, None),
        absEta=(None, None),
        nMeasurementsMin=4,
    )

    TRACKSELECTORCONFIGMS = TrackSelectorConfig(
        pt=(0 * u.MeV, None),
        absEta=(None, None),
        nMeasurementsMin=6,
    )

    AMBIGUITYRESOLUTIONCONFIG = AmbiguityResolutionConfig(
        maximumSharedHits=1,
        nMeasurementsMin=4,
        maximumIterations=1000000,
    )

    AMBIGUITYRESOLUTIONCONFIGMS = AmbiguityResolutionConfig(
        maximumSharedHits=1,
        nMeasurementsMin=6,
        maximumIterations=1000000,
    )

    SEEDFINDEROPTIONSARG= SeedFinderOptionsArg( 
        beamPos=(0 * u.mm, 0 * u.mm, 0 * u.mm),
        bFieldInZ=-1.47 * u.T,
    )

    SPACEPOINTGRIDCONFIGARG=SpacePointGridConfigArg(
        rMax=500 * u.mm,
        zBinEdges=[-150.,-75.,75.,150.],
        impactMax=1 * u.mm,
        phiBinDeflectionCoverage=1,
    )

    SEEDINGALGORITHMCONFIGARG=SeedingAlgorithmConfigArg(
        zBinNeighborsTop=[[0,1],[-1,1],[-1,0]],
        zBinNeighborsBottom=[[0,1],[-1,1],[-1,0]],
        numPhiNeighbors=1,
        allowSeparateRMax=True,
    )

    s = acts.examples.Sequencer(
        events=nEvents,
        numThreads=1,
        outputDir=str(outputDir),
        trackFpes=False,
    )

    field = acts.examples.MagneticFieldMapXyz("bfield/field_map.txt")
    
    addParticleReader(
        s,
        inputDir=inputDir,
        outputDirRoot=outputDir,
        det_suffix ="",
        outputParticles = "particles_generated_selected"

    )

    if useFATRAS:
        addFatras(
            s,
            trackingGeometry,
            field,
            rnd,
            outputDirRoot=outputDir
        )

    else:
        addSimHitsReader(
            s,
            inputDir = inputDir,
            outputSimHits = "simhits",
            outputDirRoot = outputDir
        )
        addGenParticleSelection(
            s,
            ParticleSelectorConfig(),
        )
        

    addDigitization(
        s,
        trackingGeometry,
        field,
        digiConfigFile="geometry/digismear.json",
        outputDirRoot=outputDir,
        rnd=rnd,
    )

    if useFATRAS:
        addDigiParticleSelection(
            s,
            ParticleSelectorConfig(
                pt=(None, None),
                eta=(None, None),
                measurements=(None, None),
                removeNeutral=True,
            ),
        )

    addSeeding(
        s,
        trackingGeometry,
        field,
        seedingAlgorithm=SeedingAlgorithm.Default,
        geoSelectionConfigFile=seedConfig,
        outputDirRoot=outputDir,
        logLevel=acts.logging.DEBUG,
        seedFinderOptionsArg=SEEDFINDEROPTIONSARG,
        seedFinderConfigArg=SeedFinderConfigArg(
            maxSeedsPerSpM=1,
            cotThetaMax=0.25,  # 0.25,
            deltaZMax=10,
            sigmaScattering=5.0,
            radLengthPerSeed=10,
            minPt=100 * u.MeV,
            impactMax=1 * u.mm,
            interactionPointCut=True,
            deltaRTopSP=(
                50 * u.mm,
                150 * u.mm,
            ),  # min and max R between Middle and Top SP
            deltaRBottomSP=(
                50 * u.mm,
                150 * u.mm,
            ),  # min and max R between Middle and Bottom SP
            collisionRegion=(-1 * u.mm, 1 * u.mm),
            r=(0 * u.mm, 250 * u.mm),
            rMiddle=(139 * u.mm, 210 * u.mm)
        ),
        seedFilterConfigArg=SeedFilterConfigArg(
            seedConfirmation=False,
            maxSeedsPerSpMConf=2,
            maxQualitySeedsPerSpMConf=1,
            compatSeedLimit=0,
        ),
        spacePointGridConfigArg=SPACEPOINTGRIDCONFIGARG,
        seedingAlgorithmConfigArg=SEEDINGALGORITHMCONFIGARG,
        particleHypothesis=acts.ParticleHypothesis.pion,
        initialVarInflation=[100, 100, 100, 100, 100, 100],
    )

        

    """
    addTrackletVertexing(
        s,
        outputDirRoot=outputDir,
    )
    """

    addCKFTracks(
        s,
        trackingGeometry,
        field,
        ckfConfig=CkfConfig(
            chi2CutOffMeasurement=15,
            chi2CutOffOutlier=25,
            numMeasurementsCutOff=1,
            maxSteps=None,
            pixelVolumes=None,
            stripVolumes=None,
            maxPixelHoles=None,
            maxStripHoles=None,
            trimTracks=False,
            seedDeduplication=True,
            stayOnSeed=True,
            constrainToVolumes=constrainToVolumesVS,
            endOfWorldVolumes=endOfWorldVolumesVS
        ),
        trackSelectorConfig = TRACKSELECTORCONFIG,
        twoWay = False,
        outputDirRoot=outputDir,
        writeTrackSummary=False
    )


    addAmbiguityResolution(
        s,
        AMBIGUITYRESOLUTIONCONFIG,
        outputDirRoot=outputDir,
    )

    addUsedMeasurementsFilter(
        s,
        logLevel = acts.logging.VERBOSE,
        outputDirRoot=outputDir,
        trackingGeometry=trackingGeometry,
        inputMeasurements="measurements",
        outputMeasurements="measurements_ms",
        inputSimHits = "simhits",
        outputMeasurementParticlesMap = "measurement_particles_map_ms",
        outputMeasurementSimHitsMap = "measurement_simhits_map_ms",
        outputParticleMeasurementsMap = "particle_measurements_map_ms",
        outputSimHitMeasurementsMap = "simhit_measurements_map_ms",
        inputSimHitMeasurementsMap = "simhit_measurements_map",
        fileSuffix="_ms"
    )
    

    #########
    # Muon spectrometer
    #########

    addSeeding(
        s,
        trackingGeometry,
        field,
        seedingAlgorithm=SeedingAlgorithm.Default,
        geoSelectionConfigFile="geometry/seed_config_MS_layers123.json",
        outputDirRoot=outputDir,
        logLevel=acts.logging.DEBUG,
        seedFinderConfigArg=SeedFinderConfigArg(
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
        seedFinderOptionsArg = SeedFinderOptionsArg(
            beamPos=(0 * u.mm, 0 * u.mm, 0 * u.mm),
            bFieldInZ=0.5 * u.T,
        ),
        seedFilterConfigArg=SeedFilterConfigArg(  # not used, why?
            seedConfirmation=False,
            maxSeedsPerSpMConf=2,
            maxQualitySeedsPerSpMConf=1,
            compatSeedLimit=0,  # added 4/10/23 (value not tuned)
        ),
        spacePointGridConfigArg=SpacePointGridConfigArg(
            rMax=9000 * u.mm,
            zBinEdges=[-3200.0, -1600.0, 1600.0, 3200.0],
            impactMax=0.1 * u.mm,
            phiBinDeflectionCoverage=1,
        ),
        seedingAlgorithmConfigArg=SeedingAlgorithmConfigArg(
            zBinNeighborsTop=[[0, 1], [-1, 1], [-1, 0]],
            zBinNeighborsBottom=[[0, 1], [-1, 1], [-1, 0]],
            numPhiNeighbors=1,
        ),
        particleHypothesis=acts.ParticleHypothesis.muon,
        inputMeasurements="measurements_ms",
        outputSpacePoints="spacepoints_ms",
        outputSeeds="seeds_ms",
        outputTrackParameters="estimatedparameters_ms",
        outputTrackParticleMatching="seed_particle_matching_ms",
        outputParticleTrackMatching="particle_seed_matching_ms",
        fileSuffix="_ms",
        initialVarInflation=[100, 100, 100, 100, 100, 100],
    )

    addCKFTracks(
        s,
        trackingGeometry,
        field,
        ckfConfig=CkfConfig(
            seedDeduplication=True,
            stayOnSeed=True,
            constrainToVolumes=[5],
        ),
        trackSelectorConfig = TRACKSELECTORCONFIG,
        twoWay = True,
        outputDirRoot=outputDir,
        writeTrackSummary=False,
        fileSuffix="_ms",
    )

    addAmbiguityResolution(
        s,
        AMBIGUITYRESOLUTIONCONFIGMS,
        outputDirRoot=outputDir,
        tracks="tracks_ms",
        fileSuffix="_ms",
    )

    addMatching(
        s,
        trackingGeometry,
        field,
        inputTracksMS = "ambi_tracks_ms",
        inputTracksVT = "ambi_tracks",
        outputTracksVT = "matchVT",
        outputTracksMS = "matchMS",
        outputMatchedTracks = "outputMatchedTracks",
        inputParticles = "particles",
        inputMeasurementParticlesMapVT = "measurement_particles_map",
        inputMeasurementParticlesMapMS = "measurement_particles_map",
        chi2max = 10000000,
        suffixOut = "matched",
        outputDirRoot=outputDir,
        writeTrackSummary = True,
        writeTrackStates = True,
        writePerformance = True,
        writeCovMat = False,
        logLevel = acts.logging.Level.DEBUG,
        geoIdForPropagation = 360288038909116416
    )
    return s
    #########
    # STEP 2
    #########

    addSeeding(
        s,
        trackingGeometry,
        field,
        seedingAlgorithm=SeedingAlgorithm.Default,
        geoSelectionConfigFile="geometry/seed_config_VS_layers234.json",
        outputDirRoot=outputDir,
        logLevel=acts.logging.DEBUG,
        seedFinderOptionsArg=SEEDFINDEROPTIONSARG,
        seedFinderConfigArg=SeedFinderConfigArg(
            maxSeedsPerSpM=1,
            cotThetaMax=0.25,  # 0.25,
            deltaZMax=10,
            sigmaScattering=5.0,
            radLengthPerSeed=10,
            minPt=100 * u.MeV,
            impactMax=1 * u.mm,
            interactionPointCut=True,
            deltaRTopSP=(
                50 * u.mm,
                150 * u.mm,
            ),  # min and max R between Middle and Top SP
            deltaRBottomSP=(
                50 * u.mm,
                150 * u.mm,
            ),  # min and max R between Middle and Bottom SP
            collisionRegion=(-1 * u.mm, 1 * u.mm),
            r=(0 * u.mm, 350 * u.mm),
            rMiddle=(0 * u.mm, 500 * u.mm)
        ),
        seedFilterConfigArg=SeedFilterConfigArg(
            seedConfirmation=False,
            maxSeedsPerSpMConf=2,
            maxQualitySeedsPerSpMConf=1,
            compatSeedLimit=0,
        ),
        spacePointGridConfigArg=SPACEPOINTGRIDCONFIGARG,
        seedingAlgorithmConfigArg=SEEDINGALGORITHMCONFIGARG,
        particleHypothesis=acts.ParticleHypothesis.pion,
        inputMeasurements="measurements_step2",
        outputSpacePoints="spacepoints_step2",
        outputSeeds="seeds_step2",
        outputTrackParameters="estimatedparameters_step2",
        outputTrackParticleMatching="seed_particle_matching_step2",
        outputParticleTrackMatching="particle_seed_matching_step2",
        fileSuffix="_step2"
    )

    return s



if "__main__" == __name__:
    
    options = getArgumentParser().parse_args()

    logLevel = acts.logging.INFO
    customLogLevel = acts.examples.defaultLogging(logLevel=logLevel)

    matDeco = (
        acts.IMaterialDecorator.fromFile("geometry/fullgeo/material-map.json")
    )

    geometryFile = '/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/geometry/fullgeo/geometry.root'
    seedConfig="geometry/seed_config_VS_layers123.json"

    constrainToVolumesVS = [3]
    endOfWorldVolumesVS = [5]
    if options.remove_vs:
        geometryFile = '/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/geometry/fullgeo_novs/geometry.root'
    elif options.remove_ms:
        geometryFile = '/home/giacomo/acts_for_NA60+/ACTS-Analysis-Scripts/geometry/fullgeo_noms/geometry.root'
        seedConfig="geometry/seed_config_VS_noms_layers123.json"

        matDeco = (
            acts.IMaterialDecorator.fromFile("geometry/fullgeo_noms/material-map.json")
        )

        constrainToVolumesVS = [1]
        endOfWorldVolumesVS = None
    detector = buildDICEgeometry(geometryFile = geometryFile,
                                 matDeco=matDeco,
                                 addVS=not options.remove_vs, addMS=not options.remove_ms)
    trackingGeometry = detector.trackingGeometry()
    decorators = detector.contextDecorators()
 
    inputDir = pathlib.Path.cwd() / options.inputDir
    outputDir = pathlib.Path.cwd() / options.outputDir

    rnd = acts.examples.RandomNumbers(seed=44)

    runFullChain(
        nEvents=options.nev,
        useFATRAS=options.fatras,
        inputDir = inputDir,
        outputDir = outputDir,
        seedConfig=seedConfig,
        constrainToVolumesVS = constrainToVolumesVS,
        endOfWorldVolumesVS = endOfWorldVolumesVS
    ).run()

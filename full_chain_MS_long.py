l#!/usr/bin/env python3
import argparse, pathlib, acts, acts.examples
from acts.examples.simulation import (
    addParticleGun,
    addParticleReader,
    MomentumConfig,
    EtaConfig,
    PhiConfig,
    ParticleConfig,
    addPythia8,
    addFatras,
    addGeant4,
    ParticleSelectorConfig,
    addDigitization,
)
from acts.examples.reconstruction import (
    addSeeding,
    addNewSeeding,
    TruthSeedRanges,
    SeedFinderConfigArg,
    SeedFinderOptionsArg,
    SeedFilterConfigArg,
    SpacePointGridConfigArg,
    SeedingAlgorithmConfigArg,
    SeedFinderConfigArgNA60,
    SeedFinderOptionsArgNA60,
    SeedFilterConfigArgNA60,
    TruthEstimatedSeedingAlgorithmConfigArg,
    SeedingAlgorithm,
    ParticleSmearingSigmas,
    addCKFTracks,
    TrackSelectorConfig,
    CkfConfig,
    addAmbiguityResolution,
    AmbiguityResolutionConfig,
    addVertexFitting,
    VertexFinder,
)
from acts.examples import TGeoDetector

# from acts import UnitConstants as u

parser = argparse.ArgumentParser(
    description="Full chain with the NA60+ MS setup - rotated SP"
)

parser.add_argument("--events", "-n", help="Number of events", type=int, default=100)

args = vars(parser.parse_args())

u = acts.UnitConstants
inputDir = (
    pathlib.Path.cwd() / "event_generation/events_40GeV_noBkg_jpsi_muons"
)

outputDir = pathlib.Path.cwd() / "output/fullchainMS_output"

matDeco = acts.IMaterialDecorator.fromFile(
    "geometry/geomMuonsLongSetup/material-map_muons_longsetup.json"
)
jsonFile = "geometry/geomMuonsLongSetup/tgeo-config_muons_longsetup.json"
tgeo_fileName = "geometry/geomMuonsLongSetup/geom_muons_longsetup.root"

logLevel = acts.logging.VERBOSE
customLogLevel = acts.examples.defaultLogging(logLevel=logLevel)

detector, trackingGeometry, decorators = TGeoDetector.create(
    jsonFile=str(jsonFile),
    fileName=str(tgeo_fileName),
    surfaceLogLevel=customLogLevel(),
    layerLogLevel=customLogLevel(),
    volumeLogLevel=customLogLevel(),
    mdecorator=matDeco,
)

field = acts.examples.MagneticFieldMapXyz("bfield/BFieldNA60plus_longsetup.txt")
field2 = acts.examples.MagneticFieldMapXyz(
    "bfield/BFieldNA60plus_ZRotated_longsetup.txt"
)
# parameters defined in Examples/Python/src/MagneticField.cpp

rnd = acts.examples.RandomNumbers(seed=123)

s = acts.examples.Sequencer(
    events=args["events"],
    numThreads=1,
    outputDir=str(outputDir),
    trackFpes=False,  # to remove mechanism to catch FPE during algorithm execution (see Mattermost 18/7/23 Igor)
)


addParticleReader(
    s,
    inputDir=inputDir,
    #   outputDirCsv=outputDir,
    outputDirRoot=outputDir,
    printParticles=False,
)

addFatras(
    s,
    trackingGeometry,
    field,
    rnd,
    preSelectParticles=ParticleSelectorConfig(
        pt=(0 * u.MeV, None),
        removeNeutral=True,
    ),
    outputDirRoot=outputDir,
    logLevel=acts.logging.INFO,
)


addDigitization(
    s,
    trackingGeometry,
    field,
    digiConfigFile="geometry/geomMuonsLongSetup/digismearMS.json",
    applyHole=False,
    outputDirRoot=outputDir,
    rnd=rnd,
)



TRUTHSEEDRANGE = TruthSeedRanges(
    pt=(100 * u.MeV, None),
    eta=(None, None),
    nHits=(4, None),
    keep=(True, True),
)

SPACEPOINTGRIDCONFIGARG = SpacePointGridConfigArg(  # NO NEED TO CHANGE THIS
    rMax=9000 * u.mm, zBinEdges=[-4000.0, -1500.0, 1500.0, 4000.0]
)

CKFCONFIG = CkfConfig(
    chi2CutOff=15,
    numMeasurementsCutOff=int(10),
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
    bFieldInZ=0 * u.T,
)

addSeeding(
    s,
    trackingGeometry,
    field2,
    #    TruthSeedRanges(pt=(None,None), eta=(1.35,4.5), nHits=(None, None)),
    TruthSeedRanges(pt=(None, None), eta=(None,None), nHits=(6, None)),
    SeedFinderConfigArg(
        maxSeedsPerSpM=5,
        #       cotThetaMax= 1.,
        cotThetaMax=0.25,
        deltaZMax=500,  # was 5
        sigmaScattering=5.0,
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
        collisionRegion=(-500 * u.mm, 500 * u.mm),
        r=(20 * u.mm, 9000 * u.mm),
        rMiddle=(3400 * u.mm, 4000 * u.mm),
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
    ),
    SeedFinderOptionsArg(
        beamPos=(0 * u.mm, 0 * u.mm),
        bFieldInZ=0 * u.T,  # WHY????????37?
    ),  # why should I give the b field? to compute phi bins in SpacePointGrid.ipp
    SeedFilterConfigArg(  # not used, why?
        seedConfirmation=False,
        maxSeedsPerSpMConf=2,
        maxQualitySeedsPerSpMConf=1,
        compatSeedLimit=0,  # added 4/10/23 (value not tuned)
    ),
    SpacePointGridConfigArg(
        rMax=9000 * u.mm,
        zBinEdges=[
            -3200.0,
            -1600.0,
            1600.0,
            3200.0,
        ],  # valori sensati???? forse devo stare in =- 3790/2
        impactMax=0.1 * u.mm,  # not used if bfieldZ is 0, otherwise it's used to compute number of phi bins
        phiBinDeflectionCoverage=1,  # not used if bfieldZ is 0, otherwise it's used to compute number of phi bins
    ),
    SeedingAlgorithmConfigArg(
        zBinNeighborsTop=[[0, 1], [-1, 1], [-1, 0]],
        zBinNeighborsBottom=[[0, 1], [-1, 1], [-1, 0]],
        numPhiNeighbors=100,
    ),
    TruthEstimatedSeedingAlgorithmConfigArg(
        deltaR=(0, 10000),
    ),
    #      seedingAlgorithm=SeedingAlgorithm.TruthEstimated,
    #    seedingAlgorithm=SeedingAlgorithm.Orthogonal,
    seedingAlgorithm=SeedingAlgorithm.Default,
    geoSelectionConfigFile="geometry/geomMuonsLongSetup/seed_configMS.json",
    outputDirRoot=outputDir,
    logLevel=acts.logging.DEBUG,
    noGuessing=True,
)

addCKFTracks(
    s,
    trackingGeometry,
    field,
    TrackSelectorConfig(
        pt=(0.0, None),
        absEta=(None, None),
        nMeasurementsMin=5,
    ),
    CkfConfig(
        chi2CutOff=15,
        numMeasurementsCutOff=10,  # ???
    ),
    outputDirRoot=outputDir,
    #    logLevel=acts.logging.DEBUG,
)

addAmbiguityResolution(
    s,
    AmbiguityResolutionConfig(
        maximumSharedHits=1,  # no points in common (piu' e' basso, piu' e' stringente)
        maximumIterations=1000000,
        nMeasurementsMin=5,
        #     nMeasurementsMin=6,
    ),
    outputDirRoot=outputDir,
    logLevel=acts.logging.DEBUG,
)

# addVertexFitting(
#    s,
#    field,
#    vertexFinder=VertexFinder.Iterative,
#    outputDirRoot=outputDir,
# )

s.run()

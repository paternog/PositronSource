## PositronSource
PositronSource is a [Geant4](http://www.geant4.org/geant4/) application to simulate a positron source based both on the conventional approach and on oriented crystals. In the latter case both the single-crystal case and the hybrid scheme can be 

Both Geant4 and [CMake](https://cmake.org/) need to be installed on your machine in order to compile and run the software.
Tested with:
[![Geant4](https://img.shields.io/badge/Geant4-11.02.p01-blue.svg)](https://geant4.web.cern.ch/) [![CMake](https://img.shields.io/badge/CMake-3.22.1-blue.svg)](https://cmake.org/)

The source code can be downloaded either as a ZIP archive, from the Code drop-down menu (look at the green box above), or directly from the terminal (open in your project working directory) via:
```
git clone git@github.com:paternog/PositronSource.git
```
Note: the latter requires [Git](https://git-scm.com/) installed on your machine.

### Settings
The application is thought to be used without the need of a deep understaing of how a Geant4 simulation works. Indeed, through a set of custom macro commands, the user can define the geometry and the physics he wants to use (see the macro run.mac inside the macros folder). 

### Geomtery
The application can be used to simulate a calorimeter composed of a set of oriented crystals. The calorimeter can be segmented both transversely and longitudinally. The number of layers and the number of crystals per layer can be set through specific macro commands. To simulate the OREO case, the first layer only can be composed crystals oriented along a specific axis, whereas the other layer are composed of randmoly oriented crystals. Furtherore, a pre-shower layer at a given distance from the main calorimeter can be positioned. All the crystal features can be set via macro commands. 

<p align="center">
    <img src="./readme_pics/geom.jpg" alt="" width="990" height="465">
</p>

<p align="center">
    <img src="./readme_pics/shower.png" alt="" width="990" height="465">
</p>

### Physics
There are many settings that can be used to tune the simulation of particle interactions inside oriented crystals.
A series of commands can be used to tag the particles and simulate their trajectories. In this case, a text file will be produced at the end of the run with the following columns:
```
"volume", "eventID", "trackID", "x [mm]", "tx [rad], "z [mm]", "xx [mm]", "yy [mm]""
```
In order to plot the trajectory use z and xx variables, while for tagging use x, tx and z.

Some important parameters of the Bair-Katkov (BK) model for radiation simulation can also be set, for instance the minimum photon energy the crystal can emit (`/crystal/setLEth 0.3 MeV`) or the number of small steps of the particle trajectory before radiation probability is evualeted with BK (`/crystal/setNSmallTrajectorySteps 1000`).

The FastSim channeling model can be deactivated through the macro command: `/crystal/setOCeffects false`. This is very useful to simulate the random/amorphous case without changing the angles with respect to selected orientation (with `/crystal/setCrystalAngleX angX_rad`, `/crystal/setCrystalAngleY angY_rad` commands).

For the other interactions, standard Geant4 physics model are used. In particular `FTFP_BERT` has been selected as base physics list and the electromagnetic models has been chosen among the `Livermore` ones.

### Event generation
The beam features (primary particles) have to be set through standard `G4GeneralParticleSource` (gps) commands as in the run.mac macro.

### Output
PositronSource is optimised for an ntuple-based output, in which data from sensitive detectors are written event by event. The default output file format is the [ROOT](https://root.cern/) file (`.root`), which contains the ntuples as [tree objects](https://root.cern.ch/doc/master/classTTree.html). The output file is saved in `output/` (in the build path) at the end of the program execution; its name can be set in through the command `/run/setfilenamesave output/NAME_YOU_WANT.root`. Alternatively, different file formats can be chosen, e.g. the CSV. 

The output file contains at least two ntuples. The first and the second ntuple are called `EdepPreShowerCrystals` and `EdepCalorimeterCrystals`, respectively, and contain the following variables (columns):
```
"eventID", "crystalID", "edep" 
```
Which represents:
- the event number within the run (column 0),
- the crystal ID (column 1),
- the energy deposited (GeV) in the crystal (column 2).

The two ntuples, as suggested by their name, refer to the energy deposited in the preshower and the calorimeter, respectively.

A third ntuple (`SecondariesInCrystals`) to score the features of the secondary particles produced in the calorimeter can be also recorded, provided that the related lines in RunAction.cc and SteppingAction.cc are decommented. This is not the default, since it can imply an unnecessary amount of storage and computation time. The columns of this ntuple are:
```
"eventID", "type", "energy", "theta"
```
where type is the name (gamma, electron, positron, other), energy is the kinetic energy (GeV) and theta is the transverse angle (rad) of the produced particle.

The information scored in the first two ntuples can also be accumulated by the custom `Run` class and saved in a txt file written at the end of the simulation (in `RunAction::EndOfRunAction(const G4Run* run)`). The aforementioned file will contain the following columns: "volumeID crystalID Edep(GeV) stdEdep(GeV)", where volumeID=0 for PreShower, whereas volumeID=1 for calorimeter. This modality can be turned on and off through the macro command: `/det/setBIscoring true or false`. Similarly, the three-dimensional distributions of energy deposition int PreShower and Calorimeter can be scored through the custom `VoxelScorer` class, which is similar to the native box mesh scorer, but, unlike this, it is intrinsically linked to volumes, thus it is less prone to positioning errors. This modality can be turned on and off through the macro command: `/det/setVoxelization true or false`, also specific commands allow the user to define the number of voxel in the three directions (see the general macro: macos/run.mac).

## Analysis scripts
A root script as well as a series of python notebooks and custom libraries useful for the simulation result analysis is provided. The code is commented and require only the setting of proper input parameters, which are mainly located at the beginning of the notebook, well separed from calculations. 

## Quick Start
Create a directory where you want to have both the source code and the working (build) directory (e.g. geant4-apps) and move (from terminal) to this path, then clone the repository with:
```
git clone git@github.com:paternog/PositronSource.git
```
Create the build directory and move to this path:
```
mkdir PositronSource-build
cd PositronSource-build
```
Create the makefile with Cmake and compile the application:
```
cmake ../PositronSource
make -j[Ncores]
```
Run the app GUI which will show you a default geometry. You could also launch a test particle by selecting from the menu Run>Beam and then pushing the green "play" button on the command bar or typing `/run/beanOn 1` in the UI terminal (where you see "Session:").

## Version
version: 1.0,
date: 02-10-2024

## Contact
Please, if you make some modifications, does not push them directly on the main branch.

For any question, comment, suggestion or bug report, please contact me at <paterno@fe.infn.it>. 

//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// gpaterno, September 2024
//
/// \file DetectorConstructionMessenger.hh
/// \brief Description of the DetectorConstruction messenger class
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstructionMessenger_h
#define DetectorConstructionMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADouble;
class G4UIcmdWithABool;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithAString;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// Detector construction messenger class to define custom commands
/// to control the geometry and other settings.

class DetectorConstructionMessenger: public G4UImessenger
{
public:
    DetectorConstructionMessenger(DetectorConstruction* mpga);
    ~DetectorConstructionMessenger();

    virtual void SetNewValue(G4UIcommand* command, G4String newValues);

private:
    DetectorConstruction* fDetector;
    
    G4UIcmdWithABool* fHybridSourceCmd;
    
    G4UIdirectory* fCmdDir;   
    G4UIcmdWithAString* fCrystalMaterialCmd;
    G4UIcmdWith3VectorAndUnit* fCrystalSizeCmd;
    G4UIcmdWithAString* fCrystalLatticeCmd;
    G4UIcmdWithADouble* fCrystalAngleXCmd;
    G4UIcmdWithADouble* fCrystalAngleYCmd;
    G4UIcmdWithADouble* fCrystalBendingAngleCmd;
    G4UIcmdWithABool* fRadModelCmd;
    G4UIcmdWithABool* fOCeffectsCmd;
    G4UIcmdWithABool* fRadiatorCmd;
    G4UIcmdWithAString* fPotentialPathCmd;
    G4UIcmdWithABool* fTaggingCmd;
    G4UIcmdWithAString* fTaggingFilenameCmd;
    G4UIcmdWithAnInteger* fTaggingIntervalCmd;
    
    G4UIcmdWithADoubleAndUnit* fRadiatorConverterSepDistanceCmd;
    G4UIcmdWith3VectorAndUnit* fConverterSizeCmd;
    G4UIcmdWithABool* fGranularConverterCmd;
    G4UIcmdWithADoubleAndUnit* fSphereRadiusCmd;
    G4UIcmdWithAString* fConverterMaterialCmd;
    
    G4UIcmdWithABool* fMagneticFieldCmd;
    G4UIcmdWithADoubleAndUnit* fFieldValueCmd;
    G4UIcmdWithADoubleAndUnit* fFieldRegionLengthCmd;
    
    G4UIcmdWithABool* fCollimatorCmd;
    G4UIcmdWithADoubleAndUnit* fCollimatorApertureCmd;
    G4UIcmdWithAString* fCollimatorHoleCmd;
    G4UIcmdWithADoubleAndUnit* fCollimatorThicknessCmd;
    G4UIcmdWithADoubleAndUnit* fCollimatorSideCmd;
    G4UIcmdWithADoubleAndUnit* fRadiatorCollimatorSepDistanceCmd;
    
    G4UIcmdWith3VectorAndUnit* fVirtualDetectorSizeCmd;
    
    G4UIcmdWithABool* fSetVoxelizationCmd;
    G4UIcmdWithAnInteger* fAbsorberColumnsCmd;
    G4UIcmdWithAnInteger* fAbsorberRowsCmd;
    G4UIcmdWithAnInteger* fAbsorberSlicesCmd;
    G4UIcmdWithADoubleAndUnit* fAbsorberxVoxelSpacingCmd;
    G4UIcmdWithADoubleAndUnit* fAbsorberyVoxelSpacingCmd;
    G4UIcmdWithADoubleAndUnit* fAbsorberzVoxelSpacingCmd;  
    
    G4UIcmdWithABool* fScoringCrystalExitCmd; 
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


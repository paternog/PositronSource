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
/// \file PrimaryGeneratorActionMessenger.hh
/// \brief Description of the PrimaryGeneratorActionMessenger class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef PrimaryGeneratorActionMessenger_h
#define PrimaryGeneratorActionMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class PrimaryGeneratorAction;
class G4UIcmdWithABool;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithADouble;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// messenger for PrimaryGenerator class. It is possible to set if we want to
/// simulate using FastSim model or import an external phase-space.

class PrimaryGeneratorActionMessenger: public G4UImessenger
{
public:
    PrimaryGeneratorActionMessenger(PrimaryGeneratorAction*);
    virtual ~PrimaryGeneratorActionMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
private:   
    PrimaryGeneratorAction* fPrimaryGeneratorAction;
    
    G4UIcmdWithABool* fReadFromFileCmd;
    G4UIcmdWithAString* fSetFileNameCmd;
    G4UIcmdWithABool* fUseGPSCmd;
    
    G4UIcmdWithAString* fPrimaryTypeCmd; 
    G4UIcmdWithADoubleAndUnit* fPrimaryEnergyCmd; 
    G4UIcmdWithADouble* fPrimaryRelSigmaEnergyCmd; 
    G4UIcmdWithADoubleAndUnit* fPrimaryXCmd;  
    G4UIcmdWithADoubleAndUnit* fPrimaryYCmd;
    G4UIcmdWithADoubleAndUnit* fPrimaryZCmd;
    G4UIcmdWithADoubleAndUnit* fPrimaryTCmd;
    G4UIcmdWithADoubleAndUnit* fPrimaryXpCmd;
    G4UIcmdWithADoubleAndUnit* fPrimaryYpCmd;  
    G4UIcmdWithADoubleAndUnit* fPrimarySigmaXCmd; 
    G4UIcmdWithADoubleAndUnit* fPrimarySigmaYCmd; 
    G4UIcmdWithADoubleAndUnit* fPrimarySigmaZCmd;
    G4UIcmdWithADoubleAndUnit* fPrimarySigmaTCmd; 
    G4UIcmdWithADoubleAndUnit* fPrimarySigmaXpCmd; 
    G4UIcmdWithADoubleAndUnit* fPrimarySigmaYpCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


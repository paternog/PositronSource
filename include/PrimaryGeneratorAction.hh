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
/// \file PrimaryGeneratorAction.hh
/// \brief Description of the PrimaryGeneratorAction class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "PrimaryGeneratorActionMessenger.hh"

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

#include "FileReader.hh"

class G4GeneralParticleSource;
class G4ParticleGun;
class G4Event;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// PrimaryGeneratorAction class. We have both the GPS and a Particle Gun. The 
/// latter can be used to set the primary beam according to the phase-space file
/// obtained trough an external code and read through the FileReader class. 
/// A set of custom commands based on the Particle Gun are also defined to better 
/// simulate a bunch of particles in a particle accelerator.  

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    PrimaryGeneratorAction();
    virtual ~PrimaryGeneratorAction();
    
    void GeneratePrimaries(G4Event*);
    
    void ReadFromFile(G4bool vBool) {fReadFromFile = vBool;};   
    void SetFileName(G4String);
    void SetUseGPS(G4bool vBool) {fUseGPS = vBool;};
    
    void SetType(G4String val) {fType = val;}   
    void SetEnergy(G4double val) {fEnergy = val;}
    void SetRelSigmaEnergy(G4double val) {fRelSigmaEnergy = val;}
    void SetX(G4double val) {fX = val;}
    void SetY(G4double val) {fY = val;}
    void SetZ(G4double val) {fZ = val;}
    void SetT(G4double val) {fT = val;}
    void SetXp(G4double val) {fXp = val;}
    void SetYp(G4double val) {fYp = val;}
    void SetSigmaX(G4double val) {fSigmaX = val;}
    void SetSigmaY(G4double val) {fSigmaY = val;}
    void SetSigmaZ(G4double val) {fSigmaZ = val;}
    void SetSigmaT(G4double val) {fSigmaT = val;}
    void SetSigmaXp(G4double val) {fSigmaXp = val;}
    void SetSigmaYp(G4double val) {fSigmaYp = val;}

private:
    PrimaryGeneratorActionMessenger* fMessenger;

    G4GeneralParticleSource* fGPS;
    G4ParticleGun* fGun;
    
    G4bool fReadFromFile;
    G4String fFileName;
    static FileReader* fFileReader;
    
    G4bool fUseGPS;
    G4double fCrystalZ;
    
    G4String fType;
    G4double fEnergy;
    G4double fRelSigmaEnergy;
    G4double fX;
    G4double fY;
    G4double fZ;
    G4double fT;
    G4double fXp;
    G4double fYp;
    G4double fSigmaX;
    G4double fSigmaY;
    G4double fSigmaZ;
    G4double fSigmaT;
    G4double fSigmaXp;
    G4double fSigmaYp;  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


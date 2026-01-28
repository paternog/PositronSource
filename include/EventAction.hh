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
/// \file EventAction.hh
/// \brief Definition of the EventAction class
//
// gpaterno, January 2026
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"

#include "globals.hh"
#include <iostream>
#include <fstream>
#include <map>

#include "Run.hh"
#include "G4RunManager.hh"

class RunAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// Event action class. Used to score the energy deposited in some important
/// volumes (radiator, converter and spheres of a granluar traget).
/// Also, Edep in a custom voxelization is scored through AddEdep method.

class EventAction : public G4UserEventAction
{
public:
    EventAction();
    ~EventAction() override;

    void BeginOfEventAction(const G4Event* event) override;
    void EndOfEventAction(const G4Event* event) override;

    //Custom methods 
    inline void SetVerbose(G4int val) {fVerboseLevel = val;}
    inline G4int GetVerbose() const {return fVerboseLevel;}
    
    G4int GetEventID() const {
        return G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
    }
    
    void AddEdepRad(G4double val) {fEdepRad += val;}
    void AddEdepConv(G4double val) {fEdepConv += val;}
    
    void AddEdep(G4double x, G4double y, G4double z, G4double E) {
        fRun = static_cast<Run*>
            (G4RunManager::GetRunManager()->GetNonConstCurrentRun());
        fRun->AddEdep(x, y, z, E);
    }
     
    void AddEdepInSpheres(G4int, G4double);

private:
    Run* fRun{nullptr};
    
    G4int fSensitiveDetector_ID = -1;
    G4int fVerboseLevel = 0;
    
    G4double fEdepRad = 0.;
    G4double fEdepConv = 0.;
    
    G4int fNSpheres = 0;
    std::map<G4int,G4double> fEdepSpheres;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


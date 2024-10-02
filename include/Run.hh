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
/// \file Run.hh
/// \brief Description of the Run class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Run_h
#define Run_h 1

#include "G4Run.hh"
#include "globals.hh"

#include "VoxelScorer.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// Run class. Used to define a custom VoxelScorer and accumulate Edep and  
/// Edep^2 in a given volume.

class Run : public G4Run
{
public:
    Run();
    virtual ~Run();

    virtual void RecordEvent(const G4Event*);
    virtual void Merge(const G4Run*);
    
    G4double GetEdep() const {return fEdep;}    
    G4double GetEdep2() const {return fEdep2;}    
    void AddEdep(G4double, G4double, G4double, G4double);
    
    G4double GetEdep(G4int t, G4int r, G4int u) const {
        return fEdepScorer->GetValue(t,r,u);
    }
    G4int GetVoxNx() const {return fEdepScorer->GetX_vox();}
    G4int GetVoxNy() const {return fEdepScorer->GetY_vox();}
    G4int GetVoxNz() const {return fEdepScorer->GetZ_vox();}     
    void CleanVoxel (G4int t, G4int r, G4int u) const { 
        fEdepScorer->Reset(t,r,u);
    }
     
private:
    G4int fCollID_edep;   
    G4double fEdep;
    G4double fEdep2;        
    
    G4bool fVoxelization;
    VoxelScorer* fEdepScorer;
    
    G4bool fRadiator;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
  

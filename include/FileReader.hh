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
/// \file FileReader.hh
/// \brief Description of the FileReader class
//
// gpaterno, January 2026
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef FileReader_h
#define FileReader_h 1

#include "G4ThreeVector.hh"
#include "globals.hh"
#include <vector>
#include <fstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// File Reader class used to read the phase-space obtained through another code.

class FileReader
{
public:
    FileReader(G4String);
    FileReader();
    ~FileReader();
    
public:
    void SetFileName(G4String);
    void StoreEvents();
    
    G4String GetAnEventParticle(G4int);
    G4ThreeVector GetAnEventPosition(G4int);
    G4ThreeVector GetAnEventMomentum(G4int);
    G4double GetAnEventEnergy(G4int);
    G4int GetNumberOfEvents();

private:
    G4String fFileName = "";
    std::ifstream inputFile;
    std::vector<G4String> evListPart;
    std::vector<G4ThreeVector> evListPos;
    std::vector<G4ThreeVector> evListMom;
    std::vector<G4double> evListEnergy;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


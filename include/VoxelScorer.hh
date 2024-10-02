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
/// \file VoxelScorer.hh
/// \brief Definition of the VoxelScorer class
//
// based on AGATA code, optimezed by gpaterno on 08/05/2022
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef VoxelScorer_h
#define VoxelScorer_h 1

#include <iostream>
#include "G4SystemOfUnits.hh"

#include "G4SDManager.hh"
#include <cmath>
#include <stdlib.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// It is a class that define a custom voxel scorer.

class VoxelScorer 
{
public:
    VoxelScorer(G4double, G4double, G4double, 
                G4int, G4int, G4int,
                G4double, G4double, G4double);

    virtual ~VoxelScorer();
    
    void AddValue(G4double, G4double, G4double, G4double);
    G4double GetValue(G4int, G4int, G4int);
    void Reset(G4int, G4int, G4int); 
    
    G4int GetX_vox() {return x_vox;}; 
    G4int GetY_vox() {return y_vox;}; 
    G4int GetZ_vox() {return z_vox;}; 
    
    G4double*** fScoreMatrix; //if it is private, merging does not work!
    
private :
    G4double x_start;
    G4double y_start;
    G4double z_start;
    G4double x_stop;
    G4double y_stop;
    G4double z_stop;
    G4double x_voxsize;
    G4double y_voxsize;
    G4double z_voxsize;
    G4int x_vox;
    G4int y_vox;
    G4int z_vox;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


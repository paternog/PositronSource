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
/// \file VoxelScorer.cc
/// \brief Definition of the VoxelScorer class
//
// based on AGATA code, optimezed by gpaterno on 08/05/2022
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "VoxelScorer.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

VoxelScorer::VoxelScorer(G4double x_len, G4double y_len, G4double z_len,
                         G4int nx_voxel, G4int ny_voxel, G4int nz_voxel,
                         G4double x_position, 
                         G4double y_position, 
                         G4double z_position)
                               
{
    x_start = x_position - (x_len/2);
    y_start = y_position - (y_len/2);
    z_start = z_position - (z_len/2);
    x_stop = x_position + (x_len/2);
    y_stop = y_position + (y_len/2);
    z_stop = z_position + (z_len/2);
    x_voxsize = x_len/nx_voxel;
    y_voxsize = y_len/ny_voxel;
    z_voxsize = z_len/nz_voxel;
    x_vox = nx_voxel;
    y_vox = ny_voxel;
    z_vox = nz_voxel;
        
    fScoreMatrix = new G4double**[nx_voxel];
    for (G4int i(0); i < nx_voxel; i++) {
        fScoreMatrix[i] = new G4double*[ny_voxel];
        for (G4int j(0); j < ny_voxel; j++) {
            fScoreMatrix[i][j] = new G4double[nz_voxel];
            for (G4int k(0); k < nz_voxel; k++) {
                fScoreMatrix[i][j][k] = 0.;
            }
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

VoxelScorer::~VoxelScorer()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void VoxelScorer::AddValue(G4double x_p, G4double y_p, 
                           G4double z_p, G4double quantity)
{ 
    if (x_p >= x_start && x_p < x_stop && 
        y_p >= y_start && y_p < y_stop && 
        z_p >= z_start && z_p < z_stop) {

        G4int x_pos = (G4int)((x_p - x_start)/x_voxsize);
        G4int y_pos = (G4int)((y_p - y_start)/y_voxsize);
        G4int z_pos = (G4int)((z_p - z_start)/z_voxsize);

        fScoreMatrix[x_pos][y_pos][z_pos] += quantity;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double VoxelScorer::GetValue(G4int x_n, G4int y_n, G4int z_n)
{
    return fScoreMatrix[x_n][y_n][z_n];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void VoxelScorer::Reset(G4int x_n, G4int y_n, G4int z_n)
{
    fScoreMatrix[x_n][y_n][z_n] = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


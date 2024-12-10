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
// gpaterno, December 2024
//
/// \file Run.cc
/// \brief Implementation of the Run class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Run.hh"
#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"

#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4THitsMap.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::Run(): 
G4Run(),
fCollID_edep(-1),
fEdep(0.),
fEdep2(0.),
fGoodEvents(0)
{
    //get an instance of the DetectorConstruction
    const DetectorConstruction* detectorConstruction = 
        static_cast<const DetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    
    //retrieve information about Absorber (virtual) voxelization
    fVoxelization = detectorConstruction->GetVoxelization();
        
    if (fVoxelization) {                    
        G4int totalColumns = detectorConstruction->GetTotalColumns();
        G4int totalRows = detectorConstruction->GetTotalRows(); 
        G4int totalSlices = detectorConstruction->GetTotalSlices();       
        G4double xVoxelSpacing = detectorConstruction->GetxVoxelSpacing();
        G4double yVoxelSpacing = detectorConstruction->GetyVoxelSpacing(); 
        G4double zVoxelSpacing = detectorConstruction->GetzVoxelSpacing();           
        G4double absorberZ = detectorConstruction->GetAbsorberZ();
        /*  
        G4double absorberDensity = detectorConstruction->GetAbsorberVolume()->
                                                         GetMaterial()->GetDensity();
        G4double voxelVolume = xVoxelSpacing*yVoxelSpacing*zVoxelSpacing;  
        G4double voxelMass = absorberDensity*voxelVolume;  
               
        G4cout << "AbsorberVolume: "
               << detectorConstruction->GetAbsorberVolume()->GetName() << G4endl;
        G4cout << "AbsorberZ: " << absorberZ/mm << " mm" << G4endl;                                           
        G4cout << "totalColumns: " << totalColumns << G4endl;
        G4cout << "totalRows: " << totalRows << G4endl;
        G4cout << "totalSlices: " << totalSlices << G4endl;
        G4cout << "xVoxelSpacing: " << xVoxelSpacing/mm << " mm" << G4endl;
        G4cout << "yVoxelSpacing: " << yVoxelSpacing/mm << " mm" << G4endl;
        G4cout << "zVoxelSpacing: " << zVoxelSpacing/mm << " mm" << G4endl;        
        G4cout << "voxel density: " << absorberDensity/(g/cm3) << " g/cm3" << G4endl; 
        G4cout << "voxel volume: " << voxelVolume/(cm3) << " cm3" << G4endl;   
        G4cout << "voxel mass: " << voxelMass/kg << " kg" << G4endl;        
        */            
        fEdepScorer = new VoxelScorer(totalColumns*xVoxelSpacing,
                                      totalRows*yVoxelSpacing,
                                      totalSlices*zVoxelSpacing,
                                      totalColumns,
                                      totalRows,                                        
                                      totalSlices,
                                      0.,
                                      0.,
                                      absorberZ);
    }
    
    //retrieve if we are using the Radiator Crystal
    fRadiator = detectorConstruction->GetRadiator();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::~Run() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::RecordEvent(const G4Event* event)
{
    //Retrieve total Edep in the Radiator Crystal
    if (fRadiator) {
        //Hit Detection System  
        if (fCollID_edep < 0) {
            G4String fCollName = "multisd/edep";
            //G4cout << "fCollName: " << fCollName << G4endl;
            fCollID_edep = 
                G4SDManager::GetSDMpointer()->GetCollectionID(fCollName);
        }
                    
        //get all Hits Collections of this Event
        G4HCofThisEvent* HCE = event->GetHCofThisEvent();
        if (!HCE) return;

        //fCollID_edep Hits Collection (map)
        G4THitsMap<G4double>* evtMap = 
            static_cast<G4THitsMap<G4double>*>(HCE->GetHC(fCollID_edep));

        //map iterator
        std::map<G4int,G4double*>::iterator itr;

        //scan the hits Collection with Energy deposit 
        G4double edep = 0.;
        G4double EdepEvent = 0.;   
        for (itr = evtMap->GetMap()->begin(); itr != evtMap->GetMap()->end(); itr++) {
            edep = *(itr->second);
        
            //accumulate Edep in the whole Run
            fEdep += edep;
            fEdep2 += edep*edep;
                       
            //accumulate Edep in the Event
            EdepEvent += edep;
        }
        if (EdepEvent > 0) fGoodEvents++;
       
    }


    //call base method
    G4Run::RecordEvent(event);      
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::Merge(const G4Run* aRun)
{
    
    //accumulate variables from each instance (one for thread) of Run
    const Run* localRun = static_cast<const Run*>(aRun);
    fEdep += localRun->fEdep;
    fEdep2 += localRun->fEdep2;
    fGoodEvents += localRun->fGoodEvents;

      
    //score Edep through custom VoxelScorer (fEdepScorer)
    if (fVoxelization) { 
        for (G4int b=0; b<(fEdepScorer->GetX_vox()) ; b++) {
            for (G4int c=0; c<(fEdepScorer->GetY_vox()) ; c++) {
                for (G4int v=0; v<(fEdepScorer->GetZ_vox()) ; v++) {
                    fEdepScorer->fScoreMatrix[b][c][v] += 
                        localRun->fEdepScorer->fScoreMatrix[b][c][v];
                    localRun->fEdepScorer->Reset(b,c,v);
                }
            }
        }
    }


    //call base method (to pass accumulated variables to the master thread)
    G4Run::Merge(aRun);     
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::AddEdep(G4double a, G4double b, G4double c, G4double quant)
{
    fEdepScorer->AddValue(a, b, c, quant/MeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


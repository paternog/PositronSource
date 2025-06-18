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
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"

#include "G4AnalysisManager.hh"
#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(EventAction* eventAction)
: fEventAction(eventAction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{
    //get an instance of the DetectorConstruction and retrieve some settings 
    const DetectorConstruction* detectorConstruction
        = static_cast<const DetectorConstruction*>
          (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    
    //get the Sensitive Volumes
    fCrystalVolume = detectorConstruction->GetCrystalVolume(); //Radiator
    fConverterVolume = detectorConstruction->GetConverterVolume(); //Converter/Target
    fAbsorberVolume = detectorConstruction->GetAbsorberVolume(); //Absorber (it could be both)
    if (fScoringVolume.size() == 0)
        fScoringVolume = detectorConstruction->GetScoringVolume(); //Spheres of the granular target
    G4int NScoringVolumes = fScoringVolume.size();
      
    //get if I want to score Edep in a voxelized Volume
    G4bool Voxelization = detectorConstruction->GetVoxelization();
    
    //get if I want to score the features of particles
    //exiting the raditior and /or the target (27/09/2024)
    G4bool scoreCrystalExit = detectorConstruction->GetScoringCrystalExit();
  
  
    //get pre and post step points
    G4StepPoint* preStepPoint = step->GetPreStepPoint();
    G4StepPoint* postStepPoint = step->GetPostStepPoint();
        
    //get the current and next particle postion
    G4ThreeVector preStepPos = preStepPoint->GetPosition();
    G4ThreeVector postStepPos = postStepPoint->GetPosition();
    G4ThreeVector pos = preStepPos + G4UniformRand()*(postStepPos - preStepPos);
          
    //get the volume of the current step
    G4LogicalVolume* volume = 
        preStepPoint->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
    G4String volumeName = volume->GetName();
  	//G4cout << "volume: " << volumeName << G4endl;
      
    //get track and particle name
    G4Track* track = step->GetTrack();
    //G4double trackWeight = track->GetWeight();
    G4String partName = track->GetDefinition()->GetParticleName();
    G4int trackID = track->GetTrackID();
    //G4int parentID = track->GetParentID();

    //get the energy deposited during in this step
    G4double edep = step->GetTotalEnergyDeposit();
    
    //get eventID
    G4int eventID = 
        G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
               
    
    //instantiating The Analysis Manager
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
       
    
    //declaration of variables useful for the scoring 
    //of the features of the particles exiting the crystals
    //G4ThreeVector postStepPos;
    G4ThreeVector postStepMom;
    G4double postStepTime;
    G4LogicalVolume* volumeNext;
    
    
    //score the Edep in the Crystal (the Radiator)
    if (volume == fCrystalVolume) {     
        fEventAction->AddEdepRad(edep);
        
        //score the features of the particle exiting the volume
        //if (scoreCrystalExit && postStepPoint->GetStepStatus() == fGeomBoundary) { //it doesn't work...
        //postStepPos = postStepPoint->GetPosition();
		postStepMom = postStepPoint->GetMomentum(); 
		postStepTime = postStepPoint->GetGlobalTime();
		volumeNext = 
		    postStepPoint->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
        //G4cout << "volume: " << volumeName << ", next volume: " << volumeNextName << G4endl;
        if (scoreCrystalExit && volumeNext && volume != volumeNext) {
            analysisManager->FillNtupleSColumn(4,0,partName);
            analysisManager->FillNtupleDColumn(4,1,postStepPos.x()/CLHEP::mm);
            analysisManager->FillNtupleDColumn(4,2,postStepPos.y()/CLHEP::mm);
            analysisManager->FillNtupleDColumn(4,3,postStepPos.z()/CLHEP::mm);
            analysisManager->FillNtupleDColumn(4,4,postStepMom.x()/CLHEP::MeV);
            analysisManager->FillNtupleDColumn(4,5,postStepMom.y()/CLHEP::MeV);
            analysisManager->FillNtupleDColumn(4,6,postStepMom.z()/CLHEP::MeV);
            analysisManager->FillNtupleDColumn(4,7,postStepTime/CLHEP::ns);
            analysisManager->FillNtupleIColumn(4,8,eventID);
            analysisManager->FillNtupleIColumn(4,9,trackID);
            analysisManager->AddNtupleRow(4);          
        }
        
    }
    
    
    //score the Edep in the Converter (the Target)
    if (volume == fConverterVolume) {     
        fEventAction->AddEdepConv(edep);
        
        //score the features of the particle exiting the volume
        //postStepPos = postStepPoint->GetPosition();
		postStepMom = postStepPoint->GetMomentum(); 
		postStepTime = postStepPoint->GetGlobalTime();
		volumeNext = 
		    postStepPoint->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
        if (scoreCrystalExit && volumeNext && volume != volumeNext) {
            analysisManager->FillNtupleSColumn(4,0,partName);
            analysisManager->FillNtupleDColumn(4,1,postStepPos.x()/CLHEP::mm);
            analysisManager->FillNtupleDColumn(4,2,postStepPos.y()/CLHEP::mm);
            analysisManager->FillNtupleDColumn(4,3,postStepPos.z()/CLHEP::mm);
            analysisManager->FillNtupleDColumn(4,4,postStepMom.x()/CLHEP::MeV);
            analysisManager->FillNtupleDColumn(4,5,postStepMom.y()/CLHEP::MeV);
            analysisManager->FillNtupleDColumn(4,6,postStepMom.z()/CLHEP::MeV);
            analysisManager->FillNtupleDColumn(4,7,postStepTime/CLHEP::ns);
            analysisManager->FillNtupleIColumn(4,8,eventID);
            analysisManager->FillNtupleIColumn(4,9,trackID);
            analysisManager->AddNtupleRow(4);          
        } 
        
    }
    
    
    //score the Edep in the voxelized volume (the Absorber)
    if (volume == fAbsorberVolume) {     
        if (edep > 0 && Voxelization) {
            fEventAction->AddEdep(pos.x(), 
                                  pos.y(),
                                  pos.z(),
                                  edep);
        }  
    }
      
      
    //score the Edep in the spheres of the granular target 
    //(NScoringVolumes !=0 only if the target is granular)
    for (int i = 0; i < NScoringVolumes; i++) {
        if (volume == fScoringVolume[i]) {    
            if (edep > 0) {
                fEventAction->AddEdepInSpheres(i, edep);            
            } 
            
            if (edep > 0 && Voxelization) {
                fEventAction->AddEdep(pos.x(),
                                      pos.y(),
                                      pos.z(),
                                      edep);
            }            
        }
    }
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


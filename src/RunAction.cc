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
/// \file RunAction.cc
/// \brief Implementation of the RunAction class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"
#include "RunActionMessenger.hh"
#include "DetectorConstruction.hh"
#include "Run.hh"
#include "G4FastSimulationManager.hh"
#include "G4ChannelingFastSimModel.hh"

#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4RegionStore.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include <iostream>
#include <fstream>
#include <string>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction() : G4UserRunAction()
{
    G4RunManager::GetRunManager()->SetPrintProgress(1000);

    fMessenger = new RunActionMessenger(this);

    //Create the analysis manager
    fAnalysisManager = G4AnalysisManager::Instance();
    fAnalysisManager->SetDefaultFileType("root");
    fAnalysisManager->SetFileName(fFileName);
#ifdef G4MULTITHREADED
    fAnalysisManager->SetNtupleMerging(true);
#else
    fAnalysisManager->SetNtupleMerging(false);
#endif
    fAnalysisManager->SetVerboseLevel(0);

    //Create the ntuple to score the particles impinging on the virtual screens
    fAnalysisManager->CreateNtuple("scoring_ntuple","virtual scoring screens");
    fAnalysisManager->CreateNtupleIColumn("screenID");
    fAnalysisManager->CreateNtupleSColumn("particle");
    fAnalysisManager->CreateNtupleDColumn("x");
    fAnalysisManager->CreateNtupleDColumn("y");
    fAnalysisManager->CreateNtupleDColumn("px");
    fAnalysisManager->CreateNtupleDColumn("py");
    fAnalysisManager->CreateNtupleDColumn("pz");
    fAnalysisManager->CreateNtupleDColumn("t");
    fAnalysisManager->CreateNtupleIColumn("eventID");
    fAnalysisManager->FinishNtuple();

    //Create the ntuple to score Edep in the Radiatior
    fAnalysisManager->CreateNtuple("edep_rad","radiator");
    fAnalysisManager->CreateNtupleDColumn("edep");
    fAnalysisManager->CreateNtupleIColumn("eventID");
    fAnalysisManager->FinishNtuple();

    //Create the ntuple to score Edep in the Converter
    fAnalysisManager->CreateNtuple("edep_conv","converter");
    fAnalysisManager->CreateNtupleDColumn("edep");
    fAnalysisManager->CreateNtupleIColumn("eventID");
    fAnalysisManager->FinishNtuple();

    //Create the ntuple to score Edep in the spheres of the Granular Target
    fAnalysisManager->CreateNtuple("edep_spheres","granular target");
    fAnalysisManager->CreateNtupleIColumn("volumeID");
    fAnalysisManager->CreateNtupleDColumn("edep");
    fAnalysisManager->CreateNtupleIColumn("eventID");
    fAnalysisManager->FinishNtuple();

    //Create the ntuple to score the particles exiting the crystals
    fAnalysisManager->CreateNtuple("scoring_ntuple2","particle leaving crystals");
    fAnalysisManager->CreateNtupleSColumn("particle");
    fAnalysisManager->CreateNtupleDColumn("x");
    fAnalysisManager->CreateNtupleDColumn("y");
    fAnalysisManager->CreateNtupleDColumn("z");
    fAnalysisManager->CreateNtupleDColumn("px");
    fAnalysisManager->CreateNtupleDColumn("py");
    fAnalysisManager->CreateNtupleDColumn("pz");
    fAnalysisManager->CreateNtupleDColumn("t");
    fAnalysisManager->CreateNtupleIColumn("eventID");
    fAnalysisManager->CreateNtupleIColumn("trackID");
    fAnalysisManager->FinishNtuple();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
    //Write results
    if (fIsFileOpened) {
        fAnalysisManager->Write();
        //fAnalysisManager->CloseFile(); //I do it in the main   

        G4int threadID = G4Threading::G4GetThreadId();
        if (threadID > 0) {
            G4cout << "writing results of thread " << G4Threading::G4GetThreadId()
                   << " to root file: " << fFileName << G4endl;
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Run* RunAction::GenerateRun() {return new Run;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* run)
{
    //Merging is only for root
    if (fFileName.find(".csv") != std::string::npos) {
        fAnalysisManager->SetNtupleMerging(false);
    }
 
    //Open the output file
    if (!fIsFileOpened) {
        fAnalysisManager->OpenFile(fFileName);
        fIsFileOpened = true;
    }

    //Print begin message
    if (IsMaster()) {
        G4int NumberOfEventToBeProcessed = run->GetNumberOfEventToBeProcessed();
  
        G4cout
        << G4endl
        << "--------------------Begin of Global Run-----------------------" 
        << G4endl
        << "Number of events to be processed: " << NumberOfEventToBeProcessed
        << G4endl
        << "--------------------------------------------------------------"
        << G4endl;
    }    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* run)
{
    //Get the number of events
    G4int nofEvents = run->GetNumberOfEvent();
    if (nofEvents == 0) return;

    //Get results from the Edep (in the Radiator Crystal) scorer
    const Run* aRun = static_cast<const Run*>(run);
    G4double Edep  = aRun->GetEdep();
    G4double Edep2 = aRun->GetEdep2();
    G4int nGoodEvents = aRun->GoodEvents();
    G4double stdEdep = sqrt(Edep2 - Edep*Edep/nGoodEvents); 

    //Print and write results to a text file
    if (IsMaster()) {
        //DetectorConstruction instance
        const DetectorConstruction* detectorConstruction = 
            static_cast<const DetectorConstruction*>(G4RunManager::GetRunManager()
                ->GetUserDetectorConstruction());

        G4bool bRadiator = detectorConstruction->GetRadiator();
        if (bRadiator) { 
            G4String matName = 
                detectorConstruction->GetCrystalVolume()->GetMaterial()->GetName();

            //Set the significant digits to print
            G4cout.precision(8);
  
            //Print the results on screen
            G4cout << G4endl
            << "--------------------End of Global Run-----------------------"
            << G4endl
            << " The run had " << nofEvents << " events";
            G4cout << G4endl
            << " Edep in " << "the Radiator Crystal" << " (made of "
            << matName << "): " << Edep/MeV
            << " +/- " << stdEdep/MeV << " MeV" << G4endl 
            << "------------------------------------------------------------" 
            << G4endl << G4endl;
            
            /*
            //Write the Edep in the radiator on a text file  
            G4double mass = detectorConstruction->GetCrystalVolume()->GetMass();
            std::ofstream fFileOut;
            fFileOut.open("output/SimulationSummary.dat", 
                          std::ofstream::out | std::ofstream::app);
            fFileOut.precision(8);
            fFileOut << "seed: " << G4Random::getTheSeed() << " | Edep: " 
                     << Edep/MeV << " +/- " << stdEdep/MeV << 
                     << " MeV in the Radiator made of " << matName
                     << " (dose: " << (Edep/mass)/(gray) << " +/- " 
                     << (stdEdep/mass)/(gray) << " Gy)"
                     << std::endl;
            fFileOut.close();
            */
        }
                
        //Write results of Edep in the VoxelScorer to a text file
        std::string pureFileName = fFileName.substr(0, fFileName.find(".root"));
        if (pureFileName.length() == fFileName.length())
            pureFileName = fFileName.substr(0, fFileName.find(".csv"));
        //G4cout << "pureFileName: " << pureFileName << G4endl;
        G4bool Voxelization = detectorConstruction->GetVoxelization();
        if (Voxelization) { 
            G4int x_n = aRun->GetVoxNx();
            G4int y_n = aRun->GetVoxNy();
            G4int z_n = aRun->GetVoxNz();
            
            G4cout << "Writing file for voxel scorer in the Absorber..." 
                   << G4endl;
            G4cout << "x_n: " << x_n << G4endl;
            G4cout << "y_n: " << y_n << G4endl;
            G4cout << "z_n: " << z_n << G4endl;
            
            G4String edep_file;
            std::stringstream s1;
            s1 << pureFileName << "_AbsorberEdepDistribution.txt";
            edep_file = s1.str();
            std::ofstream edep_stream;
            edep_stream.open(edep_file);
            edep_stream.precision(6);

            //first line: voxel structure [nx, ny, nz]
            //second line: values (nx) of columns (end of a column is new line),
            //then for rows, then for slices
            edep_stream << x_n << " " << y_n << " " << z_n << '\n';
  
            for (G4int k=0; k < z_n; k++) {
                for (G4int i=0; i < x_n; i++) {
                    for (G4int j=0; j < y_n; j++) {
                        edep_stream << aRun->GetEdep(i,j,k) << " ";
                        aRun->CleanVoxel(i,j,k);
                    }
                    edep_stream << '\n';
                }
            }

            edep_stream.close();
            G4cout << "...file written!" << G4endl << G4endl;
        }
        
    }
   
    //End message
    G4cout
        << G4endl
        << " The run consisted of " << nofEvents << " particles" << G4endl
        << "------------------------------------------------------------"
        << G4endl << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::SetFileName(G4String filename) 
{
    if (filename != "") fFileName = filename;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


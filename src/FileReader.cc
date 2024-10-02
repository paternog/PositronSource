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
/// \file FileReader.cc
/// \brief Implementation of the FileReader class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "FileReader.hh"

#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FileReader::FileReader(G4String fileName)
{
    fFileName = fileName;
    inputFile.open(fFileName.data());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FileReader::~FileReader()
{
    inputFile.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FileReader::SetFileName(G4String vFileName)
{
    fFileName=vFileName;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FileReader::StoreEvents()
{    
    if (evListPos.size() == 0) {     
        G4String particle = "geantino";  
        G4double x = 0.;
        G4double y = 0.;      
        G4double z = 0.;     
        G4double px = 0.;
        G4double py = 0.;
        G4double pz = 0.;
        //G4double u = 0.;
        //G4double v = 0.;
        //G4double w = 0.;
        G4double p = 0.;
        G4double mass = 0.;
        G4double E = 0.;
        G4double K = 0.;
                
        G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
        
        while (inputFile.good()) {
            if (inputFile.good()) inputFile >> particle;
            if (inputFile.good()) inputFile >> x;
            if (inputFile.good()) inputFile >> y;
            if (inputFile.good()) inputFile >> px;
            if (inputFile.good()) inputFile >> py;
            if (inputFile.good()) inputFile >> pz;
            
            //x = 3.*(G4UniformRand()-0.5)*mm; //Alexei Sytov's method
            //y = 3.*(G4UniformRand()-0.5)*mm;
            
            p = sqrt(px*px+py*py+pz*pz);
            //u = px/p; //they are not required, I can use directly the momentum components!
            //v = py/p;
            //w = pz/p;                          
            mass = particleTable->FindParticle(particle)->GetPDGMass();   
            E = sqrt(p*p + mass*mass)*MeV;
            K = E - mass; //kinetic energy
            if (K <= 0) {K = 1.*eV;}

            evListPart.push_back(particle);
            evListPos.push_back(G4ThreeVector(x*cm, y*cm, z*cm));
            evListMom.push_back(G4ThreeVector(px*MeV, py*MeV, pz*MeV));
            //evListMom.push_back(G4ThreeVector(u, v, w)); 
            evListEnergy.push_back(K);     
            
            /*
            G4cout.precision(17);
            G4cout << particle << " " << x << " " << y << " " << px << " " << py << " " << pz
                   << G4endl << "p: " << p << " MeV, mass: " << mass/MeV << " MeV, E: " << E/MeV 
                   << " MeV, K: " << K/MeV << " MeV" << G4endl;
            */
        }
        
        //for some reason the last row is read two times, so we remove it!
        evListPart.pop_back();
        evListPos.pop_back();
        evListMom.pop_back();
        evListEnergy.pop_back();
        
        G4cout << "the file contains " << evListPart.size() << " particles" << G4endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4String FileReader::GetAnEventParticle(G4int vIndex)
{
    return evListPart.at(vIndex);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector FileReader::GetAnEventPosition(G4int vIndex)
{
    return evListPos.at(vIndex);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double FileReader::GetAnEventEnergy(G4int vIndex)
{
    return evListEnergy.at(vIndex);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ThreeVector FileReader::GetAnEventMomentum(G4int vIndex)
{
    return evListMom.at(vIndex);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int FileReader::GetNumberOfEvents()
{
    return evListPos.size();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


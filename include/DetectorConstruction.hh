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
/// \file DetectorConstruction.hh
/// \brief Description of the DetectorConstruction class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4ios.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include <vector>

#include "G4Region.hh"
#include "G4PVPlacement.hh"

#include "DetectorConstructionMessenger.hh"
#include "G4ChannelingFastSimModel.hh"

#define NSpheresMax 10000

class G4VPhysicalVolume;
class G4LogicalVolume;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// Detector construction class to define materials and geometry.

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    DetectorConstruction();
    ~DetectorConstruction() override;

    G4VPhysicalVolume* Construct() override;
    virtual void ConstructSDandField();
    
    //method to get the scoring volumes
    std::vector<G4LogicalVolume*> GetScoringVolume() const {
        return fScoringVolume;}        
    
    //method to set if it is a hybrid source or not
    void SetHybridSource(G4bool val) {fHybridSource = val;}
   
    //methods to set the Crystal (Radiator) features
    void SetCrystalMaterial(G4String val) {fCrystalMaterialStr = val;}
    void SetCrystalSize(G4ThreeVector val) {fCrystalSize = val;}
    void SetCrystalBendingAngle(G4double val) {fBendingAngle = val;}    
    void SetCrystalLattice(G4String val) {fLattice = val;}
    void SetCrystalAngleX(G4double val) {fAngleX = val;}
    void SetCrystalAngleY(G4double val) {fAngleY = val;}
    G4double GetCrystalZ() const {return fCrystalZ;}
    void SetRadiationModel(G4bool val) {fActivateRadiationModel = val;}
    void SetOCeffects(G4bool val) {fActivateOCeffects = val;}
    G4bool GetOCeffects() const {return fActivateOCeffects;}
    G4LogicalVolume* GetCrystalVolume() const {return fCrystalLogic;}
    void SetRadiator(G4bool val) {fRadiator = val;}
    G4bool GetRadiator() const {return fRadiator;}
    void SetPotentialPath(const G4String path){fPotentialPath = path;}
    void SetTagging(G4bool bval) {fTagging = bval;}
    void SetTaggingFilename(G4String filename) {fTaggingFilename = filename;}
    void SetTaggingInterval(G4int printInterval) {fTaggingInterval = printInterval;}   
    
    //method to set/get the Converter (Target) features
    void SetRadiatorConverterSepDistance(G4double val) {
        fRadiatorConverterSepDistance = val;}
    G4double GetRadiatorConverterSepDistance() const {
        return fRadiatorConverterSepDistance;}
    void SetConverterSize(G4ThreeVector val) {fConverterSize = val;}
    void SetConverterMaterial(G4String val) {fConverterMaterialStr = val;}
    void SetGranularConverter(G4bool val) {fGranularConverter = val;}
    void SetSphereRadius(G4double val) {fSphereRadius = val;}
    G4int GetNSpheres() const {return fNSpheres;}
    G4LogicalVolume* GetConverterVolume() const {return fConverterLogic;}
       
    //methods to set the Magnetic field features
    void SetMagneticField(G4bool val) {fSetMagneticField = val;}
    void SetFieldValue(G4double val) {fFieldValue = val;}
    void SetFieldRegionLength(G4double val) {fFieldRegionLength = val;}
       
    //methods to set the Collimator features
    void SetCollimator(G4bool val) {fSetCollimator = val;}
    void SetCollimatorHole(G4String val) {fCollimatorHole = val;}
    void SetCollimatorAperture(G4double val) {fCollimatorAperture = val;}
    void SetCollimatorThickness(G4double val) {fCollimatorThickness = val;}
    void SetCollimatorSide(G4double val) {fCollimatorSide = val;}
    void SetRadiatorCollimatorSepDistance(G4double val) {
        fRadiatorCollimatorSepDistance = val;}
    G4double GetRadiatorCollimatorSepDistance() const {
        return fRadiatorCollimatorSepDistance;}
       
    //methods to set/Get the VirtualDetector features
    void SetVirtualDetectorSize(G4ThreeVector val) {fVirtualDetectorSize = val;}
    std::vector<G4ThreeVector> GetVirtualDetectorPositionVector() const {
        return fVirtualDetectorPositionVector;}
    
    //for the voxelization of the Absorber
    void SetVoxelization(G4bool bval) {fIWantVoxelization = bval;} 
    void SetTotalColumns(G4int cols) {ftotalColumns = cols;} 
    void SetTotalRows(G4int rows) {ftotalRows = rows;} 
    void SetTotalSlices(G4int slices) {ftotalSlices = slices;} 
    void SetxVoxelSpacing(G4double dx) {fxVoxelSpacing = dx;}
    void SetyVoxelSpacing(G4double dy) {fyVoxelSpacing = dy;}
    void SetzVoxelSpacing(G4double dz) {fzVoxelSpacing = dz;}
    G4bool GetVoxelization() const {return fIWantVoxelization;}
    G4int GetTotalColumns() const {return ftotalColumns;}
    G4int GetTotalRows() const {return ftotalRows;}
    G4int GetTotalSlices() const {return ftotalSlices;}
    G4double GetxVoxelSpacing() const {return fxVoxelSpacing;}
    G4double GetyVoxelSpacing() const {return fyVoxelSpacing;}
    G4double GetzVoxelSpacing() const {return fzVoxelSpacing;}
    G4double GetAbsorberZ() const {return fAbsorberZ;}
    G4LogicalVolume* GetAbsorberVolume() const {return fAbsorberLogic;}
    
    //methods to set and get ScoreCrystalExit (27/09/2024)
    void SetScoringCrystalExit(G4bool bval) {fScoringCrystalExit = bval;} 
    G4bool GetScoringCrystalExit() const {return fScoringCrystalExit;}
        
protected:
  std::vector<G4LogicalVolume*> fScoringVolume; //for spheres only

private:
    DetectorConstructionMessenger* fMessenger;  
    
    G4bool fHybridSource;
        
    G4Region* fCrystalRegion;
    G4LogicalVolume* fCrystalLogic;
    G4String fCrystalMaterialStr;
    G4Material* fCrystalMaterial;
    G4ThreeVector fCrystalSize;
    G4double fBendingAngle;
    G4String fLattice;  
    G4double fAngleX;
    G4double fAngleY;
    G4double fCrystalZ;
    G4bool fActivateRadiationModel;
    G4bool fActivateOCeffects;
    G4bool fRadiator;
    G4String fPotentialPath;
    G4bool fTagging = false;
    G4String fTaggingFilename = "output/output_for_tagging.txt";
    G4int fTaggingInterval = 100;
    
    G4double fRadiatorConverterSepDistance;
    G4ThreeVector fConverterSize;
    G4double fConverterZ;
    G4LogicalVolume* fConverterLogic;
    G4bool fGranularConverter;
    G4String fConverterMaterialStr;
    G4Material* fConverterMaterial;
    G4double fSphereRadius;
    G4LogicalVolume* fSphereLogic[NSpheresMax];
    G4int fNSpheres;
    G4bool fConverter;
       
    G4bool fSetMagneticField;
    G4double fFieldValue;
    G4double fFieldRegionLength;
    G4LogicalVolume* fMFlogic; 
    
    G4bool fSetCollimator;
    G4double fCollimatorAperture;
    G4String fCollimatorHole;
    G4double fCollimatorThickness;
    G4double fCollimatorSide;
    G4double fRadiatorCollimatorSepDistance;
    G4LogicalVolume* fCollimatorLogic; 
      
    G4ThreeVector fVirtualDetectorSize;
    std::vector<G4ThreeVector> fVirtualDetectorPositionVector;
    G4LogicalVolume* fVirtualDetectorLogic0;
    G4LogicalVolume* fVirtualDetectorLogic1;
    G4LogicalVolume* fVirtualDetectorLogic2;
        
    G4bool fIWantVoxelization;
    G4int ftotalColumns;
    G4int ftotalRows;
    G4int ftotalSlices;
    G4double fxVoxelSpacing; 
    G4double fyVoxelSpacing;     
    G4double fzVoxelSpacing;
    G4double fAbsorberZ;
    G4LogicalVolume* fAbsorberLogic = nullptr;
    
    G4bool fScoringCrystalExit = false;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

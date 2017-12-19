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
// $Id: B4DetectorConstruction.cc 87359 2014-12-01 16:04:27Z gcosmo $
// 
/// \file B4DetectorConstruction.cc
/// \brief Implementation of the B4DetectorConstruction class

#include "B4DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4OpBoundaryProcess.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal 
G4GlobalMagFieldMessenger* B4DetectorConstruction::fMagFieldMessenger = 0; 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4DetectorConstruction::B4DetectorConstruction()
 : G4VUserDetectorConstruction(),
   modulePV(0),
   fCheckOverlaps(true)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4DetectorConstruction::~B4DetectorConstruction()
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B4DetectorConstruction::Construct()
{
  // Define materials 
  DefineMaterials();
  
  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4DetectorConstruction::DefineMaterials()
{ 
  // Copper material defined using NIST Manager
  // I use Cu as default absorber material but you can switch to lead
  G4NistManager* CunistManager = G4NistManager::Instance();
  CunistManager->FindOrBuildMaterial("G4_Cu");
  
  // Lead material defined using NIST Manager
 // G4NistManager* PbnistManager = G4NistManager::Instance();
 // PbnistManager->FindOrBuildMaterial("G4_Pb");

  // Polystyrene material defined using NIST Manager
  // I use this material for the core of plastic scintillating fibers
  // cannot find any G4_Polystyrene, I build it later
  //G4NistManager* PynistManager = G4NistManager::Instance();
  //PynistManager->FindOrBuildMaterial("G4_Polystyrene");

  // PMMA material, there's no default G4_PMMA, I build it (C502H8)
  G4String name, symbol;    // a=mass of a mole;
  G4double a, z;            // z=mean number of protons;  
  
  // create elements
  a = 1.01*g/mole;
  G4Element* elH  = new G4Element(name="Hydrogen",symbol="H" , z= 1., a); //Hidrogen

  a = 12.01*g/mole;
  G4Element* elC  = new G4Element(name="Carbon"  ,symbol="C" , z= 6., a); //Carbon

  a = 16.00*g/mole;
  G4Element* elO  = new G4Element(name="Oxygen"  ,symbol="O" , z= 8., a); //Oxygen

  a = 28.09*g/mole;
  G4Element* elSi = new G4Element(name="Silicon", symbol="Si", z=14., a); //Silicon
  
  a = 18.9984*g/mole;
  G4Element* elF  = new G4Element("Fluorine",symbol="F" , z= 9., a); //Fluorine
  
  // create PMMA
  G4Material* PMMA = new G4Material("PMMA", 1.19*g/cm3, 3);
  PMMA -> AddElement(elC, 5);
  PMMA -> AddElement(elO, 2);
  PMMA -> AddElement(elH, 8); //PMMA building complete

  // create Polystyrene (C5H5)
  G4Material* Polystyrene = new G4Material("Polystyrene", 1.05*g/cm3, 2);
  Polystyrene -> AddElement(elC, 8);
  Polystyrene -> AddElement(elH, 8); //Polystyrene building complete

  // create Fluorinated Polymer (C2F2)
  // I use it for the cladding of the Cherenkov fibers
  G4Material* fluorinatedPolymer =
  new G4Material("Fluorinated_Polymer", 1.43*g/cm3, 2);
  fluorinatedPolymer->AddElement(elC,2);
  fluorinatedPolymer->AddElement(elF,2);
  //fluorinatedPolymer->AddElement(H,2); //Fluorinated Polymer building complete

  // create Glass (SiO2)
  G4Material* Glass = new G4Material("Glass", 2.4*g/cm3, 2);
  Glass -> AddElement(elSi, 1);
  Glass -> AddElement(elO, 2); //Glass building complete

  // Vacuum material defined using NIST Manager
  G4NistManager* VanistManager = G4NistManager::Instance();
  VanistManager->FindOrBuildMaterial("G4_Galactic");

  // Silicon material defined using NIST Manager
  G4NistManager* SinistManager = G4NistManager::Instance();
  SinistManager->FindOrBuildMaterial("G4_Si");

  // Air material defined using NIST Manager
  // You can use Air instead of vacuum
  //G4NistManager* AinistManager = G4NistManager::Instance();
  //AinistManager->FindOrBuildMaterial("G4_AIR");

  // Print materials 
  // I don't want to print materials all the times,
  // if you want uncomment it
  //G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B4DetectorConstruction::DefineVolumes()
{
  // Geometry parameters of world, module, fibers, SiPM

  // Geometry parameters of the module
  G4int Nofmodules = 71; //the actual number of modules is Nofmodules^2, choose 3,5,7,9
  G4int NofFibers = 32; // 32 of each type
  G4int NofScinFibers = NofFibers/2;
  G4int NofCherFibers = NofFibers/2;
  G4int NofFibersrow = NofFibers/4;
  G4int NofFiberscolumn = NofFibersrow;
  G4double moduleZ = 2623.*mm;
  G4double moduleX = 10.14*mm; 
  G4double moduleY = moduleX;

  // Geometry parameters of the world, world is a box
  G4double worldX = 200 * moduleX;
  G4double worldY = 200 * moduleY;
  G4double worldZ = 60 * moduleZ;

  // Geometry parameters of the fiber
  G4double fiberradius = 0.5*mm;
  G4double fiberZ = moduleZ;

  // Geometry parameters of the core
  G4double coreradius = 0.485*mm;
  G4double coreZ = moduleZ;

  // Geometry parameters of the cladding
  G4double claddingradiusmin = 0.485*mm;
  G4double claddingradiusmax = 0.50*mm;
  G4double claddingZ = moduleZ;

  // Geometry parameters of the SiPM
  G4double SiPMX = 1.*mm;
  G4double SiPMY = SiPMX;
  G4double SiPMZ = 0.36*mm;

  // Geometry parameters of the SiPM, active silicon layer
  G4double SiX = 1.*mm;
  G4double SiY = SiX;
  G4double SiZ = 0.05*mm;

  // Geometry parameters of the module equipped with SiPM
  // I build it so i can replicate the entire module + SiPM 
  G4double moduleequippedZ = moduleZ + SiPMZ;
  G4double moduleequippedX = moduleX; 
  G4double moduleequippedY = moduleY;

  // Get materials for vacuum, absorber, scintillating and cherenkov fibers, SiPM
  G4Material* defaultMaterial = G4Material::GetMaterial("G4_Galactic"); // you can switch to air 
  G4Material* absorberMaterial = G4Material::GetMaterial("G4_Cu"); // G4_Cu
  G4Material* ScinMaterial = G4Material::GetMaterial("Polystyrene");
  G4Material* CherMaterial = G4Material::GetMaterial("PMMA");
  G4Material* GlassMaterial = G4Material::GetMaterial("Glass");
  G4Material* SiMaterial = G4Material::GetMaterial("G4_Si");
  G4Material* CladCherMaterial = G4Material::GetMaterial("Fluorinated_Polymer");

  // I need to specify the optical proprieties of the scintillating fiber material,
  // optical proprieties are different from scintillating proprieties and 
  // scintillating proprieties will be defined later.
  // We don't have to add WLS proprieties to scintillating fibers

  const G4int ENTRIES = 32;
  
  G4double photonEnergy[ENTRIES] =                    // Use Energy(eV)=1.24/waevelenght(um)
            { 2.034*eV, 2.068*eV, 2.103*eV, 2.139*eV, // 2.034eV is 610nm RED  
              2.177*eV, 2.216*eV, 2.256*eV, 2.298*eV,     
              2.341*eV, 2.386*eV, 2.433*eV, 2.481*eV,
              2.532*eV, 2.585*eV, 2.640*eV, 2.697*eV,
              2.757*eV, 2.820*eV, 2.885*eV, 2.954*eV, // 2.75eV is 450nm BLUE (peak of scintillating fibers)
              3.026*eV, 3.102*eV, 3.181*eV, 3.265*eV, // 3.09eV is 400nm VIOLET (end of visible)
              3.353*eV, 3.446*eV, 3.545*eV, 3.649*eV,
              3.760*eV, 3.877*eV, 4.002*eV, 4.136*eV }; //4.1eV is 300nm UV (cherenkov peak is 310-350nm)

  G4double rindexScin[ENTRIES] =
            { 1.59, 1.59, 1.59, 1.59,
              1.59, 1.59, 1.59, 1.59,
              1.59, 1.59, 1.59, 1.59,
              1.59, 1.59, 1.59, 1.59,
              1.59, 1.59, 1.59, 1.59,
              1.59, 1.59, 1.59, 1.59,
              1.59, 1.59, 1.59, 1.59,
              1.59, 1.59, 1.59, 1.59 };

  G4double absorptionScin[ENTRIES] =
             { 400*cm, 400*cm, 400*cm, 400*cm,
               400*cm, 400*cm, 400*cm, 400*cm,
               400*cm, 400*cm, 400*cm, 400*cm,
               400*cm, 400*cm, 400*cm, 400*cm,
               400*cm, 400*cm, 400*cm, 400*cm,
               400*cm, 400*cm, 400*cm, 400*cm,
               400*cm, 400*cm, 400*cm, 400*cm,
               400*cm, 400*cm, 400*cm, 400*cm };

  // I don't want any ABSLENGTH for the scintillating and cherenkov fibers
  // I take into account in parameterization of photon transportation
  // if you want uncomment it             
  G4MaterialPropertiesTable *MPTScin = new G4MaterialPropertiesTable();
  MPTScin -> AddProperty("RINDEX", photonEnergy, rindexScin, ENTRIES)->SetSpline(true);
  //MPTScin -> AddProperty("ABSLENGTH", photonEnergy, absorptionScin, ENTRIES)->SetSpline(true);

  // I need to specify the optical proprieties of the cherenkov fiber material
  // there are no scintillating proprieties for PMMA (clear fibres)
  // we don't have to add WLS proprieties

  G4double rindexCher[ENTRIES] =
            { 1.49, 1.49, 1.49, 1.49,
              1.49, 1.49, 1.49, 1.49,
              1.49, 1.49, 1.49, 1.49,
              1.49, 1.49, 1.49, 1.49,
              1.49, 1.49, 1.49, 1.49,
              1.49, 1.49, 1.49, 1.49,
              1.49, 1.49, 1.49, 1.49,
              1.49, 1.49, 1.49, 1.49 };

 G4double absorptionCher[ENTRIES] = 
            { 890*cm, 890*cm, 890*cm, 890*cm,
              890*cm, 890*cm, 890*cm, 890*cm,
              890*cm, 890*cm, 890*cm, 890*cm,
              890*cm, 890*cm, 890*cm, 890*cm,
              890*cm, 890*cm, 890*cm, 890*cm,
              890*cm, 890*cm, 890*cm, 890*cm,
              890*cm, 890*cm, 890*cm, 890*cm,
              890*cm, 890*cm, 890*cm, 890*cm };

  G4MaterialPropertiesTable *MPTCher = new G4MaterialPropertiesTable();
  MPTCher -> AddProperty("RINDEX", photonEnergy, rindexCher, ENTRIES)->SetSpline(true);
  //MPTCher -> AddProperty("ABSLENGTH", photonEnergy, absorptionCher, ENTRIES)->SetSpline(true);
  CherMaterial -> SetMaterialPropertiesTable(MPTCher);

  // I need to specify the optical proprieties of the cherenkov cladding material

  G4double rindexCherclad[ENTRIES] =
            { 1.42, 1.42, 1.42, 1.42,
              1.42, 1.42, 1.42, 1.42,
              1.42, 1.42, 1.42, 1.42,
              1.42, 1.42, 1.42, 1.42,
              1.42, 1.42, 1.42, 1.42,
              1.42, 1.42, 1.42, 1.42,
              1.42, 1.42, 1.42, 1.42,
              1.42, 1.42, 1.42, 1.42 };

  G4MaterialPropertiesTable *MPTCherclad = new G4MaterialPropertiesTable();
  MPTCherclad -> AddProperty("RINDEX", photonEnergy, rindexCherclad, ENTRIES)->SetSpline(true);
  CladCherMaterial -> SetMaterialPropertiesTable(MPTCherclad);

  // We don't set any optical proriety for the absorber, if you want uncomment it
  // I need to specify the optical proprieties of the absorber material

  /*G4double rindexabsorber[ENTRIES] =
            { 1.1, 1.1, 1.1, 1.1,
              1.1, 1.1, 1.1, 1.1,
              1.1, 1.1, 1.1, 1.1,
              1.1, 1.1, 1.1, 1.1,
              1.1, 1.1, 1.1, 1.1,
              1.1, 1.1, 1.1, 1.1,
              1.1, 1.1, 1.1, 1.1,
              1.1, 1.1, 1.1, 1.1 };

  G4double absorptionabsorber[ENTRIES] =   // I set 1nm absorption lenght so light doesn't 
            { 1.*nm, 1.*nm, 1.*nm, 1.*nm,  // propagate inside copper
              1.*nm, 1.*nm, 1.*nm, 1.*nm,
              1.*nm, 1.*nm, 1.*nm, 1.*nm,
              1.*nm, 1.*nm, 1.*nm, 1.*nm,
              1.*nm, 1.*nm, 1.*nm, 1.*nm,
              1.*nm, 1.*nm, 1.*nm, 1.*nm,
              1.*nm, 1.*nm, 1.*nm, 1.*nm,
              1.*nm, 1.*nm, 1.*nm, 1.*nm };
             
  G4MaterialPropertiesTable *MPTabsorber = new G4MaterialPropertiesTable();
  MPTabsorber -> AddProperty("RINDEX", photonEnergy, rindexabsorber, ENTRIES)->SetSpline(true);
  MPTabsorber -> AddProperty("ABSLENGTH", photonEnergy, absorptionabsorber, ENTRIES)->SetSpline(true);
  absorberMaterial -> SetMaterialPropertiesTable(MPTabsorber);*/

  // I need to specify the optical proprieties of the air material
  // turn on if you use Air

  /*G4double rindexair[ENTRIES] =
              { 1.0003, 1.0003, 1.0003, 1.0003,
                1.0003, 1.0003, 1.0003, 1.0003,
                1.0003, 1.0003, 1.0003, 1.0003,
                1.0003, 1.0003, 1.0003, 1.0003,
                1.0003, 1.0003, 1.0003, 1.0003,
                1.0003, 1.0003, 1.0003, 1.0003,
                1.0003, 1.0003, 1.0003, 1.0003,
                1.0003, 1.0003, 1.0003, 1.0003 };

  G4MaterialPropertiesTable *MPTair = new G4MaterialPropertiesTable();
  MPTair -> AddProperty("RINDEX", photonEnergy, rindexair, ENTRIES)->SetSpline(true);
  defaultMaterial -> SetMaterialPropertiesTable(MPTair);
*/
  // I need to specify the optical proprieties of the vacuum material

  G4double rindexvacuum[ENTRIES] = // By definition
            { 1.0, 1.0, 1.0, 1.0,
              1.0, 1.0, 1.0, 1.0,
              1.0, 1.0, 1.0, 1.0,
              1.0, 1.0, 1.0, 1.0,
              1.0, 1.0, 1.0, 1.0,
              1.0, 1.0, 1.0, 1.0,
              1.0, 1.0, 1.0, 1.0,
              1.0, 1.0, 1.0, 1.0 };

  G4MaterialPropertiesTable *MPTvacuum = new G4MaterialPropertiesTable();
  MPTvacuum -> AddProperty("RINDEX", photonEnergy, rindexvacuum, ENTRIES)->SetSpline(true);
  defaultMaterial -> SetMaterialPropertiesTable(MPTvacuum);

  // I need to specify the optical proprieties of the glass material

  G4double rindexglass[ENTRIES] =
            { 1.51, 1.51, 1.51, 1.51,
              1.51, 1.51, 1.51, 1.51,
              1.51, 1.51, 1.51, 1.51,
              1.51, 1.51, 1.51, 1.51,
              1.51, 1.51, 1.51, 1.51,
              1.51, 1.51, 1.51, 1.51,
              1.51, 1.51, 1.51, 1.51,
              1.51, 1.51, 1.51, 1.51 };

  G4MaterialPropertiesTable *MPTglass = new G4MaterialPropertiesTable();
  MPTglass -> AddProperty("RINDEX", photonEnergy, rindexglass, ENTRIES)->SetSpline(true);
  GlassMaterial -> SetMaterialPropertiesTable(MPTglass);

  // I need to specify the optical proprieties of the Si material

  G4double rindexSi[ENTRIES] =
            { 3.42, 3.42, 3.42, 3.42,
              3.42, 3.42, 3.42, 3.42,
              3.42, 3.42, 3.42, 3.42,
              3.42, 3.42, 3.42, 3.42,
              3.42, 3.42, 3.42, 3.42,
              3.42, 3.42, 3.42, 3.42,
              3.42, 3.42, 3.42, 3.42,
              3.42, 3.42, 3.42, 3.42 };

  G4double absorptionSi[ENTRIES] = 
            { 0.001*mm, 0.001*mm, 0.001*mm, 0.001*mm,
              0.001*mm, 0.001*mm, 0.001*mm, 0.001*mm,
              0.001*mm, 0.001*mm, 0.001*mm, 0.001*mm,
              0.001*mm, 0.001*mm, 0.001*mm, 0.001*mm,
              0.001*mm, 0.001*mm, 0.001*mm, 0.001*mm,
              0.001*mm, 0.001*mm, 0.001*mm, 0.001*mm,
              0.001*mm, 0.001*mm, 0.001*mm, 0.001*mm,
              0.001*mm, 0.001*mm, 0.001*mm, 0.001*mm };

  G4MaterialPropertiesTable *MPTSi = new G4MaterialPropertiesTable();
  MPTSi -> AddProperty("RINDEX", photonEnergy, rindexSi, ENTRIES)->SetSpline(true);
  MPTSi -> AddProperty("ABSLENGHT", photonEnergy, absorptionSi, ENTRIES)->SetSpline(true);
  SiMaterial -> SetMaterialPropertiesTable(MPTSi); 
  
  // I need to specify the SCINTILLATING proprieties of the scintillating fiber material
  // I specify also the Birk Constant of the polystyrene

  G4double Scin_FAST[ENTRIES] = // Emission spectrum for the fast component 
            { 0., 0., 0., 0.,
              0., 0., 0., 0.,
              0., 0., 0., 0.1,
              0.2, 0.4, 0.6, 0.8,
              1., 0.8, 0.6, 0.1,
              0., 0., 0., 0.,
              0., 0., 0., 0.,
              0., 0., 0., 0. };

  G4double Scin_SLOW[ENTRIES] = // Emission spectrum for the slow component
            { 0., 0., 0., 0.,
              0., 0., 0., 0.,
              0., 0., 0., 0.,
              0., 0., 0., 0.,
              0., 0., 0., 0.,
              0., 0., 0., 0.,
              0., 0., 0., 0.,
              0., 0., 0., 0. };

  // Set Briks Constant for scintillator
  ScinMaterial->GetIonisation()->SetBirksConstant(0.126*mm/MeV);

  MPTScin -> AddProperty("FASTCOMPONENT", photonEnergy, Scin_FAST, ENTRIES);
  MPTScin -> AddProperty("SLOWCOMPONENT", photonEnergy, Scin_SLOW, ENTRIES);
  MPTScin -> AddConstProperty("SCINTILLATIONYIELD", 10000./MeV); // Typical is 10000./MeV (this is what makes full simulations long as hell)
  MPTScin -> AddConstProperty("RESOLUTIONSCALE", 1.0); // Broad the fluctuation of photons produced
  MPTScin -> AddConstProperty("FASTTIMECONSTANT", 2.8*ns);
  MPTScin -> AddConstProperty("SLOWTIMECONSTANT", 10.*ns);
  MPTScin -> AddConstProperty("YIELDRATIO", 1.0); // I don't want a slow component, if you want it must change
  ScinMaterial -> SetMaterialPropertiesTable(MPTScin);
  
  if ( ! defaultMaterial || ! absorberMaterial || ! ScinMaterial || ! CherMaterial || ! GlassMaterial || ! CladCherMaterial ) {
    G4ExceptionDescription msg;
    msg << "Cannot retrieve materials already defined."; 
    G4Exception("B4DetectorConstruction::DefineVolumes()",
      "MyCode0001", FatalException, msg);
  }
   
  // Building the calorimeter

  // Here I build the world

  G4VSolid* worldS 
    = new G4Box("World",                        // its name
                 worldX/2, worldY/2, worldZ/2); // its size
                         
  G4LogicalVolume* worldLV
    = new G4LogicalVolume(
                 worldS,           // its solid
                 defaultMaterial,  // its material (Galactic or Air)
                 "World");         // its name
  
  // I set the world as invisible
  worldLV->SetVisAttributes(G4VisAttributes::Invisible);
                                   
  G4VPhysicalVolume* worldPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 worldLV,          // its logical volume                         
                 "World",          // its name
                 0,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

   // Here I build the module equipped with SiPM

   G4VSolid* moduleequippedS
    = new G4Box("moduleequipped",                                          // its name
                 moduleequippedX/2, moduleequippedY/2, moduleequippedZ/2); // its size
                         
  G4LogicalVolume* moduleequippedLV
    = new G4LogicalVolume(
                 moduleequippedS,           // its solid
                 defaultMaterial,           // its material
                 "moduleequipped");         // its name

  // Here I place the single module equipped with SiPM in the center of the world
  // I rotate it a bit in order to avoid problems due to the perfect alinment 
  // between the beam and the module. (1.;1.51)*deg
  /*
  G4RotationMatrix rotm  = G4RotationMatrix();
  rotm.rotateY(2.0*deg); // Set the rotation angles
  rotm.rotateX(2.0*deg); // 0.*deg no rotation!     
  G4ThreeVector position;
  position.setX(0.);
  position.setY(0.);
  position.setZ(0.);
  G4Transform3D transform = G4Transform3D(rotm,position);

  G4VPhysicalVolume* moduleequippedPV = new G4PVPlacement(
                                                transform,        // its position and rotation
                                                moduleequippedLV, // its logical volume                         
                                                "moduleequipped", // its name
                                                worldLV,          // its mother  volume
                                                false,            // no boolean operation
                                                0,                // copy number
                                                fCheckOverlaps);  // checking overlaps 
  
  */
  // Here I build the calorimeter itself. As calorimeter I mean the matrix of
  // modules equipped. Uncomment it only if you want more than one module.
  
   G4VSolid* CalorimeterS 
    = new G4Box("CalorimeterS",                                                                  // its name
                 moduleequippedX*Nofmodules/2, moduleequippedY*Nofmodules/2, moduleequippedZ/2); // its size
                         
  G4LogicalVolume* CalorimeterLV
    = new G4LogicalVolume(
                 CalorimeterS,           // its solid
                 defaultMaterial,        // its material 
                 "CalorimeterLV");       // its name
  
  // Here I place the modules equipped inside the calorimeter
  // There is no rotation of the modules, I will later rotate the entire calorimeter

  G4double m_x, m_y;
  G4ThreeVector vec_m;
  G4VPhysicalVolume* physi_moduleequipped[Nofmodules][Nofmodules];
  for(int row=0; row<Nofmodules; row++){
     for(int column=0; column<Nofmodules; column++){
        m_x = ((Nofmodules-1)/2)*moduleX - moduleX*row;
        m_y = ((Nofmodules-1)/2)*moduleY - moduleY*column;
           
        vec_m.setX(m_x);
        vec_m.setY(m_y);
        vec_m.setZ(0.);

        physi_moduleequipped[row][column] = new G4PVPlacement(0,
                                                        vec_m,              
                                                        moduleequippedLV,     
                                                        "moduleequipped",                        
                                                        CalorimeterLV,                      
                                                        false,                          
                                                        0); 
      };
   };           

  // Here I place and rotate the entire calorimeter

  G4RotationMatrix rotm  = G4RotationMatrix();
  rotm.rotateY(0.75*deg);  // Set the rotation angles//0.75
  rotm.rotateX(1.0*deg);      //1.0
  G4ThreeVector position;
  position.setX(0.);
  position.setY(0.);
  position.setZ(0.);
  G4Transform3D transform = G4Transform3D(rotm,position); 

  G4VPhysicalVolume* CalorimeterPV = new G4PVPlacement(
                                                transform,        // its position and rotation
                                                CalorimeterLV,    // its logical volume                         
                                                "Calorimeter",    // its name
                                                worldLV,          // its mother  volume
                                                false,            // no boolean operation
                                                0,                // copy number
                                                fCheckOverlaps);  // checking overlaps 
              
  // Here I build the module: to do that I build the rectangular absorber
  // I will later put fibers into it  

  G4VSolid* moduleS
    = new G4Box("module",                          // its name
                 moduleX/2, moduleY/2, moduleZ/2); // its size
                         
  G4LogicalVolume* moduleLV
    = new G4LogicalVolume(
                 moduleS,           // its solid
                 absorberMaterial,  // its material
                 "module");         // its name

  G4ThreeVector pos_module;
  pos_module.setX(0.);
  pos_module.setY(0.);
  pos_module.setZ(-0.18);
                              
  G4VPhysicalVolume* modulePV = new G4PVPlacement(
                                                0,                // no rotation
                                                pos_module,       // at (0,0,-0.18)
                                                moduleLV,         // its logical volume                         
                                                "module",         // its name
                                                moduleequippedLV,          // its mother  volume
                                                false,            // no boolean operation
                                                0,                // copy number
                                                fCheckOverlaps);  // checking overlaps 


  // Here I define the Optical Surface PROPRIETIES between the glass and the Si of the SiPM

  G4OpticalSurface* OpSurfaceGlassSi = new G4OpticalSurface("OpSurfaceGlassSi");
  
  OpSurfaceGlassSi -> SetType(dielectric_metal);
  OpSurfaceGlassSi -> SetModel(glisur);
  OpSurfaceGlassSi -> SetFinish(polished);

  G4double efficiencyOpSurfaceGlassSi[ENTRIES] =     // detection efficiency 
                                    { 0.4, 0.4, 0.4, 0.4,
                                      0.4, 0.4, 0.4, 0.4,
                                      0.4, 0.4, 0.4, 0.4,
                                      0.4, 0.4, 0.4, 0.4,
                                      0.4, 0.4, 0.4, 0.4,
                                      0.4, 0.4, 0.4, 0.4,
                                      0.4, 0.4, 0.4, 0.4,
                                      0.4, 0.4, 0.4, 0.4};
                                      

   G4double reflectivityOpSurfaceGlassSi[ENTRIES] =  // 0% reflection
                                    { 0., 0., 0., 0.,
                                      0., 0., 0., 0.,
                                      0., 0., 0., 0.,
                                      0., 0., 0., 0.,
                                      0., 0., 0., 0.,
                                      0., 0., 0., 0.,
                                      0., 0., 0., 0.,
                                      0., 0., 0., 0. };

  G4MaterialPropertiesTable* MPTOpSurfaceGlassSi = new G4MaterialPropertiesTable();
  MPTOpSurfaceGlassSi -> AddProperty("EFFICIENCY", photonEnergy, efficiencyOpSurfaceGlassSi, ENTRIES)->SetSpline(true);
  MPTOpSurfaceGlassSi -> AddProperty("REFLECTIVITY", photonEnergy, reflectivityOpSurfaceGlassSi, ENTRIES)->SetSpline(true);
  OpSurfaceGlassSi -> SetMaterialPropertiesTable(MPTOpSurfaceGlassSi);

  // Here I build the SiPM

  G4VSolid* SiPMS
    = new G4Box("SiPM",                      // its name
                 SiPMX/2, SiPMY/2, SiPMZ/2); // its size
                         
  G4LogicalVolume* SiPMLV
    = new G4LogicalVolume(
                 SiPMS,             // its solid
                 GlassMaterial,     // its material
                 "SiPM");           // its name

 // Here I build the Si of the SiPM
 
 G4VSolid* SiS
   = new G4Box("Si",                     // its name
                SiX/2, SiY/2, SiZ/2);       // its size
                         
 G4LogicalVolume* SiLV
   = new G4LogicalVolume(
                 SiS,            // its solid
                 SiMaterial,     // its material
                 "Si");          // its name

 // I put the Si inside the SiPM, I will put the SiPMs next to fibers later

 G4ThreeVector vec_Si;
 vec_Si.setX(0.);
 vec_Si.setY(0.);
 vec_Si.setZ(SiPMZ/2-SiZ/2); // Si at the end of SiPM
                             
 G4VPhysicalVolume* SiPV = new G4PVPlacement(
                                             0,                 // no rotation
                                             vec_Si,  
                                             SiLV,              // its logical volume                         
                                             "Si",              // its name
                                             SiPMLV,            // its mother  volume
                                             false,             // no boolean operation
                                             0,                 // copy number
                                             fCheckOverlaps);   // checking overlaps 
 
  // I set the visualization attributes of the Si of the SiPM
  G4VisAttributes* SiVisAtt = new G4VisAttributes(G4Colour(0.0,0.8,0.0)); //green
  SiVisAtt->SetVisibility(true);
  SiVisAtt->SetForceWireframe(true);
  SiVisAtt->SetForceSolid(true);
  SiLV->SetVisAttributes(SiVisAtt); //end of visualization attributes

  // Here I place the Logical Skin Surface around the silicon of the SiPM
  G4LogicalSkinSurface* OpsurfaceSi = new G4LogicalSkinSurface("OpsurfaceSi", SiLV, OpSurfaceGlassSi);

  // Here I define the Optical Surface PROPRIETIES between the scintillating fibers and the default material
  // air or vacuum
  // I'm trying to define an optical surface completly blacked because we absorb the light at one end of fibers

  G4OpticalSurface* OpSurfacedefault = new G4OpticalSurface("OpSurfacedefault");
  
  OpSurfacedefault -> SetType(dielectric_dielectric);
  OpSurfacedefault -> SetModel(unified);
  OpSurfacedefault -> SetFinish(polishedbackpainted); // Painted from inside the fibers, light is absorbed
 
  // Here I build the Scintillating fiber with its core and cladding
  // I will put the fibers later inside the module

  G4Tubs* S_fiber = new G4Tubs("S_fiber", 0., fiberradius, fiberZ/2, 0., 2.*pi);

  G4LogicalVolume* logic_S_fiber = new G4LogicalVolume(S_fiber,          //its solid
                                                       defaultMaterial,  //its material
                                                       "S_fiber");       //its name

  G4Tubs* Core_S_fiber = new G4Tubs("Core_S_fiber", 0., coreradius, coreZ/2, 0., 2.*pi);

  G4LogicalVolume* logic_Core_S_fiber = new G4LogicalVolume(Core_S_fiber,   //its solid
                                                            ScinMaterial,   //its material
                                                            "Core_S_fiber");//its name

  // I set the visualization attributes of the scintillating core fibers
  G4VisAttributes* ScincoreVisAtt = new G4VisAttributes(G4Colour(0.0,0.0,0.8)); //blue
  ScincoreVisAtt->SetVisibility(true);
  ScincoreVisAtt->SetForceWireframe(true);
  ScincoreVisAtt->SetForceSolid(true);
  logic_Core_S_fiber->SetVisAttributes(ScincoreVisAtt); //end of visualization attributes

  G4ThreeVector vec_Core_S;
  vec_Core_S.setX(0.);
  vec_Core_S.setY(0.);
  vec_Core_S.setZ(0.); 
                             
  G4VPhysicalVolume* Core_S_PV = new G4PVPlacement(
                                             0,                        // no rotation
                                             vec_Core_S,               // its position
                                             logic_Core_S_fiber,       // its logical volume                         
                                             "Core_S_fiber",           // its name
                                             logic_S_fiber,            // its mother  volume
                                             false,                    // no boolean operation
                                             0,                        // copy number
                                             fCheckOverlaps);          // checking overlaps
 
  // Here I place the optical surface "OpSurfacedefault" between the scintillatinf core and the default material
  G4LogicalBorderSurface* logic_OpSurface_SCoredefault;
  logic_OpSurface_SCoredefault = new G4LogicalBorderSurface("logic_OpSurface_SCoredefault", Core_S_PV, worldPV, OpSurfacedefault);

  G4Tubs* Clad_S_fiber = new G4Tubs("Clad_S_fiber", claddingradiusmin, claddingradiusmax, claddingZ/2, 0., 2.*pi);

  G4LogicalVolume* logic_Clad_S_fiber = new G4LogicalVolume(Clad_S_fiber,   //its solid
                                                            CherMaterial,   //its material
                                                            "Clad_S_fiber");//its name

  // I set the visualization attributes of the scintillating clad fibers
  G4VisAttributes* ScincladVisAtt = new G4VisAttributes(G4Colour(0.0,1.0,1.0));//light blue
  ScincladVisAtt->SetVisibility(true);
  ScincladVisAtt->SetForceWireframe(true);
  ScincladVisAtt->SetForceSolid(true);
  logic_Clad_S_fiber->SetVisAttributes(ScincladVisAtt); //end of visualization attributes

 G4ThreeVector vec_Clad_S;
 vec_Clad_S.setX(0.);
 vec_Clad_S.setY(0.);
 vec_Clad_S.setZ(0.); 
                             
 G4VPhysicalVolume* Clad_S_PV = new G4PVPlacement(
                                             0,                        // no rotation
                                             vec_Clad_S,               // its position
                                             logic_Clad_S_fiber,       // its logical volume                         
                                             "Clad_S_fiber",           // its name
                                             logic_S_fiber,            // its mother  volume
                                             false,                    // no boolean operation
                                             0,                        // copy number
                                             fCheckOverlaps);          // checking overlaps

// Here I place the optical surface "OpSurfacedefault" between the scintillating clad and the default material
G4LogicalBorderSurface* logic_OpSurface_SCladdefault;
logic_OpSurface_SCladdefault = new G4LogicalBorderSurface("logic_OpSurface_SCladdefault", Clad_S_PV, worldPV, OpSurfacedefault);

// Here I build the Cherenkov fiber with its cladding
// I will put the fibers later inside the module

G4Tubs* C_fiber = new G4Tubs("C_fiber", 0., fiberradius, fiberZ/2, 0., 2.*pi);

G4LogicalVolume* logic_C_fiber = new G4LogicalVolume(C_fiber,       //it solid
                                                     CherMaterial,  //its material
                                                     "C_fiber");     //its name

G4Tubs* Core_C_fiber = new G4Tubs("Core_C_fiber", 0., coreradius, coreZ/2, 0., 2.*pi);

G4LogicalVolume* logic_Core_C_fiber = new G4LogicalVolume(Core_S_fiber,   //its solid
                                                          CherMaterial,   //its material
                                                          "Core_C_fiber");//its name

// I set the visualization attributes of the cherenkov core fibers
G4VisAttributes* ChercoreVisAtt = new G4VisAttributes(G4Colour(1.0,0.0,0.0)); //red
ChercoreVisAtt->SetVisibility(true);
ChercoreVisAtt->SetForceWireframe(true);
ChercoreVisAtt->SetForceSolid(true);
logic_Core_C_fiber->SetVisAttributes(ChercoreVisAtt); //end of visualization attributes

 G4ThreeVector vec_Core_C;
 vec_Core_C.setX(0.);
 vec_Core_C.setY(0.);
 vec_Core_C.setZ(0.); 
                             
 G4VPhysicalVolume* Core_C_PV = new G4PVPlacement(
                                             0,                        // no rotation
                                             vec_Core_C,               // its position
                                             logic_Core_C_fiber,       // its logical volume                         
                                             "Core_C_fiber",           // its name
                                             logic_C_fiber,            // its mother  volume
                                             false,                    // no boolean operation
                                             0,                        // copy number
                                             fCheckOverlaps);          // checking overlaps

// Here I place the optical surface "OpSurfacedefault" between the cherenkov core and the default material
G4LogicalBorderSurface* logic_OpSurface_CCoredefault;
logic_OpSurface_CCoredefault = new G4LogicalBorderSurface("logic_OpSurface_CCoredefault", Core_C_PV, worldPV, OpSurfacedefault);

  G4Tubs* Clad_C_fiber = new G4Tubs("Clad_C_fiber", claddingradiusmin, claddingradiusmax, claddingZ/2, 0., 2.*pi);

  G4LogicalVolume* logic_Clad_C_fiber = new G4LogicalVolume(Clad_C_fiber,   //its solid
                                                            CladCherMaterial,   //its material
                                                            "Clad_C_fiber");//its name

  // I set the visualization attributes of the cherenkov clad fibers
  G4VisAttributes* ChercladVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));//yellow 
  ChercladVisAtt->SetVisibility(true);
  ChercladVisAtt->SetForceWireframe(true);
  ChercladVisAtt->SetForceSolid(true);
  logic_Clad_C_fiber->SetVisAttributes(ChercladVisAtt); //end of visualization attributes

 G4ThreeVector vec_Clad_C;
 vec_Clad_C.setX(0.);
 vec_Clad_C.setY(0.);
 vec_Clad_C.setZ(0.); 
                             
 G4VPhysicalVolume* Clad_C_PV = new G4PVPlacement(
                                             0,                        // no rotation
                                             vec_Clad_C,               // its position
                                             logic_Clad_C_fiber,       // its logical volume                         
                                             "Clad_C_fiber",           // its name
                                             logic_C_fiber,            // its mother  volume
                                             false,                    // no boolean operation
                                             0,                        // copy number
                                             fCheckOverlaps);          // checking overlaps

 // Here I place the optical surface "OpSurfacedefault" between the cherenkov clad and the default material
G4LogicalBorderSurface* logic_OpSurface_CCladdefault;
logic_OpSurface_CCladdefault = new G4LogicalBorderSurface("logic_OpSurface_CCladdefault", Clad_C_PV, worldPV, OpSurfacedefault);

  // Here I place the Scintillating fibers and the SiPM next to them
  // Attention: I place an optical surface painted (blacked) from the moduleequippedPV 
  // to the SiPMPV, in so doing I completly avoid any cross talk between SiPMs
 
  G4VPhysicalVolume* physi_S_fiber[NofFibersrow][NofFiberscolumn];
  G4VPhysicalVolume* physi_SiPM[NofFibersrow][NofFiberscolumn];  
  G4LogicalBorderSurface* logic_OpSurface_defaultAir[NofFibersrow][NofFiberscolumn];

  for(int row=0; row<NofFibersrow; row++){
     std::stringstream S_fiber_row;
     S_fiber_row.str("");
     S_fiber_row << row;
     for(int column=0; column<NofFiberscolumn; column++){
        std::stringstream S_fiber_column;
        S_fiber_column.str("");
        S_fiber_column << column;
        std::string S_name;
        std::string SiPM_name;
        S_name = "S_row" + S_fiber_row.str() + "_column_" + S_fiber_column.str(); 
        SiPM_name = "SiPMS_row" + S_fiber_row.str() + "_column_" + S_fiber_column.str();

        // I need to specify the position of each scintillating fiber before placing them
        G4double S_x, S_y, S_z;
        G4ThreeVector vec_S_fiber;
        G4ThreeVector vec_SiPM;
        
        if(row == 0 || row == 2 || row == 4 || row == 6){
          if(column == 0 || column == 2 || column == 4 || column == 6){
            S_x = -moduleX/2 + 0.13375 + fiberradius + (0.2675+fiberradius*2)*column;
            S_y = -moduleY/2 + 0.13375 + fiberradius + (0.2675+fiberradius*2)*row;
         
            vec_S_fiber.setX(S_x);
            vec_S_fiber.setY(S_y);
            vec_S_fiber.setZ(0.);

            vec_SiPM.setX(S_x);
            vec_SiPM.setY(S_y);
            vec_SiPM.setZ(fiberZ/2+SiPMZ/2-0.18);

            // I need to place the scintillating fibers
            physi_S_fiber[row][column] = new G4PVPlacement(0,
                                                         vec_S_fiber,     //its position
                                                         logic_S_fiber,   //its logical volume
                                                         S_name,          //its name
                                                         moduleLV,        //its mother
                                                         false,           //no boulean operat
                                                         0); 

            // I need to place the SiPMs
            physi_SiPM[row][column] = new G4PVPlacement(0,
                                                        vec_SiPM,                      //its position
                                                        SiPMLV,                        //its logical volume
                                                        SiPM_name,                    //its name
                                                        moduleequippedLV,                      //its mother
                                                        false,                        //no boulean operat
                                                        0); 

          logic_OpSurface_defaultAir[NofFibersrow][NofFiberscolumn] = new G4LogicalBorderSurface("logic_OpSurface_defaultAir", CalorimeterPV, 
            physi_SiPM[row][column], OpSurfacedefault);
          }
        }
       if(row == 1 || row == 3 || row == 5 || row == 7){
         if(column == 1 || column == 3 || column == 5 || column == 7){
         S_x = -moduleX/2 + 0.13375 + fiberradius + (0.2675+fiberradius*2)*column;
         S_y = -moduleY/2 + 0.13375 + fiberradius + (0.2675+fiberradius*2)*row;

         vec_S_fiber.setX(S_x);
         vec_S_fiber.setY(S_y);
         vec_S_fiber.setZ(0.);

         vec_SiPM.setX(S_x);
         vec_SiPM.setY(S_y);
         vec_SiPM.setZ(fiberZ/2+SiPMZ/2-0.18);

         // I need to place the scintillating fibers
         physi_S_fiber[row][column] = new G4PVPlacement(0,
                                                        vec_S_fiber,     //its position
                                                        logic_S_fiber,   //its logical volume
                                                        S_name,          //its name
                                                        moduleLV,        //its mother
                                                        false,           //no boulean operat
                                                        0); 

         // I need to place the SiPMs
         physi_SiPM[row][column] = new G4PVPlacement(0,
                                                     vec_SiPM,                      //its position
                                                     SiPMLV,                        //its logical volume
                                                     SiPM_name,                    //its name
                                                     moduleequippedLV,                      //its mother
                                                     false,                        //no boulean operat
                                                     0); 
         logic_OpSurface_defaultAir[NofFibersrow][NofFiberscolumn] = new G4LogicalBorderSurface("logic_OpSurface_defaultAir", CalorimeterPV, 
           physi_SiPM[row][column], OpSurfacedefault);
         }
       }
     };
  };

  // Here I place the Cherenkov fibers
  G4VPhysicalVolume* physi_C_fiber[NofFibersrow][NofFiberscolumn];
  
  for(int row=0; row<NofFibersrow; row++){
     std::stringstream C_fiber_row;
     C_fiber_row.str("");
     C_fiber_row << row;
     for(int column=0; column<NofFiberscolumn; column++){
        std::stringstream C_fiber_column;
        C_fiber_column.str("");
        C_fiber_column << column;
        std::string C_name;
        std::string SiPM_name;
        C_name = "C_row" + C_fiber_row.str() + "_column_" + C_fiber_column.str(); 
        SiPM_name = "SiPMC_row" + C_fiber_row.str() + "_column_" + C_fiber_column.str();

        // I need to specify the position of each cherenkov fiber
        G4double C_x, C_y, C_z;
        G4ThreeVector vec_C_fiber;
        G4ThreeVector vec_SiPM;
        
        if(row == 0 || row == 2 || row == 4 || row == 6){
          if(column == 1 || column == 3 || column == 5 || column == 7){
            C_x = -moduleX/2 + 0.13375 + fiberradius + (0.2675+fiberradius*2)*(column);
            C_y = -moduleY/2 + 0.13375 + fiberradius + (0.2675+fiberradius*2)*(row);
         
            vec_C_fiber.setX(C_x);
            vec_C_fiber.setY(C_y);
            vec_C_fiber.setZ(0.);

            vec_SiPM.setX(C_x);
            vec_SiPM.setY(C_y);
            vec_SiPM.setZ(fiberZ/2+SiPMZ/2-0.18);

            // I need to place the cherenkov fibers
            physi_C_fiber[row][column] = new G4PVPlacement(0,
                                                         vec_C_fiber,      //its position
                                                         logic_C_fiber,    //its logical volume
                                                         C_name,           //its name
                                                         moduleLV,         //its mother
                                                         false,            //no boulean operat
                                                         0);

            // I need to place the SiPMs
            physi_SiPM[row][column] = new G4PVPlacement(0,
                                                        vec_SiPM,            //its position
                                                        SiPMLV,              //its logical volume
                                                        SiPM_name,           //its name
                                                        moduleequippedLV,    //its mother
                                                        false,               //no boulean operat
                                                        0); 
            logic_OpSurface_defaultAir[NofFibersrow][NofFiberscolumn] = new G4LogicalBorderSurface("logic_OpSurface_defaultAir", CalorimeterPV, 
             physi_SiPM[row][column], OpSurfacedefault);
          }
        }
       if(row == 1 || row == 3 || row == 5 || row == 7){
         if(column == 0 || column == 2 || column == 4 || column == 6){
         C_x = -moduleX/2 + 0.13375 + fiberradius + (0.2675+fiberradius*2)*(column);
         C_y = -moduleY/2 + 0.13375 + fiberradius + (0.2675+fiberradius*2)*(row);

         vec_C_fiber.setX(C_x);
         vec_C_fiber.setY(C_y);
         vec_C_fiber.setZ(0.);
          
         vec_SiPM.setX(C_x);
         vec_SiPM.setY(C_y);
         vec_SiPM.setZ(fiberZ/2+SiPMZ/2-0.18);

         // I need to place the cherenkov fibers
         physi_C_fiber[row][column] = new G4PVPlacement(0,
                                                         vec_C_fiber,                   //its position
                                                         logic_C_fiber,     //its logical volume
                                                         C_name,                        //its name
                                                         moduleLV,                      //its mother
                                                         false,                          //no boulean operat
                                                         0); 
         
          // I need to place the SiPMs
          physi_SiPM[row][column] = new G4PVPlacement(0,
                                                      vec_SiPM,                      //its position
                                                      SiPMLV,                        //its logical volume
                                                      SiPM_name,                    //its name
                                                      moduleequippedLV,                      //its mother
                                                      false,                        //no boulean operat
                                                      0); 
          logic_OpSurface_defaultAir[NofFibersrow][NofFiberscolumn] = new G4LogicalBorderSurface("logic_OpSurface_defaultAir", CalorimeterPV, 
           physi_SiPM[row][column], OpSurfacedefault);
         }
       }
     };
  };

  // I return the physical World

  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4DetectorConstruction::ConstructSDandField()
{ 
  // Create global magnetic field messenger,
  // Uniform magnetic field is then created automatically if
  // the field value is not zero
  G4ThreeVector fieldValue = G4ThreeVector();
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);
  
  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

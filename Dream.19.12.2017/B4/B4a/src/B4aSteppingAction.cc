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
// $Id: B4aSteppingAction.cc 68058 2013-03-13 14:47:43Z gcosmo $
// 
/// \file B4aSteppingAction.cc
/// \brief Implementation of the B4aSteppingAction class

#include "B4aSteppingAction.hh"
#include "B4aEventAction.hh"
#include "B4DetectorConstruction.hh"
#include "G4Material.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"

#include "G4OpBoundaryProcess.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4aSteppingAction::B4aSteppingAction(
                      const B4DetectorConstruction* detectorConstruction,
                      B4aEventAction* eventAction)
  : G4UserSteppingAction(),
    fDetConstruction(detectorConstruction),
    fEventAction(eventAction)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4aSteppingAction::~B4aSteppingAction()
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4aSteppingAction::UserSteppingAction(const G4Step* step)
{
  
  // get volume of the current step
  G4VPhysicalVolume* volume 
    = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
      
  /*if ( volume->GetName() == "module" ) {
    //Function to add up energy depoisted in absorber
    fEventAction->Addmodule(step->GetTotalEnergyDeposit());
  }*/

  if ( volume->GetName() != "World" ) {
     //Function to add up energy deposited in the whole calorimeter
     fEventAction->Addenergy(step->GetTotalEnergyDeposit());
  }
  
  //I add up the energy deposited by electrons and positrons
  if (volume->GetName() != "World") {
    if(step->GetTrack()->GetDefinition()->GetParticleName() == "e-" || 
       step->GetTrack()->GetDefinition()->GetParticleName() == "e+"){
      fEventAction->Addem(step->GetTotalEnergyDeposit());
    }
  }

  std::string Fiber;
  std::string S_fiber = "S_fiber";
  std::string C_fiber = "C_fiber";
  Fiber = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName();


  G4StepPoint* pPostStepPoint = step->GetPostStepPoint(); //post step point of each particle
  G4VPhysicalVolume* thePostPV = pPostStepPoint->GetPhysicalVolume(); //volume of post step point
  G4double distance; // will be the distance a photon travels before reaching a SiPM
  G4double p2,p3,p4,ptot; // used as probabilities for parameterization of light
  p2=G4UniformRand(); // random numeber between 0 and 1
  p3=0.4; // SiPM photon detection efficiency 40%
  G4ThreeVector Prestep;
  G4ThreeVector Postsep;
  G4ThreeVector Momentum; // will be the versor of the momentum of each photon inside fibres
  G4double costheta; // will be the angle of emission of each photon inside fibres

  if ( strstr(Fiber.c_str(),S_fiber.c_str())){
    //Function to add up energy depoisted in scintillating fibers
    double k_B = 0.126; //Birk constant
    double edep = step->GetTotalEnergyDeposit();
    double edepattenuated = 0;
    double stepl = step->GetStepLength();
    if(step->GetTrack()->GetDefinition()->GetPDGCharge() != 0.){
        if (stepl != 0)
                {
                    double newedep = (edep/stepl) / ( 1+k_B*(edep/stepl) ) * stepl;
                    edep = newedep;
                }
    }
    else if ( step->GetTrack()->GetDynamicParticle()->GetDefinition()->GetParticleName() == "neutron"
              || step->GetTrack()->GetDynamicParticle()->GetDefinition()->GetParticleName() == "anti_neutron" ) {
                edep = 0.;
            }
    fEventAction->AddScin(edep);
  }

  if ( strstr(Fiber.c_str(),C_fiber.c_str())){
    //Function to add up energy deposited in Cherenkov fibres
    fEventAction->AddCher(step->GetTotalEnergyDeposit());
  }
  
  //Print out the names of particles produced in showers
  //uncomment only if you want this information
  //G4cout<< "particella "<<step->GetTrack()->GetDefinition()->GetParticleName()<<G4endl;

  //Print out the process that produced each photon: Cherenkov or scintillation
  //uncomment only if you want this information
  //if ( step->GetTrack()->GetDefinition()->GetParticleName() == "opticalphoton"){
  //G4cout<< "optical creator process " <<step->GetTrack()->GetCreatorProcess()->GetProcessName()<<G4endl;
  //}

G4OpBoundaryProcessStatus theStatus = Undefined;

G4ProcessManager* OpManager =
                     G4OpticalPhoton::OpticalPhoton()->GetProcessManager();

 if (OpManager) {
     G4int MAXofPostStepLoops =
              OpManager->GetPostStepProcessVector()->entries();
     G4ProcessVector* fPostStepDoItVector =
              OpManager->GetPostStepProcessVector(typeDoIt);

     for ( G4int i=0; i<MAXofPostStepLoops; i++) {
         G4VProcess* fCurrentProcess = (*fPostStepDoItVector)[i];
         fOpProcess = dynamic_cast<G4OpBoundaryProcess*>(fCurrentProcess);
         if (fOpProcess) { theStatus = fOpProcess->GetStatus(); break;}
     }
  }

  std::string SiPMC = "SiPMC";
  std::string SiPMS = "SiPMS";
  std::string SiPMdetection;

  //If the particle is an optical photon...
  if(step->GetTrack()->GetDefinition()->GetParticleName() == "opticalphoton"){

     switch (theStatus){

        case TotalInternalReflection: //it's a photon reflected inside at the core cladding fibre boundary
        // Here starts the parameterization of light, if you don't want it and want to
        // have the complete full simulation with light transportation comment from here to break line!
        // Warning: if you want full simulation with both Cherenkov and scintillation on make sure to have
        // a small light yield in DetectorConstruction.cc otherwise simulation takes forever

           Fiber = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName();//get fibre name
           if(strstr(Fiber.c_str(),S_fiber.c_str())){ //it's a scintillating fibre
               Prestep = step->GetPreStepPoint()->GetPosition(); //get pre step point   
               Postsep = step->GetPostStepPoint()->GetPosition(); //get post step point
               Momentum = step->GetTrack()->GetMomentumDirection(); //get momentum direction of the photon
               if(Momentum.z()>0.){ //if the photon is going towards SiPM 
                 costheta = Momentum.z(); //cosine of the photon respect to fibre axis
                 if(costheta>0.94) { //correspond to 20.4 degree maximum acceptance angle of fibres
                   distance = (560.9-Prestep.z())/costheta; //distance the photon travels before SiPM
                   p4 = std::exp(-(distance/4000)); //exponential decay for light attenuation with attenuation lenght of 5 m
                   ptot =p3; //p4*p3probability of not being asborbed * probability of being detected by SiPM
                   if(p2<ptot){ // it that's the case
                     fEventAction->AddScintillation(); //add one photoelectron from scintillation
                     step->GetTrack()->SetTrackStatus(fStopAndKill); // kill photon in order to not track it, if not simulation becomes too long
                   }
                 }
               }
             }

             if(strstr(Fiber.c_str(),C_fiber.c_str())){ //it's a Cherenkov fibre
               Prestep = step->GetPreStepPoint()->GetPosition();   
               Postsep = step->GetPostStepPoint()->GetPosition();
               Momentum = step->GetTrack()->GetMomentumDirection();
              if(Momentum.z()>0.){
                costheta = Momentum.z();
                if(costheta>0.94){
                  distance = (1560.9-Prestep.z())/costheta;
                  p4 = std::exp(-(distance/8900));
                  ptot =p3;//p4*p3
                  if(p2<ptot){                 
                    fEventAction->AddCherenkov();// add one photoelectron from Cherenkov
                    step->GetTrack()->SetTrackStatus(fStopAndKill);
                  }
                }
              }
             }
    break;

  case Detection:
  // if you want no parameterization and complete full simulation uncomment this part
    
   /* SiPMdetection = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName();
    if (strstr(SiPMdetection.c_str(),SiPMC.c_str()))
     {
       fEventAction->AddCherenkov();
     } 
   
    if (strstr(SiPMdetection.c_str(),SiPMS.c_str()))
    {
      fEventAction->AddScintillation();
    }*/
  
  break;

  default: 
     //only for parameterization, comment for full simulation
     step->GetTrack()->SetTrackStatus(fStopAndKill);
  break;
  }
}

  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// - AUTHOR: Hector Alvarez-Pol 04/2008
/******************************************************************
 * Copyright (C) 2005-2016, Hector Alvarez-Pol                     *
 * All rights reserved.                                            *
 *                                                                 *
 * License according to GNU LESSER GPL (see lgpl-3.0.txt).         *
 * For the list of contributors see CREDITS.                       *
 ******************************************************************/
//////////////////////////////////////////////////////////////////
/// \class ActarSimSciRingDetectorConstruction
/// Scintillator detector description
/////////////////////////////////////////////////////////////////

#include "ActarSimSciRingDetectorConstruction.hh"
#include "ActarSimDetectorConstruction.hh"
#include "ActarSimGasDetectorConstruction.hh"
//#include "ActarSimSciRingDetectorMessenger.hh"
#include "ActarSimROOTAnalysis.hh"
#include "ActarSimSciRingSD.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4Cons.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4RotationMatrix.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4RunManager.hh"
#include "G4Transform3D.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "globals.hh"

//////////////////////////////////////////////////////////////////
/// Constructor. Sets the material and the pointer to the Messenger
ActarSimSciRingDetectorConstruction::
ActarSimSciRingDetectorConstruction(ActarSimDetectorConstruction* det)
  :	sciBulkMaterial(0),detConstruction(det) {
  SetSciBulkMaterial("CsI");

  //Options for Silicon and scintillator coverage:
  // 6 bits to indicate which sci wall is present (1) or absent (0)
  // order is:
  // bit1 (lsb) beam output wall
  // bit2 lower (gravity based) wall
  // bit3 upper (gravity based) wall
  // bit4 left (from beam point of view) wall
  // bit5 right (from beam point of view) wall
  // bit6 (msb) beam entrance wall
  SetSideCoverage(1);

  SetXBoxSciHalfLength(100*cm);
  SetYBoxSciHalfLength(100*cm);
  SetZBoxSciHalfLength(100*cm);

  // create commands for interactive definition of the calorimeter
  // sciMessenger = new ActarSimSciDetectorMessenger(this);
}

//////////////////////////////////////////////////////////////////
/// Destructor
ActarSimSciRingDetectorConstruction::~ActarSimSciRingDetectorConstruction(){
  // delete sciMessenger;
}

//////////////////////////////////////////////////////////////////
/// Wrap for the construction functions
G4VPhysicalVolume* ActarSimSciRingDetectorConstruction::Construct(G4LogicalVolume* worldLog) {
  //Introduce here other constructors for materials around the TOF (windows, frames...)
  //which can be controlled by the calMessenger
  //ConstructTOFWorld(worldLog);
  return ConstructSci(worldLog);
}

//////////////////////////////////////////////////////////////////
/// Real construction work is performed here.
G4VPhysicalVolume* ActarSimSciRingDetectorConstruction::ConstructSci(G4LogicalVolume* worldLog) {
  //Chamber Y,Z length
  //G4double chamberSizeY=detConstruction->GetChamberYLength();
  G4double chamberSizeZ=detConstruction->GetChamberSizeZ();

  //Gas chamber position inside the chamber
  ActarSimGasDetectorConstruction* gasDet = detConstruction->GetGasDetector();
  G4double zGasBoxPosition=gasDet->GetGasBoxCenterZ();

  //----------------------------- the Silicon and CsI disks
  G4double Rmax=48*mm;
  G4double Rmin=24*mm;
  G4double Phi_0=0*deg;
  G4double Phi_f=360*deg;
  G4double ZlengthCsI=10.00*mm;

  G4Tubs *CsIring=new G4Tubs("CsIring",Rmin,Rmax,ZlengthCsI,Phi_0,Phi_f);

  G4LogicalVolume* CsIring_log= new G4LogicalVolume(CsIring,sciBulkMaterial,"CsIring_log",0,0,0);

  G4double sectorPhi=Phi_f/16.;

  G4Tubs *CsISector=new G4Tubs("CsISector",Rmin,Rmax,ZlengthCsI,0.,sectorPhi);

  G4LogicalVolume *CsISector_log=new G4LogicalVolume(CsISector,sciBulkMaterial,"CsISector_log",0,0,0);
  G4VisAttributes* CsISectorVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,1.0));

  CsISectorVisAtt->SetVisibility(true);

  G4double CsIPos_x=0*mm;
  //G4double CsIPos_y=chamberSizeY/2-1*mm;
  G4double CsIPos_y=0;
  G4double CsIPos_z=0.0*mm;

  G4double distance[3]={320*mm,550*mm,1090*mm}; //For 10Li experiment

  //G4double distance[3]={570*mm,0*mm,0*mm}; //For 16C experiment Only one ring detector

  for(G4int k=0;k<3;k++){
    CsIPos_z=distance[k]+chamberSizeZ-zGasBoxPosition;

    G4VPhysicalVolume *CsIring_phys=new G4PVPlacement(0,
						      G4ThreeVector(CsIPos_x,CsIPos_y,CsIPos_z),
						      CsIring_log,"CsIringdet",worldLog,false,k);

    if(CsIring_phys){;}
  }

  G4VPhysicalVolume *CsISector_phys=new G4PVReplica("CsISector",CsISector_log,CsIring_log,kPhi,16,sectorPhi);

  //------------------------------------------------
  // Sensitive detectors
  //------------------------------------------------
  CsISector_log->SetSensitiveDetector( detConstruction->GetSciRingSD() );

  //------------------------------------------------------------------
  // Visualization attributes
  //------------------------------------------------------------------
  return CsISector_phys;
}

//////////////////////////////////////////////////////////////////
/// Set the material the scintillator bulk is made of
void ActarSimSciRingDetectorConstruction::SetSciBulkMaterial (G4String mat) {
  G4Material* pttoMaterial = G4Material::GetMaterial(mat);
  if (pttoMaterial) sciBulkMaterial = pttoMaterial;
}

//////////////////////////////////////////////////////////////////
/// Updates Scintillator detector
void ActarSimSciRingDetectorConstruction::UpdateGeometry() {
  Construct(detConstruction->GetWorldLogicalVolume());
  G4RunManager::GetRunManager()->
    DefineWorldVolume(detConstruction->GetWorldPhysicalVolume());
}

//////////////////////////////////////////////////////////////////
/// Prints Scintillator detector parameters. To be filled
void ActarSimSciRingDetectorConstruction::PrintDetectorParameters() {
  G4cout << "##################################################################"
	 << G4endl
	 << "####  ActarSimSciRingDetectorConstruction::PrintDetectorParameters() ####"
	 << G4endl;
  G4cout << "##################################################################"
	 << G4endl;
}

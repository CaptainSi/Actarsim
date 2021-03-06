// - AUTHOR: Hector Alvarez-Pol 04/2008
/******************************************************************
 * Copyright (C) 2005-2016, Hector Alvarez-Pol                     *
 * All rights reserved.                                            *
 *                                                                 *
 * License according to GNU LESSER GPL (see lgpl-3.0.txt).         *
 * For the list of contributors see CREDITS.                       *
 ******************************************************************/
//////////////////////////////////////////////////////////////////
/// \class ActarSimGasDetectorConstruction
/// Gas volume detector description
/////////////////////////////////////////////////////////////////

#include "ActarSimGasDetectorConstruction.hh"
#include "ActarSimDetectorConstruction.hh"
#include "ActarSimGasDetectorMessenger.hh"
//#include "ActarSimDetectorMessenger.hh"
#include "ActarSimROOTAnalysis.hh"
#include "ActarSimGasSD.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4RunManager.hh"
#include "G4Transform3D.hh"

#include "globals.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//////////////////////////////////////////////////////////////////
/// Constructor
/// Sets the material and the pointer to the Messenger
ActarSimGasDetectorConstruction::
ActarSimGasDetectorConstruction(ActarSimDetectorConstruction* det)
  :	detConstruction(det){

  SetGasPressure(1.01325*bar);
  SetGasTemperature(293.15*kelvin);
  //DefineGas();
  SetBeamShieldMaterial("Iron");
  SetGasMaterial("D2");

  //Default value for the volume is a Box
  SetDetectorGeometry("box");

  //default size of GasBox (0.5x0.5x0.5 m3)
  gasBoxSizeX = 0.5 * m;
  gasBoxSizeY = 0.5 * m;
  gasBoxSizeZ = 0.5 * m;

  //default center of the gas (if the gasChamber is a box)
  gasBoxCenterX = 0.;
  gasBoxCenterY = 0.;
  gasBoxCenterZ = 0.;

  //default size of GasTub (pi x .5 x 0.5 m3 )
  radiusGasTub = 0.5 * m;
  lengthGasTub = 0.5 * m;

  //default size of BeamShieldTub (pi x 0.1 x 0.5 m3 )
  innerRadiusBeamShieldTub = 0.1*m;
  outerRadiusBeamShieldTub = 0.1001*m;
  lengthBeamShieldTub = 0.5 * m;

  // create commands for interactive definition of the calorimeter
  gasMessenger = new ActarSimGasDetectorMessenger(det,this);
}

//////////////////////////////////////////////////////////////////
/// Destructor
ActarSimGasDetectorConstruction::~ActarSimGasDetectorConstruction(){
  delete gasMessenger;
}

//////////////////////////////////////////////////////////////////
/// Wrap for the construction functions within the TOF
G4VPhysicalVolume* ActarSimGasDetectorConstruction::Construct(G4LogicalVolume* chamberLog) {
  return ConstructGas(chamberLog);
}

//////////////////////////////////////////////////////////////////
/// Constructs the Gas volume detector elements
G4VPhysicalVolume* ActarSimGasDetectorConstruction::ConstructGas(G4LogicalVolume* chamberLog) {
  //////////////////////////////////////////////////////////////////////
  //      GAS VOLUME
  // Several geometries are possible. Select the different options using
  // the messenger commands
  //////////////////////////////////////////////////////////////////////

  G4LogicalVolume* gasLog(0);                   //pointer to logic gas
  G4VPhysicalVolume* gasPhys(0);                //pointer to physic gas

  if(detectorGeometry == "box"){
    G4cout << "##################################################################" << G4endl
	   << "######  ActarSimGasDetectorConstruction::ConstructActarTPC()  #######" << G4endl
	   << " Box-like gas geometry." << G4endl;
    G4cout << " Box Parameters: " << G4endl
	   << " gasBoxSizeX = " <<  gasBoxSizeX/mm
	   << " mm,  gasBoxSizeY = " <<  gasBoxSizeY/mm
	   << " mm,  gasBoxSizeZ = " <<  gasBoxSizeZ/mm << " mm" << G4endl
	   << " gasBoxCenterX = " <<  gasBoxCenterX/mm
	   << " mm,  gasBoxCenterY = " <<  gasBoxCenterY/mm
	   << " mm,  gasBoxCenterZ = " <<  gasBoxCenterZ << " mm" << G4endl
	   << " gasMaterial: " <<  gasMaterial << G4endl;
    G4cout << "##################################################################"<< G4endl;


    if(detConstruction->GetACTARTPCGeoIncludedFlag() == "on"){
      //gas Box size: (266*170*266)mm
      gasBoxSizeX = 133.*mm;
      gasBoxSizeY = 85.*mm;
      gasBoxSizeZ = 133.*mm;
      //Pad Size : GasBox height from chamber floor = 4.54mm
      gasBoxCenterX = 0.*mm;
      gasBoxCenterY = -105.0+gasBoxSizeY+4.54*mm; // gasBox shifted to be at the bottom of chamber and above the pads
      gasBoxCenterZ = 0.*mm;
    }
    else if(detConstruction->GetACTARTPCDEMOGeoIncludedFlag() == "on"){
      //gas Box size: (74*170*138)mm
      gasBoxSizeX = 37.*mm;
      gasBoxSizeY = 85.*mm;
      gasBoxSizeZ = 69.*mm;
      //Pad Size : GasBox height from chamber floor = 4.54mm
      gasBoxCenterX = 0.*mm;
      gasBoxCenterY = -105.0+gasBoxSizeY+4.54*mm; // gasBox shifted to be at the bottom of chamber and above the pads
      gasBoxCenterZ = 0.*mm;
    }
    else {
      gasBoxSizeX = GetGasBoxSizeX();
      gasBoxSizeY = GetGasBoxSizeY();
      gasBoxSizeZ = GetGasBoxSizeZ();

      gasBoxCenterX = GetGasBoxCenterX();
      gasBoxCenterY = GetGasBoxCenterY();
      gasBoxCenterZ = GetGasBoxCenterZ();
    }

    G4Box* gasBox;
    gasBox = new G4Box("gasBox",gasBoxSizeX,gasBoxSizeY,gasBoxSizeZ);

    gasLog = new G4LogicalVolume(gasBox,gasMaterial,"gasLog");

    gasPhys = new G4PVPlacement(0,
				G4ThreeVector(gasBoxCenterX,gasBoxCenterY,gasBoxCenterZ),
				gasLog,"gasPhys",chamberLog,false,0);

    // //--------------------------
    // // Field Cage wire replaced by a copper foil around the GasBox
    // //--------------------------
    // G4double wireFoilSizeX = gasBoxSizeX;
    // G4double wireFoilSizeY = gasBoxSizeY;
    // G4double wireFoilSizeZ = 0.00125/2*mm;//That's for 20 um wire with a pitch of 1 mm

    // G4Box *wireFoil=new G4Box("wireFoilBox",wireFoilSizeX,wireFoilSizeY,wireFoilSizeZ);
    // wireFoilLog=new G4LogicalVolume(wireFoil,G4Material::GetMaterial("Copper"),"wireFoilBox");

    // G4double wireFoilPosX = 0.*mm;
    // G4double wireFoilPosY =  -105.0+gasBoxSizeY+4.54*mm;
    // G4double wireFoilPosZ = gasBoxSizeZ+wireFoilSizeZ;

    // wireFoilPhys=new G4PVPlacement(0,G4ThreeVector(wireFoilPosX,wireFoilPosY,wireFoilPosZ),
    // 				   wireFoilLog,"wireFoilBox",chamberLog,false,0);

    // G4VisAttributes* wireFoilVisAtt= new G4VisAttributes(G4Colour(1.0,0.5,0.));
    // wireFoilVisAtt->SetVisibility(true);
    // wireFoilLog->SetVisAttributes(wireFoilVisAtt);
  }
  else if(detectorGeometry == "tube"){
    G4cout << "##################################################################" << G4endl
	   << "########  ActarSimGasDetectorConstruction::ConstructActarTPC()  ########" << G4endl
	   << " Tube-like gas geometry." << G4endl;
    G4cout << " Tube Parameters: " << G4endl
	   << " radiusGasTub = " <<  radiusGasTub/mm
	   << " mm,  lengthGasTub = " <<  lengthGasTub/mm << " mm" << G4endl
	   << " gasMaterial: " <<  gasMaterial << G4endl;
    G4cout << "##################################################################" << G4endl;

    //centered in (0,0,lengthGasTub) to have origin in the detector entrance
    //gasBoxCenterZ = lengthGasTub;
    gasBoxCenterZ = 0.*mm;

    G4Tubs* gasTub;
    gasTub = new G4Tubs("gasTub",0*mm,radiusGasTub,lengthGasTub,0,twopi);

    gasLog = new G4LogicalVolume(gasTub,gasMaterial,"gasLog");

    gasPhys = new G4PVPlacement(0,
				G4ThreeVector(gasBoxCenterX,gasBoxCenterY,gasBoxCenterZ),
				gasLog,"gasPhys",chamberLog,false,0);
  }
  else {
    G4cout << G4endl
	   << " ERROR in ActarSimGasDetectorConstruction::ConstructActarTPC(). No valid volume type defined "
	   << G4endl;
  }

  G4LogicalVolume* beamShieldLog(0);                //pointer to logic
  G4VPhysicalVolume* beamShieldPhys;                //pointer to physic

  //if( beamShieldPhys){;}

  if( beamShieldGeometry == "tube"){
    G4cout << "##################################################################" << G4endl
	   << "########  ActarSimGasDetectorConstruction::ConstructActarTPC()  ########" << G4endl
	   << " Beam shielding geometry." << G4endl;
    G4cout << " Tube Parameters: " << G4endl
	   << " innerRadiusBeamShieldTub = " <<  innerRadiusBeamShieldTub/mm
	   << " mm, outerRadiusBeamShieldTub = " <<  outerRadiusBeamShieldTub/mm << " mm" << G4endl
	   << " lengthBeamShieldTub = " <<  lengthBeamShieldTub/mm
	   << "mm, beamShieldMaterial: " <<  beamShieldMaterial << G4endl;
    G4cout << "##################################################################" << G4endl;

    G4Tubs* beamShieldTub;
    beamShieldTub = new G4Tubs("beamShieldTub",innerRadiusBeamShieldTub,
			       outerRadiusBeamShieldTub,lengthBeamShieldTub,0,twopi);

    beamShieldLog = new G4LogicalVolume(beamShieldTub,beamShieldMaterial,"beamShieldLog");

    beamShieldPhys = new G4PVPlacement(0,
				       G4ThreeVector(gasBoxCenterX,gasBoxCenterY,gasBoxCenterZ),
				       beamShieldLog,"beamShieldPhys",gasLog,false,0);
  }

  //------------------------------------------------
  // Sensitive detectors
  //------------------------------------------------
  gasLog->SetSensitiveDetector( detConstruction->GetGasSD() );

  //------------------------------------------------------------------
  // Visualization attributes
  //------------------------------------------------------------------
  //worldLog->SetVisAttributes (G4VisAttributes::Invisible);
  G4VisAttributes* gasVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
  G4VisAttributes* beamShieldVisAtt = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
  gasVisAtt->SetVisibility(true);
  gasLog->SetVisAttributes(gasVisAtt);
  if( beamShieldGeometry == "tube")  beamShieldLog->SetVisAttributes(beamShieldVisAtt);

  return gasPhys;
}

//////////////////////////////////////////////////////////////////
/// Sets the material the gas is made of
///
/// STP used are P = 1atm and T = 20ºC
void ActarSimGasDetectorConstruction::SetGasMaterial (G4String mat) {
  //Gas Pressure & Temperature
  G4double pressure=GetGasPressure();
  G4double temperature=GetGasTemperature();

  G4double density;
  G4double a;  // atomic mass
  G4double z;  // atomic number
  G4double n;

  //Gas Mix
  //const G4int NGasMix=NumberOfGasMix;

  //HAPOL NOTE: REMOVE REPETITIONS
  G4Element* ele_H  = new G4Element("Hydrogen" ,"H" , z= 1., a=   1.00794*g/mole);
  //G4Element* ele_D  = new G4Element("Deuterium","D" , z= 1., a=    2.0140*g/mole);
  //G4Element* ele_He = new G4Element("Helium"   ,"He", z= 2., a=    4.0026*g/mole);
  G4Element* ele_C  = new G4Element("Carbon"   ,"C",  z=6.,  a=   12.0107*g/mole);
  //G4Element* ele_N  = new G4Element("Nitrogen" ,"N" , z= 7., a=  14.00674*g/mole);
  //G4Element* ele_O  = new G4Element("Oxygen"   ,"O" , z= 8., a=   15.9994*g/mole);
  G4Element* ele_F  = new G4Element("Fluorine" ,"F",  z=9.,  a=18.9984032*g/mole);
  //G4Element* ele_Na = new G4Element("Sodium"   ,"Na", z=11., a=  22.98977*g/mole);
  //G4Element* ele_S  = new G4Element("Sulphur"  ,"S",  z=16., a=    32.066*g/mole);
  //G4Element* ele_Ar = new G4Element("Argon"    ,"Ar", z=18., a=   39.9481*g/mole);
  /*
  G4Element* ele_Zn = new G4Element("Zinc",     "Zn", z=30., a=     65.39*g/mole);
  G4Element* ele_Ge = new G4Element("Germanium","Ge", z=32., a=     72.61*g/mole);
  G4Element* ele_Br = new G4Element("Bromine"  ,"Br", z=35., a=    79.904*g/mole);
  G4Element* ele_Cd = new G4Element("Cadmium"  ,"Cd", z=48., a=   112.411*g/mole);
  G4Element* ele_Te = new G4Element("Tellurium","Te", z=52., a=    127.60*g/mole);
  G4Element* ele_I  = new G4Element("Iodine"   ,"I",  z=53., a= 126.90447*g/mole);
  G4Element* ele_Cs = new G4Element("Cesium"   ,"Cs", z=55., a= 132.90545*g/mole);
  G4Element* ele_Ba = new G4Element("Barium"   ,"Ba", z=56., a=   137.327*g/mole);
  G4Element* ele_La = new G4Element("Lanthanum","La", z=57., a=  138.9055*g/mole);
  G4Element* ele_Ce = new G4Element("Cerium"   ,"Ce", z=58., a=   140.116*g/mole);
  G4Element* ele_Lu = new G4Element("Lutecium" ,"Lu", z=71., a=   174.967*g/mole);
  G4Element* ele_W  = new G4Element("Tungsten" ,"W" , z=74., a=    183.84*g/mole);
  G4Element* ele_Pb = new G4Element("Lead"     ,"Pb", z=82., a=    207.20*g/mole);
  G4Element* ele_Bi = new G4Element("Bismuth"  ,"Bi", z=83., a= 208.98038*g/mole);
  */

  G4int ncomponents, natoms;
  G4double fractionmass, abundance;

  G4double Vm=0.08206*temperature*atmosphere/(pressure*kelvin);

  G4Isotope* iso_H2= new G4Isotope("iso_H2",z=1,n=2, a=2.0140*g/mole);
  G4Element* ele_D= new G4Element("Deuterium","D" , ncomponents=1);
  ele_D->AddIsotope(iso_H2, abundance = 100.*perCent);

  //material definition by user's T and P

  //H2 (default  0.083812*mg/cm3 STP)
  //density	= (0.083812*293.15*kelvin*pressure)/(1.01325*bar*temperature)*mg/cm3;
  density	= (2*1.00794/Vm)*mg/cm3;
  G4Material* H2 =
    new G4Material("H2", density, ncomponents=2, kStateGas, temperature, pressure);
  H2->AddElement(ele_H, natoms=1);
  H2->AddElement(ele_H, natoms=1);

  //D2 (default  0.16746*mg/cm3 STP)
  //density	=(0.16746*293.15*kelvin*pressure)/(1.01325*bar*temperature)*mg/cm3;
  density	= (2*2.0140/Vm)*mg/cm3;
  G4Material* D2 =
    new G4Material("D2", density, ncomponents=2, kStateGas, temperature, pressure);
  D2->AddElement(ele_D, natoms=1);
  D2->AddElement(ele_D, natoms=1);

  //He (default  0.16642*mg/cm3 STP)
  //density	=(0.16642*293.15*kelvin*pressure)/(1.01325*bar*temperature)*mg/cm3;
  density	= (4.0026/Vm)*mg/cm3;
  G4Material* He =
    new G4Material("He", z=2, a=4.0026*g/mole, density, kStateGas, temperature, pressure);

  //Ar (default  0.16642*mg/cm3 STP)
  //density	=(0.16642*293.15*kelvin*pressure)/(1.01325*bar*temperature)*mg/cm3;
  density	= (39.9481/Vm)*mg/cm3;
  G4Material* Ar =
    new G4Material("Ar", z=2, a=39.9481*g/mole, density, kStateGas, temperature, pressure);

  //CF4 (default  3.6586*mg/cm3 STP)
  //density	=(3.6586*293.15*kelvin*pressure)/(1.01325*bar*temperature)*mg/cm3;
  density	= ((12.0107+4*18.9984032)/Vm)*mg/cm3;
  G4Material* CF4 =
    new G4Material("CF4", density, ncomponents=2, kStateGas, temperature, pressure);
  CF4->AddElement(ele_C, natoms=1);
  CF4->AddElement(ele_F, natoms=4);

  //Methane (default  0.66697*mg/cm3 STP)
  //density = (0.6669*293.15*kelvin*pressure)/(1.01325*bar*temperature)*mg/cm3;
  density	= ((12.0107+4*1.00794)/Vm)*mg/cm3;
  G4Material* CH4 =
    new G4Material("CH4", density, ncomponents=2, kStateGas, temperature, pressure) ;
  CH4->AddElement(ele_C,1);
  CH4->AddElement(ele_H,4);

  //Isobutane (default  2.41464*mg/cm3 STP)
  //density = (2.41464*293.15*kelvin*pressure)/(1.01325*bar*temperature)*mg/cm3;
  density	= ((4*12.0107+10*1.00794)/Vm)*mg/cm3;
  G4Material* iC4H10 =
    new G4Material("iC4H10", density, ncomponents=2, kStateGas, temperature, pressure) ;
  iC4H10->AddElement(ele_C,4);
  iC4H10->AddElement(ele_H,10);

  if(mat=="H2"){
    gasMaterial = H2;
    detConstruction->SetUpdateChamberMaterial(H2);
  }
  else if(mat=="D2"){
    gasMaterial = D2;
    detConstruction->SetUpdateChamberMaterial(D2);
  }
  else if(mat=="He"){
    gasMaterial = He;
    detConstruction->SetUpdateChamberMaterial(He);
  }
  else if(mat=="Ar"){
    gasMaterial = Ar;
    detConstruction->SetUpdateChamberMaterial(Ar);
  }
  else if(mat=="CF4"){
    gasMaterial = CF4;
    detConstruction->SetUpdateChamberMaterial(CF4);
  }
  else if(mat=="CH4"){
    gasMaterial = CH4;
    detConstruction->SetUpdateChamberMaterial(CH4);
  }
  else if(mat=="iC4H10"){
    gasMaterial = iC4H10;
    detConstruction->SetUpdateChamberMaterial(iC4H10);
  }
  else if(mat=="GasMix"){

    density = 0*mg/cm3;
    G4double DensitySum=0;
    //G4double FractionMass[NGasMix];
    G4double FractionMass[10];
    //G4Material pttoMaterial[NGasMix];
    G4Material *pttoMaterial[10];

    for(G4int i=0;i<NumberOfGasMix;i++) {
      //pttoMaterial[i] = G4Material::GetMaterial(GasMat[i]);
      if(gasMixMaterial[i]=="H2")pttoMaterial[i]=H2;
      else if(gasMixMaterial[i]=="D2")pttoMaterial[i]=D2;
      else if(gasMixMaterial[i]=="He")pttoMaterial[i]=He;
      else if(gasMixMaterial[i]=="Ar")pttoMaterial[i]=Ar;
      else if(gasMixMaterial[i]=="CF4")pttoMaterial[i]=CF4;
      else if(gasMixMaterial[i]=="CH4")pttoMaterial[i]=CH4;
      else if(gasMixMaterial[i]=="iC4H10")pttoMaterial[i]=iC4H10;

      density+= ((gasMixRatio[i]*pttoMaterial[i]->GetDensity()));
      //G4cout <<" Gas Mat "<<gasMixMaterial[i]<<" Gas Ratio "<<gasMixRatio[i]<<" Mat Density "<<gasMixRatio[i]*pttoMaterial[i]->GetDensity()*cm3/mg<< G4endl;
      DensitySum+=pttoMaterial[i]->GetDensity();
    }

    for(G4int i=0;i<NumberOfGasMix;i++) {
      FractionMass[i]=pttoMaterial[i]->GetDensity()/DensitySum;
    }

    G4Material* GasMix =
      new G4Material("GasMix", density, ncomponents=NumberOfGasMix, kStateGas, temperature, pressure);

    for(G4int i=0;i<NumberOfGasMix;i++) {
      GasMix->AddMaterial( pttoMaterial[i], fractionmass = FractionMass[i] ) ;
    }

    gasMaterial = GasMix;
    detConstruction->SetUpdateChamberMaterial(GasMix);
  }
}

//////////////////////////////////////////////////////////////////
///  Sets the material the medium is made of
void ActarSimGasDetectorConstruction::SetBeamShieldMaterial(G4String mat) {
  G4Material* pttoMaterial = G4Material::GetMaterial(mat);
  if (pttoMaterial) beamShieldMaterial = pttoMaterial;
}

//////////////////////////////////////////////////////////////////
///  Sets the geometry of the detector (box or tube)
void ActarSimGasDetectorConstruction::SetDetectorGeometry(G4String type) {
  detectorGeometry = type;
}

//////////////////////////////////////////////////////////////////
///  Sets the geometry of the detector (box or tube)
void ActarSimGasDetectorConstruction::SetBeamShieldGeometry(G4String type) {
  beamShieldGeometry = type;
}

//////////////////////////////////////////////////////////////////
///  Updates Gas detector
void ActarSimGasDetectorConstruction::UpdateGeometry() {
  // Construct(detConstruction->GetWorldLogicalVolume());
  // G4RunManager::GetRunManager()->
  //   DefineWorldVolume(detConstruction->GetWorldPhysicalVolume());
  Construct(detConstruction->GetChamberLogicalVolume());
  G4RunManager::GetRunManager()->
    DefineWorldVolume(detConstruction->GetChamberPhysicalVolume());
}

//////////////////////////////////////////////////////////////////
///  Prints Gas volume detector parameters
void ActarSimGasDetectorConstruction::PrintDetectorParameters() {
  G4cout << "##################################################################" << G4endl
	 << "##  ActarSimGasDetectorConstruction::PrintDetectorParameters() ###" << G4endl
	 << " The gas volume is a " ;
  if(detectorGeometry == "box")
    G4cout << "box; its parameters are:" << G4endl;
  if(detectorGeometry == "tube")
    G4cout << "tube; its parameters are:" << G4endl;
  G4cout << " The gas material is: " << gasMaterial  << G4endl;
  if(detectorGeometry == "box")
    G4cout << " The gasBox size is : " << gasBoxSizeX/mm << "x" << gasBoxSizeY/mm
	   << "x" << gasBoxSizeZ/mm << " mm3 " << G4endl << G4endl ;
  if(detectorGeometry == "tube")
    G4cout << " The gasTube parameters are: " << G4endl
	   << " radiusGasTub = " <<  radiusGasTub/mm
	   << "mm,  lengthGasTub = " <<  lengthGasTub/mm << " mm" << G4endl ;
  if( beamShieldGeometry == "tube"){
    G4cout << " The beam shielding parameters are:"  << G4endl
	   << " innerRadiusBeamShieldTub = " <<  innerRadiusBeamShieldTub/mm
	   << " mm, outerRadiusBeamShieldTub = " <<  outerRadiusBeamShieldTub/mm << " mm" << G4endl
	   << " lengthBeamShieldTub = " <<  lengthBeamShieldTub/mm
	   << " mm, beamShieldMaterial: " <<  beamShieldMaterial << G4endl;
  }
  G4cout << "##################################################################"
	 << G4endl;
}

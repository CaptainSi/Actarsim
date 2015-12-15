 ///////////////////////////////////////////////////////////////
//*-- AUTHOR : Hector Alvarez Pol
//*-- Date: 04/2008
//*-- Last Update: 15/12/14
// --------------------------------------------------------------
// Description:
//   Messenger of the Gas ACTAR detector
//
// --------------------------------------------------------------
// Comments:
//
//
// --------------------------------------------------------------
/////////////////////////////////////////////////////////////////

#include "ActarSimGasDetectorMessenger.hh"

#include "G4RunManager.hh"

#include "ActarSimDetectorConstruction.hh"
#include "ActarSimGasDetectorConstruction.hh"
#include "ActarSimPrimaryGeneratorAction.hh"

#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4Tokenizer.hh"


ActarSimGasDetectorMessenger::
ActarSimGasDetectorMessenger(ActarSimDetectorConstruction* ActarSimDet,ActarSimGasDetectorConstruction* ActarSimGasDet)
  :ActarSimDetector(ActarSimDet), ActarSimGasDetector(ActarSimGasDet){
  //
  // Constructor with fully functionality
  //

  detDir = new G4UIdirectory("/ActarSim/det/gas/");
  detDir->SetGuidance("gas detector control");

  gasMaterCmd = new G4UIcmdWithAString("/ActarSim/det/gas/setGasMat",this);
  gasMaterCmd->SetGuidance("Select Material of the Gas (for the Gas box and the Chamber).");
  gasMaterCmd->SetParameterName("gasMat",false);
  gasMaterCmd->SetDefaultValue("D2");
  gasMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  gasPresCmd = new G4UIcmdWithADoubleAndUnit("/ActarSim/det/gas/setGasPressure",this);
  gasPresCmd->SetGuidance("Select the Gas Pressure (for the Gas box and the Chamber).");
  gasPresCmd->SetParameterName("gasPressure",false);
  gasPresCmd->SetRange("gasPressure>=0.");
  gasPresCmd->SetUnitCategory("Pressure");
  gasPresCmd->SetDefaultUnit("bar");
  gasPresCmd->SetDefaultValue(1.01325);
  gasPresCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  gasTempCmd = new G4UIcmdWithADoubleAndUnit("/ActarSim/det/gas/setGasTemperature",this);
  gasTempCmd->SetGuidance("Select the Gas Temperature (for the Gas box and the Chamber).");
  gasTempCmd->SetParameterName("gasTemperature",false);
  gasTempCmd->SetRange("gasTemperature>=0.");
  gasTempCmd->SetUnitCategory("Temperature");
  gasTempCmd->SetDefaultUnit("kelvin");
  gasTempCmd->SetDefaultValue(293.15);
  gasTempCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  beamShieldMaterCmd = new G4UIcmdWithAString("/ActarSim/det/gas/setBeamShieldMat",this);
  beamShieldMaterCmd->SetGuidance("Select Material of the beam shield.");
  beamShieldMaterCmd->SetParameterName("gasMat",false);
  beamShieldMaterCmd->SetDefaultValue("iC4H10");
  beamShieldMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  detectorGeometryCmd = new G4UIcmdWithAString("/ActarSim/det/gas/setDetectorGeometry",this);
  detectorGeometryCmd->SetGuidance("Select the geometry of the detector.");
  detectorGeometryCmd->SetGuidance("Choice : box(default), tube.");
  detectorGeometryCmd->SetParameterName("choice",false);
  detectorGeometryCmd->SetDefaultValue("box");
  detectorGeometryCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  setBeamShieldCmd = new G4UIcmdWithAString("/ActarSim/det/gas/setBeamShield",this);
  setBeamShieldCmd->SetGuidance("Sets a beam shield and selects the geometry.");
  setBeamShieldCmd->SetGuidance("Choice : tube, off.");
  setBeamShieldCmd->SetParameterName("choice",false);
  setBeamShieldCmd->SetDefaultValue("tube");
  setBeamShieldCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  gasBoxSizeXCmd = new G4UIcmdWithADoubleAndUnit("/ActarSim/det/gas/setXLengthGasBox",this);
  gasBoxSizeXCmd->SetGuidance("Select the half-length X dimension of the Gas Box.");
  gasBoxSizeXCmd->SetParameterName("gasBoxSizeX",false);
  gasBoxSizeXCmd->SetRange("gasBoxSizeX>=0.");
  gasBoxSizeXCmd->SetUnitCategory("Length");
  gasBoxSizeXCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  gasBoxSizeYCmd = new G4UIcmdWithADoubleAndUnit("/ActarSim/det/gas/setYLengthGasBox",this);
  gasBoxSizeYCmd->SetGuidance("Select the half-length Y dimension of the Gas Box.");
  gasBoxSizeYCmd->SetParameterName("gasBoxSizeY",false);
  gasBoxSizeYCmd->SetRange("gasBoxSizeY>=0.");
  gasBoxSizeYCmd->SetUnitCategory("Length");
  gasBoxSizeYCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  gasBoxSizeZCmd = new G4UIcmdWithADoubleAndUnit("/ActarSim/det/gas/setZLengthGasBox",this);
  gasBoxSizeZCmd->SetGuidance("Select the half-length Z dimension of the Gas Box.");
  gasBoxSizeZCmd->SetParameterName("gasBoxSizeZ",false);
  gasBoxSizeZCmd->SetRange("gasBoxSizeZ>=0.");
  gasBoxSizeZCmd->SetUnitCategory("Length");
  gasBoxSizeZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  radiusGasTubCmd = new G4UIcmdWithADoubleAndUnit("/ActarSim/det/gas/setRadiusGasTub",this);
  radiusGasTubCmd->SetGuidance("Select the external radius of the Gas Tube.");
  radiusGasTubCmd->SetParameterName("radiusGasTub",false);
  radiusGasTubCmd->SetRange("radiusGasTub>=0.");
  radiusGasTubCmd->SetUnitCategory("Length");
  radiusGasTubCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  lengthGasTubCmd = new G4UIcmdWithADoubleAndUnit("/ActarSim/det/gas/setLengthGasTub",this);
  lengthGasTubCmd->SetGuidance("Select the half-length of the Gas Tube.");
  lengthGasTubCmd->SetParameterName("lengthGasTub",false);
  lengthGasTubCmd->SetRange("lengthGasTub>=0.");
  lengthGasTubCmd->SetUnitCategory("Length");
  lengthGasTubCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  innerRadiusBeamShieldTubCmd = new G4UIcmdWithADoubleAndUnit("/ActarSim/det/gas/setInnerRadiusBeamShieldTub",this);
  innerRadiusBeamShieldTubCmd->SetGuidance("Select the external radius of the Gas Tube.");
  innerRadiusBeamShieldTubCmd->SetParameterName("innerRadiusBeamShieldTub",false);
  innerRadiusBeamShieldTubCmd->SetRange("innerRadiusBeamShieldTub>=0.");
  innerRadiusBeamShieldTubCmd->SetUnitCategory("Length");
  innerRadiusBeamShieldTubCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  outerRadiusBeamShieldTubCmd = new G4UIcmdWithADoubleAndUnit("/ActarSim/det/gas/setRadiusBeamShieldTub",this);
  outerRadiusBeamShieldTubCmd->SetGuidance("Select the internal radius of the Gas Tube.");
  outerRadiusBeamShieldTubCmd->SetParameterName("outerRadiusBeamShieldTub",false);
  outerRadiusBeamShieldTubCmd->SetRange("outerRadiusBeamShieldTub>=0.");
  outerRadiusBeamShieldTubCmd->SetUnitCategory("Length");
  outerRadiusBeamShieldTubCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  lengthBeamShieldTubCmd = new G4UIcmdWithADoubleAndUnit("/ActarSim/det/gas/setLengthBeamShieldTub",this);
  lengthBeamShieldTubCmd->SetGuidance("Select the half-length of the Gas Tube.");
  lengthBeamShieldTubCmd->SetParameterName("lengthBeamShieldTub",false);
  lengthBeamShieldTubCmd->SetRange("lengthBeamShieldTub>=0.");
  lengthBeamShieldTubCmd->SetUnitCategory("Length");
  lengthBeamShieldTubCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  printCmd = new G4UIcmdWithoutParameter("/ActarSim/det/gas/print",this);
  printCmd->SetGuidance("Prints geometry.");
  printCmd->AvailableForStates(G4State_Idle);

  detDirMix = new G4UIdirectory("/ActarSim/det/gas/mixture/");
  detDirMix->SetGuidance("gas mixture control");

  gasMixtureFlagCmd = new G4UIcmdWithAnInteger("/ActarSim/det/gas/mixture/GasMixture",this);
  gasMixtureFlagCmd->SetGuidance("Set a Gas Mixture with N number of elements (default 0).");
  gasMixtureFlagCmd->SetGuidance("  Choice : 0 (default) to 10");
  //gasMixtureFlagCmd->SetParameterName("GasNumber",false);
  gasMixtureFlagCmd->SetParameterName("N",true,true);
  gasMixtureFlagCmd->SetRange("N>0");
  gasMixtureFlagCmd->SetDefaultValue(2);
  //gasMixtureFlagCmd->SetCandidates("2 3 4 5 6 7 8 9 10");
  //gasMixtureFlagCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  gasMixtureCmd = new G4UIcommand("/ActarSim/det/gas/mixture/setGasMix",this);
  gasMixtureCmd->SetGuidance("Set a Gas Mixture (for the Gas box and the Chamber).");
  gasMixtureCmd->SetGuidance("[usage] /ActarSim/det/gas/setGasMix GasNum GasMat GasRatio");
  gasMixtureCmd->SetGuidance("        GasNum:(int) GasNumber (from 1 to 7)");
  gasMixtureCmd->SetGuidance("        GasMat:(string) Gas Material from the list");
  gasMixtureCmd->SetGuidance("        GasRatio:(double) Gas Ratio in Mixture (from 0 to 1)");

  G4UIparameter* gasMixtureParam;
  gasMixtureParam = new G4UIparameter("GasNum",'i',false);
  gasMixtureParam->SetDefaultValue("1");
  gasMixtureParam->SetParameterRange("GasNum>0");
  gasMixtureParam->SetParameterRange("GasNum<10");
  gasMixtureCmd->SetParameter(gasMixtureParam);
  gasMixtureParam = new G4UIparameter("GasMat",'s',false);
  gasMixtureParam->SetDefaultValue("D2");
  gasMixtureParam->SetParameterCandidates("H2 D2 He Ar CF4 iC4H10");
  gasMixtureCmd->SetParameter(gasMixtureParam);
  gasMixtureParam = new G4UIparameter("GasRatio",'d',false);
  gasMixtureParam->SetDefaultValue("0");
  gasMixtureParam->SetParameterRange("GasRatio>=0.");
  gasMixtureParam->SetParameterRange("GasRatio<=1.");
  gasMixtureCmd->SetParameter(gasMixtureParam);
}


ActarSimGasDetectorMessenger::~ActarSimGasDetectorMessenger() {
  //
  // Destructor
  //
  delete detDir;
  delete detDirMix;
  delete gasMaterCmd;
  delete gasMixtureFlagCmd;
  delete gasMixtureCmd;
  delete gasPresCmd;
  delete gasTempCmd;
  delete beamShieldMaterCmd;
  delete detectorGeometryCmd;
  delete setBeamShieldCmd;
  delete gasBoxSizeXCmd;
  delete gasBoxSizeYCmd;
  delete gasBoxSizeZCmd;
  delete radiusGasTubCmd;
  delete lengthGasTubCmd;
  delete innerRadiusBeamShieldTubCmd;
  delete outerRadiusBeamShieldTubCmd;
  delete lengthBeamShieldTubCmd;
  delete printCmd;
}


void ActarSimGasDetectorMessenger::SetNewValue(G4UIcommand* command,
					    G4String newValue) {
  //
  // Setting the new values and connecting to actions
  //

  if(command == gasMaterCmd)
    {
      //ActarSimGasDetector->DefineGas();
      ActarSimGasDetector->SetGasMaterial(newValue);
      //ActarSimDetector->SetChamberMaterial(newValue);
      ActarSimDetector->UpdateGeometry();
      //ActarSimDetector->PrintDetectorParameters();
    }

  if(command == gasMixtureFlagCmd)
    {
      ActarSimGasDetector->SetGasMixture(gasMixtureFlagCmd->GetNewIntValue(newValue));
    }

  if(command == gasMixtureCmd)
    {
      GasMixtureCommand(newValue);
    }

  if(command == gasTempCmd)
  {
    ActarSimGasDetector->SetGasTemperature(gasTempCmd->GetNewDoubleValue(newValue));
  }

  if(command == beamShieldMaterCmd)
  {
    ActarSimGasDetector->SetBeamShieldMaterial(newValue);
  }

  if(command == detectorGeometryCmd)
  {
    ActarSimGasDetector->SetDetectorGeometry(newValue);
  }

  if(command == setBeamShieldCmd)
  {
    ActarSimGasDetector->SetBeamShieldGeometry(newValue);
  }

  if(command == gasBoxSizeXCmd)
  {
    ActarSimGasDetector->SetGasBoxSizeX(gasBoxSizeXCmd->GetNewDoubleValue(newValue));
  }

  if(command == gasBoxSizeYCmd)  {
    //HAPOL NOTE CHECK THIS DOUBLE COMMAND!!!!
    ActarSimGasDetector->SetGasBoxSizeY(gasBoxSizeYCmd->GetNewDoubleValue(newValue));
    ActarSimGasDetector->SetGasBoxCenterY(gasBoxSizeYCmd->GetNewDoubleValue(newValue));
  }

  if(command == gasBoxSizeZCmd)  {
    //HAPOL NOTE CHECK THIS DOUBLE COMMAND!!!!
    ActarSimGasDetector->SetGasBoxSizeZ(gasBoxSizeZCmd->GetNewDoubleValue(newValue));
    ActarSimGasDetector->SetGasBoxCenterZ(gasBoxSizeZCmd->GetNewDoubleValue(newValue));
  }

  if(command == radiusGasTubCmd)
  {
    ActarSimGasDetector->SetRadiusGasTub(radiusGasTubCmd->GetNewDoubleValue(newValue));
  }

  if(command == lengthGasTubCmd)
  {
    ActarSimGasDetector->SetLengthGasTub(lengthGasTubCmd->GetNewDoubleValue(newValue));
  }

  if(command == innerRadiusBeamShieldTubCmd)
  {
    ActarSimGasDetector->SetInnerRadiusBeamShieldTub(innerRadiusBeamShieldTubCmd->GetNewDoubleValue(newValue));
  }

  if(command == outerRadiusBeamShieldTubCmd)
  {
    ActarSimGasDetector->SetOuterRadiusBeamShieldTub(outerRadiusBeamShieldTubCmd->GetNewDoubleValue(newValue));
  }

  if(command == lengthBeamShieldTubCmd)
  {
    ActarSimGasDetector->SetLengthBeamShieldTub(lengthBeamShieldTubCmd->GetNewDoubleValue(newValue));
  }


  if( command == printCmd )
    { ActarSimGasDetector->PrintDetectorParameters(); }

}

void ActarSimGasDetectorMessenger::GasMixtureCommand(G4String newValues){
  //
  // Selection of the mixture gases.
  //

  G4Tokenizer next( newValues );
  // check argument
  fGasNumber = StoI(next());
  fGasMaterial = next();
  fGasRatio = StoD(next());

  //G4cout << " fGasNumber "<<fGasNumber<< " fGasMaterial "<<fGasMaterial<< " fGasRatio "<<fGasRatio<< G4endl;

  ActarSimGasDetector->SetGasMixMaterial(fGasNumber,fGasMaterial,fGasRatio);

}

// - AUTHOR: Hector Alvarez-Pol 11/2004
/******************************************************************
 * Copyright (C) 2005-2016, Hector Alvarez-Pol                     *
 * All rights reserved.                                            *
 *                                                                 *
 * License according to GNU LESSER GPL (see lgpl-3.0.txt).         *
 * For the list of contributors see CREDITS.                       *
 ******************************************************************/
//////////////////////////////////////////////////////////////////
/// \class ActarSimPrimaryGeneratorMessenger
/// Messenger for the primary event generator.
/////////////////////////////////////////////////////////////////

#include "ActarSimPrimaryGeneratorMessenger.hh"

#include "ActarSimPrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4Ions.hh"
#include "G4ios.hh"
#include "G4Tokenizer.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//////////////////////////////////////////////////////////////////
/// Constructor
/// command included in this ActarSimPrimaryGeneratorMessenger:
/// - /ActarSim/gun/List
/// - /ActarSim/gun/particle
/// - /ActarSim/gun/realisticBeam
/// - /ActarSim/gun/beamInteraction
/// - /ActarSim/gun/emittance
/// - /ActarSim/gun/beamDirection
/// - /ActarSim/gun/beamTheta
/// - /ActarSim/gun/beamPhi
/// - /ActarSim/gun/beamPosition
/// - /ActarSim/gun/beamRadiusAtEntrance
/// - /ActarSim/gun/reactionFromEvGen
/// - /ActarSim/gun/reactionFromFile
/// - /ActarSim/gun/reactionFromCrossSection
/// - /ActarSim/gun/reactionFile
/// - /ActarSim/gun/reactionFromCine
/// - /ActarSim/gun/Cine/randomTheta
/// - /ActarSim/gun/randomTheta
/// - /ActarSim/gun/randomPhi
/// - /ActarSim/gun/alphaSource
/// - /ActarSim/gun/randomThetaVal
/// - /ActarSim/gun/randomPhiVal
/// - /ActarSim/gun/Cine/randomThetaVal
/// - /ActarSim/gun/Cine/incidentIon
/// - /ActarSim/gun/Cine/targetIon
/// - /ActarSim/gun/Cine/scatteredIon
/// - /ActarSim/gun/Cine/recoilIon
/// - /ActarSim/gun/Cine/reactionQ
/// - /ActarSim/gun/Cine/labEnergy
/// - /ActarSim/gun/Cine/thetaLabAngle
/// - /ActarSim/gun/reactionFromKine
/// - /ActarSim/gun/Kine/randomThetaCM
/// - /ActarSim/gun/Kine/randomPhiAngle
/// - /ActarSim/gun/Kine/randomThetaRange
/// - /ActarSim/gun/Kine/incidentIon
/// - /ActarSim/gun/Kine/targetIon
/// - /ActarSim/gun/Kine/scatteredIon
/// - /ActarSim/gun/Kine/recoilIon
/// - /ActarSim/gun/Kine/labEnergy
/// - /ActarSim/gun/Kine/userThetaCM
/// - /ActarSim/gun/Kine/userPhiAngle
/// - /ActarSim/gun/Kine/vertexPosition
/// - /ActarSim/gun/energy
/// - /ActarSim/gun/direction
/// - /ActarSim/gun/position
/// - /ActarSim/gun/time
/// - /ActarSim/gun/randomVertexZPosition
/// - /ActarSim/gun/randomVertexZRange
/// - /ActarSim/gun/vertexZPosition
/// - /ActarSim/gun/polarization
/// - /ActarSim/gun/number
/// - /ActarSim/gun/ion
ActarSimPrimaryGeneratorMessenger::ActarSimPrimaryGeneratorMessenger(ActarSimPrimaryGeneratorAction* actarSimGun)
  : actarSimActionGun(actarSimGun) {

  particleTable = G4ParticleTable::GetParticleTable();
  ionTable = G4IonTable::GetIonTable();

  G4bool omitable;
  G4UIparameter* parameter;

  gunDir = new G4UIdirectory("/ActarSim/gun/");
  gunDir->SetGuidance("PrimaryGenerator control");

  listCmd = new G4UIcmdWithoutParameter("/ActarSim/gun/List",this);
  listCmd->SetGuidance("List available particles.");
  listCmd->SetGuidance(" Invoke G4ParticleTable.");

  particleCmd = new G4UIcmdWithAString("/ActarSim/gun/particle",this);
  particleCmd->SetGuidance("Select the incident particle.");
  particleCmd->SetGuidance(" (proton is default)");
  particleCmd->SetGuidance(" (ion can be specified for shooting ions)");
  particleCmd->SetParameterName("particle",false);
  particleCmd->SetDefaultValue("proton");
  G4String candidateList;
  G4int nPtcl = particleTable->entries();
  for(G4int i=0;i<nPtcl;i++) {
    if(!(particleTable->GetParticle(i)->IsShortLived())) {
      candidateList += particleTable->GetParticleName(i);
      candidateList += " ";
    }
  }
  candidateList += "ion ";
  particleCmd->SetCandidates(candidateList);

  realisticBeamCmd = new G4UIcmdWithAString("/ActarSim/gun/realisticBeam",this);
  realisticBeamCmd->SetGuidance("Simulates beam emittance according to emittance parameters.");
  realisticBeamCmd->SetGuidance("  Choice : on, off(default)");
  realisticBeamCmd->SetParameterName("choice",true);
  realisticBeamCmd->SetDefaultValue("off");
  realisticBeamCmd->SetCandidates("on off");
  realisticBeamCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  beamInteractionCmd = new G4UIcmdWithAString("/ActarSim/gun/beamInteraction",this);
  beamInteractionCmd->SetGuidance("Simulates the beam energy loss in gas.");
  beamInteractionCmd->SetGuidance("  Choice : on, off(default)");
  beamInteractionCmd->SetParameterName("choice",true);
  beamInteractionCmd->SetDefaultValue("off");
  beamInteractionCmd->SetCandidates("on off");
  beamInteractionCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  emittanceCmd = new G4UIcmdWithADouble("/ActarSim/gun/emittance",this);
  emittanceCmd->SetGuidance("Selects the value of the emittance [in mm mrad].");
  emittanceCmd->SetGuidance(" Default value is 1 mm mrad. ");
  emittanceCmd->SetParameterName("emittance",false);
  emittanceCmd->SetDefaultValue(1.);
  emittanceCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  beamDirectionCmd = new G4UIcmdWith3Vector("/ActarSim/gun/beamDirection",this);
  beamDirectionCmd->SetGuidance("Set beam momentum direction.");
  beamDirectionCmd->SetGuidance("Direction needs not to be a unit vector.");
  beamDirectionCmd->SetParameterName("Px","Py","Pz",true,true);
  beamDirectionCmd->SetRange("Px != 0 || Py != 0 || Pz != 0");

  beamThetaCmd = new G4UIcmdWithADoubleAndUnit("/ActarSim/gun/beamTheta",this);
  beamThetaCmd->SetGuidance("Sets theta angle for beam (in degrees)");
  beamThetaCmd->SetParameterName("beamTheta",false);
  //beamThetaCmd->SetRange("userThetaCM>=0.");
  beamThetaCmd->SetUnitCategory("Angle");
  beamThetaCmd->SetDefaultValue(0.);
  beamThetaCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  beamPhiCmd = new G4UIcmdWithADoubleAndUnit("/ActarSim/gun/beamPhi",this);
  beamPhiCmd->SetGuidance("Sets phi angle for beam (in degrees)");
  beamPhiCmd->SetParameterName("beamPhi",false);
  //beamPhiCmd->SetRange("userPhiCM>=0.");
  beamPhiCmd->SetUnitCategory("Angle");
  beamPhiCmd->SetDefaultValue(0.);
  beamPhiCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  beamPositionCmd = new G4UIcmdWith3VectorAndUnit("/ActarSim/gun/beamPosition",this);
  beamPositionCmd->SetGuidance("Set beam starting position.");
  beamPositionCmd->SetParameterName("X","Y","Z",true,true);
  beamPositionCmd->SetDefaultUnit("mm");
  //beamPositionCmd->SetUnitCategory("Length");
  //beamPositionCmd->SetUnitCandidates("microm mm cm m km");

  beamRadiusAtEntranceCmd = new G4UIcmdWithADoubleAndUnit("/ActarSim/gun/beamRadiusAtEntrance",this);
  beamRadiusAtEntranceCmd->SetGuidance("Selects the beam radius at entrance of ACTAR.");
  beamRadiusAtEntranceCmd->SetGuidance("Used with the emittance to calculate the position and angle");
  beamRadiusAtEntranceCmd->SetGuidance("distributions of the beam when a realisticBeam option is set.");
  beamRadiusAtEntranceCmd->SetParameterName("beamRadius",false);
  beamRadiusAtEntranceCmd->SetRange("beamRadius>0.");
  beamRadiusAtEntranceCmd->SetUnitCategory("Length");
  beamRadiusAtEntranceCmd->SetDefaultValue(1.);
  beamRadiusAtEntranceCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  //HAPOL NOTE: REMOVE OR REBUILT COMPLETELY
  reactionFromEvGenCmd = new G4UIcmdWithAString("/ActarSim/gun/reactionFromEvGen",this);
  reactionFromEvGenCmd->SetGuidance("DO NOT USE. Simulates beam/target from event generator. DO NOT USE.");
  reactionFromEvGenCmd->SetGuidance("  Choice : on, off(default)");
  reactionFromEvGenCmd->SetParameterName("choice",true);
  reactionFromEvGenCmd->SetDefaultValue("off");
  reactionFromEvGenCmd->SetCandidates("on off");
  reactionFromEvGenCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  //commands affecting the input file selection for the reaction
  reactionFromFileCmd = new G4UIcmdWithAString("/ActarSim/gun/reactionFromFile",this);
  reactionFromFileCmd->SetGuidance("Select a reaction from an input file");
  reactionFromFileCmd->SetGuidance("  Choice : on, off(default)");
  reactionFromFileCmd->SetParameterName("choice",true);
  reactionFromFileCmd->SetDefaultValue("off");
  reactionFromFileCmd->SetCandidates("on off");
  reactionFromFileCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  //HAPOL NOTE: REMOVE OR REBUILT COMPLETELY
  reactionFromCrossSectionCmd = new G4UIcmdWithAString("/ActarSim/gun/reactionFromCrossSection",this);
  reactionFromCrossSectionCmd->SetGuidance("DO NOT USE. Simulates beam/target from the cross-sections. DO NOT USE.");
  reactionFromCrossSectionCmd->SetGuidance("  Choice : on, off(default)");
  reactionFromCrossSectionCmd->SetParameterName("choice",true);
  reactionFromCrossSectionCmd->SetDefaultValue("off");
  reactionFromCrossSectionCmd->SetCandidates("on off");
  reactionFromCrossSectionCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  reactionFileCmd = new G4UIcmdWithAString("/ActarSim/gun/reactionFile",this);
  reactionFileCmd->SetGuidance("Select the reaction definition file.");
  reactionFileCmd->SetParameterName("reactionFile",false);
  reactionFileCmd->SetDefaultValue("He8onC12Elastic.dat");

  //commands affecting the Cine kinematic reaction generator
  reactionFromCineCmd = new G4UIcmdWithAString("/ActarSim/gun/reactionFromCine",this);
  reactionFromCineCmd->SetGuidance("Select a reaction using Cine");
  reactionFromCineCmd->SetGuidance("  Choice : on, off(default)");
  reactionFromCineCmd->SetParameterName("choice",true);
  reactionFromCineCmd->SetDefaultValue("off");
  reactionFromCineCmd->SetCandidates("on off");
  reactionFromCineCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  CineDir = new G4UIdirectory("/ActarSim/gun/Cine/");
  CineDir->SetGuidance("Cine generator control");

  CinerandomThetaCmd = new G4UIcmdWithAString("/ActarSim/gun/Cine/randomTheta",this);
  CinerandomThetaCmd->SetGuidance("Select a random Theta angle for the scattered particle.");
  CinerandomThetaCmd->SetGuidance("  Choice : on(default), off");
  CinerandomThetaCmd->SetParameterName("choice",true);
  CinerandomThetaCmd->SetDefaultValue("on");
  CinerandomThetaCmd->SetCandidates("on off");
  CinerandomThetaCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  randomThetaCmd = new G4UIcmdWithAString("/ActarSim/gun/randomTheta",this);
  randomThetaCmd->SetGuidance("Select a random Theta angle for the scattered particle.");
  randomThetaCmd->SetGuidance("  Choice : on(default), off, ext");
  randomThetaCmd->SetParameterName("choice",true);
  randomThetaCmd->SetDefaultValue("on");
  randomThetaCmd->SetCandidates("on off ext");
  randomThetaCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  randomPhiCmd = new G4UIcmdWithAString("/ActarSim/gun/randomPhi",this);
  randomPhiCmd->SetGuidance("Select a random Phi angle for the scattered particle.");
  randomPhiCmd->SetGuidance("  Choice : on(default), off");
  randomPhiCmd->SetParameterName("choice",true);
  randomPhiCmd->SetDefaultValue("on");
  randomPhiCmd->SetCandidates("on off");
  randomPhiCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  alphaSourceCmd = new G4UIcmdWithAString("/ActarSim/gun/alphaSource",this);
  alphaSourceCmd->SetGuidance("CHECK THIS COMMAND!");
  alphaSourceCmd->SetGuidance("  Choice : on(default), off");
  alphaSourceCmd->SetParameterName("choice",true);
  alphaSourceCmd->SetDefaultValue("off");
  alphaSourceCmd->SetCandidates("on off");
  alphaSourceCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  randomThetaValCmd = new G4UIcommand("/ActarSim/gun/randomThetaVal", this);
  randomThetaValCmd->SetGuidance("Sets the limits in the Theta angle for the scattered particle.");
  randomThetaValCmd->SetGuidance("The value is randomly chosen between the limits.");
  parameter = new G4UIparameter("thetaMin", 'd', omitable = true);
  parameter->SetDefaultValue(0.);
  randomThetaValCmd->SetParameter(parameter);
  parameter = new G4UIparameter("thetaMax", 'd', omitable = true);
  parameter->SetDefaultValue(180.);
  randomThetaValCmd->SetParameter(parameter);
  parameter = new G4UIparameter("unit", 's', omitable = true);
  parameter->SetDefaultValue("deg");
  randomThetaValCmd->SetParameter(parameter);

  randomPhiValCmd = new G4UIcommand("/ActarSim/gun/randomPhiVal", this);
  randomPhiValCmd->SetGuidance("Sets the limits in the Phi angle for the scattered particle.");
  randomPhiValCmd->SetGuidance("The value is randomly chosen between the limits.");
  parameter = new G4UIparameter("phiMin", 'd', omitable = true);
  parameter->SetDefaultValue(0.);
  randomPhiValCmd->SetParameter(parameter);
  parameter = new G4UIparameter("phiMax", 'd', omitable = true);
  parameter->SetDefaultValue(180.);
  randomPhiValCmd->SetParameter(parameter);
  parameter = new G4UIparameter("unit", 's', omitable = true);
  parameter->SetDefaultValue("deg");
  randomPhiValCmd->SetParameter(parameter);

  CinerandomThetaValCmd = new G4UIcommand("/ActarSim/gun/Cine/randomThetaVal", this);
  CinerandomThetaValCmd->SetGuidance("Sets the limist in the Theta angle for the scattered particle.");
  CinerandomThetaValCmd->SetGuidance("The value is randomly chosen between the limits.");
  parameter = new G4UIparameter("thetaMin", 'd', omitable = true);
  parameter->SetDefaultValue(0.);
  CinerandomThetaValCmd->SetParameter(parameter);
  parameter = new G4UIparameter("thetaMax", 'd', omitable = true);
  parameter->SetDefaultValue(180.);
  CinerandomThetaValCmd->SetParameter(parameter);
  parameter = new G4UIparameter("unit", 's', omitable = true);
  parameter->SetDefaultValue("deg");
  CinerandomThetaValCmd->SetParameter(parameter);

  incidentIonCmd = new G4UIcommand("/ActarSim/gun/Cine/incidentIon",this);
  incidentIonCmd->SetGuidance("Set properties of incident ion to be generated.");
  incidentIonCmd->SetGuidance("[usage] /ActarSim/gun/Cine/incidentIon Z A Q E");
  incidentIonCmd->SetGuidance("        Z:(int) AtomicNumber");
  incidentIonCmd->SetGuidance("        A:(int) AtomicMass");
  incidentIonCmd->SetGuidance("        Q:(int) Charge of ion (in unit of e)");
  incidentIonCmd->SetGuidance("        E:(double) Excitation energy (in keV)");

  G4UIparameter* incidentParam;
  incidentParam = new G4UIparameter("Z",'i',false);
  incidentParam->SetDefaultValue("1");
  incidentIonCmd->SetParameter(incidentParam);
  incidentParam = new G4UIparameter("A",'i',false);
  incidentParam->SetDefaultValue("1");
  incidentIonCmd->SetParameter(incidentParam);
  incidentParam = new G4UIparameter("Q",'i',false);
  incidentParam->SetDefaultValue("0");
  incidentIonCmd->SetParameter(incidentParam);
  incidentParam = new G4UIparameter("E",'d',true);
  incidentParam->SetDefaultValue("0.0");
  incidentIonCmd->SetParameter(incidentParam);

  targetIonCmd = new G4UIcommand("/ActarSim/gun/Cine/targetIon",this);
  targetIonCmd->SetGuidance("Set properties of target ion to be generated.");
  targetIonCmd->SetGuidance("[usage] /ActarSim/gun/Cine/targetIon Z A Q E");
  targetIonCmd->SetGuidance("        Z:(int) AtomicNumber");
  targetIonCmd->SetGuidance("        A:(int) AtomicMass");
  targetIonCmd->SetGuidance("        Q:(int) Charge of ion (in unit of e)");
  targetIonCmd->SetGuidance("        E:(double) Excitation energy (in keV)");

  G4UIparameter* targetParam;
  targetParam = new G4UIparameter("Z",'i',false);
  targetParam->SetDefaultValue("1");
  targetIonCmd->SetParameter(targetParam);
  targetParam = new G4UIparameter("A",'i',false);
  targetParam->SetDefaultValue("1");
  targetIonCmd->SetParameter(targetParam);
  targetParam = new G4UIparameter("Q",'i',false);
  targetParam->SetDefaultValue("0");
  targetIonCmd->SetParameter(targetParam);
  targetParam = new G4UIparameter("E",'d',true);
  targetParam->SetDefaultValue("0.0");
  targetIonCmd->SetParameter(targetParam);

  scatteredIonCmd = new G4UIcommand("/ActarSim/gun/Cine/scatteredIon",this);
  scatteredIonCmd->SetGuidance("Set properties of scattered ion to be generated.");
  scatteredIonCmd->SetGuidance("[usage] /ActarSim/gun/Cine/scatteredIon Z A Q E");
  scatteredIonCmd->SetGuidance("        Z:(int) AtomicNumber");
  scatteredIonCmd->SetGuidance("        A:(int) AtomicMass");
  scatteredIonCmd->SetGuidance("        Q:(int) Charge of ion (in unit of e)");
  scatteredIonCmd->SetGuidance("        E:(double) Excitation energy (in keV)");

  G4UIparameter* scatteredParam;
  scatteredParam = new G4UIparameter("Z",'i',false);
  scatteredParam->SetDefaultValue("1");
  scatteredIonCmd->SetParameter(scatteredParam);
  scatteredParam = new G4UIparameter("A",'i',false);
  scatteredParam->SetDefaultValue("1");
  scatteredIonCmd->SetParameter(scatteredParam);
  scatteredParam = new G4UIparameter("Q",'i',false);
  scatteredParam->SetDefaultValue("0");
  scatteredIonCmd->SetParameter(scatteredParam);
  scatteredParam = new G4UIparameter("E",'d',true);
  scatteredParam->SetDefaultValue("0.0");
  scatteredIonCmd->SetParameter(scatteredParam);

  recoilIonCmd = new G4UIcommand("/ActarSim/gun/Cine/recoilIon",this);
  recoilIonCmd->SetGuidance("Set properties of recoil ion to be generated.");
  recoilIonCmd->SetGuidance("[usage] /ActarSim/gun/Cine/recoilIon Z A Q E");
  recoilIonCmd->SetGuidance("        Z:(int) AtomicNumber");
  recoilIonCmd->SetGuidance("        A:(int) AtomicMass");
  recoilIonCmd->SetGuidance("        Q:(int) Charge of ion (in unit of e)");
  recoilIonCmd->SetGuidance("        E:(double) Excitation energy (in keV)");

  G4UIparameter* recoilParam;
  recoilParam = new G4UIparameter("Z",'i',false);
  recoilParam->SetDefaultValue("1");
  recoilIonCmd->SetParameter(recoilParam);
  recoilParam = new G4UIparameter("A",'i',false);
  recoilParam->SetDefaultValue("1");
  recoilIonCmd->SetParameter(recoilParam);
  recoilParam = new G4UIparameter("Q",'i',false);
  recoilParam->SetDefaultValue("0");
  recoilIonCmd->SetParameter(recoilParam);
  recoilParam = new G4UIparameter("E",'d',true);
  recoilParam->SetDefaultValue("0.0");
  recoilIonCmd->SetParameter(recoilParam);

  reactionQCmd = new G4UIcmdWithADoubleAndUnit("/ActarSim/gun/Cine/reactionQ",this);
  reactionQCmd->SetGuidance("Sets the reaction Q ");
  reactionQCmd->SetParameterName("reactionQ",false);
  //reactionQCmd->SetRange("reactionQ>=0.");
  reactionQCmd->SetUnitCategory("Energy");
  reactionQCmd->SetDefaultValue(12.);
  reactionQCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  labEnergyCmd = new G4UIcmdWithADoubleAndUnit("/ActarSim/gun/Cine/labEnergy",this);
  labEnergyCmd->SetGuidance("Sets the laboratory energy ");
  labEnergyCmd->SetParameterName("labEnergy",false);
  labEnergyCmd->SetRange("labEnergy>=0.");
  labEnergyCmd->SetUnitCategory("Energy");
  labEnergyCmd->SetDefaultValue(100.);
  labEnergyCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  thetaLabAngleCmd = new G4UIcmdWithADoubleAndUnit("/ActarSim/gun/Cine/thetaLabAngle",this);
  thetaLabAngleCmd->SetGuidance("Sets theta lab angle for the scattered particle (degrees)");
  thetaLabAngleCmd->SetParameterName("thetaLabAngle",false);
  thetaLabAngleCmd->SetRange("thetaLabAngle>=0.");
  thetaLabAngleCmd->SetUnitCategory("Angle");
  thetaLabAngleCmd->SetDefaultValue(0.5);
  thetaLabAngleCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  //commands affecting the Cine kinematic reaction generator
  reactionFromKineCmd = new G4UIcmdWithAString("/ActarSim/gun/reactionFromKine",this);
  reactionFromKineCmd->SetGuidance("Select a reaction using Kine");
  reactionFromKineCmd->SetGuidance("  Choice : on(default), off");
  reactionFromKineCmd->SetParameterName("choice",true);
  reactionFromKineCmd->SetDefaultValue("on");
  reactionFromKineCmd->SetCandidates("on off");
  reactionFromKineCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  KineDir = new G4UIdirectory("/ActarSim/gun/Kine/");
  KineDir->SetGuidance("Kine generator control");

  KineRandomThetaCmd = new G4UIcmdWithAString("/ActarSim/gun/Kine/randomThetaCM",this);
  KineRandomThetaCmd->SetGuidance("Randomize Theta_CM of outgoing particles");
  KineRandomThetaCmd->SetGuidance("  Choice : on(default), off");
  KineRandomThetaCmd->SetParameterName("choice",true);
  KineRandomThetaCmd->SetDefaultValue("on");
  KineRandomThetaCmd->SetCandidates("on off");
  KineRandomThetaCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  KineRandomPhiAngleCmd = new G4UIcmdWithAString("/ActarSim/gun/Kine/randomPhiAngle",this);
  KineRandomPhiAngleCmd->SetGuidance("Randomize Lab Phi angles of out-going particles");
  KineRandomPhiAngleCmd->SetGuidance("  Choice : on (default), off");
  KineRandomPhiAngleCmd->SetParameterName("choice",true);
  KineRandomPhiAngleCmd->SetDefaultValue("on");
  KineRandomPhiAngleCmd->SetCandidates("on off");
  KineRandomPhiAngleCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  KineRandomThetaRangeCmd = new G4UIcommand("/ActarSim/gun/Kine/randomThetaRange", this);
  KineRandomThetaRangeCmd->SetGuidance("Sets the limits in the Theta angle for the scattered particle.");
  KineRandomThetaRangeCmd->SetGuidance("The value is randomly chosen between the limits.");
  parameter = new G4UIparameter("thetaMin", 'd', omitable = true);
  parameter->SetDefaultValue(0.);
  KineRandomThetaRangeCmd->SetParameter(parameter);
  parameter = new G4UIparameter("thetaMax", 'd', omitable = true);
  parameter->SetDefaultValue(180.);
  KineRandomThetaRangeCmd->SetParameter(parameter);
  parameter = new G4UIparameter("unit", 's', omitable = true);
  parameter->SetDefaultValue("deg");
  KineRandomThetaRangeCmd->SetParameter(parameter);

  KineIncidentIonCmd = new G4UIcommand("/ActarSim/gun/Kine/incidentIon",this);
  KineIncidentIonCmd->SetGuidance("Set properties of incident ion to be generated.");
  KineIncidentIonCmd->SetGuidance("[usage] /ActarSim/gun/Kine/incidentIon Z A Q E Mass");
  KineIncidentIonCmd->SetGuidance("        Z:(int) AtomicNumber");
  KineIncidentIonCmd->SetGuidance("        A:(int) AtomicMass (in Atomic mass unit u)");
  KineIncidentIonCmd->SetGuidance("        Q:(int) Charge of ion (in unit of e)");
  KineIncidentIonCmd->SetGuidance("        E:(double) Excitation energy (in MeV)");
  KineIncidentIonCmd->SetGuidance("     Mass:(double) mass in u");

  G4UIparameter* KineIncidentParam;
  KineIncidentParam = new G4UIparameter("Z",'i',false);
  KineIncidentParam->SetDefaultValue("1");
  KineIncidentIonCmd->SetParameter(KineIncidentParam);
  KineIncidentParam = new G4UIparameter("A",'i',false);
  KineIncidentParam->SetDefaultValue("1");
  KineIncidentIonCmd->SetParameter(KineIncidentParam);
  KineIncidentParam = new G4UIparameter("Q",'i',false);
  KineIncidentParam->SetDefaultValue("0");
  KineIncidentIonCmd->SetParameter(KineIncidentParam);
  KineIncidentParam = new G4UIparameter("E",'d',true);
  KineIncidentParam->SetDefaultValue("0.0");
  KineIncidentIonCmd->SetParameter(KineIncidentParam);
  KineIncidentParam = new G4UIparameter("Mass",'d',true);
  KineIncidentParam->SetDefaultValue("1.0");
  KineIncidentIonCmd->SetParameter(KineIncidentParam);

  KineTargetIonCmd = new G4UIcommand("/ActarSim/gun/Kine/targetIon",this);
  KineTargetIonCmd->SetGuidance("Set properties of target ion to be generated.");
  KineTargetIonCmd->SetGuidance("[usage] /ActarSim/gun/Cine/targetIon Z A Q E Mass");
  KineTargetIonCmd->SetGuidance("        Z:(int) AtomicNumber");
  KineTargetIonCmd->SetGuidance("        A:(int) AtomicMass in (in u)");
  KineTargetIonCmd->SetGuidance("        Q:(int) Charge of ion (in unit of e)");
  KineTargetIonCmd->SetGuidance("        E:(double) Excitation energy (in MeV)");
  KineTargetIonCmd->SetGuidance("     Mass:(double) mass in u");

  G4UIparameter* KineTargetParam;
  KineTargetParam = new G4UIparameter("Z",'i',false);
  KineTargetParam->SetDefaultValue("1");
  KineTargetIonCmd->SetParameter(KineTargetParam);
  KineTargetParam = new G4UIparameter("A",'i',false);
  KineTargetParam->SetDefaultValue("1");
  KineTargetIonCmd->SetParameter(KineTargetParam);
  KineTargetParam = new G4UIparameter("Q",'i',false);
  KineTargetParam->SetDefaultValue("0");
  KineTargetIonCmd->SetParameter(KineTargetParam);
  KineTargetParam = new G4UIparameter("E",'d',true);
  KineTargetParam->SetDefaultValue("0.0");
  KineTargetIonCmd->SetParameter(KineTargetParam);
  KineTargetParam = new G4UIparameter("Mass",'d',true);
  KineTargetParam->SetDefaultValue("1.0");
  KineTargetIonCmd->SetParameter(KineTargetParam);

  KineScatteredIonCmd = new G4UIcommand("/ActarSim/gun/Kine/scatteredIon",this);
  KineScatteredIonCmd->SetGuidance("Set properties of scattered ion to be generated.");
  KineScatteredIonCmd->SetGuidance("[usage] /ActarSim/gun/Cine/scatteredIon Z A Q E Mass");
  KineScatteredIonCmd->SetGuidance("        Z:(int) AtomicNumber");
  KineScatteredIonCmd->SetGuidance("        A:(int) AtomicMass");
  KineScatteredIonCmd->SetGuidance("        Q:(int) Charge of ion (in unit of e)");
  KineScatteredIonCmd->SetGuidance("        E:(double) Excitation energy (in MeV)");
  KineScatteredIonCmd->SetGuidance("     Mass:(double) mass in u");

  G4UIparameter* KineScatteredParam;
  KineScatteredParam = new G4UIparameter("Z",'i',false);
  KineScatteredParam->SetDefaultValue("1");
  KineScatteredIonCmd->SetParameter(KineScatteredParam);
  KineScatteredParam = new G4UIparameter("A",'i',false);
  KineScatteredParam->SetDefaultValue("1");
  KineScatteredIonCmd->SetParameter(KineScatteredParam);
  KineScatteredParam = new G4UIparameter("Q",'i',false);
  KineScatteredParam->SetDefaultValue("0");
  KineScatteredIonCmd->SetParameter(KineScatteredParam);
  KineScatteredParam = new G4UIparameter("E",'d',true);
  KineScatteredParam->SetDefaultValue("0.0");
  KineScatteredIonCmd->SetParameter(KineScatteredParam);
  KineScatteredParam = new G4UIparameter("Mass",'d',true);
  KineScatteredParam->SetDefaultValue("1.0");
  KineScatteredIonCmd->SetParameter(KineScatteredParam);

  KineRecoilIonCmd = new G4UIcommand("/ActarSim/gun/Kine/recoilIon",this);
  KineRecoilIonCmd->SetGuidance("Set properties of recoil ion to be generated.");
  KineRecoilIonCmd->SetGuidance("[usage] /ActarSim/gun/Cine/recoilIon Z A Q E Mass");
  KineRecoilIonCmd->SetGuidance("        Z:(int) AtomicNumber");
  KineRecoilIonCmd->SetGuidance("        A:(int) AtomicMass");
  KineRecoilIonCmd->SetGuidance("        Q:(int) Charge of ion (in unit of e)");
  KineRecoilIonCmd->SetGuidance("        E:(double) Excitation energy (in MeV)");
  KineRecoilIonCmd->SetGuidance("     Mass:(double) mass in u");

  G4UIparameter* KineRecoilParam;
  KineRecoilParam = new G4UIparameter("Z",'i',false);
  KineRecoilParam->SetDefaultValue("1");
  KineRecoilIonCmd->SetParameter(KineRecoilParam);
  KineRecoilParam = new G4UIparameter("A",'i',false);
  KineRecoilParam->SetDefaultValue("1.");
  KineRecoilIonCmd->SetParameter(KineRecoilParam);
  KineRecoilParam = new G4UIparameter("Q",'i',false);
  KineRecoilParam->SetDefaultValue("0");
  KineRecoilIonCmd->SetParameter(KineRecoilParam);
  KineRecoilParam = new G4UIparameter("E",'d',true);
  KineRecoilParam->SetDefaultValue("0.0");
  KineRecoilIonCmd->SetParameter(KineRecoilParam);
  KineRecoilParam = new G4UIparameter("Mass",'d',true);
  KineRecoilParam->SetDefaultValue("1.0");
  KineRecoilIonCmd->SetParameter(KineRecoilParam);

  KineLabEnergyCmd = new G4UIcmdWithADoubleAndUnit("/ActarSim/gun/Kine/labEnergy",this);
  KineLabEnergyCmd->SetGuidance("Sets the laboratory energy.");
  KineLabEnergyCmd->SetParameterName("labEnergy",false);
  KineLabEnergyCmd->SetRange("labEnergy>=0.");
  KineLabEnergyCmd->SetUnitCategory("Energy");
  KineLabEnergyCmd->SetDefaultValue(100.);
  KineLabEnergyCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  KineUserThetaCMCmd = new G4UIcmdWithADoubleAndUnit("/ActarSim/gun/Kine/userThetaCM",this);
  KineUserThetaCMCmd->SetGuidance("Sets theta CM angle for scattered particle (in degrees)");
  KineUserThetaCMCmd->SetParameterName("userThetaCM",false);
  KineUserThetaCMCmd->SetRange("userThetaCM>=0.");
  KineUserThetaCMCmd->SetUnitCategory("Angle");
  KineUserThetaCMCmd->SetDefaultValue(0.);
  KineUserThetaCMCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  // user set the phi angle of particles, useful when testing the kinematics reconstruction methods
  KineUserPhiAngleCmd = new G4UIcmdWithADoubleAndUnit("/ActarSim/gun/Kine/userPhiAngle",this);
  KineUserPhiAngleCmd->SetGuidance("User set phi angle for outgoing particle in the Lab system (in degrees)");
  KineUserPhiAngleCmd->SetParameterName("userPhiAngle",false);
  KineUserPhiAngleCmd->SetRange("userPhiAngle>=0.");
  KineUserPhiAngleCmd->SetUnitCategory("Angle");
  KineUserPhiAngleCmd->SetDefaultValue(0.0);
  KineUserPhiAngleCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  vertexPositionCmd = new G4UIcmdWith3VectorAndUnit("/ActarSim/gun/Kine/vertexPosition",this);
  vertexPositionCmd->SetGuidance("Set the position of the vertex.");
  vertexPositionCmd->SetParameterName("X","Y","Z",true,true);
  vertexPositionCmd->SetDefaultUnit("cm");

  //commands affecting individual particles
  energyCmd  = new G4UIcmdWithADoubleAndUnit("/ActarSim/gun/energy",this);
  energyCmd->SetGuidance("Sets the kinetic energy of the primary particle");
  energyCmd->SetParameterName("energy",false);
  energyCmd->SetRange("energy>=0.");
  energyCmd->SetUnitCategory("Energy");
  energyCmd->SetDefaultValue(1.);
  energyCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  directionCmd = new G4UIcmdWith3Vector("/ActarSim/gun/direction",this);
  directionCmd->SetGuidance("Set momentum direction.");
  directionCmd->SetGuidance("Direction does not need to be a unit vector.");
  directionCmd->SetParameterName("Px","Py","Pz",true,true);
  directionCmd->SetRange("Px != 0 || Py != 0 || Pz != 0");

  positionCmd = new G4UIcmdWith3VectorAndUnit("/ActarSim/gun/position",this);
  positionCmd->SetGuidance("Set starting position of the particle.");
  positionCmd->SetParameterName("X","Y","Z",true,true);
  positionCmd->SetDefaultUnit("cm");
  //positionCmd->SetUnitCategory("Length");
  //positionCmd->SetUnitCandidates("microm mm cm m km");

  timeCmd = new G4UIcmdWithADoubleAndUnit("/ActarSim/gun/time",this);
  timeCmd->SetGuidance("Set initial time of the particle.");
  timeCmd->SetParameterName("t0",true,true);
  timeCmd->SetDefaultUnit("ns");
  //timeCmd->SetUnitCategory("Time");
  //timeCmd->SetUnitCandidates("ns ms s");

  randomVertexZPositionCmd = new G4UIcmdWithAString("/ActarSim/gun/randomVertexZPosition",this);
  randomVertexZPositionCmd->SetGuidance("Randomize the reaction vertex Z position");
  randomVertexZPositionCmd->SetGuidance("Choice : on(default), off , ext");
  randomVertexZPositionCmd->SetParameterName("choice",true);
  randomVertexZPositionCmd->SetDefaultValue("on");
  randomVertexZPositionCmd->SetCandidates("on off ext");
  randomVertexZPositionCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  randomVertexZPositionRangeCmd = new G4UIcommand("/ActarSim/gun/randomVertexZRange", this);
  randomVertexZPositionRangeCmd->SetGuidance("Set the min and max Z-value of random vertex position");
  randomVertexZPositionRangeCmd->SetGuidance("The value is randomly chosen between the limits.");
  parameter = new G4UIparameter("randomVertexZMin", 'd', omitable = true);
  parameter->SetDefaultValue(0.);
  randomVertexZPositionRangeCmd->SetParameter(parameter);
  parameter = new G4UIparameter("randomVertexZMax", 'd', omitable = true);
  parameter->SetDefaultValue(300.);
  randomVertexZPositionRangeCmd->SetParameter(parameter);
  parameter = new G4UIparameter("unit", 's', omitable = true);
  parameter->SetDefaultValue("mm");
  randomVertexZPositionRangeCmd->SetParameter(parameter);

  vertexZPositionCmd = new G4UIcmdWithADoubleAndUnit("/ActarSim/gun/vertexZPosition",this);
  vertexZPositionCmd->SetGuidance("Set the Z-value of the reaction vertex.");
  vertexZPositionCmd->SetParameterName("Z0",true,true);
  vertexZPositionCmd->SetDefaultUnit("mm");
  vertexZPositionCmd->SetUnitCategory("Length");
  vertexZPositionCmd->SetDefaultValue(0.0);
  vertexZPositionCmd->SetUnitCandidates("mm cm m");

  polCmd = new G4UIcmdWith3Vector("/ActarSim/gun/polarization",this);
  polCmd->SetGuidance("Set polarization.");
  polCmd->SetParameterName("Px","Py","Pz",true,true);
  polCmd->SetRange("Px>=-1.&&Px<=1.&&Py>=-1.&&Py<=1.&&Pz>=-1.&&Pz<=1.");

  numberCmd = new G4UIcmdWithAnInteger("/ActarSim/gun/number",this);
  numberCmd->SetGuidance("Set number of particles to be generated in a single event.");
  numberCmd->SetParameterName("N",true,true);
  numberCmd->SetRange("N>0");

  ionCmd = new G4UIcommand("/ActarSim/gun/ion",this);
  ionCmd->SetGuidance("Set properties of ion to be generated.");
  ionCmd->SetGuidance("[usage] /ActarSim/gun/ion Z A Q E");
  ionCmd->SetGuidance("        Z:(int) AtomicNumber");
  ionCmd->SetGuidance("        A:(int) AtomicMass");
  ionCmd->SetGuidance("        Q:(int) Charge of Ion (in unit of e)");
  ionCmd->SetGuidance("        E:(double) Excitation energy (in keV)");

  G4UIparameter* param;
  param = new G4UIparameter("Z",'i',false);
  param->SetDefaultValue("1");
  ionCmd->SetParameter(param);
  param = new G4UIparameter("A",'i',false);
  param->SetDefaultValue("1");
  ionCmd->SetParameter(param);
  param = new G4UIparameter("Q",'i',true);
  param->SetDefaultValue("0");
  ionCmd->SetParameter(param);
  param = new G4UIparameter("E",'d',true);
  param->SetDefaultValue("0.0");
  ionCmd->SetParameter(param);
}

//////////////////////////////////////////////////////////////////
/// Destructor
ActarSimPrimaryGeneratorMessenger::~ActarSimPrimaryGeneratorMessenger() {
  delete gunDir;
  delete listCmd;
  delete particleCmd;
  delete energyCmd;
  delete directionCmd;
  delete beamThetaCmd;
  delete beamPhiCmd;
  delete positionCmd;
  delete timeCmd;
  delete randomVertexZPositionCmd;
  delete randomVertexZPositionRangeCmd;
  delete vertexZPositionCmd;
  delete polCmd;
  delete numberCmd;
  delete ionCmd;
  delete realisticBeamCmd;
  delete beamInteractionCmd;
  delete emittanceCmd;
  delete beamDirectionCmd;
  delete beamPositionCmd;
  delete beamRadiusAtEntranceCmd;
  delete reactionFromEvGenCmd;
  delete reactionFromCrossSectionCmd;
  delete reactionFromFileCmd;
  delete reactionFileCmd;
  delete randomThetaCmd;
  delete randomPhiCmd;
  delete alphaSourceCmd;
  delete CinerandomThetaCmd;
  delete randomThetaValCmd;
  delete randomPhiValCmd;
  delete CinerandomThetaValCmd;
  delete reactionFromCineCmd;
  delete CineDir;
  delete incidentIonCmd;
  delete targetIonCmd;
  delete scatteredIonCmd;
  delete recoilIonCmd;
  delete reactionQCmd;
  delete labEnergyCmd;
  delete thetaLabAngleCmd;
  delete reactionFromKineCmd;
  delete KineDir;
  delete KineRandomThetaCmd;
  delete KineRandomThetaRangeCmd;
  delete KineRandomPhiAngleCmd;
  delete KineIncidentIonCmd;
  delete KineTargetIonCmd;
  delete KineScatteredIonCmd;
  delete KineRecoilIonCmd;
  delete KineLabEnergyCmd;           // in MeV
  delete KineUserThetaCMCmd;         // in degrees
  delete KineUserPhiAngleCmd;        // in degrees
  delete vertexPositionCmd;
}

//////////////////////////////////////////////////////////////////
/// Setting the values using the interface
void ActarSimPrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command,
						    G4String newValues) {
  if( command==listCmd )
    particleTable->DumpTable();

  if( command == particleCmd ) {
    if (newValues =="ion") {
      fShootIon = true;
    } else {
      fShootIon = false;
      G4ParticleDefinition* pd = particleTable->FindParticle(newValues);
      if(pd != 0)
	{ actarSimActionGun->SetParticleDefinition( pd ); }
    }
  }

  if( command == realisticBeamCmd )
    actarSimActionGun->SetRealisticBeamFlag(newValues);

  if( command == beamInteractionCmd )
    actarSimActionGun->SetBeamInteractionFlag(newValues);

  if( command == emittanceCmd)
    actarSimActionGun->SetEmittance(emittanceCmd->GetNewDoubleValue(newValues));

  if( command==beamDirectionCmd )
    actarSimActionGun->SetBeamMomentumDirection(beamDirectionCmd->GetNew3VectorValue(newValues));

  if( command == beamThetaCmd )
    actarSimActionGun->
      SetUserThetaAngle(beamThetaCmd->GetNewDoubleValue(newValues));

  if( command == beamPhiCmd )
    actarSimActionGun->
      SetUserPhiAngle(beamPhiCmd->GetNewDoubleValue(newValues));

  if( command==beamPositionCmd )
    actarSimActionGun->SetBeamPosition(beamPositionCmd->GetNew3VectorValue(newValues));

  if(command == beamRadiusAtEntranceCmd)
    actarSimActionGun->SetBeamRadiusAtEntrance(beamRadiusAtEntranceCmd->GetNewDoubleValue(newValues));

  if( command == reactionFromEvGenCmd )
    actarSimActionGun->SetReactionFromEvGenFlag(newValues);

  if( command == reactionFromCrossSectionCmd )
    actarSimActionGun->SetReactionFromCrossSectionFlag(newValues);

  if( command == reactionFromFileCmd )
    actarSimActionGun->SetReactionFromFileFlag(newValues);

  if( command == reactionFromCineCmd )
    actarSimActionGun->SetReactionFromCineFlag(newValues);

  if( command == reactionFileCmd )
    actarSimActionGun->SetReactionFile(newValues);

  if( command == randomThetaCmd )
    actarSimActionGun->SetRandomThetaFlag(newValues);

  if( command == randomPhiCmd )
    actarSimActionGun->SetRandomPhiFlag(newValues);

  if( command == alphaSourceCmd )
    actarSimActionGun->SetAlphaSourceFlag(newValues);

  if( command == CinerandomThetaCmd )
    actarSimActionGun->SetRandomThetaFlag(newValues);

  if( command == randomThetaValCmd ){
    G4double thetaMax, thetaMin;
    //    ConvertToDoublePair(newValues, thetaMin, thetaMax);

    G4double x, y;
    char unts[30];
    std::istringstream is(newValues);
    is >> x >> y >> unts;
    G4String unt = unts;

    thetaMin = x*G4UIcommand::ValueOf(unt);
    thetaMax = y*G4UIcommand::ValueOf(unt);

    actarSimActionGun->SetRandomThetaVal(thetaMin,thetaMax);
  }

  if( command == randomPhiValCmd ){
    G4double phiMax, phiMin;
    //    ConvertToDoublePair(newValues, thetaMin, thetaMax);

    G4double x, y;
    char unts[30];
    std::istringstream is(newValues);
    is >> x >> y >> unts;
    G4String unt = unts;

    phiMin = x*G4UIcommand::ValueOf(unt);
    phiMax = y*G4UIcommand::ValueOf(unt);

    actarSimActionGun->SetRandomPhiVal(phiMin,phiMax);
  }

  if( command == CinerandomThetaValCmd ){
    G4double thetaMax, thetaMin;
    //    ConvertToDoublePair(newValues, thetaMin, thetaMax);

    G4double x, y;
    char unts[30];
    std::istringstream is(newValues);
    is >> x >> y >> unts;
    G4String unt = unts;

    thetaMin = x*G4UIcommand::ValueOf(unt);
    thetaMax = y*G4UIcommand::ValueOf(unt);

    actarSimActionGun->SetRandomThetaVal(thetaMin,thetaMax);
  }

  if( command == incidentIonCmd )
    incidentIonCommand(newValues);

  if( command == targetIonCmd )
    targetIonCommand(newValues);

  if( command == scatteredIonCmd )
    scatteredIonCommand(newValues);

  if( command == recoilIonCmd )
    recoilIonCommand(newValues);

  if( command == reactionFromKineCmd )
    actarSimActionGun->SetReactionFromKineFlag(newValues);

  if( command == KineRandomThetaCmd )
    actarSimActionGun->SetRandomThetaFlag(newValues);

  if( command == KineRandomPhiAngleCmd )
    actarSimActionGun->SetRandomPhiAngleFlag(newValues);

  if( command == KineRandomThetaRangeCmd ){
    G4double thetaMax, thetaMin;

    G4double x, y;
    char unts[30];
    std::istringstream is(newValues);
    is >> x >> y >> unts;
    G4String unt = unts;

    thetaMin = x*G4UIcommand::ValueOf(unt);
    thetaMax = y*G4UIcommand::ValueOf(unt);

    actarSimActionGun->SetRandomThetaVal(thetaMin,thetaMax);
  }

  if( command == KineIncidentIonCmd )
    KineIncidentIonCommand(newValues);

  if( command == KineTargetIonCmd )
    KineTargetIonCommand(newValues);

  if( command == KineScatteredIonCmd )
    KineScatteredIonCommand(newValues);

  if( command == KineRecoilIonCmd )
    KineRecoilIonCommand(newValues);

  if( command == KineLabEnergyCmd ){
    G4double incidentEnergyTmp;
    incidentEnergyTmp=KineLabEnergyCmd->GetNewDoubleValue(newValues);
    actarSimActionGun->SetLabEnergy(incidentEnergyTmp);
    actarSimActionGun->SetIncidentEnergy(incidentEnergyTmp);//Piotr: why do we have incident energy and lab energy?
  }

  if( command == KineUserThetaCMCmd )
    actarSimActionGun->
      SetThetaCMAngle(KineUserThetaCMCmd->GetNewDoubleValue(newValues));

  if( command == KineUserPhiAngleCmd )
    actarSimActionGun->
      SetUserPhiAngle(KineUserPhiAngleCmd->GetNewDoubleValue(newValues));

  if( command==randomVertexZPositionCmd )
    actarSimActionGun->SetRandomVertexZPositionFlag(newValues);

  if( command==vertexPositionCmd )
    actarSimActionGun->SetVertexPosition(vertexPositionCmd->GetNew3VectorValue(newValues));


  if( command == reactionQCmd )
    actarSimActionGun->
      SetReactionQ(reactionQCmd->GetNewDoubleValue(newValues));

  if( command == labEnergyCmd )
    actarSimActionGun->
      SetLabEnergy(labEnergyCmd->GetNewDoubleValue(newValues));

  if( command == thetaLabAngleCmd )
    actarSimActionGun->
      SetThetaLabAngle(thetaLabAngleCmd->GetNewDoubleValue(newValues));

  if( command == energyCmd ){
    actarSimActionGun->SetParticleEnergy(energyCmd->GetNewDoubleValue(newValues));//Piotr: Doesn't set the energy when modify
    actarSimActionGun->SetIncidentEnergy(energyCmd->GetNewDoubleValue(newValues));
  }

  if( command==directionCmd )
    actarSimActionGun->SetParticleMomentumDirection(directionCmd->GetNew3VectorValue(newValues));

  if( command==positionCmd )
    actarSimActionGun->SetParticlePosition(positionCmd->GetNew3VectorValue(newValues));

  if( command==timeCmd )
    actarSimActionGun->SetParticleTime(timeCmd->GetNewDoubleValue(newValues));

  if( command==vertexZPositionCmd )
    actarSimActionGun->SetVertexZPosition(vertexZPositionCmd->GetNewDoubleValue(newValues));

  if( command == randomVertexZPositionRangeCmd ){
    G4double vertexZMin, vertexZMax;

    G4double x, y;
    char unts[30];
    std::istringstream is(newValues);
    is >> x >> y >> unts;
    G4String unt = unts;

    vertexZMin = x*G4UIcommand::ValueOf(unt);
    vertexZMax = y*G4UIcommand::ValueOf(unt);

    actarSimActionGun->SetRandomVertexZPositionVal(vertexZMin,vertexZMax);
  }

  if( command==polCmd )
    actarSimActionGun->SetParticlePolarization(polCmd->GetNew3VectorValue(newValues));

  if( command==numberCmd )
    actarSimActionGun->SetNumberOfParticles(numberCmd->GetNewIntValue(newValues));

  if( command==ionCmd )
    IonCommand(newValues);
}

//////////////////////////////////////////////////////////////////
/// Get current value from commands
G4String ActarSimPrimaryGeneratorMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4String cv;

  if( command==directionCmd) {
    cv = directionCmd->ConvertToString(actarSimActionGun->GetParticleMomentumDirection());
  }
  else if( command==particleCmd) {
    cv = actarSimActionGun->GetParticleDefinition()->GetParticleName();
  }
  else if( command==energyCmd) {
    cv = energyCmd->ConvertToString(actarSimActionGun->GetParticleEnergy(),"GeV");
  }
  else if( command==positionCmd) {
    cv = positionCmd->ConvertToString(actarSimActionGun->GetParticlePosition(),"cm");
  }
  else if( command==timeCmd) {
    cv = timeCmd->ConvertToString(actarSimActionGun->GetParticleTime(),"ns");
  }
  else if( command==polCmd) {
    cv = polCmd->ConvertToString(actarSimActionGun->GetParticlePolarization());
  }
  else if( command==numberCmd) {
    cv = numberCmd->ConvertToString(actarSimActionGun->GetNumberOfParticles());
  }
  else if( command==ionCmd) {
    if (fShootIon) {
      cv = ItoS(fAtomicNumber) + " " + ItoS(fAtomicMass) + " ";
      cv += ItoS(fIonCharge);
    } else {
      cv = "";
    }
  }
  return cv;
}

//////////////////////////////////////////////////////////////////
/// Particular behavior of the ion command. Ion state should be selected.
void ActarSimPrimaryGeneratorMessenger::IonCommand(G4String newValues) {
  //if (fShootIon) {
  G4Tokenizer next( newValues );
  // check argument
  fAtomicNumber = StoI(next());
  fAtomicMass = StoI(next());
  G4String sQ = next();
  if (sQ.isNull()) {
    fIonCharge = fAtomicNumber;
  } else {
    fIonCharge = StoI(sQ);
    sQ = next();
    if (sQ.isNull()) {
      fIonExciteEnergy = 0.0;
    } else {
      fIonExciteEnergy = StoD(sQ) * keV;
    }
  }

  G4ParticleDefinition* ion;
  ion =  ionTable->GetIon( fAtomicNumber, fAtomicMass, fIonExciteEnergy);
  if (ion==0) {
    G4cout << "##################################################################"
	   << G4endl
	   << "#######   ActarSimPrimaryGeneratorMessenger::IonCommand()  #######"
	   << "Ion with Z=" << fAtomicNumber
	   << " A=" << fAtomicMass << "can not be defined" << G4endl;
    G4cout << "##################################################################"
	   << G4endl;
  } else {
    actarSimActionGun->SetParticleDefinition(ion);
    actarSimActionGun->SetParticleCharge(fIonCharge*eplus);
  }
  /*
    } else {
    G4cout << "##################################################################"
	   << G4endl
	   << "#######   ActarSimPrimaryGeneratorMessenger::IonCommand()  #######"
	   << "Set /gun/particle to ion before using /gun/ion command" << G4endl;
    G4cout << "##################################################################"
	   << G4endl;
	   }
  */
}

//////////////////////////////////////////////////////////////////
/// Particular behavior of the incident ion command. Ion state should be selected.
void ActarSimPrimaryGeneratorMessenger::incidentIonCommand(G4String newValues){
  G4Tokenizer next( newValues );
  // check argument
  fAtomicNumber = StoI(next());
  fAtomicMass = StoI(next());
  G4String sQ = next();
  if (sQ.isNull()) {
    fIonCharge = fAtomicNumber;
  } else {
    fIonCharge = StoI(sQ);
    sQ = next();
    if (sQ.isNull()) {
      fIonExciteEnergy = 0.0;
    } else {
      fIonExciteEnergy = StoD(sQ) * keV;
    }
  }

  G4Ions* ion;
  ion = (G4Ions*) ionTable->GetIon(fAtomicNumber,
					fAtomicMass,
					fIonExciteEnergy);
  if (ion==0) {
    G4cout << "##################################################################"
	   << G4endl
	   << "####  ActarSimPrimaryGeneratorMessenger::incidentIonCommand() ####"
	   << "Ion with Z=" << fAtomicNumber
	   << " A=" << fAtomicMass << "can not be defined" << G4endl;
    G4cout << "##################################################################"
	   << G4endl;
  } else {
    actarSimActionGun->SetIncidentIon(ion);
    actarSimActionGun->SetIncidentIonCharge(fIonCharge*eplus);
    //actarSimActionGun->SetIncidentIonExcEnergy(fIonExciteEnergy);
  }
}

//////////////////////////////////////////////////////////////////
/// Particular behavior of the (KINE) incident ion command. Ion state should be selected.
void ActarSimPrimaryGeneratorMessenger::KineIncidentIonCommand(G4String newValues){
  G4Tokenizer next( newValues );
  // check argument
  fAtomicNumber = StoI(next());
  fAtomicMass = StoI(next());
  fIonCharge = StoI(next());
  fIonExciteEnergy = StoD(next());
  fIonMass=StoD(next());

  G4Ions* ion;
  ion = (G4Ions*) ionTable->GetIon(fAtomicNumber,
					fAtomicMass,
					fIonExciteEnergy);
  if (ion==0) {
    G4cout << "##################################################################"
	   << G4endl
	   << "####  ActarSimPrimaryGeneratorMessenger::KineIncidentIonCommand() ####"
	   << "Ion with Z=" << fAtomicNumber
	   << " A=" << fAtomicMass << "can not be defined" << G4endl;
    G4cout << "##################################################################"
	   << G4endl;
  } else {
    actarSimActionGun->SetIncidentIon(ion);
    actarSimActionGun->SetIncidentIonCharge(fIonCharge*eplus);
    actarSimActionGun->SetMassOfProjectile(fIonMass);
    actarSimActionGun->SetExEnergyOfProjectile(fIonExciteEnergy);
    //actarSimActionGun->SetIncidentIonExcEnergy(fIonExciteEnergy);
  }
}

//////////////////////////////////////////////////////////////////
/// Particular behavior of the target ion command. Ion state should be selected.
void ActarSimPrimaryGeneratorMessenger::targetIonCommand(G4String newValues){
  G4Tokenizer next( newValues );
  // check argument
  fAtomicNumber = StoI(next());
  fAtomicMass = StoI(next());
  G4String sQ = next();
  if (sQ.isNull()) {
    fIonCharge = fAtomicNumber;
  } else {
    fIonCharge = StoI(sQ);
    sQ = next();
    if (sQ.isNull()) {
      fIonExciteEnergy = 0.0;
    } else {
      fIonExciteEnergy = StoD(sQ) * keV;
    }
  }

  G4Ions* ion;
  ion = (G4Ions*) ionTable->GetIon(fAtomicNumber,
					fAtomicMass,
					fIonExciteEnergy);
  if (ion==0) {
    G4cout << "##################################################################"
	   << G4endl
	   << "#####  ActarSimPrimaryGeneratorMessenger::targetIonCommand() ####"
	   << "Ion with Z=" << fAtomicNumber
	   << " A=" << fAtomicMass << "can not be defined" << G4endl;
    G4cout << "##################################################################"
	   << G4endl;
  } else {
    actarSimActionGun->SetTargetIon(ion);
    actarSimActionGun->SetTargetIonCharge(fIonCharge*eplus);
    //actarSimActionGun->SetTargetIonExcEnergy(fIonExciteEnergy);
  }
}

//////////////////////////////////////////////////////////////////
/// Particular behavior of the (KINE) target ion command. Ion state should be selected.
void ActarSimPrimaryGeneratorMessenger::KineTargetIonCommand(G4String newValues){
  G4Tokenizer next( newValues );
  // check argument
  fAtomicNumber = StoI(next());
  fAtomicMass = StoI(next());
  fIonCharge = StoI(next());
  fIonExciteEnergy = StoD(next());
  fIonMass=StoD(next());

  G4Ions* ion;
  ion = (G4Ions*) ionTable->GetIon(fAtomicNumber,
					fAtomicMass,
					fIonExciteEnergy);
  if (ion==0) {
    G4cout << "##################################################################"
	   << G4endl
	   << "#####  ActarSimPrimaryGeneratorMessenger::KineTargetIonCommand() ####"
	   << "Ion with Z=" << fAtomicNumber
	   << " A=" << fAtomicMass << "can not be defined" << G4endl;
    G4cout << "##################################################################"
	   << G4endl;
  } else {
    actarSimActionGun->SetTargetIon(ion);
    actarSimActionGun->SetTargetIonCharge(fIonCharge*eplus);
    actarSimActionGun->SetMassOfTarget(fIonMass);
    actarSimActionGun->SetExEnergyOfTarget(fIonExciteEnergy);

    //actarSimActionGun->SetTargetIonExcEnergy(fIonExciteEnergy);
  }
}

//////////////////////////////////////////////////////////////////
/// Particular behavior of the scattered ion command. Ion state should be selected.
void ActarSimPrimaryGeneratorMessenger::scatteredIonCommand(G4String newValues){
  G4Tokenizer next( newValues );
  // check argument
  fAtomicNumber = StoI(next());
  fAtomicMass = StoI(next());
  G4String sQ = next();
  if (sQ.isNull()) {
    fIonCharge = fAtomicNumber;
  } else {
    fIonCharge = StoI(sQ);
    sQ = next();
    if (sQ.isNull()) {
      fIonExciteEnergy = 0.0;
    } else {
      fIonExciteEnergy = StoD(sQ) * keV;
    }
  }

  G4Ions* ion;
  ion = (G4Ions*) ionTable->GetIon(fAtomicNumber,
					fAtomicMass,
					fIonExciteEnergy);
  if (ion==0) {
    G4cout << "##################################################################"
	   << G4endl
	   << "###  ActarSimPrimaryGeneratorMessenger::scatteredIonCommand()  ###"
	   << "Ion with Z=" << fAtomicNumber
	   << " A=" << fAtomicMass << "can not be defined" << G4endl;
    G4cout << "##################################################################"
	   << G4endl;
  } else {
    actarSimActionGun->SetScatteredIon(ion);
    actarSimActionGun->SetScatteredIonCharge(fIonCharge*eplus);
    //actarSimActionGun->SetScatteredIonExcEnergy(fIonExciteEnergy);
  }
}

//////////////////////////////////////////////////////////////////
/// Particular behavior of the (KINE) scattered ion command. Ion state should be selected.
void ActarSimPrimaryGeneratorMessenger::KineScatteredIonCommand(G4String newValues){
  // Selection of scattered ion command. Ion state should be selected.
  G4Tokenizer next( newValues );
  // check argument
  fAtomicNumber = StoI(next());
  fAtomicMass = StoI(next());
  fIonCharge = StoI(next());
  fIonExciteEnergy = StoD(next());
  fIonMass=StoD(next());

  G4Ions* ion;
  ion = (G4Ions*) ionTable->GetIon(fAtomicNumber,
					fAtomicMass,
					fIonExciteEnergy);
  if (ion==0) {
    G4cout << "##################################################################"
	   << G4endl
	   << "###  ActarSimPrimaryGeneratorMessenger::KineScatteredIonCommand()  ###"
	   << "Ion with Z=" << fAtomicNumber
	   << " A=" << fAtomicMass << "can not be defined" << G4endl;
    G4cout << "##################################################################"
	   << G4endl;
  } else {
    actarSimActionGun->SetScatteredIon(ion);
    actarSimActionGun->SetScatteredIonCharge(fIonCharge*eplus);
    actarSimActionGun->SetMassOfScattered(fIonMass);
    actarSimActionGun->SetExEnergyOfScattered(fIonExciteEnergy);
    //    G4cout << "ActarSimPrimaryGeneratorMessenger::KineScatteredIonCommand(): excitation energy=" << fIonExciteEnergy << G4endl;
    //actarSimActionGun->SetScatteredIonExcEnergy(fIonExciteEnergy);
  }
}

//////////////////////////////////////////////////////////////////
/// Particular behavior of the recoil ion command. Ion state should be selected.
void ActarSimPrimaryGeneratorMessenger::recoilIonCommand(G4String newValues){
  G4Tokenizer next( newValues );
  // check argument
  fAtomicNumber = StoI(next());
  fAtomicMass = StoI(next());
  G4String sQ = next();
  if (sQ.isNull()) {
    fIonCharge = fAtomicNumber;
  } else {
    fIonCharge = StoI(sQ);
    sQ = next();
    if (sQ.isNull()) {
      fIonExciteEnergy = 0.0;
    } else {
      fIonExciteEnergy = StoD(sQ) * keV;
    }
  }

  G4Ions* ion;
  ion =  (G4Ions*) ionTable->GetIon(fAtomicNumber,
					 fAtomicMass,
					 fIonExciteEnergy);
  if (ion==0) {
    G4cout << "##################################################################"
	   << G4endl
	   << "#####  ActarSimPrimaryGeneratorMessenger::recoilIonCommand()  ####"
	   << "Ion with Z=" << fAtomicNumber
	   << " A=" << fAtomicMass << "can not be defined" << G4endl;
    G4cout << "##################################################################"
	   << G4endl;
  } else {
    actarSimActionGun->SetRecoilIon(ion);
    actarSimActionGun->SetRecoilIonCharge(fIonCharge*eplus);
    //actarSimActionGun->SetRecoilIonExcEnergy(fIonExciteEnergy);
  }
}

//////////////////////////////////////////////////////////////////
/// Particular behavior of the (KINE) recoil ion command. Ion state should be selected.
void ActarSimPrimaryGeneratorMessenger::KineRecoilIonCommand(G4String newValues){
  G4Tokenizer next( newValues );
  // check argument
  fAtomicNumber = StoI(next());
  fAtomicMass = StoI(next());
  fIonCharge = StoI(next());
  fIonExciteEnergy = StoD(next());
  fIonMass=StoD(next());

  G4Ions* ion;
  ion =  (G4Ions*) ionTable->GetIon(fAtomicNumber,
					 fAtomicMass,
					 fIonExciteEnergy);
  if (ion==0) {
    G4cout << "##################################################################"
	   << G4endl
	   << "#####  ActarSimPrimaryGeneratorMessenger::KineRecoilIonCommand()  ####"
	   << "Ion with Z=" << fAtomicNumber
	   << " A=" << fAtomicMass << "can not be defined" << G4endl;
    G4cout << "##################################################################"
	   << G4endl;
  } else {
    actarSimActionGun->SetRecoilIon(ion);
    actarSimActionGun->SetRecoilIonCharge(fIonCharge*eplus);
    actarSimActionGun->SetMassOfRecoiled(fIonMass);
    actarSimActionGun->SetExEnergyOfRecoiled(fIonExciteEnergy);
  }
}

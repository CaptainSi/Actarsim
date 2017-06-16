// - AUTHOR: Hector Alvarez-Pol 11/2004
/******************************************************************
 * Copyright (C) 2005-2016, Hector Alvarez-Pol                     *
 * All rights reserved.                                            *
 *                                                                 *
 * License according to GNU LESSER GPL (see lgpl-3.0.txt).         *
 * For the list of contributors see CREDITS.                       *
 ******************************************************************/
//////////////////////////////////////////////////////////////////
/// \class ActarSimPhysicsList
/// Physics List
/////////////////////////////////////////////////////////////////

#include "ActarSimPhysicsList.hh"
#include "ActarSimPhysicsListMessenger.hh"

#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4EmStandardPhysics_option2.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronInelasticQBBC.hh"
#include "G4IonBinaryCascadePhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4StoppingPhysics.hh"

#include "PhysListEmStandard.hh"
#include "PhysListEmStandardWVI.hh"
#include "PhysListEmStandardSS.hh"
#include "PhysListEmStandardGS.hh"
#include "HadrontherapyIonStandard.hh"

//#include "ActarSimParticlesBuilder.hh"
#include "ActarSimStepLimiterBuilder.hh"

#include "G4UnitsTable.hh"
#include "G4LossTableManager.hh"
#include "G4EmProcessOptions.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4EmConfigurator.hh"

#include "G4IonFluctuations.hh"
#include "G4IonParametrisedLossModel.hh"
#include "G4UniversalFluctuation.hh"

#include "G4BraggIonGasModel.hh"
#include "G4BetheBlochIonGasModel.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4IonPhysics.hh"


//#include "/Applications/ROOT6/root-6.08.06/core/base/inc/TString.h"
#include "TMath.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include "TRandom.h"
#include <cstdio>
#include "TRandom1.h"


//////////////////////////////////////////////////////////////////
/// Constructor. Initializing values
ActarSimPhysicsList::ActarSimPhysicsList():  G4VModularPhysicsList(){

  emBuilderIsRegisted = false;
  stepLimiterIsRegisted = false;
  helIsRegisted = false;
  bicIsRegisted = false;
  ionIsRegisted = false;
  gnucIsRegisted = false;
  gasIsRegisted = false;
  stopIsRegisted = false;
  verbose = 0;
  G4LossTableManager::Instance()->SetVerbose(verbose);
  //G4LossTableManager::Instance();
  defaultCutValue = 1.*mm;
  cutForGamma     = defaultCutValue;
  cutForElectron  = defaultCutValue;
  cutForPositron  = defaultCutValue;

  pMessenger = new ActarSimPhysicsListMessenger(this);

  // EM physics
  //emPhysicsList = new PhysListEmStandard("local");
  //emPhysicsList = new G4EmStandardPhysics(1);

  // Add Physics builders
  RegisterPhysics(new G4DecayPhysics());
  //RegisterPhysics(new ActarSimParticlesBuilder());
  steplimiter = new ActarSimStepLimiterBuilder();
}

//////////////////////////////////////////////////////////////////
/// Destructor. Nothing to do
ActarSimPhysicsList::~ActarSimPhysicsList() {
  delete emPhysicsList;
  delete pMessenger;
}

//////////////////////////////////////////////////////////////////
/// Registering the physics processes
void ActarSimPhysicsList::AddPhysicsList(const G4String& name){
    
    
    //-----Read the file now and save that in a CrossSectionVariable
    //should put here a flag but it doesn't work yet - MBabo
    if(std::ifstream("cross_section.dat")){
       // _CrossSectionINTER_ = ReadCrossSectionFile("cross_section.dat");
    }
    
    
    
  if(verbose > 0) {
    G4cout << "Add Physics <" << name
           << "> emBuilderIsRegisted= " << emBuilderIsRegisted
           << G4endl;
  }
  if ((name == "emstandard") && !emBuilderIsRegisted) {
    RegisterPhysics(new G4EmStandardPhysics(1));
    emBuilderIsRegisted = true;
    G4cout << "ActarSimPhysicsList::AddPhysicsList <" << name << ">" << G4endl;
  } else if (name == "local" && !emBuilderIsRegisted) {
    emPhysicsList = new PhysListEmStandard(name);
    emBuilderIsRegisted = true;
    G4cout << "ActarSimPhysicsList::AddPhysicsList <" << name << ">" << G4endl;
  } else if (name == "emstandard_opt1" && !emBuilderIsRegisted) {
    RegisterPhysics(new G4EmStandardPhysics_option1());
    emBuilderIsRegisted = true;
    G4cout << "ActarSimPhysicsList::AddPhysicsList <" << name << ">" << G4endl;
  } else if (name == "emstandard_opt2" && !emBuilderIsRegisted) {
    RegisterPhysics(new G4EmStandardPhysics_option2());
    emBuilderIsRegisted = true;
    G4cout << "ActarSimPhysicsList::AddPhysicsList <" << name << ">" << G4endl;
  } else if (name == "emstandard_opt3" && !emBuilderIsRegisted) {
    RegisterPhysics(new G4EmStandardPhysics_option3());
    emBuilderIsRegisted = true;
    G4cout << "ActarSimPhysicsList::AddPhysicsList <" << name << ">" << G4endl;
  } else if (name == "emstandard_opt4" && !emBuilderIsRegisted) {
    RegisterPhysics(new G4EmStandardPhysics_option4());
    emBuilderIsRegisted = true;
    G4cout << "ActarSimPhysicsList::AddPhysicsList <" << name << ">" << G4endl;
  } else if (name == "emlivermore" && !emBuilderIsRegisted) {
    RegisterPhysics(new G4EmLivermorePhysics());
    emBuilderIsRegisted = true;
    G4cout << "ActarSimPhysicsList::AddPhysicsList <" << name << ">" << G4endl;
  } else if (name == "empenelope" && !emBuilderIsRegisted) {
    RegisterPhysics(new G4EmPenelopePhysics());
    emBuilderIsRegisted = true;
    G4cout << "ActarSimPhysicsList::AddPhysicsList <" << name << ">" << G4endl;
  } else if (name == "standardSS" && !emBuilderIsRegisted) {
    emPhysicsList = new PhysListEmStandardSS(name);
    emBuilderIsRegisted = true;
    G4cout << "ActarSimPhysicsList::AddPhysicsList <" << name << ">" << G4endl;
  } else if (name == "standardWVI" && !emBuilderIsRegisted) {
    emPhysicsList = new PhysListEmStandardWVI(name);
    emBuilderIsRegisted = true;
    G4cout << "ActarSimPhysicsList::AddPhysicsList <" << name << ">" << G4endl;
  } else if (name == "standardGS" && !emBuilderIsRegisted) {
    emPhysicsList = new PhysListEmStandardGS(name);
    emBuilderIsRegisted = true;
    G4cout << "ActarSimPhysicsList::AddPhysicsList <" << name << ">" << G4endl;
  }

  // Register Low Energy  processes for protons and ions
  // Stopping power parameterisation: ICRU49 (default model)
  // Register Standard processes for protons and ions

  else if (name == "ion-standard") {
    if (ionIsRegisted)
      G4cout << "ActarSimPhysicsList::AddPhysicsList: " << name
	     << " cannot be registered ---- ion List already existing"
	     << G4endl;
    else {
      RegisterPhysics( new HadrontherapyIonStandard(name) );
      ionIsRegisted = true;
      G4cout << "ActarSimPhysicsList::AddPhysicsList <" << name << ">" << G4endl;
    }
  } else if (name == "ionGasModels" && !gasIsRegisted && emBuilderIsRegisted) {
    //AddPhysicsList("emstandard");
    AddIonGasModels();
    gasIsRegisted = true;
    G4cout << "ActarSimPhysicsList::AddPhysicsList <" << name << ">" << G4endl;
  } else if (name == "elastic" && !helIsRegisted && emBuilderIsRegisted) {
    RegisterPhysics(new G4HadronElasticPhysics());
    helIsRegisted = true;
    G4cout << "ActarSimPhysicsList::AddPhysicsList <" << name << ">" << G4endl;
  } else if (name == "binary" && !bicIsRegisted && emBuilderIsRegisted) {
    RegisterPhysics(new G4HadronInelasticQBBC());
    bicIsRegisted = true;
    G4cout << "ActarSimPhysicsList::AddPhysicsList <" << name << ">" << G4endl;
  } else if (name == "binary_ion" && !ionIsRegisted && emBuilderIsRegisted) {
    RegisterPhysics(new G4IonBinaryCascadePhysics());
    ionIsRegisted = true;
    G4cout << "ActarSimPhysicsList::AddPhysicsList <" << name << ">" << G4endl;
  } else if (name == "gamma_nuc" && !gnucIsRegisted && emBuilderIsRegisted) {
    RegisterPhysics(new G4EmExtraPhysics());
    gnucIsRegisted = true;
    G4cout << "ActarSimPhysicsList::AddPhysicsList <" << name << ">" << G4endl;
  } else if (name == "stopping" && !stopIsRegisted && emBuilderIsRegisted) {
    RegisterPhysics(new G4StoppingPhysics());
    stopIsRegisted = true;
    G4cout << "ActarSimPhysicsList::AddPhysicsList <" << name << ">" << G4endl;
  } else if(!emBuilderIsRegisted) {
    G4cout << "PhysicsList::AddPhysicsList <" << name << ">"
           << " fail - EM physics should be registered first " << G4endl;
  } else {
    G4cout << "ActarSimPhysicsList::AddPhysicsList <" << name << ">"
           << " fail - module is already regitered or is unknown " << G4endl;
  }
}

// Bosons
#include "G4ChargedGeantino.hh"
#include "G4Geantino.hh"
#include "G4Gamma.hh"
#include "G4OpticalPhoton.hh"

// leptons
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4NeutrinoE.hh"
#include "G4AntiNeutrinoE.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4NeutrinoMu.hh"
#include "G4AntiNeutrinoMu.hh"

// Hadrons
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"

//////////////////////////////////////////////////////////////////
/// Construct Particles
void ActarSimPhysicsList::ConstructParticle() {

  if(verbose > 0)
    G4cout << "Construct Particles" << G4endl;

  // pseudo-particles
  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();

  // gamma
  G4Gamma::GammaDefinition();

  // optical photon
  G4OpticalPhoton::OpticalPhotonDefinition();

  // leptons
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  G4MuonPlus::MuonPlusDefinition();
  G4MuonMinus::MuonMinusDefinition();

  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();
  G4NeutrinoMu::NeutrinoMuDefinition();
  G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();

  // mesons
  G4MesonConstructor mConstructor;
  mConstructor.ConstructParticle();

  // barions
  G4BaryonConstructor bConstructor;
  bConstructor.ConstructParticle();

  // ions
  G4IonConstructor iConstructor;
  iConstructor.ConstructParticle();

  //G4VModularPhysicsList::ConstructParticle();
}

//////////////////////////////////////////////////////////////////
/// Construct Processes
void ActarSimPhysicsList::ConstructProcess() {

  if(verbose > 0)
    G4cout << "Construct Processes" << G4endl;

  if(!emBuilderIsRegisted) { AddPhysicsList("standard"); }
  if(!emPhysicsList) { G4VModularPhysicsList::ConstructProcess(); }
  else{
    AddTransportation();
    emPhysicsList->ConstructProcess();
  }
  // Define energy interval for loss processes
  G4EmProcessOptions emOptions;
  emOptions.SetMinEnergy(0.1*keV);
  emOptions.SetMaxEnergy(100.*GeV);
  emOptions.SetDEDXBinning(90);
  emOptions.SetLambdaBinning(90);
  //emOptions.SetBuildPreciseRange(false);
  //emOptions.SetApplyCuts(true);
  //emOptions.SetVerbose(0);
}

//////////////////////////////////////////////////////////////////
/// Sets the cut on the physics interaction calculations.
///  "G4VUserPhysicsList::SetCutsWithDefault" method sets
///  the default cut value for all particle types
void ActarSimPhysicsList::SetCuts() {
  SetCutValue(cutForGamma, "gamma");
  SetCutValue(cutForElectron, "e-");
  SetCutValue(cutForPositron, "e+");

  if (verbose>0) DumpCutValuesTable();
}

//////////////////////////////////////////////////////////////////
/// Selecting verbosity
void ActarSimPhysicsList::SetVerbose(G4int val){
  verbose = val;
}

//////////////////////////////////////////////////////////////////
/// Setting cut value for the gammas
void ActarSimPhysicsList::SetCutForGamma(G4double cut){
  cutForGamma = cut;
  SetParticleCuts(cutForGamma, G4Gamma::Gamma());
}

//////////////////////////////////////////////////////////////////
/// Setting cut value for the electron
void ActarSimPhysicsList::SetCutForElectron(G4double cut){
  cutForElectron = cut;
  SetParticleCuts(cutForElectron, G4Electron::Electron());
}

//////////////////////////////////////////////////////////////////
/// Setting cut value for the positron
void ActarSimPhysicsList::SetCutForPositron(G4double cut) {
  cutForPositron = cut;
  SetParticleCuts(cutForPositron, G4Positron::Positron());
}

//////////////////////////////////////////////////////////////////
/// Adds the ion gas model
void ActarSimPhysicsList::AddIonGasModels() {
  G4EmConfigurator* em_config = G4LossTableManager::Instance()->EmConfigurator();

  G4ParticleTable::G4PTblDicIterator* aaParticleIterator = G4ParticleTable::GetParticleTable()->GetIterator();
  aaParticleIterator->reset();
  while ((*aaParticleIterator)()) {
    G4ParticleDefinition* particle = aaParticleIterator->value();
    G4String partname = particle->GetParticleName();
    if(partname == "alpha" || partname == "He3" || partname == "GenericIon") {
      G4BraggIonGasModel* mod1 = new G4BraggIonGasModel();
      G4BetheBlochIonGasModel* mod2 = new G4BetheBlochIonGasModel();
      G4double eth = 2.*MeV*particle->GetPDGMass()/proton_mass_c2;
      em_config->SetExtraEmModel(partname,"ionIoni",mod1,"",0.0,eth,
                                 new G4IonFluctuations());
      em_config->SetExtraEmModel(partname,"ionIoni",mod2,"",eth,100*TeV,
                                 new G4UniversalFluctuation());
    }
  }
}


//////////////////////////////////////////////////////////////////
/// Add angular distribution in transfer reactions as well as energy dependency of total cross section (MB)
//void ActarSimPhysicsList::AddAngularCrossSection(){}

/// Constructor. Initializing values
CrossSectionVariable::CrossSectionVariable(){
    
    int Ei = 3; // 1, 2, X ??
    int Ai = 3;
    
    dXS =(float**) malloc(sizeof(float*) * Ei);
    for (int e=0; e<Ei; e++) dXS[e]=(float*) malloc(sizeof(float) * Ai );
    XAngle =(float*) malloc(sizeof(float*) * Ai);
    XEnergy =(float*) malloc(sizeof(float*) * Ei);
    XSIntegral =(float*) malloc(sizeof(float*) * Ei);
    
    //Initialization
    for (int e=0; e<Ei ; e++)
    {
        XSIntegral[e]=0;
        for (int a=0; a<Ai; a++)
            dXS[e][a]=0;
    }
    for (int a=0; a<Ai; a++)
        XAngle[a]=-1;
    
    //dXS[1][1]={{-1.}};
    //XEnergy[1]={-1.};
    //XAngle[1]={-1.};
    
}

/// Destructor. Nothing to do
CrossSectionVariable::~CrossSectionVariable() {
}

/*

CrossSectionVariable GetCrossSectionVariable(CrossSectionVariable XSV){
    CrossSectionVariable CrossSectionCopy;
    CrossSectionCopy = XSV;
    return CrossSectionCopy;
}
*/

CrossSectionVariable ReadCrossSectionFile(TString name_of_xsfile)
    {
        CrossSectionVariable crosssection_var;
 
        //float** dXS;
        //float* XEnergy;
        //float* XAngle;
        //float *XSIntegral;
        
        //TH1F* histo_d = new TH1F ("histo_d","histo_d", 360, 0, 180);
        //TH1F* histo_r = new TH1F ("histo_r","histo_r", 360, 0, 180);
 
        
        //float var[1]={0};
        //float dumb[1]={0};
        float ienergy[2]={0,0}, iangle[2]={0,0}, isigma[2]={0.0};
        float varI=0.;
        //float theta_dwba= 180.;
        
        //int N_iter=100000;
        int li=0;
        int Ei=1, Ai=-1;
        
        TString name_file = name_of_xsfile;
        
        // format of the file has to be < float_energy float_angle float_sigma >
        // read the file to find the total cross section and the number of lignes and columns
        FILE* finput_i = fopen(name_file, "r");
        
        while (!feof(finput_i))
        {
            fscanf(finput_i, "%f %f %f", &ienergy[0], &iangle[0], &isigma[0]);
            fscanf(finput_i, "%f %f %f", &ienergy[1], &iangle[1], &isigma[1]);
            //G4cout <<iangle[0] << " then : " << iangle[1] << G4endl;
            //G4cout <<isigma[0] <<G4endl;
            if (ienergy[0]!=ienergy[1]) {Ei++;}
            if (iangle[0]<iangle[1]) {Ai++;Ai++;}
            //else {iangle[0]=-1;Ai=0;}
            //cout << var[0] << endl ;
            varI=varI+isigma[0]+isigma[1];
            li++; li++;
            //G4cout << Ai << G4endl;
            
            // well, I don't know why, but I fixed it here
            Ai=180;
        }
        Ei++;
        //G4cout << "Ei=, Ai= " << Ei << "  " << Ai << G4endl;
        FILE* finput_f = fopen(name_file, "r");
        //dXS=new *float[Ei][Ai]
        
        crosssection_var.dXS =(float**) malloc(sizeof(float*) * (Ei+1));
        for (int e=0; e<Ei; e++) crosssection_var.dXS[e]=(float*) malloc(sizeof(float) * (Ai+1) );
        
        crosssection_var.XAngle =(float*) malloc(sizeof(float) * (Ai));
        //G4cout << "dimensions : " << sizeof(crosssection_var.XAngle)/sizeof(crosssection_var.XAngle[0]) << G4endl;
        
        crosssection_var.XEnergy =(float*) malloc(sizeof(float*) * (Ei+1));
        //G4cout << "dimensions : " << sizeof(crosssection_var.XEnergy)/sizeof(crosssection_var.XEnergy[0]) << G4endl;

        
        crosssection_var.XSIntegral =(float*) malloc(sizeof(float*) * (Ei+1));
        
        //Initialization
        for (int e=0; e<Ei ; e++)
        {
            crosssection_var.XSIntegral[e]=0;
            for (int a=0; a<Ai; a++)
            {crosssection_var.dXS[e][a]=0;
            crosssection_var.XAngle[a]=0;}
        }
        
        
        crosssection_var.XEnergy[0]=crosssection_var.XAngle[0]=0;
        
        float varenergy[2], varangle[2], varsigma[2];
        
        //fill the dXS pointer at 2 dimensions
        for (int e=1; e<Ei+1 ; e++)
            for (int a=1; a<Ai+1; a++)
            {
                fscanf(finput_f, "%f %f %f", &varenergy[0], &varangle[0], &varsigma[0]);
                crosssection_var.dXS[e-1][a-1]=varsigma[0];
                crosssection_var.XSIntegral[e-1]=crosssection_var.XSIntegral[e-1]+2.*TMath::Pi()*varsigma[0]*TMath::Sin(varangle[0]/180.*TMath::Pi())*1./180*TMath::Pi();
                if(crosssection_var.XEnergy[e-1]!=varenergy[0]){crosssection_var.XEnergy[e-1]=varenergy[0];}
                if(crosssection_var.XAngle[a-1]<=crosssection_var.XAngle[a]) {crosssection_var.XAngle[a-1]=varangle[0];}
                
            }
        crosssection_var.number_angles = Ai;
        crosssection_var.number_energies = Ei;
        
        
        //  printf ("Random string: %s\n",buffer);
        //free (buffer);
        
        
        //for (int e=0; e<Ei ; e++) G4cout << " Cross sections " << crosssection_var.XSIntegral[e] << G4endl;
       // cout << Ei << " Energies over "  << Ai << " degrees " << endl ;
return crosssection_var;
    
}
    
float SetAngleFromCrossSection(CrossSectionVariable crosssection_var, int Energycase)
{
    //---------------Give an angle of the emitted particle according to a random---------------//
    //------------distribution weigthed by the calculated DWBA total cross section--------------//
    
    //float** idXS
    
    float randangle=-1.;
    int number_of_angles = 180;//sizeof dXS[0] /sizeof(float);

    
    //cout << "In file : " << n << " lines and total integral is " << integral << endl ;
    
    float angle[number_of_angles];
    float sigma[number_of_angles];
    //int i=0;
    
    sigma[0]=crosssection_var.dXS[Energycase][0];
    float integral=sigma[0];
    for(int a=1; a<number_of_angles; a++)
    {
        sigma[a]=sigma[a-1]+crosssection_var.dXS[Energycase][a];
        integral=integral+crosssection_var.dXS[Energycase][a];
    }
    
    // randomization through a weighted distribution
    TRandom1 *start = new TRandom1;
    start->SetSeed(rand());
    double ran = start->TRandom1::Rndm();
    double r = ran*integral; //normalization over the range of angles
    
    for(int k=0; k<number_of_angles; k++)
        if(r<sigma[k+1] && r>sigma[k])
            randangle=crosssection_var.XAngle[k];
    
return (float)randangle; //deg ??
}


float SetEnergybeamFromCrossSection(CrossSectionVariable crosssection_var)
{
    
    //---------------Give an energy beam according to a random distribution---------------//
    //----------------weigthed by the calculated DWBA total cross section-----------------//
    
    float randenergy=-1.;
    int number_of_energies = crosssection_var.number_energies;
    //int number_of_energies = ((sizeof crosssection_var.XEnergy) / (sizeof crosssection_var.XEnergy));
   // G4cout << "number_of_energies " << number_of_energies << G4endl;
    
    float Ttsigma[number_of_energies];
   // int i=0;
    
    Ttsigma[0]=crosssection_var.XSIntegral[0];
    float integral=Ttsigma[0];
    
    G4cout <<"Ttsigma[0] " << Ttsigma[0] <<G4endl;

    for(int e=1; e<number_of_energies; e++)
    {
        Ttsigma[e]=Ttsigma[e-1]+crosssection_var.XSIntegral[e];
        integral=integral+crosssection_var.XSIntegral[e];
        
        // G4cout <<"Ttsigma[" << e << "] " << Ttsigma[e] <<G4endl;
         //G4cout <<"integral " << integral <<G4endl;
        
    }
    
    // randomization through a weighted distribution
    TRandom1 *start = new TRandom1;
    start->SetSeed(rand());
    double ran = start->TRandom1::Rndm();
    double r = ran*integral; //normalization over the range of angles
    
    
    
    randenergy=crosssection_var.XEnergy[0];
    for(int k=0; k<number_of_energies; k++)
        if(r>=Ttsigma[k])
            randenergy=crosssection_var.XEnergy[k+1];

   // G4cout << "r = " << r << G4endl ;
    //G4cout <<"randenergy " << randenergy <<G4endl;
    
return (float)randenergy;
}

float GetIntegralCrossSection(CrossSectionVariable crosssection_var, float Beam_energy)
{
    int e=0;
    while (Beam_energy<crosssection_var.XEnergy[e])
    {e++;}
    
    return crosssection_var.XSIntegral[e];
}



 float GetZpositionVertex(float Beam_energy)
 {
     //--------------this function should be changed for---------------------//
     //--------------every gas mixture and/or incident beam !!!--------------//
     
     float Zpos = (-0.24961*Beam_energy*Beam_energy-3.07658*Beam_energy+31.02802); // Beam Energy in MeV/A
     return Zpos; // en cm
 }




/*

float CrossSectionVariable::DrawAngularCrossSection()
{
        int N=100000000;
        srand((unsigned)time(0));
        TH1F* histo_d = new TH1F ("histo_d","histo_d", 180, 0, 180);
        float angle;
        float** test=read_multi_distribs("cross_section_2.txt");
        
        for (int i=0; i<N; i++)
        {
            angle =  give_angle(test, 1);
            histo_d->Fill(angle);
        }
        new TCanvas;
        histo_d->Draw();
        
        return angle ;
        
}
*/

    //----------------------------------------------------------------------//



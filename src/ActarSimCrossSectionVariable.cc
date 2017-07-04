//
//  ActarSimCrossSectionVariable.cc
//  
//
//  Created by Mathieu BABO on 6/19/17.
//
//

#include "ActarSimCrossSectionVariable.hh"
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
    
    
    for (int e=0; e<Ei ; e++) std::cout << " Cross sections " << crosssection_var.XSIntegral[e] << std::endl;
    // cout << Ei << " Energies over "  << Ai << " degrees " << endl ;
    return crosssection_var;
    
    
    fclose (finput_i);
    fclose (finput_f);
}

float CrossSectionVariable::SetAngleFromCrossSection(int Energycase)
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
    
    sigma[0]=dXS[Energycase][0];
    float integral=sigma[0];
    for(int a=1; a<number_of_angles; a++)
    {
        sigma[a]=sigma[a-1]+dXS[Energycase][a];
        integral=integral+dXS[Energycase][a];
    }
    
    // randomization through a weighted distribution
    TRandom1 *start = new TRandom1;
    start->SetSeed(rand());
    double ran = start->TRandom1::Rndm();
    double r = ran*integral; //normalization over the range of angles
    
    for(int k=0; k<number_of_angles; k++)
        if(r<sigma[k+1] && r>sigma[k])
            randangle=XAngle[k];
    
    return (float)randangle; //deg ??
}


float CrossSectionVariable::SetEnergybeamFromCrossSection()
{
    
    //---------------Give an energy beam according to a random distribution---------------//
    //----------------weigthed by the calculated DWBA total cross section-----------------//
    
    float randenergy=-1.;
    int number_of_energies = number_energies;
    //int number_of_energies = ((sizeof crosssection_var.XEnergy) / (sizeof crosssection_var.XEnergy));
    // G4cout << "number_of_energies " << number_of_energies << G4endl;
    
    float Ttsigma[number_of_energies];
    // int i=0;
    
    Ttsigma[0]=XSIntegral[0];
    float integral=Ttsigma[0];
    
   // G4cout <<"Ttsigma[0] " << Ttsigma[0] <<G4endl;
    
    for(int e=1; e<number_of_energies; e++)
    {
        Ttsigma[e]=Ttsigma[e-1]+XSIntegral[e];
        integral=integral+XSIntegral[e];
        
        // G4cout <<"Ttsigma[" << e << "] " << Ttsigma[e] <<G4endl;
        //G4cout <<"integral " << integral <<G4endl;
        
    }
    
    // randomization through a weighted distribution
    TRandom1 *start = new TRandom1;
    start->SetSeed(rand());
    double ran = start->TRandom1::Rndm();
    double r = ran*integral; //normalization over the range of angles
    
    
    
    randenergy=XEnergy[0];
    for(int k=0; k<number_of_energies; k++)
        if(r>=Ttsigma[k])
            randenergy=XEnergy[k+1];
    
    // G4cout << "r = " << r << G4endl ;
    //G4cout <<"randenergy " << randenergy <<G4endl;
    
    return (float)randenergy;
}

float CrossSectionVariable::GetIntegralCrossSection(float Beam_energy)
{
    int e=0;
    while (Beam_energy<XEnergy[e])
    {e++;}
    
    return XSIntegral[e];
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

//
//  ActarSimCrossSectionVariable.hh
//  
//
//  Created by Mathieu BABO on 6/19/17.
//
//

#ifndef ActarSimCrossSectionVariable_h
#define ActarSimCrossSectionVariable_h 1

#include <stdio.h>
//#include "globals.hh"


#include "/Applications/ROOT6/root-6.08.06/core/base/inc/TString.h"


// MBabo variables

class CrossSectionVariable
{
public:
    float** dXS; // array of the cross section XS[beam_energy][angle]
    float *XEnergy; // array of the energy of the beam
    float *XAngle; // array of the angle with respect to beam direction
    float *XSIntegral; // total cross section as a function of the energy
    int number_energies;
    int number_angles;
    
public:
    CrossSectionVariable();
    ~CrossSectionVariable();
    float SetAngleFromCrossSection(int);
    float DrawAngularCrossSection();
    float GetIntegralCrossSection(float);
    float SetEnergybeamFromCrossSection();
   
};

CrossSectionVariable ReadCrossSectionFile(TString);
float GetZpositionVertex(float);
/*CrossSectionVariable ReadCrossSectionFile(TString);
//void AddAngularCrossSection();
float SetAngleFromCrossSection(CrossSectionVariable, int);
//float** ReadCrossSectionFile(TString);
float DrawAngularCrossSection();
float GetIntegralCrossSection(CrossSectionVariable, float);
float SetEnergybeamFromCrossSection(CrossSectionVariable);
float GetZpositionVertex(float);*/



#endif /* ActarSimCrossSectionVariable_h */

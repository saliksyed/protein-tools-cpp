//
//  EnergyEvaluator.cpp
//  proteintools
//
//  Created by Salik Syed on 10/8/17.
//  Copyright © 2017 N/A. All rights reserved.
//

#include "EnergyEvaluator.hpp"
#include <math.h>
#include <iostream>


float EnergyEvaluator::getEnergy(Eigen::Matrix4Xf * atomParams, Eigen::Matrix4Xf* atoms) {
    return 0.0;
}

float CPUEvaluator::getEnergy(Eigen::Matrix4Xf * atomParams, Eigen::Matrix4Xf* atoms) {
    float energy = 0.0;
    for(size_t i = 0; i < _atomCount - 1; i++) {
        float charge1 = atomParams->col(i).y();
        float sigma1 = atomParams->col(i).z();
        float epsilon1 = atomParams->col(i).w();
        Eigen::Vector3f pos1 = atoms->col(i).head<3>();
        for(size_t j = i + 1; j < _atomCount; j++) {
            float charge2 = atomParams->col(j).y();
            float sigma2 = atomParams->col(j).z();
            float epsilon2 = atomParams->col(j).w();
            Eigen::Vector3f pos2 = atoms->col(j).head<3>();
            Eigen::Vector3f l = pos1 - pos2;
            float rij = l.norm();
            float r0ij = pow(2.0, 1/6.0) * (sigma1 + sigma2) * 0.5;
            float eij = _lj14scale * epsilon1 * epsilon2;
            float e0 = _c14scale;
            float electrostatic = (charge1 * charge2) / (4.0 * M_PI  * e0 * rij);
            float r0ij_rij = (r0ij/rij);
            float lennardJones =  eij * pow(r0ij_rij, 12.0f) - 2.0f * pow(r0ij_rij, 6.0f);
            energy += lennardJones + electrostatic;
        }
    }
    return energy;
}

OpenGLEvaluator::OpenGLEvaluator(unsigned int atomCount, float lj14scale, float c14scale) : EnergyEvaluator(atomCount, lj14scale, c14scale) {
    
}


float OpenGLEvaluator::getEnergy(Eigen::Matrix4Xf * atomParamsIn, Eigen::Matrix4Xf* atomsIn) {
    return 0.0;
}

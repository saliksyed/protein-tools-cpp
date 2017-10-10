//
//  EnergyEvaluator.hpp
//  proteintools
//
//  Created by Salik Syed on 10/8/17.
//  Copyright © 2017 N/A. All rights reserved.
//

#ifndef EnergyEvaluator_hpp
#define EnergyEvaluator_hpp

#include <Eigen/Core>
#include <stdio.h>

class EnergyEvaluator {
public:
    EnergyEvaluator(float lj14scale, float c14scale) : _c14scale(c14scale), _lj14scale(lj14scale) {}
    virtual float getEnergy(Eigen::Matrix4Xf* atomParams, Eigen::Matrix4Xf* atoms, size_t numAtoms);
protected:
    float _c14scale;
    float _lj14scale;
};

class OpenCLEvaluator : EnergyEvaluator {
    float getEnergy(Eigen::Matrix4Xf* atomParams, Eigen::Matrix4Xf* atoms, size_t numAtoms) override;
};

class CPUEvaluator : EnergyEvaluator {
    float getEnergy(Eigen::Matrix4Xf* atomParams, Eigen::Matrix4Xf* atoms, size_t numAtoms) override;
};

#endif /* EnergyEvaluator_hpp */

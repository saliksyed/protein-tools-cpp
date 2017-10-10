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
#include <OpenCL/opencl.h>
#include <stdio.h>

class EnergyEvaluator {
public:
    EnergyEvaluator(unsigned int atomCount, float lj14scale, float c14scale) : _c14scale(c14scale), _lj14scale(lj14scale), _atomCount(atomCount) {}
    virtual float getEnergy(Eigen::Matrix4Xf* atomParams, Eigen::Matrix4Xf* atoms);
protected:
    float _c14scale;
    float _lj14scale;
    unsigned int _atomCount;
};


class CPUEvaluator : EnergyEvaluator {
public:
    CPUEvaluator(unsigned int atomCount, float lj14scale, float c14scale) : EnergyEvaluator(atomCount, lj14scale, c14scale) {}
    float getEnergy(Eigen::Matrix4Xf* atomParams, Eigen::Matrix4Xf* atoms) override;
};

class OpenGLEvaluator : EnergyEvaluator {
public:
    OpenGLEvaluator(unsigned int atomCount, float lj14scale, float c14scale);
    float getEnergy(Eigen::Matrix4Xf* atomParams, Eigen::Matrix4Xf* atoms) override;

};


#endif /* EnergyEvaluator_hpp */

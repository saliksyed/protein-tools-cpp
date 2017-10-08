//
//  Simulator.hpp
//  proteintools
//
//  Created by Salik Syed on 10/5/17.
//  Copyright © 2017 N/A. All rights reserved.
//

#ifndef Simulator_hpp
#define Simulator_hpp

#include <stdio.h>
#include "Forcefield.hpp"
#include "Chain.hpp"
#include "Conformation.hpp"
#include <vector>
#include <Eigen/Core>
#include <assert.h>
using namespace std;

class Simulator {
public:
    Simulator(Chain& chain, ForceField& forcefield);
    ~Simulator();
    float getEnergy();
    void setConformation(Conformation& conformation);
    void getConformation(Conformation& conformation);
    
protected:
    Eigen::Matrix4Xf* _atoms;
    Eigen::Matrix4Xf* _atomsTransformed;
    Eigen::Matrix4Xf* _atomParams;
    
    vector<const Residue*> _residues;
    vector<Eigen::Matrix4f> _matrixStack;
    map<size_t, Eigen::Matrix4f> _torsionTransforms;
    vector<Eigen::Matrix4f> _translationTransforms;
    Conformation* _conformation;
    
};

#endif /* Simulator_hpp */

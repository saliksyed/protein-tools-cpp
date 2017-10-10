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
#include "EnergyEvaluator.hpp"
#include <vector>
#include <Eigen/Core>
#include <assert.h>
using namespace std;

struct AtomInfo {
    AtomType atom;
    ResidueType residue;
    bool geometryValid;
    bool isBackboneParent;
    bool isBackboneChild;
    Eigen::Vector4f position;
};

class Simulator {
public:
    Simulator(Chain& chain, ForceField& forcefield);
    ~Simulator();
    float getEnergy() const;
    void setConformation(Conformation& conformation, bool forceUpdate=false);
    Conformation getConformation() const;
    void getAtoms(vector<AtomInfo> & atoms) const;
    void getBonds(vector<pair<int, int>> & bonds) const;
protected:
    Eigen::Matrix4Xf* _atoms;
    Eigen::Matrix4Xf* _atomsTransformed;
    Eigen::Matrix4Xf* _atomParams;
    size_t _numAtoms;
    EnergyEvaluator* _evaluator;
    vector<const Residue*> _residues;
    vector<Eigen::Matrix4f> _matrixStack;
    vector<Eigen::Matrix4f> _translationTransforms;
    Conformation* _conformation;
    
};

#endif /* Simulator_hpp */

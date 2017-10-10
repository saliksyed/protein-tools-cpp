//
//  Simulator.cpp
//  proteintools
//
//  Created by Salik Syed on 10/5/17.
//  Copyright © 2017 N/A. All rights reserved.
//

#include "Simulator.hpp"
#include "PDBGeometry.hpp"
#include "Residue.hpp"
#include <iostream>
#include <Eigen/Geometry>
using namespace std;

#define TORSION_EPSILON 0.00001f

Simulator::Simulator(Chain& chain, ForceField& forcefield) {
    
    size_t atom_count = 0;
    ResidueTypeIterator it = chain.first_residue();

    while (it != chain.last_residue()) {
        const Residue * r = forcefield.getResidue(*it,
                                            it==chain.first_residue(),
                                            it==chain.last_residue() - 1);
        _residues.push_back(r);
        atom_count += r->numAtoms();
        it++;
    }
    
    // initialize the conformation
    _conformation = new Conformation(_residues.size());
    
    // init torsion params
    for(size_t i = 0; i < _conformation->numTorsionParameters(); i++) {
        _conformation->setTorsion(i, 0.0f);
    }
    
    // Allocate space for original and transformed atom positions
    _atoms = new Eigen::Matrix4Xf(4, atom_count);
    _atomsTransformed = new Eigen::Matrix4Xf(4, atom_count);
    
    // Allocate space for the charge parameters:
    _atomParams = new Eigen::Matrix4Xf(4, atom_count);
    
    it = chain.first_residue();
    size_t curr_atom = 0;

    for(size_t i = 0; i < _residues.size(); i++) {
        const Residue * r = _residues[i];
        ResidueType type = r->geometry_name();
        PDBGeometry & geom = PDBGeometry::get(type);
        for(size_t j = 0; j < r->numAtoms(); j++) {
            pair<AtomName,AtomType> a = r->getAtom(j);
            _atomParams->col(curr_atom)<<
            -1.0f,
            a.second.charge,
            a.second.sigma,
            a.second.epsilon;
            AtomName name(a.first);
            if (geom.hasGeometry(name)) {
                _atomParams->col(curr_atom)<<
                1.0f,
                a.second.charge,
                a.second.sigma,
                a.second.epsilon;
                Eigen::Vector4f pos = geom.position(name);
                _atoms->col(curr_atom)<<pos.x(), pos.y(), pos.z(), 1.0;
                _atomsTransformed->col(curr_atom)<<pos.x(), pos.y(), pos.z(), 1.0;
            } else {
                _atoms->col(curr_atom)<<0.0, 0.0, 0.0, 0.0;
                _atomsTransformed->col(curr_atom)<<0.0, 0.0, 0.0, 0.0;
            }
            curr_atom++;
        }
        if (i != 0) {
            Eigen::Matrix4f m = _residues[i - 1]->getTransformForChild(r);
            _translationTransforms.push_back(m);
            _matrixStack.push_back(m*_matrixStack[_matrixStack.size() -1]);
        } else {
            _translationTransforms.push_back(Eigen::Matrix4f::Identity(4, 4));
            _matrixStack.push_back(Eigen::Matrix4f::Identity(4, 4));
        }
    }
    setConformation(*_conformation, true);
}

void Simulator::setConformation(Conformation &conformation, bool forceUpdate) {
    // Check which parameters have changed
    size_t firstUpdateIdx = numeric_limits<size_t>::max();
    for(size_t i = 0; i < _conformation->numTorsionParameters(); i++) {
        float currentParamValue = _conformation->getTorsion(i);
        float newParamValue = conformation.getTorsion(i);
        
        if (abs(currentParamValue - newParamValue) > TORSION_EPSILON) {
            // update the rotation matrix if non-zero:
            if (abs(newParamValue) > TORSION_EPSILON) {
                Eigen::Matrix4f rot;
                rot.setIdentity();
                rot.block<3,3>(0,0) = Eigen::AngleAxisf(newParamValue,
                                                          _residues[i]->getBondAxis()).matrix();
                rot.rightCols(1) = Eigen::Vector4f(0.0, 0.0, 0.0, 1.0);
                _torsionTransforms[i] = rot;
            } else {
                _torsionTransforms.erase(i);
            }
            firstUpdateIdx = min(firstUpdateIdx, i);
        }
    }
    
    // update the transformed atom data
    size_t curr = 0;
    for(size_t i = 1; i < _residues.size(); i++) {
        size_t numAtoms = _residues[i]->numAtoms();
        if (i >= firstUpdateIdx || forceUpdate) {
            if (_torsionTransforms.find(i) != _torsionTransforms.end()) {
                _matrixStack[i] = _translationTransforms[i] *
                                _torsionTransforms[i] *
                                _matrixStack[i - 1];
            } else {
                _matrixStack[i] = _translationTransforms[i] *
                                _matrixStack[i - 1];
            }
            for(size_t j = 0; j < numAtoms; j++) {
                _atomsTransformed->col(curr + j) << _matrixStack[i] * _atoms->col(curr + j);
            }
        }
        curr += numAtoms;
    }
}

void Simulator::getConformation(Conformation& conformation) const {
    assert(conformation.numTorsionParameters() == _conformation->numTorsionParameters());
    for(size_t i = 0; i < _conformation->numTorsionParameters(); i++) {
        conformation.setTorsion(i, _conformation->getTorsion(i));
    }
}

void Simulator::getAtoms(vector<AtomInfo>& atomList) const{
    size_t idx = 0;
    for(size_t i = 0; i < _residues.size(); i++) {
        const Residue * r = _residues[i];
        ResidueType type = r->name();
        for(size_t j = 0; j < r->numAtoms(); j++){
            pair<AtomName, AtomType> a = r->getAtom(j);
            AtomInfo info;
            AtomName name(a.first);
            PDBGeometry & geom = PDBGeometry::get(r->geometry_name().c_str());
            info.geometryValid = geom.hasGeometry(name);
            info.atom = a.second;
            info.isBackboneChild = j == r->parentAtomIndex();
            info.isBackboneChild = j == r->childAtomIndex();
            info.residue = type;
            info.position = _atomsTransformed->col(idx);
            atomList.push_back(info);
            idx++;
        }
    }
}
void Simulator::getBonds(vector<pair<int, int>> & bondList) const{
    size_t offset = 0;
    for(size_t i = 0; i < _residues.size(); i++) {
        const Residue * r = _residues[i];
        BondIterator a_it = r->first_bond();
        while(a_it != r->last_bond()) {
            size_t a1 = a_it->atom1_idx + offset;
            size_t a2 = a_it->atom2_idx + offset;
            bondList.push_back(pair<int, int>((int)a1, (int)a2));
            a_it++;
        }
        offset += r->numAtoms();
    }
}

float Simulator::getEnergy() const {
    // TODO: write OpenCL kernel using _atomsTransformed.data() and _atomParams.data()
    return 0.0;
}

Simulator::~Simulator() {
    delete _atoms;
    delete _atomsTransformed;
    delete _atomParams;
}
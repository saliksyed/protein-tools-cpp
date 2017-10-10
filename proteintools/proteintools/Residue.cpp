//
//  Residue.cpp
//  proteintools
//
//  Created by Salik Syed on 10/4/17.
//  Copyright © 2017 N/A. All rights reserved.
//

#include "Residue.hpp"
#include "PDBGeometry.hpp"
#include <Eigen/Geometry>

#define DEFAULT_BACKBONE_BOND_LENGTH 1.47f

Residue::Residue() {
    _bondAxis = Eigen::Vector4f(-1.0 * DEFAULT_BACKBONE_BOND_LENGTH, 0.0, 0.0, 1.0);
}

void Residue::setName(ResidueType name) {
    _name = name;
}

ResidueType Residue::geometry_name() const {
    if (_name.length() == 3 && _name != "HIE") {
        return _name;
    } else if (_name == "HIE") {
        return "HIS";
    } else {
        return _name.substr(1, 3);
    }
}

void Residue::addAtom(AtomName& name, AtomType& type) {
    _atoms[name] = type;
    _atomIdxToName[(int)_atomsOrdered.size()] = name;
    _atomsOrdered.push_back(type);
    
}

pair<AtomName, AtomType> Residue::getAtom(size_t idx) const {
    return pair<AtomName, AtomType>(((map<int, AtomName>) _atomIdxToName)[(int)idx], _atomsOrdered[idx]);
}

void Residue::addBond(Bond& b) {
    _bonds.push_back(b);
}

void Residue::addExternalBond(Bond& b) {

    if (_atomsOrdered[b.atom1_idx].element[0] == 'N') {
        _parentAtomIdx = b.atom1_idx;
    } else if (_atomsOrdered[b.atom1_idx].element[0] == 'C') {
        _childAtomIdx = b.atom1_idx;
    }
    _externalBonds.push_back(b);
}


Eigen::Matrix4f Residue::getTransformForChild(const Residue *child) const {
    // find parent (Nitrogen atom) of current residues external bonds
    ResidueType name = geometry_name();
    ResidueType child_name = child->geometry_name();
    PDBGeometry & parent_geom = PDBGeometry::get(name);
    PDBGeometry & child_geom = PDBGeometry::get(child_name);

    float bondLength = DEFAULT_BACKBONE_BOND_LENGTH;
    Eigen::Vector4f bondAxis(-1.0 * bondLength, 0.0, 0.0, 0.0);
    
    AtomName parentAtomName = ((map<int, AtomName>) _atomIdxToName)[_parentAtomIdx];
    AtomName childAtomName = ((map<int, AtomName>)child->_atomIdxToName)[child->_childAtomIdx];
    if (parent_geom.hasGeometry(parentAtomName) &&
        child_geom.hasGeometry(childAtomName)) {
        Eigen::Vector4f parent_pos = parent_geom.position(parentAtomName);
        Eigen::Vector4f child_pos = child_geom.position(childAtomName);
        Eigen::Vector4f t = (parent_pos - child_pos) + bondAxis;
        Eigen::Affine3f transform(Eigen::Translation3f(t.x(), t.y(), t.z()));
        return transform.matrix();
    }

    return Eigen::Matrix4f();
}
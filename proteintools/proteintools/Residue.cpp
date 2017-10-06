//
//  Residue.cpp
//  proteintools
//
//  Created by Salik Syed on 10/4/17.
//  Copyright © 2017 N/A. All rights reserved.
//

#include "Residue.hpp"


Residue::Residue() {}

void Residue::setName(ResidueType name) {
    _name = name;
}

void Residue::addAtom(AtomName& name, AtomType& type) {
    _atoms[name] = type;
}

void Residue::addBond(Bond& b) {
    _bonds.push_back(b);
}


void Residue::addExternalBond(Bond& b) {
    _externalBonds.push_back(b);
}
//
//  Residue.hpp
//  proteintools
//
//  Created by Salik Syed on 10/4/17.
//  Copyright © 2017 N/A. All rights reserved.
//

#ifndef Residue_hpp
#define Residue_hpp

#include <stdio.h>
#include <string>
#include <Eigen/Core>
#include <map>
#include <vector>
#include "Atom.hpp"
#include "PDBGeometry.hpp"
using namespace std;


struct Bond {
    Bond(int a, int b) : atom1_idx(a), atom2_idx(b) {};
    unsigned int atom1_idx;
    unsigned int atom2_idx;
};

typedef string ResidueType;
typedef string AtomName;

class Residue {
public:
    Residue(string name, map<AtomName, AtomType>& atoms, vector<Bond>& bonds, vector<Bond>& externalBonds);
    void setGeometry(PDBGeometry & geometry);
};

#endif /* Residue_hpp */

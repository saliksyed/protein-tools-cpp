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
#include "Atom.hpp"

using namespace std;

enum ResidueType {
    ALA,
    ARG,
    ASN,
    ASP,
    CYS,
    GLN,
    GLU,
    GLY,
    HIS,
    HIS_DEFAULT_ISOMER,
    ILE,
    LEU,
    LYS,
    MET,
    PRO,
    PHE,
    SER,
    THR,
    TRP,
    TYR,
    VAL
};

class ResidueGeometry {
public:
    ResidueGeometry(string & path);
    map<AtomType, Eigen::Vector4d>& getGeometry(ResidueType & residue);
protected:
    map<ResidueType,map<AtomType, Eigen::Vector4d>> _geometry;
};


#endif /* Residue_hpp */

//
//  ForceField.h
//  proteintools
//
//  Created by Salik Syed on 10/4/17.
//  Copyright © 2017 N/A. All rights reserved.
//

#ifndef ForceField_h
#define ForceField_h

#include <string>
#include <map>
#include <Eigen/Core>
#include "Residue.hpp"
#include "Atom.hpp"
#include "Bond.hpp"
using namespace std;

class ForceField {
public:
    ForceField(const char* path);
    
protected:
    double _lj14scale;
    double _c14scale;
    map<string,Eigen::Vector4f> _geometry;
    map<int, AtomType> _atomTypes;
    map<AtomType,Eigen::Vector4f> _chargeParams;
    map<ResidueType,AtomType> _residueAtoms;
    map<ResidueType,Bond> _residueBonds;
    map<ResidueType,Bond> _residueExternalBonds;
    
    void parse(const char * path);

};

#endif /* ForceField_h */

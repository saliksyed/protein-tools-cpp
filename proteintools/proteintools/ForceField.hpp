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
using namespace std;

class ForceField {
public:
    ForceField(const char* path);
    Residue* getResidue(ResidueType r, bool isChainStart=false, bool isChainEnd=false);
    Residue* getResidue(const char* r, bool isChainStart=false, bool isChainEnd=false);
protected:
    double _lj14scale;
    double _c14scale;
    map<int, AtomType> _atomTypes;
    map<ResidueType,Residue*> _residues;
    
    void parse(const char * path);

};

#endif /* ForceField_h */

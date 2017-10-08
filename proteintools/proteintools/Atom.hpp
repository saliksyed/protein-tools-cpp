//
//  Atom.hpp
//  proteintools
//
//  Created by Salik Syed on 10/4/17.
//  Copyright © 2017 N/A. All rights reserved.
//

#ifndef Atom_hpp
#define Atom_hpp

#include <stdio.h>
#include <Eigen/Core>
using namespace std;

#define ATOM_NAME_MAX_CHAR 8

typedef string AtomName;

struct AtomType {
    AtomType() : charge(0.0), sigma(0.0), epsilon(1.0) {}
    int ID;
    char element[ATOM_NAME_MAX_CHAR];
    char atomClass[ATOM_NAME_MAX_CHAR];
    float mass;
    float charge;
    float sigma;
    float epsilon;
};


#endif /* Atom_hpp */

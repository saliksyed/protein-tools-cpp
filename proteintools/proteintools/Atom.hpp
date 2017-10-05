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

struct AtomType {
    int name;
    char element[ATOM_NAME_MAX_CHAR];
    char atomClass[ATOM_NAME_MAX_CHAR];
    double mass;
    
};

class Atom {
public:
    Atom(AtomType & definition);
    Eigen::Vector4f getPosition();
    AtomType getType();
protected:
    AtomType& _type;
    Eigen::Vector4f _position;
};

#endif /* Atom_hpp */

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

using namespace std;

class PDBGeometry;

struct Bond {
    Bond(int a, int b) : atom1_idx(a), atom2_idx(b) {};
    unsigned int atom1_idx;
    unsigned int atom2_idx;
};

typedef string ResidueType;
typedef typename map<AtomName,AtomType>::iterator AtomIterator;
typedef typename vector<Bond>::iterator BondIterator;


class Residue {
public:
    Residue();
    ResidueType name() { return _name; }
    ResidueType geometry_name();
    void setName(ResidueType name);
    void addAtom(AtomName& name, AtomType& type, int idx);
    void addBond(Bond& b);
    void addExternalBond(Bond& b);

    int childAtomIndex() { return _childAtomIdx; };
    int parentAtomIndex() { return _parentAtomIdx; };
    Eigen::Vector3f getBondAxis() { return _bondAxis; };
    Eigen::Matrix4f getTransformForChild(Residue *child);

    size_t numAtoms() { return _atoms.size(); }
    AtomIterator first_atom() { return _atoms.begin(); }
    AtomIterator last_atom() { return _atoms.end(); }
    BondIterator first_bond() { return _bonds.begin(); }
    BondIterator last_bond() { return _bonds.end(); }
    BondIterator first_external_bond() { return _externalBonds.begin(); }
    BondIterator last_external_bond() { return _externalBonds.end(); }

protected:
    int _childAtomIdx;
    int _parentAtomIdx;
    ResidueType _name;
    Eigen::Vector3f _bondAxis;
    map<AtomName, AtomType> _atoms;
    map<int, AtomName> _atomIdxToName;
    vector<Bond> _bonds;
    vector<Bond> _externalBonds;
    
};

#endif /* Residue_hpp */

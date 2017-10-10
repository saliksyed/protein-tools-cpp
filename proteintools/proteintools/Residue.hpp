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
    Bond(unsigned int a, unsigned int b) : atom1_idx(a), atom2_idx(b) {};
    unsigned int atom1_idx;
    unsigned int atom2_idx;
};


typedef string ResidueType;
typedef typename vector<Bond>::const_iterator BondIterator;


class Residue {
public:
    Residue();
    ResidueType name() const { return _name; }
    ResidueType geometry_name() const;
    void setName(ResidueType name);
    void addAtom(AtomName& name, AtomType& type);
    void addBond(Bond& b);
    void addExternalBond(Bond& b);

    int childAtomIndex() const{ return _childAtomIdx; };
    int parentAtomIndex() const{ return _parentAtomIdx; };
    Eigen::Vector4f getBondAxis() const{ return _bondAxis; };
    Eigen::Matrix4f getTransformForChild(const Residue *child) const;

    size_t numAtoms() const { return _atoms.size(); }
    pair<AtomName, AtomType> getAtom(size_t idx) const;
    BondIterator first_bond() const { return _bonds.cbegin(); }
    BondIterator last_bond() const { return _bonds.cend(); }
    BondIterator first_external_bond() const { return _externalBonds.cbegin(); }
    BondIterator last_external_bond() const { return _externalBonds.cend(); }

protected:
    int _childAtomIdx;
    int _parentAtomIdx;
    ResidueType _name;
    Eigen::Vector4f _bondAxis;
    map<AtomName, AtomType> _atoms;
    map<int, AtomName> _atomIdxToName;
    vector<AtomType> _atomsOrdered;
    vector<Bond> _bonds;
    vector<Bond> _externalBonds;
    
};

#endif /* Residue_hpp */

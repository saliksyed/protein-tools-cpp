//
//  Chain.hpp
//  proteintools
//
//  Created by Salik Syed on 10/4/17.
//  Copyright © 2017 N/A. All rights reserved.
//

#ifndef Chain_hpp
#define Chain_hpp

#include <stdio.h>


#include <Eigen/Core>
#include "Residue.hpp"
#include "Atom.hpp"
#include <vector>
using namespace std;

typedef typename vector<ResidueType>::iterator ResidueTypeIterator;

class Chain {
public:
    Chain(const char* sequence);
    ResidueTypeIterator first_residue() { return _residues.begin(); }
    ResidueTypeIterator last_residue() { return _residues.end(); }
protected:
    static ResidueType symbolToResidue(char symbol);
    vector<ResidueType> _residues;
};

#endif /* Chain_hpp */

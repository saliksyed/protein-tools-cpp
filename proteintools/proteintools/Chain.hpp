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

class Chain {
public:
    Chain(string& sequence);
    void writeToPDB(string& path);
    void setConformation();
    void getConformation();
    double getEnergy();
protected:
    vector<ResidueType> _residues;
    vector<Eigen::Matrix4Xd> transforms;
    float * _atoms;
};

#endif /* Chain_hpp */

//
//  PDBGeometry.hpp
//  proteintools
//
//  Created by Salik Syed on 10/4/17.
//  Copyright © 2017 N/A. All rights reserved.
//

#ifndef PDBGeometry_hpp
#define PDBGeometry_hpp

#include <stdio.h>
#include <Eigen/Core>
#include "Atom.hpp"
#include "Residue.hpp"
#include <map>
using namespace std;

class PDBGeometry {
public:
    static void load(const char * folder);
    static PDBGeometry& get(ResidueType& res);
    static PDBGeometry& get(const char* res);
    PDBGeometry(ResidueType r, const char * path);
    ResidueType name();
    bool hasGeometry(AtomName& name);
    Eigen::Vector4d position(AtomName& name);
protected:
    ResidueType _residueName;
    map<AtomName, Eigen::Vector4d> _geometry;
    void parse(const char * path);
private:
    static map<ResidueType, PDBGeometry*> _data;
    static const char* residueTypes[];
};

#endif /* PDBGeometry_hpp */

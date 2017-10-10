//
//  PDBGeometry.cpp
//  proteintools
//
//  Created by Salik Syed on 10/4/17.
//  Copyright © 2017 N/A. All rights reserved.
//

#include "PDBGeometry.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <iterator>
#include <vector>
#include <string>
#include <stdio.h>
#include <stdlib.h>
using namespace std;

#define ATOM_POS_DELIMITER "ATOM"
#define NUM_RESIDUES 20

const char* PDBGeometry::residueTypes[] = {
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL"
};

map<ResidueType, PDBGeometry*> PDBGeometry::_data = map<ResidueType, PDBGeometry*>();

void PDBGeometry::load(const char * folder) {
    if(PDBGeometry::_data.size() > 0) {
        return;
    } else {
        for(int i = 0; i <NUM_RESIDUES; i++) {
            const char * residueType = PDBGeometry::residueTypes[i];
            char path[1024];
            sprintf(path, "%s/%s.pdb", folder, residueType);
            PDBGeometry::_data[ResidueType(residueType)] = new PDBGeometry(residueType, path);
        }
    }
}

PDBGeometry& PDBGeometry::get(ResidueType& res) {
    return *PDBGeometry::_data[res];
}


PDBGeometry& PDBGeometry::get(const char* res) {
    return *PDBGeometry::_data[ResidueType(res)];
}

PDBGeometry::PDBGeometry(ResidueType name, const char * path): _residueName(name) {
    parse(path);
}

string PDBGeometry::name() {
    return _residueName;
}

void PDBGeometry::parse(const char *path) {
    ifstream f;
    f.open(path);
    if (f.is_open()) {
        string line;
        while (!f.eof()) {
            getline(f, line);
            if (line.length() == 0) break;
            istringstream iss(line);
            vector<string> tokens{istream_iterator<string>{iss},
                istream_iterator<string>{}};
            if (tokens[0] == ATOM_POS_DELIMITER) {
                AtomName name = tokens[2];
                double x = strtod(tokens[5].c_str(), NULL);
                double y = strtod(tokens[6].c_str(), NULL);
                double z = strtod(tokens[7].c_str(), NULL);
                Eigen::Vector4f vec(x, y, z, 1.0);
                _geometry[name] = vec;
            }
        }
    }
    f.close();
}

bool PDBGeometry::hasGeometry(AtomName& name) {
    return _geometry.find(name) != _geometry.end();
}

Eigen::Vector4f PDBGeometry::position(AtomName& name) {
    return _geometry[name];
}


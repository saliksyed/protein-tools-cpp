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
#include <string>
using namespace std;

PDBGeometry::PDBGeometry(const char * path) {
    parse(path);
}

void PDBGeometry::parse(const char *path) {
    ifstream f;
    f.open(path);
    if (f.is_open()) {
        string line;
        while (!f.eof()) {
            getline(f, line);
            cout<<line<<endl;
        }
    }
    f.close();
}


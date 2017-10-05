//
//  ForceField.cpp
//  proteintools
//
//  Created by Salik Syed on 10/4/17.
//  Copyright © 2017 N/A. All rights reserved.
//

#include <stdio.h>
#include <iostream>
#include "ForceField.hpp"
#include "tinyxml2.h"
using namespace std;
using namespace tinyxml2;

ForceField::ForceField(const char * path) {
    parse(path);
}

void ForceField::parse(const char * path) {
    XMLDocument doc;
    doc.LoadFile(path);
    
    // parse atom types
    XMLElement * curr = doc.FirstChildElement("ForceField")->FirstChildElement("AtomTypes")->FirstChildElement("Type");
    while (curr) {
        AtomType a;
        strcpy(a.atomClass, curr->Attribute("class"));
        strcpy(a.element,curr->Attribute("element"));
        curr->QueryIntAttribute("name", &a.name);
        curr->QueryDoubleAttribute("mass", &a.mass);
        curr = curr->NextSiblingElement();
        _atomTypes.insert(pair<int, AtomType>(a.name, a));
        
    }
}
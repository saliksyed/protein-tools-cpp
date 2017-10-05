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
        _atomTypes.insert(pair<int, AtomType>(a.name, a));
        curr = curr->NextSiblingElement();
    }
    
    // parse charge params
    curr = doc.FirstChildElement("ForceField")->FirstChildElement("NonbondedForce");
    curr->QueryDoubleAttribute("coulomb14scale", &this->_c14scale);
    curr->QueryDoubleAttribute("lj14scale", &this->_lj14scale);
    curr = curr->FirstChildElement("Atom");
    while (curr) {
        int key;
        curr->QueryIntAttribute("type", &key);
        map<int,AtomType>::iterator t = _atomTypes.find(key);
        
        if (t != _atomTypes.end()) {
            curr->QueryDoubleAttribute("charge", &t->second.charge);
            curr->QueryDoubleAttribute("sigma", &t->second.sigma);
            curr->QueryDoubleAttribute("epsilon", &t->second.epsilon);
        }
        curr = curr->NextSiblingElement();
    }
    
    // parse the residues
    curr = doc.FirstChildElement("ForceField")->FirstChildElement("Residues")->FirstChildElement("Residue");;
    while (curr) {
        XMLElement* atomsElem = curr->FirstChildElement("Atom");
        XMLElement* bondsElem = curr->FirstChildElement("Bond");
        XMLElement* externalBondsElem = curr->FirstChildElement("ExternalBond");
        map<AtomName, AtomType> atoms;
        vector<Bond> bonds;
        vector<Bond> externalBonds;
        while (atomsElem) {
            int type;
            AtomName name = AtomName(curr->Attribute("name"));
            curr->QueryIntAttribute("type", &type);
            map<int,AtomType>::iterator t = _atomTypes.find(type);
            if (t!= _atomTypes.end()) {
                atoms[name] = _atomTypes[type];
            }
            atomsElem = atomsElem->NextSiblingElement();
        }
        while (bondsElem) {
            int a, b;
            bondsElem->QueryIntAttribute("from", &a);
            bondsElem->QueryIntAttribute("to", &b);
            bonds.push_back(Bond(a,b));
            bondsElem = bondsElem->NextSiblingElement();
        }
        while (externalBondsElem) {
            int a;
            externalBondsElem->QueryIntAttribute("from", &a);
            externalBonds.push_back(Bond(a,-1));
            externalBondsElem = externalBondsElem->NextSiblingElement();
        }
        curr = curr->NextSiblingElement();
    }
}

Residue * ForceField::createResidue(ResidueType r, bool isChainStart=false, bool isChainEnd=false) {
    return NULL;
}

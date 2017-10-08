//
//  Chain.cpp
//  proteintools
//
//  Created by Salik Syed on 10/4/17.
//  Copyright © 2017 N/A. All rights reserved.
//

#include "Chain.hpp"

ResidueType Chain::symbolToResidue(char symbol) {
    switch(symbol) {
        case 'A':
            return ResidueType("ALA");
            break;
        case 'R':
            return ResidueType("ARG");
            break;
        case 'N':
            return ResidueType("ASN");
            break;
        case 'D':
            return ResidueType("ASP");
            break;
        case 'C':
            return ResidueType("CYS");
            break;
        case 'Q':
            return ResidueType("GLN");
            break;
        case 'E':
            return ResidueType("GLU");
            break;
        case 'G':
            return ResidueType("GLY");
            break;
        case 'H':
            return ResidueType("HIS");
            break;
        case 'I':
            return ResidueType("ILE");
            break;
        case 'L':
            return ResidueType("LEU");
            break;
        case 'K':
            return ResidueType("LYS");
            break;
        case 'M':
            return ResidueType("MET");
            break;
        case 'P':
            return ResidueType("PRO");
            break;
        case 'F':
            return ResidueType("PHE");
            break;
        case 'S':
            return ResidueType("SER");
            break;
        case 'T':
            return ResidueType("THR");
            break;
        case 'W':
            return ResidueType("TRP");
            break;
        case 'Y':
            return ResidueType("TYR"); 
            break;
        case 'V':
            return ResidueType("VAL"); 
            break;
        default:
            return NULL;
    }
}

Chain::Chain(const char* sequence) {
    while(*sequence) {
        _residues.push_back(Chain::symbolToResidue(*sequence));
        sequence++;
    }
}

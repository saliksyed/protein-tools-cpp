//
//  Conformation.cpp
//  proteintools
//
//  Created by Salik Syed on 10/8/17.
//  Copyright © 2017 N/A. All rights reserved.
//

#include "Conformation.hpp"


Conformation::Conformation(size_t sz) : _torsionAngles(sz) {}

void Conformation::setTorsion(size_t idx, float angle) {
    _torsionAngles[idx] = angle;
}

float Conformation::getTorsion(size_t idx) {
    return _torsionAngles[idx];
}

size_t Conformation::numTorsionParameters() {
    return _torsionAngles.size();
}
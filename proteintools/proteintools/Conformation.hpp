//
//  Conformation.hpp
//  proteintools
//
//  Created by Salik Syed on 10/8/17.
//  Copyright © 2017 N/A. All rights reserved.
//

#ifndef Conformation_hpp
#define Conformation_hpp

#include <stdio.h>
#include <vector>
using namespace std;

class Conformation {
public:
    Conformation(size_t sz);
    void setTorsion(size_t idx, float angle);
    float getTorsion(size_t idx);
    size_t numTorsionParameters();
private:
    vector<float> _torsionAngles;
};

#endif /* Conformation_hpp */

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

class PDBGeometry {
public:
    PDBGeometry(const char * path);
protected:
    parse(const char * path);
};

#endif /* PDBGeometry_hpp */

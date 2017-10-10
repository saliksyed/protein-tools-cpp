//
//  main.cpp
//  proteintools
//
//  Created by Salik Syed on 10/4/17.
//  Copyright © 2017 N/A. All rights reserved.
//

#include "Server.hpp"
#include "ProteinTools.h"
#include <iostream>
using namespace std;


int main(int argc, char** argv)
{
    ForceField f = ForceField("/Users/saliksyed/Desktop/amber99sb.xml");
    PDBGeometry::load("/Users/saliksyed/src/protein-tools-cpp/proteintools/proteintools/data/v3PDB");
    Chain c("CQCQQQ");
    Simulator sim(c, f);
    Server* server = new Server(sim);
    server->run(1283);
    return 0;
}

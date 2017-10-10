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
#include <time.h>
using namespace std;


int main(int argc, char** argv)
{
    ForceField f = ForceField("/Users/saliksyed/Desktop/amber99sb.xml");
    PDBGeometry::load("/Users/saliksyed/src/protein-tools-cpp/proteintools/proteintools/data/v3PDB");
    Chain c("LSDEDFKAVFGMTRSAFANLPLWKQQHLKKEKGLF");
    Simulator sim(c, f);
//    printf ("Starting timer!\n");
//    clock_t t = clock();
//    sim.getEnergy();
//    t = clock() - t;
//    printf ("It took me %f seconds.\n",((float)t)/CLOCKS_PER_SEC);

    Conformation conf = sim.getConformation();
    conf.setTorsion(3, 0.2);
    sim.setConformation(conf);
    Server* server = new Server(sim);
    server->run(1283);
    return 0;
}

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
#include <stdlib.h>
using namespace std;

Simulator * currSim;
Conformation * conformation;

void sampleConformation() {
    while(true) {
        int param = rand() % conformation->numTorsionParameters();
        float prev = conformation->getTorsion(param);
        conformation->setTorsion(param, prev + 0.0001);
        currSim->setConformation(*conformation);
    }
}

void runServer() {
    // start server to visualize results
    Server* server = new Server(*currSim);
    server->run(1283);
}

int main(int argc, char** argv)
{
    
    // load forcefield + geometry
    ForceField f = ForceField("/Users/saliksyed/Desktop/amber99sb.xml");
    PDBGeometry::load("/Users/saliksyed/src/protein-tools-cpp/proteintools/proteintools/data/v3PDB");

    // create chain + sim
    Chain c("LSDEDFKAVFGMTRSAFANLPLWKQQHLKKEKGLF");
    Simulator sim(c, f);

    Conformation conf = sim.getConformation();
    
    currSim = &sim;
    conformation = &conf;

    
    // run thread to sample conformations
    thread sample(sampleConformation);
    
    // run thread to serve results
    thread server(runServer);
    
    sample.join();
    server.join();

    return 0;
}

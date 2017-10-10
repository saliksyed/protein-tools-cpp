//
//  Server.hpp
//  proteintools
//
//  Created by Salik Syed on 10/8/17.
//  Copyright © 2017 N/A. All rights reserved.
//

#ifndef Server_hpp
#define Server_hpp
#include <uWS/uWS.h>
#include <stdio.h>
#include "Simulator.hpp"
#include <string>
#include <thread>
using namespace std;

class Server {
public:
    Server(Simulator& sim);
    void run(int port);
protected:
    string serializeTopology(vector<AtomInfo> & atomInfo, vector<pair<int, int>> & bonds) const;
    string serializeState(vector<AtomInfo> & atomInfo) const;
    Simulator& _sim;
};

#endif /* Server_hpp */

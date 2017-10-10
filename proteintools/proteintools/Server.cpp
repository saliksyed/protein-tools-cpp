//
//  Server.cpp
//  proteintools
//
//  Created by Salik Syed on 10/8/17.
//  Copyright © 2017 N/A. All rights reserved.
//

#include "Server.hpp"
#include "json11.hpp"
#include <iostream>
using namespace json11;


string Server::serializeTopology(vector<AtomInfo> & atomInfo, vector<pair<int, int>> & bonds) const {
    vector<Json> items;
    for(AtomInfo & atom : atomInfo) {
        Json item  = Json::object {
            { "element", atom.atom.element },
            { "atomClass", atom.atom.atomClass },
            { "residueId", atom.residue },
            { "isParent", atom.isBackboneParent },
            { "isChild", atom.isBackboneChild },
            { "geometryValid", atom.geometryValid },
            { "x", atom.position.x() },
            { "y", atom.position.y() },
            { "z", atom.position.z() }
        };
        items.push_back(item);
    }
    
    vector<Json> bondItems;
    
    for(pair<int, int> & bond : bonds) {
        bondItems.push_back(Json::array{bond.first, bond.second});
    }
    
    Json data = Json::object {
        {"atoms", items},
        {"bonds", bondItems}
    };

    return Json(data).dump();
}

string Server::serializeState(vector<AtomInfo> & atomInfo) const {
    vector<Json> items;
    for(AtomInfo & atom : atomInfo) {
        Json item  = Json::object {
            { "x", atom.position.x() },
            { "y", atom.position.y() },
            { "z", atom.position.z() }
        };
        items.push_back(item);
    }
    Json data = Json::object {
        {"atoms", items},
    };
    return Json(data).dump();
}

Server::Server(Simulator& sim): _sim(sim) {}

void Server::run(int port) {
    uWS::Hub _uwsHub;
    _uwsHub.onMessage([this](uWS::WebSocket<uWS::SERVER> *ws, char *message, size_t length, uWS::OpCode opCode) {
        if (message[0] == 'T') {
            vector<AtomInfo> info;
            this->_sim.getAtoms(info);
            vector<pair<int, int>> bondInfo;
            this->_sim.getBonds(bondInfo);
            string ret = serializeTopology(info, bondInfo);
            ws->send(ret.c_str(), ret.length(), uWS::OpCode::TEXT);
        } else if (message[0] == 'S') {
            vector<AtomInfo> info;
            this->_sim.getAtoms(info);
            string ret = serializeState(info);
            ws->send(ret.c_str(), ret.length(), uWS::OpCode::TEXT);
        } else {
            vector<AtomInfo> info;
            this->_sim.getAtoms(info);
            vector<pair<int, int>> bonds;
            this->_sim.getBonds(bonds);
            string ret = serializeTopology(info, bonds);
            ws->send(ret.c_str(), ret.length(), uWS::OpCode::TEXT);
        }
    });
    _uwsHub.listen(port);
    _uwsHub.run();
}
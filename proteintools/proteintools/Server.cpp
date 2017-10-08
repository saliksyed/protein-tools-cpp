//
//  Server.cpp
//  proteintools
//
//  Created by Salik Syed on 10/8/17.
//  Copyright © 2017 N/A. All rights reserved.
//

#include "Server.hpp"


Server::Server(int port) {

    _uwsHub.onMessage([](uWS::WebSocket<uWS::SERVER> *ws, char *message, size_t length, uWS::OpCode opCode) {
        if (strcmp(message, "topology")==0) {
            // TODO: read topology from simulator, serialize it, send it
            ws->send(message, length, uWS::OpCode::BINARY);
        } else if (strcmp(message, "state")==0) {
            // TODO: read the state for each atom from the simulator, serialize it, send it.
            ws->send(message, length, uWS::OpCode::BINARY);
        }
    });
    
    if (_uwsHub.listen(port)) {
        _uwsHub.run();
    }
}
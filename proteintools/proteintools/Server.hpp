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

class Server {
public:
    Server(int port);
    
protected:
    uWS::Hub _uwsHub;
};

#endif /* Server_hpp */

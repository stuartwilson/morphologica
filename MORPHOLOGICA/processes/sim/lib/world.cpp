//
//  world.cpp
//

#include "world.h"

world::world(const char* processName,
             const char* logfileLocation,
             int seed,
             int portID,
             double dt){
    
    this->processName=processName;      // process name
    std::stringstream ss;       // logfile location
    ss<<logfileLocation;                // logfile location
    srand(seed);               // random seed
    this->portID=portID;           // tcpip port ID
    master.init(portID);
    vector<Client> ports;       // Remember this can be used for multi inputs
    
    TIME = 0;
    this->dt = dt;
    
    logfile.open(ss.str().c_str(),ios::out|ios::app);
    ss.clear();
    time_t timer = time(NULL);
    logfile<<"*********"<<endl;
    logfile<<"   HI!"<<endl;
    logfile<<"*********"<<endl;
    logfile<<"Time now: "<<timer<<endl;
    logfile<<"Sim name: "<<processName<<endl;
    logfile<<"**********"<<endl<<flush;
};

//world::world(const world &obj){

//};


vector<string> world::getCommand(vector<double*>msgOut){
    
    stringstream out;
    out.clear();
    out.setf(ios::fixed,ios::floatfield);
    
    for(int i=0;i<msgOut.size();i++){
        out<<*msgOut[i]<<",";
    }
    
    vector <string> command;
    string messageI=master.exchange(out.str().c_str());
    stringstream ss(messageI);
    while (ss.good()){
        string substr;
        getline(ss,substr,',');
        command.push_back(substr);
    } ss.clear();
    
    return command;
};

const char* world::timeStamp(void){
    const char* TIMEcs;
    std::stringstream TIMEss;
    TIMEss<<setw(10)<<setfill('0')<<TIME;
    TIMEcs = TIMEss.str().c_str();
    return TIMEcs;
}

//const char* TIMEcs;
//std::stringstream TIMEss;
//TIMEss<<setw(10)<<setfill('0')<<TIME;
//TIMEcs = TIMEss.str().c_str();


/*
 vector<string> world::getCommand(void){
 
 //const char* TIMEcs;
 //std::stringstream TIMEss;
 //TIMEss<<setw(10)<<setfill('0')<<TIME;
 //TIMEcs = TIMEss.str().c_str();
 
 stringstream out;
 out.clear();
 out.setf(ios::fixed,ios::floatfield);
 
 // DEFINE OUTPUT MESSAGE
 //for(int i=0;i<f.size();i++){out<<*f[i]<<",";}
 
 const char * outputMesg;
 outputMesg = "empty";//out.str().c_str();
 
 vector <string> command;
 //command.clear();
 string messageI=master.exchange(outputMesg);
 stringstream ss(messageI);
 while (ss.good()){
 string substr;
 getline(ss,substr,',');
 command.push_back(substr);
 } ss.clear();
 
 
 logfile<<"input: ";
 for(int i=0;i<command.size();i++){
 logfile<<stoi(command[i])<<',';
 }logfile<<endl;
 
 logfile<<"output: "<<outputMesg<<endl;
 
 
 return command;
 };
 */

/*
 vector<vector<string> > world::getPortComms(void){
 
 vector<vector<string> > commands(ports.size());
 
 for(int i=0;i<ports.size();i++){
 string messageI=ports[i].exchange(outputMessage);
 stringstream ss(messageI);
 while (ss.good()){
 string substr;
 getline(ss,substr,',');
 commands[i].push_back(substr);
 } ss.clear();
 }
 return commands;
 };
 */

world::~world(){
    
    logfile<<"*********"<<endl;
    logfile<<"   FIN"<<endl;
    logfile<<"*********"<<endl<<flush;
    
    logfile.close();
    //    master.closeSocket();
    master.~Client();
    
    for(int i=0;i<ports.size();i++){
        //        ports[i].closeSocket();
        ports[i].~Client();
    }
    
};


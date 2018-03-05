#include "lib/world.h"
#include "lib/sockserve.h"
#include "lib/display.h"
#include "lib/tools.h"
#include "animats/basic.h"

using namespace basic;

#include <iostream>
#include <fstream>
#include <vector>
#include <armadillo>
#include <numeric>
#include <iomanip>

#include "animats/animatTools.h"

using namespace std;
using namespace arma;
using namespace tools;

vector <vector <int> > collisionDetection(std::vector<animat>&,double);
void boundary(animat&, double, int, double, int);
void boundaryCircle(animat&, double, bool, double, int);
void drawAnimat(Gdisplay&, animat, bool);
void calculateCollision(vector<vec>&, vector<double>&, vector<bool>&, double);
void drawEnvironment(Gdisplay&);
void shuffle(std::vector<animat>&);
vector<vector<int> > graph(vector<vector<int> >, int);

int main(int argc, char **argv) {
    
    srand(atoi(argv[3]));
    
    
    double dt = .01;            // integration timestep
    
    // INITIALIZATION
    world W(argv[1],argv[2],atoi(argv[3]),atoi(argv[4]),dt);
    
    // DISPLAYS
    vector<Gdisplay> displays;
    vector<double>fix(3,0.);
    vector<double>eye(3,0.);
    vector<double>rot(3,0.);
    displays.push_back(Gdisplay(600,"arena view",2.,0.,0.));
    displays[0].resetDisplay(fix,eye,rot);
    displays[0].redrawDisplay();
    
    displays.push_back(Gdisplay(600,"pup view",7.,0.,0.));
    displays[1].resetDisplay(fix,eye,rot);
    displays[1].redrawDisplay();
    
    
    vector<animat> A(12,basic::getAnimat(7,&W));
    shuffle(A);
    
    
    double Ta = 0.5, Ta_init=Ta;                // ambient temperature
    double boundX = 2.,boundX_init=boundX;
    double boundY = 2.,boundY_init=boundY;
    
    double floorH = -0.4;
    double mint = 0.;
    double maxt = 1.;
    double oneOverMaxtMinusMint = 1./(maxt-mint);
    
    unsigned int frameN = 0;
    unsigned int F = 0; // focal pup
    
    vector <double> red(3,0.); red[0]=1.;
    vector <double> gre(3,0.); gre[1]=1.;
    vector <double> blu(3,0.); blu[2]=1.;
    
    double TIME=0.;
    
    vector<vector<int> > groups;
    vector<bool> largest;
    
    double exposed=0.;
    double litterSize =0.;
    double macrohuddle=0.;
    double microhuddle=0.;
    double cogX=0.;
    double cogY=0.;
    vector <double*> f;
    f.push_back(&exposed);
    f.push_back(&macrohuddle);
    f.push_back(&microhuddle);
    f.push_back(&litterSize);
    f.push_back(&cogX);
    f.push_back(&cogY);
    
    double Tp = 0.6;        // Preferred temperature
    double v2 = 5.0;        // Forwards drive
    double slope = 8.0;     // Homeothermotaxis constant (smoothes phase transition)
    
    
    ////////////////// SAVING/LOADING
    vector< vector <double*> > S;
    for(int I=0;I<A.size();I++){
        vector <double*> s;
        for(int i=0;i<A[I].nSOM;i++){
            for (int j=0;j<A[I].nin;j++){
                s.push_back(&A[I].win[i][j]);
            }
            for (int j=0;j<A[I].nout;j++) {
                s.push_back(&A[I].wout[i][j]);
            }
        }
        S.push_back(s);
    }
    ////////////////// SAVING/LOADING
    
    
    bool doing = true;
    
    while(doing){
        
        std::stringstream TIMEss;
        TIMEss<<setw(10)<<setfill('0')<<TIME;
        const char* TIMEcs = TIMEss.str().c_str();
        
        std::stringstream out;
        out.clear();
        out.setf(ios::fixed,ios::floatfield);
        
        /*
         *********************
         DEFINE OUTPUT MESSAGE
         */
        
        for(int I=0;I<A.size();I++){
            out<<A[I].Tb<<",";
        }
        for(int I=0;I<A.size();I++){
            int k = 0;
            for(int i=0;i<A[I].nP;i++){
                if(A[I].contacts[i]==1){
                    k++;
                }
            }
            out<<1-((double)k/(double)A[I].Npoint)<<",";
        }
        for(int I=0;I<A.size();I++){
            out<<A[I].mapTime<<",";
        }
        out<<macrohuddle<<","<<microhuddle<<",";
        
        /*
         DEFINE OUTPUT MESSAGE
         *********************
         */
        
        vector <string> command;
        string messageI=W.master.exchange(out.str().c_str());
        stringstream ss(messageI);
        while (ss.good()){
            string substr;
            getline(ss,substr,',');
            command.push_back(substr);
        } ss.clear();
        
        F%=A.size();
        
        switch(stoi(command[0])){
                
                // *** QUIT ***
            case 0:{W.logfile<<W.processName<<"@"<<TIMEcs<<": 0=QUIT"<<endl<<flush;
                
                W.logfile.close();
                for(int i=0;i<displays.size();i++){
                    displays[i].closeDisplay();
                }
                W.master.closeSocket();
                
                for(int i=0;i<W.ports.size();i++){
                    W.ports[i].closeSocket();
                }
                doing=false;
            } break;
                
                
                
                // *** STEP ***
            case 1:{//W.logfile<<W.processName<<"@"<<TIMEcs<<": 1=STEP"<<endl<<flush;
                
                // Exchange comms
                vector< vector <string> > commands(W.ports.size());
                for(int i=0;i<W.ports.size();i++){
                    string messageI=W.ports[i].exchange(out.str().c_str());
                    stringstream ss(messageI);
                    while (ss.good()){
                        string substr;
                        getline(ss,substr,',');
                        commands[i].push_back(substr);
                    } ss.clear();
                }
                
                
                // Homeothermotaxis
                for(int I=0;I<A.size();I++){
                    
                    // External co-ordinate system
                    vector<double> cog = A[I].getCoG();
                    A[I].orientation = atan2(cog[1]-A[I].Xnow[0](1),cog[0]-A[I].Xnow[0](0));
                    double cs = cos(-A[I].orientation),  ss = sin(-A[I].orientation);
                    double TlExternal = 0., TrExternal = 0.;
                    int countExternal = 0;
                    double cL = 0., cR = 0.;
                    for(int i=0;i<A[I].nP;i++){
                        double dx = cog[0]-A[I].Xnow[i](0);
                        double dy = cog[1]-A[I].Xnow[i](1);
                        if(dx*ss+dy*cs > 0.0){
                            TrExternal += A[I].Tc[i];
                            countExternal++;
                            cR += 1.0;
                        } else {
                            TlExternal += A[I].Tc[i];
                            cL += 1.0;
                        }
                    }
                    double cRnormE = (double)countExternal/(double)A[I].nP;
                    TrExternal /= cRnormE;
                    TlExternal /= (1.-cRnormE);
                    
                    // Change in orientation
                    A[I].dOrientation = 0.;
                    if(countExternal>0){
                        if(fabs(TrExternal-TlExternal)>0.){
                            A[I].dOrientation = atan(slope*(A[I].Tb-Tp)*(TrExternal-TlExternal));
                        }
                    }
                    
                    // Sensory input to network
                    A[I].input[0] = TlExternal;
                    A[I].input[1] = TrExternal;
                    
                }
                
                // APPLY MOTOR DRIVE
                double alpha = stod(command[1]); // Teacher/network control
                for(int I=0;I<A.size();I++){
                    switch (A[I].floating){
                        case(0):{ // on the ground
                            
                            vector<double> r = A[I].computeResponse(A[I].input);
                            int maxi = A[I].findWinner(r);
                            
                            vector<double> driveA (3,0.);
                            driveA[0] = cos(A[I].orientation + A[I].dOrientation);
                            driveA[1] = sin(A[I].orientation + A[I].dOrientation);
                            driveA[2] = 0.0;
                            
                            vector<double> driveB (3,0.);
                            driveB[0] = cos(A[I].orientation + A[I].wout[maxi][0]);
                            driveB[1] = sin(A[I].orientation + A[I].wout[maxi][0]);
                            driveB[2] = 0.;
                            
                            double driveX = driveA[0]*alpha+driveB[0]*(1.-alpha);
                            double driveY = driveA[1]*alpha+driveB[1]*(1.-alpha);
                            double driveZ = driveA[2]*alpha+driveB[2]*(1.-alpha);
                            
                            // Apply motor drive
                            A[I].externalForces(v2*driveX,v2*driveY,v2*driveZ-10.0);
                            
                        } break;
                        case(1):{ // in the air
                            A[I].externalForces(0.0,0.0,-10.0);
                        } break;
                        case(2):{ // in the air as part of a group
                            A[I].externalForces(0.0,0.0,-150.0);
                        } break;
                    }
                }
                
                // UNCOMMENT FOR LEARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                // UPDATE CORTICAL NETWORK
                if(alpha>0.){
                    for(int I=0;I<A.size();I++){
                        bool dynamic = stoi(command[1])==1;
                        vector<double> target(1,A[I].dOrientation);
                        A[I].iterateSOMsuper(A[I].input,target,dynamic);
                    }
                 
                }
                 //
                
                
                // Set contacts to zero
                for(int I=0;I<A.size();I++){
                    for(int i=0;i<A[I].nP;i++){
                        A[I].contacts[i]=0;
                    }
                }
                
                
                // SW: NOT SURE WHY RESETTING BEFORE USEAGE IN NEXT PART???????????????????????????????????????
                for(int I=0;I<A.size();I++){
                    for(int i=0;i<A[I].Npoint;i++){
                        A[I].Tc[i] = Ta;
                    }
                }

                
                // Enforce arena boundary (circle)
                for(int I=0;I<A.size();I++){
                    boundaryCircle(A[I], boundX, true, 0.1, 3);
                }
                
                // Enforce floor boundary (plane)
                for(int I=0;I<A.size();I++){
                    boundary(A[I],+floorH, 2, 0., 2);
                }
                
                
                // Collision detection / resolution
                vector<vector<int> > agg=collisionDetection(A,dt);
                
                
                // Temperature update (terms 1 and 3 from Glancy et al., (2015))
                for (int I=0;I<A.size();I++){
                    double TC = 0.;
                    double Eta = 1.;
                    double TCnorm = 0.;
                    for(int i=0;i<A[I].nP;i++){
                        if(A[I].contacts[i]==1){ // pup-pup
                            TC += A[I].Tc[i];
                            TCnorm += 1.;
                        }
                    }
                    if(TCnorm>0){
                        TC /= TCnorm;
                        Eta = 1.-(TCnorm/(double)A[I].nP);
                    }
                    A[I].Tb += ( -A[I].k1*Eta*(A[I].Tb-Ta) -A[I].k2*(1.-Eta)*(A[I].Tb-TC) + A[I].G) *dt;
                }
                
                // ANALYSE GROUPINGS, RESOLVE 'FLOATERS'
                vector<vector<int> > conns;
                vector<int> conn (2,0);
                for(int i=0;i<A.size();i++){
                    for(int j=0;j<A.size();j++){
                        if(agg[i][j]){
                            conn[0]=i;
                            conn[1]=j;
                            conns.push_back(conn);
                        }
                    }
                }
                groups=graph(conns, A.size());
                
                for(int i=0;i<A.size();i++){
                    A[i].floating = 0;
                }
                vector<bool> floaters(groups.size(),true);
                for(int i=0;i<groups.size();i++){
                    for(int j=0;j<groups[i].size();j++){
                        for(int k=0;k<A[groups[i][j]].nP;k++){
                            if(A[groups[i][j]].contacts[k]==2){
                                floaters[i] = false;
                                break;
                            }
                        }
                    }
                }
                for(int i=0;i<floaters.size();i++){
                    if(floaters[i]){
                        int groupsize = groups[i].size();
                        if(groupsize==1){
                            A[groups[i][0]].floating = 1;
                        } else {
                            for(int j=0;j<groupsize;j++){
                                A[groups[i][j]].floating = 2;
                            }
                        }
                    }
                }
                
                
                
                largest.resize(A.size());
                for(int i=0;i<A.size();i++){
                    largest[i]=false;
                }
                for(int i=0;i<groups[0].size();i++){
                    largest[groups[0][i]]=true;
                }
                litterSize =(double)A.size();                   // number of pups
                macrohuddle=(double)groups[0].size();           // size of largest group
                microhuddle=(double)groups.size();              // number of groups
                cogX=0.;
                cogY=0.;
                for(int i=0;i<A[F].Npoint;i++){
                    cogX += A[F].Xnow[i](0);
                    cogY += A[F].Xnow[i](1);
                }
                cogX /= (double)A[F].Npoint;
                cogY /= (double)A[F].Npoint;
                
                
                // END STEP
                TIME++;
                
                
            } break;
                
                /*
                
                // ADD A NEW ANIMAT
            case 2:{W.logfile<<W.processName<<"@"<<TIMEcs<<": 2=ADDA"<<endl<<flush;
                
                A.push_back(basic::getAnimat(pupDensity,&W));
                
                vector<double> Porg(3,0.);
                Porg[0] = stod(command[1]);
                Porg[1] = stod(command[2]);
                Porg[2] = stod(command[3]);
                vector<double> Rorg(3,0.);
                Rorg[0] = stod(command[4]);
                Rorg[1] = stod(command[5]);
                Rorg[2] = stod(command[6]);
                A.back().reposition(Porg,Rorg);
                A.back().setOrientation(stod(command[6]));
                
            } break;
                
                // DELETE THE LAST ANIMAT
            case 3:{W.logfile<<W.processName<<"@"<<TIMEcs<<": 3=DELA"<<endl<<flush;
                if(A.size()>0){
                    A.pop_back();
                }
                if(A.size()<=0){
                    A.push_back(basic::getAnimat(pupDensity,&W));
                    vector<double> Porg(3,0.);
                    Porg[0] = stod(command[1]);
                    Porg[1] = stod(command[2]);
                    Porg[2] = stod(command[3]);
                    vector<double> Rorg(3,0.);
                    Rorg[0] = stod(command[4]);
                    Rorg[1] = stod(command[5]);
                    Rorg[2] = stod(command[6]);
                    A.back().reposition(Porg,Rorg);
                    A.back().setOrientation(stod(command[5]));
                }
                
            } break;
                
                */
                
                
            case 4:{W.logfile<<W.processName<<"@"<<TIMEcs<<": 4=PARS"<<endl<<flush;
                
                switch(stoi(command[1])){
                    case 0:{A[stoi(command[2])].k1 =stod(command[3]);}break;
                    case 1:{A[stoi(command[2])].k2 =stod(command[3]);}break;
                    case 2:{A[stoi(command[2])].G  =stod(command[3]);}break;
                }
                
            } break;
                
            case 5:{W.logfile<<W.processName<<"@"<<TIMEcs<<": 5=ENVP"<<endl<<flush;
                
                switch(stoi(command[1])){
                        
                    case 0:{Ta=stod(command[2]);}break;
                    case 1:{boundX=stod(command[2]);}break;
                    case 2:{boundY=stod(command[2]);}break;
                }
                
            } break;
                
                
            case 6:{W.logfile<<W.processName<<"@"<<TIMEcs<<": 6=SOMP"<<endl<<flush;
                for(int i=0;i<A.size();i++){
                    A[i].setParam(stoi(command[1]),stod(command[2]));
                }
            } break;

                
                /* RE-SHUFFLE ANIMATS */
            case 7:{W.logfile<<W.processName<<"@"<<TIMEcs<<": 7=SHUF"<<endl<<flush;
                
                shuffle(A);
            
            } break;
                
                
                /* ARENA PLOT */
            case 8:{//W.logfile<<W.processName<<"@"<<TIMEcs<<": 8=DISA"<<endl<<flush;
                
                displays[0].resetDisplay(fix,eye,rot);
                for(int I=0;I<A.size();I++){
                    vector<vector<double> > Xsp;
                    vector<double> xsp(3);
                    for(int i=0;i<A[I].nP;i++){
                        xsp[0]=A[I].Xnow[i](0);
                        xsp[1]=A[I].Xnow[i](1);
                        xsp[2]=A[I].Xnow[i](2);
                        Xsp.push_back(xsp);
                    }
                    vector<vector<double> > cl;
                    for(int i=0;i<A[I].nP;i++){
                        if(A[I].contacts[i]==1){
                            cl.push_back(getJetColor(fmin(fmax((A[I].Tc[i]-mint)*oneOverMaxtMinusMint,0.),1.)));
                        } else {
                            cl.push_back(getJetColor(fmin(fmax((A[I].Tb-mint)*oneOverMaxtMinusMint,0.),1.)));
                        }
                    }
                    displays[0].drawSphereFromMesh(Xsp,A[I].sphereM,cl);
                }
                // indicate focal pup
                displays[0].drawCylinder(A[F].Xnow[0](0),A[F].Xnow[0](1),A[F].Xnow[0](2),A[F].Xnow[0](0),A[F].Xnow[0](1),A[F].Xnow[0](2)+0.5,1e-6,0.01,12,vector<double>(3,0.95));
                
                // draw arena
                vector<double> p1(3,0.),p2(3,0.),p3(3,0.),p4(3,0.),nm(3,0.),cl(3,0.);
                nm[0]=0.;nm[1]=0.;nm[2]=1.;
                double floorHminusR = floorH;
                displays[0].drawCylinder(0.,0.,floorHminusR,0.,0.,floorHminusR-0.1,boundX,boundX,60,vector<double>(3,0.5));
                displays[0].redrawDisplay();
                
            } break;
                
                
            case 9:{//W.logfile<<W.processName<<"@"<<TIMEcs<<": 9=DISB"<<endl<<flush;
                
                double space = 1.25;
                displays[1].resetDisplay(fix,eye,rot);
                
                // PLOT ACTIVITY
                vector<double> r = A[F].computeResponse(A[F].input);
                for(int i=0;i<A[F].H[0].size();i++){
                    vector<double> c=getJetColor(r[i]);
                    displays[1].drawHex(A[F].H[0][i]+space,A[F].H[1][i]+2*space,0.0,0.016,c[0],c[1],c[2]);
                }

                
                {
                    // PLOT OUTPUT WEIGHTS
                    double pltMax = -1e9,  pltMin = +1e9;
                    for(int i=0;i<A[F].nSOM;i++){
                        double d = A[F].win[i][0]-A[F].win[i][1];
                        if(d>pltMax){ pltMax = d; }
                        if(d<pltMin){ pltMin = d; }
                    }
                    double pltRange = pltMax-pltMin;
                    for(int i=0;i<A[F].nSOM;i++){
                        vector<double> c = getJetColor((A[F].win[i][0]-A[F].win[i][1]-pltMin)/(pltRange));
                        displays[1].drawHex(A[F].H[0][i]+space,A[F].H[1][i]-0.*space,0.0,0.016,c[0],c[1],c[2]);
                    }
                }
                {
                    // PLOT OUTPUT WEIGHTS
                    double pltMax = -1e9,  pltMin = +1e9;
                    for(int i=0;i<A[F].nSOM;i++){
                        if(A[F].wout[i][0]>pltMax){ pltMax = A[F].wout[i][0]; }
                        if(A[F].wout[i][0]<pltMin){ pltMin = A[F].wout[i][0]; }
                    }
                    double pltRange = pltMax-pltMin;
                    for(int i=0;i<A[F].nSOM;i++){
                        vector<double> c = getJetColor((A[F].wout[i][0]-pltMin)/(pltRange));
                        displays[1].drawHex(A[F].H[0][i]-space,A[F].H[1][i]-2.*space,0.0,0.016,c[0],c[1],c[2]);
                    }
                }
                
                // PLOT FOCAL PUP
                vector<double> cog = A[F].getCoG();
                for(int k=0;k<A[F].clusterIDs.size();k++){
                    for(int i=0;i<A[F].clusterIDs[k].size();i++){
                        int j= A[F].clusterIDs[k][i];
                        if(k==5){
                            displays[1].drawSphere(A[F].Xnow[j](0)-cog[0],A[F].Xnow[j](1)-cog[1],A[F].Xnow[j](2)-cog[2],A[F].Ranimat[j],red,9);
                        } else {
                            displays[1].drawSphere(A[F].Xnow[j](0)-cog[0],A[F].Xnow[j](1)-cog[1],A[F].Xnow[j](2)-cog[2],A[F].Ranimat[j],vector<double>(3,0.9),9);
                        }
                    }
                }
                
                displays[1].redrawDisplay();
                
                
            } break;
                
            case 10:{W.logfile<<W.processName<<"@"<<TIMEcs<<": 10=RECD"<<endl<<flush;
                
                std::stringstream frameFile1;
                frameFile1<<"logs/"<<W.processName<<"frameA";
                frameFile1<<setw(5)<<setfill('0')<<frameN;
                frameFile1<<".png";
                displays[0].saveImage(frameFile1.str());
                
                std::stringstream frameFile2;
                frameFile2<<"logs/"<<W.processName<<"frameB";
                frameFile2<<setw(5)<<setfill('0')<<frameN;
                frameFile2<<".png";
                displays[1].saveImage(frameFile2.str());
                frameN++;
                
            } break;
                
                
                // *** SAVE ***
            case 11:{W.logfile<<W.processName<<"@"<<TIMEcs<<": 11=SAVE"<<endl<<flush;
                
                if(command.size()==2){
                    ofstream outFile;
                    std::stringstream oFile; oFile<<command[1];
                    outFile.open(oFile.str().c_str(),ios::out|ios::binary);
                    double n=(double)S.size();double* N=&n;
                    outFile.write((char*)N,sizeof(double));
                    for(int i=0;i<S.size();i++){
                        double m=(double)S[i].size();double* M=&m;
                        outFile.write((char*)M,sizeof(double));
                        for(int j=0;j<S[i].size();j++){
                            outFile.write((char*)S[i][j],sizeof(double));
                        }
                    }
                    outFile.close();
                }else{W.logfile<<"No output filename."<<endl<<flush;}
                
            } break;
                
                
                // *** LOAD ***
            case 12:{W.logfile<<W.processName<<"@"<<TIMEcs<<": 12=LOAD"<<endl<<flush;
                
                if(command.size()==2){
                    double dummy;
                    ifstream inFile;
                    std::stringstream iFile;
                    iFile<<command[1];
                    inFile.open(iFile.str().c_str(),ios::in|ios::binary);
                    inFile.read((char*)&dummy,sizeof(double));
                    int I=(int)dummy;
                    if(I==S.size()){
                        for(int i=0;i<I;i++){
                            inFile.read((char*)&dummy,sizeof(double));
                            int J=(int)dummy;
                            if(J==S[i].size()){
                                for(int j=0;j<J;j++){
                                    inFile.read((char*)S[i][j],sizeof(double));
                                }
                            }else{W.logfile<<"Wrong dims I."<<endl<<flush;}
                        }
                    }else{W.logfile<<"Wrong dims J."<<endl<<flush;}
                    inFile.close();
                }else{W.logfile<<"No input filename."<<endl<<flush;}
            } break;
                
                
        }
    }
    
    return 0;
};


// COLLISION DETECTION / RESOLUTION
vector<vector<int> > collisionDetection(std::vector<animat>& a,double dt){
    
    // contact pattern analysis
    vector<vector<int> > agg;
    vector<int> ag (a.size(),0.);
    for(int i=0;i<a.size();i++){
        agg.push_back(ag);
    }
    
    vec J2I = zeros<vec>(3);
    vec I2J = zeros<vec>(3);
    
    double dx, dy, dz, d, rAB, halfOverlap;
    
    // Randomize order of contact resolution
    vector<double> order(a.size(),0.);
    for(int i=0;i<a.size();i++){
        order[i]=tools::randFloat();
    }
    vector<int> Ordered = tools::sort(order);
    
    
    for(int Iorg=0;Iorg<a.size();Iorg++){
        int I = Ordered[Iorg];
        a[I].projectPositions();
        int nI = a[I].Npoint;
        for(int Jorg=0;Jorg<a.size();Jorg++){
            int J = Ordered[Jorg];
            if(Iorg<Jorg){
                int nJ = a[J].Npoint;
                for(int i=0;i<nI;i++){
                    for(int j=0;j<nJ;j++){
                        dx = a[I].Xnew[i][0]-a[J].Xnew[j][0];
                        dy = a[I].Xnew[i][1]-a[J].Xnew[j][1];
                        dz = a[I].Xnew[i][2]-a[J].Xnew[j][2];
                        d = sqrt(dx*dx+dy*dy+dz*dz);
                        rAB = a[I].Ranimat[i]+a[J].Ranimat[j];
                        if(d<rAB){
                            J2I(0) = dx/d;
                            J2I(1) = dy/d;
                            J2I(2) = dz/d;
                            if(dz>0){ // Prioritize by height
                                a[I].Xnew[i]=a[I].Xnow[i]+J2I*(rAB-d);
                                a[J].Vel[j]=-J2I;
                            } else {
                                a[J].Xnew[j]=a[J].Xnow[j]-J2I*(rAB-d);
                                a[I].Vel[i]=J2I;
                            }
                            a[I].contacts[i]=1;
                            a[J].contacts[j]=1;
                            agg[I][J]=1;
                            
                            // Just use this as a display thingy
                            a[I].Tc[i] = a[J].Tb;
                            a[J].Tc[j] = a[I].Tb;
                        }
                    }
                }
            }
        }
        a[I].integrate();
    }
    return agg;
}


// PLANE BOUNDARY, E.G., FOR FLOOR
void boundary(animat& a,double loc,int axis,double bounce, int contactFlag){
    
    if (loc<0.){
        for(int j=0;j<a.Npoint;j++){
            if (a.Xtemp[j][axis]<loc+a.Ranimat[j]){
                vec normal = zeros<vec>(3);
                normal(axis) = 1.0;
                vec vNew = -2.*dot(a.Vel[j],normal)*normal+a.Vel[j];
                a.Xnew[j]=a.Xnow[j]+vNew*a.params.dt*bounce;
                if(a.Xnew[j][axis]<loc+a.Ranimat[j]){
                    a.Xnew[j][axis]=loc+a.Ranimat[j];
                    a.contacts[j]=contactFlag;
                }
            }
        }
    }else{
        for(int j=0;j<a.Npoint;j++){
            if (a.Xtemp[j][axis]>loc-a.Ranimat[j]){
                vec normal = zeros<vec>(3);
                normal(axis) = -1.0;
                vec vNew = -2.*dot(a.Vel[j],normal)*normal+a.Vel[j];
                a.Xnew[j]=a.Xnow[j]+vNew*a.params.dt*bounce;
                if(a.Xnew[j][axis]>loc-a.Ranimat[j]){
                    a.Xnew[j][axis]=loc-a.Ranimat[j];
                    a.contacts[j]=contactFlag;
                }
            }
        }
    }
}

// ARENA BOUNDARY
void boundaryCircle(animat& a,double radius,bool slippy,double bounce, int contactFlag){
    double r2 = radius*radius;
    for(int j=0;j<a.Npoint;j++){
        double x = fabs(a.Xtemp[j][0])+a.Ranimat[j];
        double y = fabs(a.Xtemp[j][1])+a.Ranimat[j];
        if(x*x+y*y>r2){
            double t = atan2(a.Xtemp[j][1],a.Xtemp[j][0]);
            vec normal = zeros<vec>(3);
            normal(0) = -cos(t);
            normal(1) = -sin(t);
            vec vNew;
            if(slippy){ // combine point velocities with cylinder normal
                vNew=-2.*dot(a.Vel[j],normal)*normal+a.Vel[j];
            } else { // bounce
                vNew = normal;
            }
            a.Xnew[j]=a.Xnow[j]+vNew*a.params.dt*bounce;
            double x2 = fabs(a.Xnew[j][0])+a.Ranimat[j];
            double y2 = fabs(a.Xnew[j][1])+a.Ranimat[j];
            if(x2*x2+y2*y2>r2){
                double t = atan2(a.Xnew[j][1],a.Xnew[j][0]);
                
                a.Xnew[j][0]=cos(t)*(radius-a.Ranimat[j]);
                a.Xnew[j][1]=sin(t)*(radius-a.Ranimat[j]);
                a.contacts[j]=contactFlag;
            }
        }
    }
}

// QUANTIFY GROUPINGS BASED ON CONTACTS
vector<vector<int> > graph(vector<vector<int> > conn, int N){
    vector<vector<int> > graph;
    int a, b, Ain, Bin;
    bool nov, ain, bin;
    for(int k=0;k<conn.size();k++){
        a = conn[k][0];
        b = conn[k][1];
        if(a != b){
            // reset flags
            Ain = -1;
            Bin = -1;
            nov = true;
            // check if conn k is in the graph
            int I = graph.size();
            for(int i=0;i<I;i++){
                ain = false;
                bin = false;
                int J = graph[i].size();
                for(int j=0;j<J;j++){
                    if(a == graph[i][j]){
                        ain = true;
                        Ain = i;
                    }
                    if(b == graph[i][j]){
                        bin = true;
                        Bin = i;
                    }
                }
                if(ain && !bin){
                    graph[i].push_back(b);
                }
                if(!ain && bin){
                    graph[i].push_back(a);
                }
                if(ain||bin){
                    nov=false;
                }
            }
            // Add a new group
            if(nov){
                graph.push_back(conn[k]);
            }
            // Join two existing groups
            if(Ain>-1 && Bin>-1 && Ain!=Bin){
                graph[Ain].pop_back();
                graph[Bin].pop_back();
                for(int l=0;l<graph[Bin].size();l++){
                    graph[Ain].push_back(graph[Bin][l]);
                }
                graph.erase(graph.begin()+Bin);
            }
        }
    }
    
    for(int k=0;k<N;k++){
        bool isolated = true;
        int I = graph.size();
        for(int i=0;i<I;i++){
            int J = graph[i].size();
            for(int j=0;j<J;j++){
                if(k==graph[i][j]){
                    isolated = false;
                    //break;
                }
            }
        }
        if(isolated){
            vector<int> isolate(1,k);
            graph.push_back(isolate);
        }
    }
    
    // sort by descending group size
    vector<vector<int> > graphSorted;
    while(graph.size()){
        int maxVal = 0;
        int maxInd = 0;
        for(int i=0;i<graph.size();i++){
            if(graph[i].size()>maxVal){
                maxVal = graph[i].size();
                maxInd = i;
            }
        }
        graphSorted.push_back(graph[maxInd]);
        graph.erase(graph.begin()+maxInd);
    }
    
    return graphSorted;
}


void shuffle(std::vector<animat>& a){

    vector<double> order(a.size(),0.);
    for(int i=0;i<a.size();i++){
        order[i]=tools::randFloat();
    }
    vector<int> Ordered = tools::sort(order);

    for(int I=0;I<a.size();I++){
        vector<double> Porg(3,0.);
        Porg[0] = (tools::randFloat()-0.5)*0.5;
        Porg[1] = (tools::randFloat()-0.5)*0.5;
        Porg[2] = (double)Ordered[I];
        vector<double> Rorg(3,0.);
        a[I].reposition(Porg,Rorg);
        a[I].setOrientation(tools::randFloat()*M_PI*2.); // z-rotation
    }
}
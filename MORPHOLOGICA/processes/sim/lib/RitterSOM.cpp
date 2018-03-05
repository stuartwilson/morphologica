//
//  RitterSOM.cpp
//  Emerge
//

#include "RitterSOM.h"
#include "tools.h"
#include <stdlib.h>
#include <math.h>

#include <iostream>
#include <fstream>
#include <numeric>

using namespace std;

vector<vector<double> > hex2(int scale, int offset){
    vector <vector <double> > H(5);
    int S = pow(2.0, scale-1) - offset;
    double s = S - offset;
    //double rad = 0.5 / s;
    //int n = 1+3*S*(S+1);
    for (int d=0; d<=S; d++){
        for (int r=-S; r<=S; r++){
            for (int g=-S; g<=S; g++){
                for (int b=-S;b<=S;b++){
                    if (abs(r)+abs(g)+abs(b)==d*2 && (r+g+b==0)){
                        H[0].push_back((g/2.+r)/s);             //X
                        H[1].push_back(g*(sqrt(3.)/2.)/s);      //Y
                        H[2].push_back(r);                      //R
                        H[3].push_back(g);                      //G
                        H[4].push_back(b);                      //B
                    }
                }
            }
        }
    }
    return H;
};

RitterSOM::RitterSOM():
H(),
epsilonINP(0.),
epsilonOUT(0.),
epsilonVAR(0.),
epsilonTHR(0.),
nhoodSigma(0.),
radius(0.),
win(),
wout(),
noise(),
threshold(),
nSOM(0),
nin(0),
nout(0),
mapTime(0),
tF(0),
x(),
y(),
reward(0.)
{
}

RitterSOM::~RitterSOM(){

};

void RitterSOM::initSOM(int nin, int nout) {
    this->nin = nin;
    this->nout = nout;

    mapTime = 0;
    //tF = 10000.;
    //tF = 100000.;
    tF = 10000.;
    
    H = hex2(6,0);
    nSOM = static_cast<int>(H[0].size());
    x = H[0];
    y = H[1];
    

    // Assign weights
    win.resize(nSOM);
    wout.resize(nSOM);
    for(int i=0;i<nSOM;i++){
        win[i].resize(nin);
        wout[i].resize(nout);
    }
    
    randomizeWeights();
    
    responseSigma = 1.;
    nhoodSigma = 0.0;
    
    noise.resize(nSOM, 1.0);//M_PI/2.0);   // SW: WHY PI/2. ??
    threshold.resize(nSOM, 0.0);
    
}

void RitterSOM::randomizeWeights(){
    for(int i=0;i<nSOM;i++){
        for (int j=0;j<nin;j++){
            win[i][j] = tools::randFloat();
        }
        double norm = 0.;
        for (int j=0;j<nout;j++) {
            wout[i][j] = tools::randFloat()-0.5;
            norm += wout[i][j]*wout[i][j];
        }
        norm = sqrt(norm);
        for (int j=0;j<nout;j++) {
            //wout[i][j] /= norm;
        }
    }
}

vector<int> RitterSOM::maxInWeight(void){
    vector<int> maxInds(nSOM);
    for (int i=0; i<nSOM; i++){
        double maxVal = -1e9;
        int maxInd = 0;
        for (int j=0; j<nin; j++){
            if(maxVal<win[i][j]){
                maxVal = win[i][j];
                maxInd = j;
            }
        }
        maxInds[i] = maxInd;
    }
    return maxInds;
}


vector<int> RitterSOM::maxOutWeight(void){
    vector<int> maxInds(nSOM);
    for (int i=0; i<nSOM; i++){
        double maxVal = -1e9;
        int maxInd = 0;
        for (int j=0; j<nout; j++){
            if(maxVal<wout[i][j]){
                maxVal = wout[i][j];
                maxInd = j;
            }
        }
        maxInds[i] = maxInd;
    }
    return maxInds;
}



void RitterSOM::setParam(int paramInd, double paramVar){
    
    switch(paramInd){
        case(0):{
            // Learning rate for input weights
            epsilonINP = paramVar;
        } break;
        case(1):{
            // Learning rate for output weights
            epsilonOUT = paramVar;
        } break;
        case(2):{
            // smoothing of recent updates of reward function (gamma in Ritter)
            epsilonVAR = paramVar;
        } break;
        case(3):{
            //
            epsilonTHR = paramVar;
        } break;
        case(4):{
            nhoodSigma = paramVar;
            nhoodSigmaPRE = 1./(nhoodSigma*nhoodSigma);
        } break;
        case(5):{
            radius = paramVar;
            radiusPRE = radius*radius;
        } break;
        case(6):{
            mapTime = 0;
            tF = paramVar;
        } break;
    }
    
}

/*
void RitterSOM::updateDynamicVariables(void){
    
    epsilonINP = 0.5 * pow(0.04, (double)t/(double)tF);
    epsilonOUT = 0.2;
    epsilonVAR = 0.005;
    epsilonTHR = 0.05;
    
    double sigmaActual = 0.5 * pow(0.5, (double)t/(double)tF);
    sigma = 1./pow(sigmaActual,2.0); // SW: Just speeds up a bit
    
    double nradiActual = 1. * fmax(exp(-1. * t / tF / ndecay), 1./(double)nSOM);
    nradi = nradiActual*nradiActual; // SW: Just speeds up a bit
    
    t++;
}
*/


void RitterSOM::updateDynamicVariables(void){
    
    double T = (double)mapTime / tF;


    /* ORIGINALS
 
    epsilonINP = 0.0025;//fmax(0.005*pow(0.01,T),0.0005);
    epsilonOUT = 0.0025;//fmax(0.005*pow(0.01,T),0.0005);
    //epsilonVAR = 0.05*pow(0.04,T);
    epsilonVAR = 0.01; // smoothing of recent updates of reward function (gamma in Ritter)
    epsilonTHR = 0.01;
    //epsilonTHR = 0.1*pow(0.04,T);

    radius =     10.;//fmax(2.0*pow(0.04,T),2./(double)nSOM);
    nhoodSigma = fmax(5.0*pow(0.02,T),0.5);
     */

    /*
    epsilonINP = 0.0025;//fmax(0.005*pow(0.01,T),0.0005);
    epsilonOUT = 0.0025;//fmax(0.005*pow(0.01,T),0.0005);
    //epsilonVAR = 0.05*pow(0.04,T);
    epsilonVAR = 0.01; // smoothing of recent updates of reward function (gamma in Ritter)
    epsilonTHR = 0.01;
    //epsilonTHR = 0.1*pow(0.04,T);
    
    radius =     10.;//fmax(2.0*pow(0.04,T),2./(double)nSOM);
    nhoodSigma = fmax(1.0*pow(0.02,T),0.5);
     nhoodSigmaPRE = 1./(nhoodSigma*nhoodSigma);
     radiusPRE = radius*radius;
     */
    
    /*
    epsilonINP = paramVar;
    epsilonOUT = paramVar;
    epsilonVAR = paramVar;
    epsilonTHR = paramVar;
    nhoodSigma = paramVar;
    nhoodSigmaPRE = 1./(nhoodSigma*nhoodSigma);
    radius = paramVar;
    radiusPRE = radius*radius;
    */
    
    

    epsilonINP *= 0.9997; // default 0.9997
    epsilonOUT *= 0.9997; // default 0.9997
    radiusPRE *= 0.9999; // default 0.9999
    
     
    //mapTime++;
}

vector<double> RitterSOM::computeResponse(vector<double> in){
    
    vector<double> r(nSOM,0);
    for(int i=0;i<nSOM;i++){
        for(int j=0;j<nin;j++){
            r[i] += (win[i][j]-in[j])*(win[i][j]-in[j]);
            //r[i] += (win[i][j]*in[j]);//*(win[i][j]-in[j]);
        }
        r[i]=exp(-r[i]*responseSigma);
    }
    return r;
    
}

vector<double> RitterSOM::calculateOutput(vector<double> in){
    vector<double> output(nout,0);
    int maxi = findWinner(computeResponse(in));
    for(int i=0;i<nout;i++){
        output[i] = wout[maxi][i] + noise[maxi] * tools::normalDistributionValue();
    }
    return output;
}

int RitterSOM::findWinner(vector<double> r) {
    double maxv=0.;
    int maxi=0;
    for(int i=0;i<nSOM;i++){
        if(maxv<r[i]){
            maxv=r[i];
            maxi=i;
        }
    }
    return maxi;
}



void RitterSOM::iterateSOM(vector<double> in,double rewardNow,bool dynamic){
    
    
    int maxi = findWinner(computeResponse(in));
    
    double rewardDiff = rewardNow-reward-threshold[maxi];
    if(rewardDiff>0){//threshold[maxi]){
        mapTime++;
        vector<double> out = calculateOutput(in);
        for(int i=0;i<nSOM;i++){
            double dist = (pow(H[0][i]-x[maxi],2.)+pow(H[1][i]-y[maxi],2.));
            double h = exp(-dist*nhoodSigmaPRE);
            if(dist <= radiusPRE){
                for(int j=0;j<nin;j++){
                    win[i][j] += epsilonINP*(in[j]-win[i][j])*h;
                }
                for(int j=0;j<nout;j++){
                    wout[i][j] += epsilonOUT*(out[j]-wout[i][j])*h;
                }
            }
            noise[i] -= epsilonVAR*noise[i]*h;
        }
        threshold[maxi] += epsilonTHR*rewardDiff;
        
        if(dynamic){
            updateDynamicVariables();
            // SW: This was outside if statement before
            // - makes more sense but may depart from Ritter
        }
    }
    reward = rewardNow;
}


void RitterSOM::iterateSOMsuper(vector<double> in, vector<double> out, bool dynamic){
    
    
    int maxi = findWinner(computeResponse(in));

    for(int i=0;i<nSOM;i++){
    double dist = (pow(H[0][i]-x[maxi],2.)+pow(H[1][i]-y[maxi],2.));
    double h = exp(-dist*nhoodSigmaPRE);
    if(dist <= radiusPRE){
        for(int j=0;j<nin;j++){
            win[i][j] += epsilonINP*(in[j]-win[i][j])*h;;//epsilonINP*(in[j]-win[i][j])*h;
        }
        double norm = 0.;
        for(int j=0;j<nout;j++){
            wout[i][j] += epsilonOUT*(out[j]-wout[i][j])*h;
            norm += wout[i][j]*wout[i][j];
        }
        norm = sqrt(norm);
        for (int j=0;j<nout;j++) {
            //wout[i][j] /= norm; //SW: THIS IS NEEDED !!!!!!!!!!!!!!
        }

        
        }
    }
    
        mapTime++;
    
    if(dynamic){
       updateDynamicVariables();
    }

}





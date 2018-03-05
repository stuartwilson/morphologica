//
//  animat.cpp
//

#include "animat.h"

using namespace arma;

animat::animat(){
    
};


animat::animat(world* W, int n, int somNout, double sphereRadius, double taxelRadius, double TbInit){
    
    this->W=W;
    
    time = 0;
    Tb = TbInit;
    theta = 0.;
    orientation = 0.0;
    dOrientation = 0.0;
    floating = 0;
    
    k1 = 1.0;
    k2 = 1.0;
    G = 0.0;
    
    initSphere(n,sphereRadius);
    
    this->nP = sphereN*sphereN*6;
    this->N = nP;
    
    
    //initSOM(nP,somNout);
    initSOM(2,somNout);


    Ranimat.resize(nP,taxelRadius);
    
    for(int i=0;i<nP;i++){
        Xanimat.push_back(sphereX[i]);
    }
    
    Manimat.resize(nP,1./(double)nP);
    
    clusters.resize(1);
    for(int i=0;i<nP;i++){
        clusters[0].push_back(i);
    }

    Tc.resize(nP,Tb);
    input.resize(2,Tb);
    contacts.resize(nP,0);
    labels.resize(nP,0);
}

//REMEMBER TO CALL AFTER POSING BODY
void animat::init(void){
    initShape(Xanimat,Manimat,clusters,boneIDs);
}

void animat::init(vector<double>Tr,vector<double>Rot){
    init();
    reposition(Tr,Rot);
}


animat::~animat(void){
    //W->logfile<<"animat deconstructor called"<<endl;
}

void animat::attachBoneEnd(int boneIndex, bool isBoneOrigin, int clusterIndex, double boneMass, double rad){
    
    if(isBoneOrigin){
        Xanimat.push_back(boneX.col(axisIn[boneIndex]));
    }else{
        Xanimat.push_back(boneX.col(axisOut[boneIndex]));
    }
    clusters[clusterIndex].push_back(N);
    boneIDs.push_back(N);
    N++;
    
    Ranimat.push_back(rad);
    Manimat.push_back(boneMass);
    contacts.push_back(0);
    Tc.push_back(Tb);
    labels.push_back(0);
    attachedBones.push_back(boneIndex);
    boneIsOrigins.push_back(isBoneOrigin);
    
}


void animat::reposition(vector<double>Tr,vector<double>Rot){
    
    std::vector<arma::vec> X = Xanimat;
    
    // Rotate cloud about origin
    mat RX, RY, RZ;
    RX <<1.0<<0.0<<0.0<<endr<<0.0<<cos(Rot[0])<<-sin(Rot[0])<<endr<<0.0<<sin(Rot[0])<<cos(Rot[0])<<endr;
    RY <<cos(Rot[1])<<0.0<<sin(Rot[1])<<endr<<0.0<<1.0<<0.0<<endr<<-sin(Rot[1])<<0.0<<cos(Rot[1])<<endr;
    RZ <<cos(Rot[2])<<-sin(Rot[2])<<0.0<<endr<<sin(Rot[2])<<cos(Rot[2])<<0.0<<endr<<0.0<<0.0<<1.0<<endr;
    for(int i=0;i<X.size();i++){
        X[i] = RX*X[i];
        X[i] = RY*X[i];
        X[i] = RZ*X[i];
    }
    
    // Translate cloud
    for(int i=0;i<X.size();i++){
        X[i](0)+=Tr[0];
        X[i](1)+=Tr[1];
        X[i](2)+=Tr[2];
    }
    
    relocate(X);
}

void animat::setOrientation(double phi){
    orientation = phi;
}

void animat::applyRotations(void){
    if (time>=0){
        vector<vec> bonePos;
        for(int i=0;i<attachedBones.size();i++){
            if(boneIsOrigins[i]){
                bonePos.push_back(boneX.col(axisIn[attachedBones[i]]));
            }else{
                bonePos.push_back(boneX.col(axisOut[attachedBones[i]]));
            }
        }
        updateBonePosition(bonePos);
    }
    time++;
}

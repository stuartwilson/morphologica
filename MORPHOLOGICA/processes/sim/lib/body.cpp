#include "body.h"


body::body(void){
    
    this->n = 1;
    nBone = 1;
    zero = 0.;
    In.push_back(0);
    Len.push_back(0.);
    boneLength.push_back(0.);
    vector <int> ax(3);
    xyz.push_back(ax);
    
    // Axis0=dummy, Axis1=+x, Axis2=+y, Axis3=+z
	arma::mat dummy;
	dummy << 1. << 0. << 0. << endr;
	Qx =dummy.t();
	Qy =dummy.t();
	Qz =dummy.t();
	
    pathI.resize(1);
    pathI[0].push_back(0);
    
    visible = true;
};

body::~body(void){
    
};

void body::addAxis(int in, double Rx, double Ry, double Rz, double len){
    
    /*
     Iteratively construct the path into the new Axis i:
     PathI[i][ordered indices of parents from origin to i]
     */
    
    if (in < n){
        vector <int> j;
        for (int i=0; i<pathI[in].size(); i++){
            j.push_back(pathI[in][i]);
        }
        j.push_back(in);
        pathI.push_back(j);
        
		In.push_back(in);
        
		arma::mat RX, RY, RZ;
        RX <<1.0<<0.0<<0.0<<endr<<0.0<<cos(Rx)<<-sin(Rx)<<endr<<0.0<<sin(Rx)<<cos(Rx)<<endr;
        RY <<cos(Ry)<<0.0<<sin(Ry)<<endr<<0.0<<1.0<<0.0<<endr<<-sin(Ry)<<0.0<<cos(Ry)<<endr;
        RZ <<cos(Rz)<<-sin(Rz)<<0.0<<endr<<sin(Rz)<<cos(Rz)<<0.0<<endr<<0.0<<0.0<<1.0<<endr;
		
        arma::mat fx, fy, fz;
        fx<<1.<<endr<<0.<<endr<<0.<<endr;
        fy<<0.<<endr<<1.<<endr<<0.<<endr;
        fz<<0.<<endr<<0.<<endr<<1.<<endr;
        
        fx = RX*RY*RZ*fx;
        fy = RX*RY*RZ*fy;
        fz = RX*RY*RZ*fz;
        
		Qx = join_rows(Qx,fx);
		Qy = join_rows(Qy,fy);
		Qz = join_rows(Qz,fz);
		
        rQx = Qx;
        rQy = Qy;
        rQz = Qz;
        
        Len.push_back(len);
        
		n++;
    
		
    } else {
        cout << "invalid Axis" << endl;
    }
};


void body::addBone(int in, double length){
    
    vector <int> ax;
    body::addAxis(xyz[in][2],0.,0.,0.,0.);
    
    ax.push_back(n+0);                          //rot rel to parent
    ax.push_back(n+1);                          //rot orth to parent
    ax.push_back(n+2);                          //rot rel to self
    xyz.push_back(ax);
    
    body::addAxis(n-1,0.,.5*M_PI,0.,0.);
    body::addAxis(n-1,0.,0.,0.,0.);
    body::addAxis(n-1,0.,0.,0.,length);
    
    nBone++;
    boneLength.push_back(length);
};




void body::buildPath(){
	
    /*
     Iteratively construct the path out from each Axis i:
     PathO[i][ordered indices from Axis i to all of its children]
    */
    
    pathO.resize(n);
    for (int i=0; i<pathI.size(); i++){
        for (int j=0; j<pathI[i].size(); j++){
            pathO[pathI[i][j]].push_back(i);
        }
        pathO[i].push_back(i);
    }
    
    boneX.resize(3,n);
    
    rDecay.resize(3,n);
    rDecayPrev.resize(3,n);
    
    for(int j=0;j<xyz.size();j++){
        int i=xyz[j][2];
        axisIn.push_back(i);
        axisOut.push_back(pathI[i][pathI[i].size()-1]);
    }

};

void body::rotateDecay(int k, int axis, double theta, double decayRate){
    rDecayPrev(axis,k) = rDecay(axis,k);
    rDecay(axis,k) += theta;
    rDecay(axis,k) -= decayRate*rDecay(axis,k);
    rotate(k,axis,rDecay(axis,k)-rDecayPrev(axis,k));
}

void body::rotate(int k, int axis, double theta){
	
    /*
     Rotate Axis k and each of its children in turn about axis 0, 1, 2
     */
    
	double ax, ay, az;
    double c = cos(theta);
    double s = sin(theta);
    double C = 1.-c;
    
	switch (axis){
		case 0:
			ax = rQx(0,k);
			ay = rQx(1,k);
			az = rQx(2,k);
			break;
		case 1:
			ax = rQy(0,k);
			ay = rQy(1,k);
			az = rQy(2,k);
			break;
		case 2:
			ax = rQz(0,k);
			ay = rQz(1,k);
			az = rQz(2,k);
			break;
	}
    
     // Consruct Quaternion required for local rotation (about the axes of the parent)
	arma::mat Qu;
	Qu << c+ax*ax*C << ax*ay*C-az*s << ax*az*C+ay*s << endr
	<< ay*ax*C+az*s << c+ay*ay*C << ay*az*C-ax*s << endr
	<< az*ax*C-ay*s << az*ay*C+ax*s << c+az*az*C << endr;
	
    arma::mat XP  = boneX;
    arma::mat QxP = rQx;
    
    // Rotate each of k's children in turn
	for(int l=1; l<pathO[k].size(); l++){
		int j = pathO[k][l];
		QxP.col(j) = Qu * QxP.col(j);
        XP.col(j) = XP.col(In[j]) + Len[j] * QxP.col(j);
	}
    
    Qx = rQx;
    Qy = rQy;
    Qz = rQz;
    
    // Rotate the x,y,z axes of each Axis too
    for(int l=1; l<pathO[k].size(); l++){
        int j = pathO[k][l];
        Qx.col(j) = Qu * Qx.col(j);
        Qy.col(j) = Qu * Qy.col(j);
        Qz.col(j) = Qu * Qz.col(j);
    }
    
    // Update Axis positions
    for (int j=0; j<n; j++){
        boneX.col(j) = boneX.col(In[j]) + Len[j] * Qx.col(j);
    }
    
    rQx = Qx;
    rQy = Qy;
    rQz = Qz;
};

void body::rotateBone(int l, int axis, double theta){
    int i = xyz[l][axis];
    int k = pathI[i][pathI[i].size()-1];
    rotate(k,0,theta);
}

void body::rotateBoneDecay(int l, int axis, double theta, double decayRate){
    int i = xyz[l][axis];
    int k = pathI[i][pathI[i].size()-1];
    rDecayPrev(0,k) = rDecay(0,k);
    rDecay(0,k) += theta;
    rDecay(0,k) -= decayRate*rDecay(0,k);
    rotate(k,0,rDecay(0,k)-rDecayPrev(0,k));
}

void body::translate(double x, double y, double z){
    arma::mat T;
    T<<x<<endr<<y<<endr<<z<<endr;
    for (int j=0; j<n; j++){
        boneX.col(j) += T;
    }
};

void body::setOriginalPose(void){
    rQxOrg = rQx;
    rQyOrg = rQy;
    rQzOrg = rQz;
    rDecayOrg = rDecay;
    rDecayPrevOrg = rDecayPrev;
}

void body::rotateAbs(int axis, double theta){
    
	arma::mat rQxCp = rQx;
	arma::mat rQyCp = rQy;
	arma::mat rQzCp = rQz;
	
	rQx = rQxOrg;
	rQy = rQyOrg;
	rQz = rQzOrg;
		
	rotate(1,axis,theta);
	
	rQxOrg = rQx;
	rQyOrg = rQy;
	rQzOrg = rQz;
	
	rQx = rQxCp;
	rQy = rQyCp;
	rQz = rQzCp;
		
}

void body::setLocation(double x, double y, double z){
    arma::mat T;
    T<<x<<endr<<y<<endr<<z<<endr;
    arma::mat off = boneX.col(0);
    for (int j=0; j<n; j++){
        boneX.col(j) += T - off;
    }
};

void body::copyState(void){
    XCopy = boneX;
    rQxCopy = rQx;
    rQyCopy = rQy;
    rQzCopy = rQz;
    rDecayCopy = rDecay;
    rDecayPrevCopy = rDecayPrev;
};

void body::undoPose(void){
    boneX = XCopy;
    rQx = rQxCopy;
    rQy = rQyCopy;
    rQz = rQzCopy;
    rDecay = rDecayCopy;
    rDecayPrev = rDecayPrevCopy;

};

void body::resetRot(int k){
	for(int l=1; l<pathO[k].size(); l++){
		int j = pathO[k][l];
		rQx.col(j) = rQxOrg.col(j);    		
		rQy.col(j) = rQyOrg.col(j);    		
		rQz.col(j) = rQzOrg.col(j);    		
						
	}
};

void body::resetPose(void){
	rQx = rQxOrg;
    rQy = rQyOrg;
    rQz = rQzOrg;
    rDecay = rDecayOrg;
    rDecayPrev = rDecayPrevOrg;
};


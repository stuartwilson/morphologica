#include "animatTools.h"
#include <math.h>

namespace animatTools{
    
    using namespace std;
    
    vector<vector<double> > getContig(int nT,int nX,double z0, int style, double param){

        int N = nX*nT;
        
        vector<vector<double> > X(N);
        
        for(int i=0;i<N;i++){
            X[i].resize(4);
        }
        
        double z1=z0;
        double tolerance=1e-6;
        int k=0;
        for(int i=0;i<nX;i++){
            double Z = 1e9;
            double R0,R1,dR,dZ,r0,r1;
            while(Z*Z>tolerance){
                /*
                 should be possible to use any sensible function
                 by changing only the next two lines,
                 e.g., R0 = sqrt(z0) for parabola, or R0 = 1. for cylinder.
                 */
                //R0=0.05;
                //R1=0.05;
                switch(style){ // cylinder
                    case(0):{
                        R0 = 1.;
                        R1 = 1.;
                    } break;
                    case(1):{ // step (continuous)
                        R0 = tanh(z0*param);
                        R1 = tanh(z1*param);
                    } break;
                    case(2):{ // parabola
                        R0 = sqrt(z0);
                        R1 = sqrt(z1);
                    } break;
                }
                
                dR = R0-R1;
                dZ = z0-z1;
                r0 = R0*M_PI/(double)nT;
                r1 = R1*M_PI/(double)nT;
                Z = sqrt(dR*dR+dZ*dZ)-(r0+r1);
                z1+=tolerance;
            }
            
            for(int j=0;j<nT;j++){
                double theta = -2.*M_PI*((double)j+0.5)/(double)nT;
                X[k][0]=R0*cos(theta);
                X[k][1]=R0*sin(theta);
                X[k][2]=z0;
                X[k][3]=r0;
                k++;
            }
            z0 = z1;
        }
        
        return X;
        
    };
    
   
    std::vector<std::vector<double> > align(std::vector<std::vector<double> > X, double xA,double yA,double zA,double xB,double yB, double zB){
        

        // Rotate cylinder
        double dx=xB-xA;
        double dy=yB-yA;
        double dz=zB-zA;
        
        // Hack cos cylinder is undefined in XY-plane
        if (fabs(dz)<1e-6){
            dz=1e-6;
        }
        double d=sqrt(dx*dx+dy*dy+dz*dz);
        double th=acos(dz/d);
        if(dz<0.){
            th*=-1.;
        }
        double ax = -dy*dz;
        double ay = dx*dz;
        double az = 0.;
        double v=sqrt(ax*ax + ay*ay + az*az);
        if(v==0.){
            if(dy>0.){
                ax=-1.;
                ay=0.;
            }else if(dy<0.){
                ax=1.;
                ay=0.;
            }else if(dx<0.){
                ax=0.;
                ay=-1.;
            }else{
                ax=0.;
                ay=1.;
            }
            az=0.;
        }else{
            ax/=v;
            ay/=v;
            az/=v;
        }
        double c=cos(th);
        double s=sin(th);
        double C=1.-c;
        
        arma::mat Qu;
        Qu << c+ax*ax*C << ax*ay*C-az*s << ax*az*C+ay*s << endr
        << ay*ax*C+az*s << c+ay*ay*C << ay*az*C-ax*s << endr
        << az*ax*C-ay*s << az*ay*C+ax*s << c+az*az*C << endr;
        
        int N = X.size();
        
        double minz = 1e7;//X[0][2];
        double maxz = -1e7;
        for(int i=0;i<N;i++){
            if(X[i][2]<minz){
                minz = X[i][2];
            }
            if(X[i][2]>maxz){
                maxz = X[i][2];
            }
        }
        for(int i=0;i<N;i++){
            X[i][2]-=minz;
        }
        double nm = d/(maxz-minz);
        
        arma::mat Cloud;
        Cloud.resize(3,N);
        for(int i=0;i<N;i++){
            Cloud(0,i)=X[i][0]*nm;
            Cloud(1,i)=X[i][1]*nm;
            Cloud(2,i)=X[i][2]*nm;
        }
        
        arma::mat CloudR = Qu*Cloud;
        
        std::vector<std::vector<double> > Xp = X;
        for(int i=0;i<N;i++){
            Xp[i].resize(4);
            Xp[i][0]=xA+CloudR(0,i);
            Xp[i][1]=yA+CloudR(1,i);
            Xp[i][2]=zA+CloudR(2,i);
            //Xp[i][3]=X[i][3];
        }
        return Xp;
        
    };

    vector<vector<double> > torus(int nI, int nJ, double aspectRatio, double radiusOuter, double radiusInner){
        // aspect ratio - usually >1
        // radius outer, e.g., 0.2
        // radius inner, use 0 for closed
        vector<vector<double> > tor;
        vector<double> x(3,0);
        double scI = 2.0*M_PI/(double)nI;
        double scJ = 2.0*M_PI/(double)nJ;
        for(int i=0;i<nI;i++){
            x[2]=radiusOuter*sin((double)i*scI)*aspectRatio;
            double r=radiusOuter*cos((double)i*scI)+radiusInner;
            for(int j=0;j<nJ;j++){
                x[0] = r*cos((double)j*scJ);
                x[1] = r*sin((double)j*scJ);
                tor.push_back(x);
            }
        }
        return tor;
    };
    
    

    double wrapAngle(double a){
        return a-6.283185307179586*floor(a/6.283185307179586);
    };
}




    
    
    
    
    
    
    /*
    vector<vector<vector<double> > > contig(int nT,int nX,double z0){
        
        vector<vector<vector<double> > > X;
        X.resize(nX);
        for(int i=0;i<nX;i++){
            X[i].resize(nT);
        }
        for(int i=0;i<nX;i++){
            for(int j=0;j<nT;j++){
                X[i][j].resize(4);
            }
        }
        
        double z1 = z0;
        double tolerance = 1e-6;
        
        for(int i=0;i<nX;i++){
            double Z = 1e9;
            double R0,R1,dR,dZ,r0,r1;
            while(Z*Z>tolerance){
                
                 //should be possible to use any sensible function
                 //by changing only the next two lines,
                 //e.g., R0 = sqrt(z0) for parabola, or R0 = 1. for cylinder.
                 
                R0 = tanh(z0*2.);
                R1 = tanh(z1*2.);
                
                dR = R0-R1;
                dZ = z0-z1;
                r0 = R0*M_PI/(double)nT;
                r1 = R1*M_PI/(double)nT;
                Z = sqrt(dR*dR+dZ*dZ)-(r0+r1);
                z1+=tolerance;
            }
            
            for(int j=0;j<nT;j++){
                double theta = -2.*M_PI*((double)j+0.5)/(double)nT;
                X[i][j][0]=R0*cos(theta);
                X[i][j][1]=R0*sin(theta);
                X[i][j][2]=z0;
                X[i][j][3]=r0;
            }
            z0 = z1;
        }
        return X;
    } 
     */

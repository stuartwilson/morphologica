#include "basic.h"
#include "animatTools.h"
#include  <math.h>

namespace basic{
    
    animat getAnimat(int n, world * W){
        
        animat A = animat(W,n,1,0.25,0.05,0.5);
        
        /* SHAPE-MATCH PARAMETERS */
        A.params.alpha = 0.2; // was 0.2
        A.params.beta = 0.5;
        
        A.params.type = Quadratic;
        A.params.dt = 0.01;
        A.params.allowFlip = true;
        A.params.volumeConservation = true;
        
        // BASIC SPREAD-EAGLE LAYOUT
        int limbs = 6;
        int pL = 1; // bones per limb
        double scale = 0.3;
        
        double boneMass = 1.;
        vector<vector<double> > lens(limbs);
        for(int i=0;i<limbs;i++){
            lens[i].resize(pL,0.0);        }
        
        lens[0][0]  =0.1;
        lens[1][0]  =-0.1;
        lens[2][0]  =0.1;
        lens[3][0]  =-0.1;
        lens[4][0]  =0.1;
        lens[5][0]  =-0.1;
        
        {
            int k=1;
            for(int i=0;i<limbs;i++){
                A.addBone(0,lens[i][0]*scale);
                k++;
                for(int j=0;j<pL-1;j++){
                    A.addBone(k-1,lens[i][j]*scale);
                    k++;
                }
            }
            
        }
        A.buildPath();
        
        // ALIGN TO XY PLANE
        A.rotate(0,1,M_PI*0.5);
        
        A.rotateBone(1,1,M_PI*0.5);
        A.rotateBone(2,1,M_PI*0.5);
        A.rotateBone(3,0,M_PI*0.5);
        A.rotateBone(4,0,M_PI*0.5);
        A.rotateBone(3,1,M_PI*0.5);
        A.rotateBone(4,1,M_PI*0.5);
        A.setOriginalPose();
        
        // ATTACH ALL BONES TO CLUSTER 0
        for(int i=0;i<limbs;i++){
            for(int j=0;j<pL;j++){
                A.attachBoneEnd(1+pL*i+j,true ,0,boneMass,0.01);
                A.attachBoneEnd(1+pL*i+j,false,0,boneMass,0.01);
            }
        }
        A.init();
        A.precompute();
        return A;
        
    }
    
    
    
    
    
    void moveAnimat(animat& a,int cycleLength,double amplitude,vector<double> phaseDiffs,vector<double> offset){

/*
        double Tp = phaseDiffs[0];
        double sig = phaseDiffs[1];//
        double v1 = phaseDiffs[2];
        
        vector<double> cog = a.getCoG();
        
        a.orientation = atan2(cog[1]-a.Xnow[0](1),cog[0]-a.Xnow[0](0));
        
        double cs = cos(-a.orientation);
        double ss = sin(-a.orientation);
        
        double TlExternal = 0.;
        double TrExternal = 0.;
        int countExternal = 0;
        double cL = 0.;
        double cR = 0.;
        for(int i=0;i<a.nP;i++){
            
            double dx = cog[0]-a.Xnow[i](0);
            double dy = cog[1]-a.Xnow[i](1);
            
            if(dx*ss+dy*cs > 0.0){
                TrExternal += a.Tc[i];
                countExternal++;
                cR += 1.0;
                //a.labels[i] +=1 ; // FLAG FOR PLOTTING
            } else {
                TlExternal += a.Tc[i];
                a.labels[i] += 1;
                cL += 1.0;
            }
            
        }
        
        double cRnormE = (double)countExternal/(double)a.nP;
        TrExternal /= cRnormE;
        TlExternal /= (1.-cRnormE);
        
        a.Text[0] = TlExternal;
        a.Text[1] = TrExternal;
        
        // Do control based on External co-ordinates
        a.dOrientation = 0.;
        if(countExternal>0){
            if(fabs(TrExternal-TlExternal)>0.){
                a.dOrientation = atan(8.0*(a.Tb-Tp)*(TrExternal-TlExternal));
            }
        }
        */
        
        /*
         double OR = a.orientation + a.dOrientation;
         
         // TRANSLATE ALLOCENTRIC ORIENTATION INTO EGOCENTRIC QUATERNION
         
         // Use one of the 'bones' as a reference vector (i.e, which way is up for the animat!)
         double VrefX = a.Xnow[a.boneIDs[0]](0)-a.Xnow[a.boneIDs[1]](0);
         double VrefY = a.Xnow[a.boneIDs[0]](1)-a.Xnow[a.boneIDs[1]](1);
         double VrefZ = a.Xnow[a.boneIDs[0]](2)-a.Xnow[a.boneIDs[1]](2);
         double VrefN = 1./sqrt(VrefX*VrefX+VrefY*VrefY+VrefZ*VrefZ);
         VrefX *= VrefN;
         VrefY *= VrefN;
         VrefZ *= VrefN;
         
         // Use the allocentrically computed direction as a target vector
         double VtarX = cos(OR);//cos(a.orientation);
         double VtarY = sin(OR); //sin(a.orientation);
         double VtarZ = 0.;
         
         // Do the cross-product to get axis normal to plane between reference and target
         double nx = VrefY*VtarZ-VrefZ*VtarY;
         double ny = VrefZ*VtarX-VrefX*VtarZ;
         double nz = VrefX*VtarY-VrefY*VtarX;
         double nl = 1./sqrt(nx*nx+ny*ny+nz*nz);
         
         // Do dot product to get angle between reference and target
         double rotAng = acos(VrefX*VtarX+VrefY*VtarY+VrefZ*VtarZ);
         
         
         // TARGET WEIGHTS
         
         //Maintain the axis and angle of rotation as the desired transformation, i.e., the teaching signal
         
         a.VrefX = VrefX;
         a.VrefY = VrefY;
         a.VrefZ = VrefZ;
         a.targetX = ((nx*nl)+1)/2.;
         a.targetY = ((ny*nl)+1)/2.;
         a.targetZ = ((nz*nl)+1)/2.;
         a.targetT = rotAng / M_PI;
         
         double tX = a.targetX*2.0-1;
         double tY = a.targetY*2.0-1;
         double tZ = a.targetZ*2.0-1;
         double tT = a.targetT*M_PI;
         
         // Weights then define a quaternion for computing force vector from egocentric reference
         double c = cos(tT);
         double s = sin(tT);
         double C = 1.-c;
         arma::mat Q, X;
         Q << c+tX*tX*C << tX*tY*C-tZ*s << tX*tZ*C+tY*s << endr
         << tY*tX*C+tZ*s << c+tY*tY*C << tY*tZ*C-tX*s << endr
         << tZ*tX*C-tY*s << tZ*tY*C+tX*s << c+tZ*tZ*C << endr;
         X <<a.VrefX<<endr<<a.VrefY<<endr<<a.VrefZ<<endr;
         X = Q * X;
         
         // just for plotting really
         a.targetVecX = X(0);
         a.targetVecY = X(1);
         a.targetVecZ = X(2);
         */
        
        
    }
    
}
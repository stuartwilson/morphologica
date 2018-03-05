#include "shrinkWrap.h"
#include "animatTools.h"
#include  <math.h>

namespace shrinkWrap{
    
    animat getAnimat(int n, world * W){
        
        //animat A = animat(n, W);
        animat A = animat(W,n,1,1.0,0.1,37.0);
        
        /* SHAPE-MATCH PARAMETERS */
        A.params.alpha = 0.1;
        A.params.beta = 0.5;
        
        A.params.type = Quadratic;
        A.params.dt = 0.01;
        A.params.allowFlip = true;
        A.params.volumeConservation = true;
        
        // BASIC SPREAD-EAGLE LAYOUT
        int limbs = 6;
        int pL = 6; // bones per limb
        double scale = 0.3;
        
        vector<vector<double> > lens(limbs), aspect(limbs);
        vector<vector<int> > nI(limbs), nJ(limbs);
        for(int i=0;i<limbs;i++){
            lens[i].resize(pL,0.0);
            nI[i].resize(pL,0);
            nJ[i].resize(pL,0);
            aspect[i].resize(pL,0);
        }
        
        // HEAD
        lens[0][0]  =0.3;//=0.5;
        lens[0][1]  =1.0;//=0.0;
        lens[0][2]  =0.0;//=0.4;
        lens[0][3]  =0.0;//=0.3;
        lens[0][4]  =0.0;//=0.1;
        lens[0][5]  =0.0;//=0.0;
        aspect[0][0]=1.0;//=1.5;
        aspect[0][1]=2.0;//=1.5;
        aspect[0][2]=1.0;//=1.5;
        aspect[0][3]=1.0;//=1.5;
        aspect[0][4]=1.0;//=0.5;
        aspect[0][5]=1.0;//=1.0;
        nI[0][0]=0;  nJ[0][0]=0;
        nI[0][1]=16; nJ[0][1]=16;
        nI[0][2]=0; nJ[0][2]=0;
        nI[0][3]=0; nJ[0][3]=0;
        nI[0][4]=0; nJ[0][4]=0;
        nI[0][5]=0; nJ[0][5]=0;

        // TAIL
        lens[3][0]=1.0;//=0.5;
        lens[3][1]=0.0;//=0.0;
        lens[3][2]=0.0;//=0.0;
        lens[3][3]=0.0;//=0.0;
        lens[3][4]=0.0;//=0.0;
        lens[3][5]=0.0;//=0.0;
        aspect[3][0]=2.0;//=5;
        aspect[3][1]=1.0;//=5;
        aspect[3][2]=1.0;//=5;
        aspect[3][3]=1.0;//=5;
        aspect[3][4]=1.0;//=5;
        aspect[3][5]=1.0;//=5;
        nI[3][0]=64;  nJ[3][0]=32;
        nI[3][1]=0; nJ[3][1]=0;
        nI[3][2]=0; nJ[3][2]=0;
        nI[3][3]=0; nJ[3][3]=0;
        nI[3][4]=0; nJ[3][4]=0;
        nI[3][5]=0; nJ[3][5]=0;

        // FRONT LEGS
        lens[1][0]=0.5;//=0.5;
        lens[1][1]=0.1;//=0.1;
        lens[1][2]=0.1;//=0.4;
        lens[1][3]=0.1;//=0.3;
        lens[1][4]=0.1;//=0.2;
        lens[1][5]=0.1;//=0.2;
        aspect[1][0]=1.0;//=2.0;
        aspect[1][1]=1.0;//=2.0;
        aspect[1][2]=1.0;//=2.0;
        aspect[1][3]=1.0;//=2.0;
        aspect[1][4]=1.0;//=2.0;
        aspect[1][5]=1.0;//=2.0;
        nI[1][0]=0;  nJ[1][0]=0;
        nI[1][1]=0; nJ[1][1]=0;
        nI[1][2]=0; nJ[1][2]=0;
        nI[1][3]=0; nJ[1][3]=0;
        nI[1][4]=0; nJ[1][4]=0;
        nI[1][5]=0; nJ[1][5]=0;

        // BACK LEGS
        lens[2][0]=0.5;//=0.5;
        lens[2][1]=0.1;//=0.1;
        lens[2][2]=0.1;//=0.4;
        lens[2][3]=0.1;//=0.3;
        lens[2][4]=0.1;//=0.2;
        lens[2][5]=0.1;//=0.2;
        aspect[2][0]=1.0;
        aspect[2][1]=1.0;
        aspect[2][2]=1.0;
        aspect[2][3]=1.0;
        aspect[2][4]=1.0;
        aspect[2][5]=1.0;
        nI[2][0]=0;  nJ[2][0]=0;
        nI[2][1]=0; nJ[2][1]=0;
        nI[2][2]=0; nJ[2][2]=0;
        nI[2][3]=0; nJ[2][3]=0;
        nI[2][4]=0; nJ[2][4]=0;
        nI[2][5]=0; nJ[2][5]=0;

        // MAKE LEGS SYMMETRICAL
        lens[5]=lens[1];
        lens[4]=lens[2];
        nI[5]=nI[1];
        nI[4]=nI[2];
        nJ[5]=nJ[1];
        nJ[4]=nJ[2];
        aspect[5]=aspect[1];
        aspect[4]=aspect[2];


        
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

        /*
            SPREAD EAGLED
         */
        
        // ALIGN TO XY PLANE
        A.rotate(0,1,M_PI*0.5);

        // SPREAD EAGLE
        A.rotateBone(1+pL*0,1,2.*M_PI*0./6);
        A.rotateBone(1+pL*1,1,2.*M_PI*1./6);
        A.rotateBone(1+pL*2,1,2.*M_PI*2./6);
        A.rotateBone(1+pL*3,1,2.*M_PI*3./6);
        A.rotateBone(1+pL*4,1,2.*M_PI*4./6);
        A.rotateBone(1+pL*5,1,2.*M_PI*5./6);
        
        // RE-ALIGNMENTS OF AXES (no visible rotations)
        
        A.rotateBone(1+pL*0+1,2,M_PI*-0.5);    // head to rearing axis
        A.rotateBone(1+pL*3+1,2,M_PI*-0.5);    // tail to rearing axis

        A.rotateBone(1+pL*1+4,0,M_PI*+0.5);    // pre-rotate 1/2
        A.rotateBone(1+pL*1+2,2,M_PI*+0.5);    // pre-rotate 2/2
        
        A.rotateBone(1+pL*2+4,0,M_PI*+0.5);    // pre-rotate 1/2
        A.rotateBone(1+pL*2+2,2,M_PI*+0.5);    // pre-rotate 2/2
    
        A.rotateBone(1+pL*4+4,0,M_PI*+0.5);    // pre-rotate 1/2
        A.rotateBone(1+pL*4+2,2,M_PI*+0.5);    // pre-rotate 2/2
   
        A.rotateBone(1+pL*5+4,0,M_PI*+0.5);    // pre-rotate 1/2
        A.rotateBone(1+pL*5+2,2,M_PI*+0.5);    // pre-rotate 2/2
        
        A.setOriginalPose();
        
        // ATTACH ALL BONES TO CLUSTER 0
        for(int i=0;i<limbs;i++){
            for(int j=0;j<pL;j++){
                A.attachBoneEnd(1+pL*i+j,true ,0,1.,0.01);
                A.attachBoneEnd(1+pL*i+j,false,0,1.,0.01);
            }
        }
        A.init();
        
        int k=0;
        for(int i=0;i<limbs;i++){
            for(int j=0;j<pL;j++){
                vector<vector<double> > boneMesh = animatTools::torus(nI[i][j],nJ[i][j],aspect[i][j],1.,0.);
                for(int m=0;m<boneMesh.size();m++){
                    boneMesh[m].push_back(0.05*scale);
                }
                vec org=A.boneX.col(A.axisIn[A.attachedBones[k*2]]);
                vec end=A.boneX.col(A.axisOut[A.attachedBones[k*2+1]]);
                boneMesh=animatTools::align(boneMesh,org(0),org(1),org(2),end(0),end(1),end(2));
                
                for(int m=0;m<boneMesh.size();m++){
                    A.obstacles.push_back(boneMesh[m]);
                }
                k++;
            }
        }
        
        /*
        vector<vector<double> > tor = animatTools::torus(32,32,1.2,1.0,0.);
        for(int i=0;i<tor.size();i++){
            tor[i].push_back(0.05*scale);
        }
        tor = animatTools::align(tor,-lens[3][0]*scale,0,0,lens[0][0]*scale,0,0);
        for(int j=0;j<tor.size();j++){
            A.obstacles.push_back(tor[j]);
        }
         */
        
        /*
        for(int i=0;i<boneStart.size();i++){
            vector<vector<double> > boneMesh = animatTools::torus(16,16,1.5,0.2,0.);
            for(int i=0;i<boneMesh.size();i++){boneMesh[i].push_back(0.05);}
            vector<vector<double> > q =
            animatTools::align(boneMesh,boneStart[i](0),boneStart[i](1),boneStart[i](2),
                               boneEnd[i](0),boneEnd[i](1),boneEnd[i](2));
            
            for(int j=0;j<q.size();j++){
                A.obstacles.push_back(q[j]);
            }
        }
        */
        /*
        vector<vec> boneStart, boneEnd, boneLens;
        for(int i=0;i<A.attachedBones.size();i++){
            if(A.boneIsOrigins[i]){
                boneStart.push_back(A.boneX.col(A.axisIn[A.attachedBones[i]]));
                //boneLens.push_back(A.boneLen[A.attachedBones[i]]);
            }else{
                boneEnd.push_back(A.boneX.col(A.axisOut[A.attachedBones[i]]));
            }
        }
        
        for(int i=0;i<boneStart.size();i++){
            vector<vector<double> > boneMesh = animatTools::torus(16,16,1.5,0.2,0.);
            for(int i=0;i<boneMesh.size();i++){boneMesh[i].push_back(0.05);}
            vector<vector<double> > q =
            animatTools::align(boneMesh,boneStart[i](0),boneStart[i](1),boneStart[i](2),
                               boneEnd[i](0),boneEnd[i](1),boneEnd[i](2));

            for(int j=0;j<q.size();j++){
                A.obstacles.push_back(q[j]);
            }
        }
        */

        
        /*
         DO SHRINK-WRAP AROUND OBSTACLES
         */
        
        for(int i=0;i<A.nP;i++){
            A.Xorg[i] *= 1e-6; // need to maintain shape so can't just set to 0
        }
        A.precompute();
        
        // Collision detection / resolution
        //for(int t=0;t<250;t++){
        for(int t=0;t<250;t++){
            A.inflate();
            A.projectPositions();
            for(int i=0;i<A.nP;i++){
                for(int j=0;j<A.obstacles.size();j++){
                    double dx = A.Xnow[i](0)-A.obstacles[j][0];
                    double dy = A.Xnow[i](1)-A.obstacles[j][1];
                    double dz = A.Xnow[i](2)-A.obstacles[j][2];
                    double rAB = A.Ranimat[i]+A.obstacles[j][3];
                    double d = sqrt(dx*dx+dy*dy+dz*dz);
                    if(d<=rAB){
                        double halfOverlap=(rAB-d)*0.5;
                        A.Xnew[i][0]=A.Xnow[i][0]+dx*halfOverlap;
                        A.Xnew[i][1]=A.Xnow[i][1]+dy*halfOverlap;
                        A.Xnew[i][2]=A.Xnow[i][2]+dz*halfOverlap;
                    }
                }
 
            }
            A.integrate();
        }
        A.Xorg=A.Xnow;
        A.precompute();
        
        /*
         RE-ASSIGN CLUSTERS
         */
        vector<vec> XpreCluster = A.Xnow;
        A.clusters.clear();
        int J=6;
        A.clusters.resize(J);
        for(int j=0;j<J;j++){
            A.clusters[j].clear();
            for(int i=0;i<A.nP;i++){
                double a = A.sphereS[i][0];
                double seg = 2.*M_PI/(double)J;
                double b = seg*(double)j;
                double aw = a-2.*M_PI*floor(a/(2.*M_PI));
                double bw = b-2.*M_PI*floor(b/(2.*M_PI));
                double angdiff = M_PI-fabs(M_PI-fabs(aw-bw));
                
                if(angdiff<=seg*0.51){
                    A.clusters[j].push_back(i);
                }
            }
        }
        vector<vector<int> > clInds;
        vector<int> clInd(2);
        for(int i=0;i<6;i++){
            clInd[0]=i;
            clInd[1]=i;
            clInds.push_back(clInd);
        }
        
        for(int i=0;i<clInds.size();i++){
            A.attachBoneEnd(1+(pL)*clInds[i][1],true ,clInds[i][0],1.,0.);
            A.attachBoneEnd(1+(pL)*clInds[i][1],false,clInds[i][0],1.,0.);
            for(int j=0;j<pL-1;j++){
                A.attachBoneEnd(1+pL*clInds[i][1]+j+1,true ,clInds[i][0],1.,0.0); // bone point size
                A.attachBoneEnd(1+pL*clInds[i][1]+j+1,false,clInds[i][0],1.,0);
            }
        }
        
        A.init();
        for(int i=0;i<A.nP;i++){
            A.Xorg[i] = XpreCluster[i];
        }
        A.precompute();
        
        
        /* RETURN ANIMAT */
        return A;
        
    }
    
    void moveAnimat(animat& a,int cycleLength,double amplitude,vector<double> phaseDiffs,vector<double> offset){
        a.resetPose();
        
        /* DO ROATATIONS HERE */
        //a.rotateBone(2,1,sin(a.time/(double)cycleLength)*amplitude);
        double pL = 4;
        //a.rotateBone(1+pL*4+2,1,0.1);
        /* DO ROATATIONS HERE */
        
        a.applyRotations();
    }
    
    
}




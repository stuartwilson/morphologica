#ifndef ____body__
#define ____body__

#include <iostream>
#include <math.h>
#include <vector>
#include <armadillo>

using namespace std;
using namespace arma;


class body {
    
public:
    int n, nBone;
    bool visible;
    double zero;
    std::vector <double> scale, scaleOrg;
	std::vector < std::vector <int> > pathI, pathO;
    std::vector < int > In, axisIn, axisOut;
    std::vector < std::vector < int > > xyz;

    arma::mat boneX, rQx, rQy, rQz, Qx, Qy, Qz, rDecay, rDecayPrev;
    arma::mat XCopy, rQxCopy, rQyCopy, rQzCopy, rDecayCopy, rDecayPrevCopy;
	arma::mat rQxOrg, rQyOrg, rQzOrg, rDecayOrg, rDecayPrevOrg;
    
    body(void);
    virtual ~body();
    void addAxis(int, double, double, double, double);
    void addBone(int, double);
    void rotate(int,int,double);
    void rotateBone(int, int, double);
    void rotateDecay(int,int,double,double);
    void rotateBoneDecay(int,int,double,double);
    void translate(double, double, double);
    void setLocation(double, double, double);
    void buildPath(void);
    void undoPose(void); //
    void resetRot(int);
    void rotateAbs(int,double);
	void resetPose(void); //
	void setOriginalPose(void); //
	void copyState(void);
    std::vector < double> Len, boneLength;

};

#endif /* defined(____body__) */

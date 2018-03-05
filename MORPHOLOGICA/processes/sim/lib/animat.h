//
//  animat.h
//

#ifndef ____animat__
#define ____animat__

#include <vector>
#include <armadillo>
#include "shapeMatch.h"
#include "body.h"
#include "world.h"
#include "sphere.h"
#include "RitterSOM.h"

class animat: public RitterSOM, public body, public ShapeMatch, public sphere{
    
public:
    animat();
    animat(world *, int,int,double,double,double);
    virtual ~animat();
    void init(void);
    void init(vector<double>,vector<double>);
    void attachBoneEnd(int, bool, int, double, double);
    void reposition(vector<double>,vector<double>);
    void setOrientation(double);
    void applyRotations(void);

    world * W;
    
    double Tb, theta, orientation, dOrientation, k1, k2, G;
    
    int floating;
    std::vector<double> Tc, input;
    std::vector<double> Ranimat;
    std::vector<int> boneIDs;
    std::vector<int> attachedBones;
    std::vector<bool>boneIsOrigins;
    std::vector<unsigned int> contacts, labels;
    std::vector<std::vector<double> > obstacles;
    int time;
    int N, nP;
    std::vector<std::vector<int> > clusters;
private:
    std::vector<arma::vec> Xanimat;
    std::vector<double> Manimat;
    
};

#endif /* defined(____animat__) */

#ifndef ____mesh__
#define ____mesh__

#include <math.h>
#include <vector>

using namespace std;


class mesh{
public:
    std::vector< std::vector<double> > M;
    std::vector< std::vector< std::vector <double> > > X;
    int nI, nJ, form1, form2;
    double A1, B1, A2, B2;
    std::vector<double>rgb;
    std::vector<double>I;
    std::vector<double>J;
    std::vector<double>meshParams, meshParamsOrg, meshParamsSet;

    mesh(int,int,int,int,double,double,double,double,double,std::vector<double>);
    virtual ~mesh();
    void update(void);
    double polarEllipseFast(double,double,double);
    double polarEllipse(double,double,double,double,double,double);
    double monoFunction(int,double,double,double);
    
    std::vector< std::vector< std::vector <double> > > getAlignedMesh(double,double,double,double,double,double,double,double,double,double,double,double);
    
    
};

#endif /* defined(____mesh__) */

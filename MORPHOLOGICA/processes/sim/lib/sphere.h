#ifndef ____sphere__
#define ____sphere__

#include <math.h>
#include <vector>

using namespace std;

class sphere {
    
public:
    vector<vector<int> > sphereC, sphereM;
    vector<vector<double> > sphereX, sphereS;
    int sphereN;
    double radius;
    sphere(void);
    sphere(int);
    void initSphere(int,double);
    virtual ~sphere();

};

#endif /* defined(____sphere__) */

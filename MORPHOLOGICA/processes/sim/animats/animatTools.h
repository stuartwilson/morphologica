#include <vector>
#include <armadillo>

namespace animatTools{
    
    using namespace arma;
    
    std::vector<std::vector<double> > getContig(int,int,double,int,double);
    std::vector<std::vector<double> > align(std::vector<std::vector<double> >,double,double,double,double,double,double);
    std::vector<std::vector<double> > torus(int,int,double,double,double);
    double wrapAngle(double);
}
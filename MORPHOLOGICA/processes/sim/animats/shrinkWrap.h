#include "../lib/animat.h"
#include "../lib/world.h"

namespace shrinkWrap{
    
    animat getAnimat(int,world *);
    void moveAnimat(animat&,int,double,vector<double>,vector<double>);
}
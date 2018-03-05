///
///  @file RitterSOM.h
///  @brief Extended self-organising map according to Ritter's 1992 paper
///  @author Hiroki Urashima, Stuart P Wilson
///  @date 20/11/2013
///

#ifndef ____RitterSOM__
#define ____RitterSOM__

#include <vector>


using std::vector;

class RitterSOM{
        
public:
    RitterSOM(void);
    virtual ~RitterSOM();
    void initSOM(int,int);
    vector<double> computeResponse(vector<double>);
    int findWinner(vector<double>);
    void iterateSOM(vector<double>,double, bool);
    void iterateSOMsuper(vector<double>,vector<double>,bool);
    vector<double> calculateOutput(vector<double>);
    vector<int> maxInWeight(void);
    vector<int> maxOutWeight(void);
    void setParam(int,double);
    void updateDynamicVariables(void);
    void randomizeWeights(void);
    
public:

    vector <vector <double> > H;

    // DYNAMIC VARIABLES
    double epsilonINP;      // learning rate for INPUT sheet
    double epsilonOUT;      // learning rate for OUTPUT sheet
    double epsilonVAR;      // learning rate for range of noise
    double epsilonTHR;      // constant for threshold function
    
    double responseSigma;   // decay constant for response
    double nhoodSigma;      // decay constant for neighbourhood
    double radiusSigma;     // decay constant for cutoff

    double radius;          // cutoff radius

    vector<vector<double> > win;    // weights for input
    vector<vector<double> > wout;   // weights for output
    vector<double> noise;           // range of noise
    vector<double> threshold;       // mean increase of reward function

    int mapTime;                  // current iteration number
    double reward;          // reward value
    int nSOM;               // the number of neurons
    
    int nin;                // the number of input data
    int nout;               // the number of output data

private:
    
    double tF;              // total number of iterations
    vector<double> x;       // x coordinates of neurons
    vector<double> y;       // y coordinates of neurons
    
    
    double radiusPRE;
    double nhoodSigmaPRE;
    double radiusDecayPRE;
};


#endif /* defined(____RitterSOM__) */

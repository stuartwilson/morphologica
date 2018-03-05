#include "mesh.h"


mesh::mesh(int nI, int nJ, int form1, int form2, double len, double A1, double B1, double A2, double B2, std::vector<double> rgb){
    this->nI = nI;
    this->nJ = nJ;
    this->form1 = form1;
    this->form2 = form2;
    this->rgb = rgb;
    I.resize(nI);
    J.resize(nJ);
    M.resize(nI);
    X.resize(nI);
    for(int i=0;i<nI;i++){
        I[i] = (double)i/((double)nI-1.);
        M[i].resize(nJ);
    }
    for(int j=0;j<nJ;j++){
        J[j] = (double)j*2.*M_PI/(double)nJ;
    }
    X.resize(nI);
    for(int i=0;i<nI;i++){
        X[i].resize(nJ);
        for(int j=0;j<nJ;j++){
            X[i][j].resize(3);
        }
    }
    meshParams.resize(5);
    meshParams[0]=len;
    meshParams[1]=A1;
    meshParams[2]=B1;
    meshParams[3]=A2;
    meshParams[4]=B2;
    meshParamsOrg=meshParams;
    meshParamsSet=meshParams;
    
    mesh::update();
};

mesh::~mesh(void){

}

void mesh::update(void){
    double ax1, ax2;
    for(int i=0;i<nI;i++){
        ax1 = monoFunction(form1, I[i], meshParams[1], meshParams[2]);
        ax2 = monoFunction(form2, I[i], meshParams[3], meshParams[4]);
        for(int j=0;j<nJ;j++){M[i][j]=polarEllipseFast(J[j],ax1,ax2);}
    }
};

double mesh::polarEllipseFast(double theta, double a, double b){
    a=fmax(a,1e-7);
    b=fmax(b,1e-7);
    double r = (b*b-a*a)*cos(2.*theta)+a*a+b*b;
    return sqrt(2.)*a*b*sqrt(r)/(r);
};

double mesh::polarEllipse(double theta, double a, double b, double c, double r0, double t0){
    a=fmax(a,1e-7);
    b=fmax(b,1e-7);
    double p = r0*((b*b-a*a)*cos(theta+t0-2.*c)+(a*a+b*b)*cos(theta-t0));
    double r = (b*b-a*a)*cos(2.*theta-2.*c)+a*a+b*b;
    double q = sqrt(2.)*a*b*sqrt(r-2.*r0*r0*(sin(theta-t0))*(sin(theta-t0)));
    return((p+q)/r);
};

double mesh::monoFunction(int form, double x, double A, double B){
    switch(form){
        case 0:{ // line
            return A+(B-A)*x;
        } break;
        case 1:{ // step
            if (A<B){ return A+0.5*(B-A)*(tanh((2.*x-1.)*M_PI)+1.);} // rise
            else    { return B+0.5*(A-B)*(tanh((1.-x*2.)*M_PI)+1.);} // fall
        } break;
        case 2:{ // bowl
            if (A<B){return A+(B-A)*tanh(x*2.*M_PI);}              // rise
            else    {return B+(A-B)*tanh((1.-x)*2.*M_PI);}         // fall
        } break;
        default: //
            return A;
    }
};

std::vector< std::vector< std::vector <double> > > mesh::getAlignedMesh(double x, double y, double z,double rQx0, double rQy0, double rQz0, double rQx1, double rQy1, double rQz1, double rQx2, double rQy2, double rQz2){
    double PI2overnJ = 2.*M_PI/(double)nJ;
    double LenOvernIminus1 = meshParams[0]/(double)(nI-1);
    for(int j=0;j<nJ;j++){
        double J = (double)j*PI2overnJ;
        double cJ = cos(J);
        double sJ = sin(J);
        for(int i=0;i<nI;i++){
            double Ilen = (double)i*LenOvernIminus1;
            X[i][j][0]=x+M[i][j]*(rQz0*cJ+rQy0*sJ)+rQx0*Ilen;
            X[i][j][1]=y+M[i][j]*(rQz1*cJ+rQy1*sJ)+rQx1*Ilen;
            X[i][j][2]=z+M[i][j]*(rQz2*cJ+rQy2*sJ)+rQx2*Ilen;
        }
    }
    return X;
};


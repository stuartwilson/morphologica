#include "display.h"
#include <armadillo>

Gdisplay::Gdisplay(int myWindowSize, const char* title, double rhoInit, double thetaInit, double phiInit){
    
    GLint att[] = { GLX_RGBA, GLX_DEPTH_SIZE, 24, GLX_DOUBLEBUFFER, None };
    //GLint att[] = { GLX_RGBA, GLX_DEPTH_SIZE, 24, GLX_DOUBLEBUFFER, None };
    disp = XOpenDisplay(NULL);
    if(disp == NULL){
        printf("\n\tcannot connect to X server\n\n");
        exit(0);
    }
    root = DefaultRootWindow(disp);
    vi = glXChooseVisual(disp, 0, att);
    if(vi == NULL){
        printf("\n\tno appropriate visual found\n\n");
        exit(0);
    } else {
        printf("\n\tvisual %p selected\n", (void *)vi->visualid);
    }
    cmap = XCreateColormap(disp, root, vi->visual, AllocNone);
    swa.colormap = cmap;
    swa.event_mask = ExposureMask | KeyPressMask;
    win = XCreateWindow(disp, root, 0, 0, myWindowSize, myWindowSize, 0, vi->depth, InputOutput, vi->visual, CWColormap | CWEventMask, &swa);
    
    
    glc = glXCreateContext(disp, vi, NULL, GL_TRUE);
    XMapWindow(disp, win);
    XStoreName(disp, win, (char*) title);
    
    this->speed = 5.*acos(-1.)/180.; // in degrees
    this->rho = 2.5 + rhoInit;
    this->theta = (thetaInit + 0.5)*acos(-1.);
    this->phi   = (phiInit + 0.00000001)*acos(-1.);
    alpha = 0.;
    Z = 1.;
    
};

void Gdisplay::setTitle(char* title){
    XStoreName(disp, win, (char*) title);
};


void Gdisplay::closeDisplay(void){
    glXDestroyContext(disp, glc);
    XDestroyWindow(disp, win);
    XCloseDisplay(disp);
};


void Gdisplay::resetDisplay(std::vector <double> fix, std::vector <double> eye, std::vector <double> rot){
    glXMakeCurrent(disp, win, glc);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    
    gluPerspective(45.0,	//"Specifies the field of view angle, in degrees, in the y direction."
                   1.0,	    //"Aspect ratio of x to y"
                   0.1,		//"Specifies the distance from the viewer to the near clipping plane (always positive)."
                   20.0);   //"Specifies the distance from the viewer to the far clipping plane (always positive)."
    
    glClearColor(1.0, 1.0, 1.0, 1.0);
    
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    IK =0.;JL=0.;WS=0.;AD=0.;TG=0.;FH=0.;UO=0.;QE=0.;RY=0.;
    
    // event processing
    if(XCheckWindowEvent(disp, win, KeyPressMask, &xev)){
        
        if (XLookupString(&xev.xkey,text,1,&key,0)==1){
            switch(text[0]){
                    
                    // primary keyboard controls
                case 'i': IK = +1.;   break; //
                case 'k': IK = -1.;   break; //
                case 'j': JL = -1.;   break; //
                case 'l': JL = +1.;   break; //
                case 'u': UO = -1.;   break; //
                case 'o': UO = +1.;   break; //
                    
                    // secondary keyboard controls
                case 'w': WS = +1.;   break; //
                case 's': WS = -1.;   break; //
                case 'a': AD = -1.;	  break; //
                case 'd': AD = +1.;	  break; //
                case 'q': QE = -1.;   break; //
                case 'e': QE = +1.;   break; //
                    
                    // tertiary keyboard controls
                case 't': TG = +1.;   break; //
                case 'g': TG = -1.;   break; //
                case 'f': FH = -1.;	  break; //
                case 'h': FH = +1.;	  break; //
                case 'r': RY = -1.;   break; //
                case 'y': RY = +1.;   break; //
                    
            }
        }
        
    }
    
    rho     = fmax(rho-WS*speed,0.);
    phi     += IK*speed;
    theta   += JL*speed;
    
    // spherical coordinates (rho=distance,theta=azimuth, i.e., horizontal,phi=polar)
    gluLookAt(rho*sin(phi)*cos(theta), // eyex (co-ordinate of camera)
              rho*sin(phi)*sin(theta), // eyey
              rho*cos(phi),            // eyez
              0.,				       // centerx (eye fixation point, i.e., 0,0,0)
              0.,				       // centery
              0.,				       // centerz
              0.,				       // upx     (which way is up)
              0.,				       // upy
              -sin(phi));			   // upz
    
    
    /*
     if (XLookupString(&xev.xkey,text,1,&key,0)==1){
     switch(text[0]){
     case 's': rho += speed;	  break; // out
     case 'w': rho = fmax(0.,rho-speed);   break; // in
     case 'a': alpha -= speed;	  break; // forward
     case 'd': alpha += speed;	  break; // forward
     case 'i': phi += speed;	  break; // forward
     case 'k': phi -= speed;   break; // backwards
     case 'j': theta -= speed; break; // left
     case 'l': theta += speed; break; // right
     
     }
     }
     
     if (XLookupString(&xev.xkey,text,1,&key,0)==1){
     switch(text[0]){
     case 'i': X += speed*cos(theta); Y += speed*sin(theta); break; // forward
     case 'k': X -= speed*cos(theta); Y -= speed*sin(theta); break; // backward
     case 'j': theta += speed;	  break; // leftward
     case 'l': theta -= speed;	  break; // rightward
     case 'a': phi += speed;	  break; // look down
     case 'z': phi -= speed;   break; // look up
     // case 'q': closeDisplay(); break;
     }
     }
     */
    //}
    
    /*
     // spherical coordinates (rho=distance,theta=azimuth, i.e., horizontal,phi=polar)
     gluLookAt(eye[0],//+rho*sin(phi)*cos(theta), // eyex (co-ordinate of camera)
     eye[1],//+rho*sin(phi)*sin(theta), // eyey
     eye[2],//+rho*cos(phi),            // eyez
     fix[0],							// centerx (eye fixation point, i.e., 0,0,0)
     fix[1],							// centery
     fix[2],							// centerz
     rot[0],							// upx     (which way is up)
     rot[1],							// upy
     rot[2]);//-sin(phi));					// upz
     */
    
    /*
     // spherical coordinates (rho=distance,theta=azimuth, i.e., horizontal,phi=polar)
     gluLookAt(X, // eyex (co-ordinate of camera)
			  Y, // eyey
			  Z,            // eyez
			  X+cos(theta),                       // centerx (eye fixation point, i.e., 0,0,0)
			  Y+sin(theta),							// centery
			  Z+sin(phi),							// centerz
			  0.,							// upx     (which way is up)
			  0.,							// upy
			  1.);					// upz
     */
    glEnable(GL_DEPTH_TEST);
    
    glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);
    glEnable(GL_LIGHTING);

    GLfloat qaAmbientLight[] = {.2,.2,.2,1.};
    GLfloat qaDiffuseLight[] = {.8,.8,.8,1.};
    GLfloat qaSpecularLight[] = {1.,1.,1.,1.};
    GLfloat col[] = {0.,0.,0.,1.f};
    GLfloat blk[] = {0.,0.,0.,1.f};

    // TOP LIGHT
    glEnable(GL_LIGHT0);
    GLfloat qaLightPosition0[] = {.0,.0,5.,1.};
    glLightfv(GL_LIGHT0, GL_AMBIENT, qaAmbientLight);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, qaDiffuseLight);
    glLightfv(GL_LIGHT0, GL_SPECULAR, qaSpecularLight);
    glLightfv(GL_LIGHT0, GL_POSITION, qaLightPosition0);

    // TOP LIGHT
    glEnable(GL_LIGHT1);
    GLfloat qaLightPosition1[] = {.0,.0,-5.,1.};
    glLightfv(GL_LIGHT1, GL_AMBIENT, qaAmbientLight);
    glLightfv(GL_LIGHT1, GL_DIFFUSE, qaDiffuseLight);
    glLightfv(GL_LIGHT1, GL_SPECULAR, qaSpecularLight);
    glLightfv(GL_LIGHT1, GL_POSITION, qaLightPosition1);

    
    
    glMaterialfv(GL_FRONT, GL_DIFFUSE, col);
    glMaterialfv(GL_FRONT, GL_SPECULAR, blk);
    glMaterialf(GL_FRONT, GL_SHININESS, 0.);
    
};

void Gdisplay::redrawDisplay(){
    glXMakeCurrent(disp, win, glc);
    XGetWindowAttributes(disp, win, &gwa);
    glViewport(0, 0, gwa.width, gwa.height);
    glXSwapBuffers(disp, win);
};


void Gdisplay::drawHex(double x,double y,double z,double r,double red,double green,double blue){
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    GLfloat col[] = {red, green, blue, 1.f};
    //GLfloat wht[] = {1.f, 1.f, 1.f, 1.f};
    glMaterialfv(GL_FRONT, GL_DIFFUSE, col);
    
    //glColor3f(red,green,blue);
    double ry = r * 1.154700538379252; // r * 1.0/sin(pi/3.0)
    double hrx = r * 0.5;
    double hry = ry * 0.5;
    glBegin(GL_POLYGON);
    glVertex3f(x,y+ry,z);
    glVertex3f(x+r,y+hry,z);
    glVertex3f(x+r,y-hry,z);
    glVertex3f(x,y-ry,z);
    glVertex3f(x-r,y-hry,z);
    glVertex3f(x-r,y+hry,z);
    glVertex3f(x,y+ry,z);
    glEnd();
};

void Gdisplay::drawTri(std::vector <double> p1,std::vector <double> p2,std::vector <double> p3,double red,double green,double blue){
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); // fill
    glBegin(GL_TRIANGLES);
    glColor3f(red,green,blue);
    glVertex3f(p1[0], p1[1], p1[2]);
    glVertex3f(p2[0], p2[1], p2[2]);
    glVertex3f(p3[0], p3[1], p3[2]);
    glEnd();
}

void Gdisplay::drawSphere(double x,double y,double z,double r,std::vector<double> C,int res){
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); // fill
    //glColor3f(red,green,blue);
    GLfloat col[] = {C[0], C[1], C[2], 1.f};
    GLfloat wht[] = {1.f, 1.f, 1.f, 1.f};
    glMaterialfv(GL_FRONT, GL_DIFFUSE, col);
    glMaterialfv(GL_FRONT, GL_SPECULAR, wht);
    glMaterialf(GL_FRONT, GL_SHININESS, 60.);
    static GLUquadricObj* Sphere = gluNewQuadric();
    glPushMatrix();
    glTranslatef(x,y,z);
    gluSphere(Sphere, r, res, res);
    glPopMatrix();
};

void Gdisplay::drawLine(double ax,double ay,double az,double bx,double by,double bz,double red,double green,double blue,double width) {
    glColor3f(red,green,blue);
    glPointSize(width);
    glBegin(GL_LINES);
    glVertex3d(ax,ay,az);
    glVertex3d(bx,by,bz);
    glEnd();
};

void Gdisplay::addFloor(double x, double y){
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    GLfloat col[] = {0.92,0.94,0.96, 1.f};
    GLfloat wht[] = {1.f, 1.f, 1.f, 1.f};
    glMaterialfv(GL_FRONT, GL_DIFFUSE, col);
    glMaterialfv(GL_FRONT, GL_SPECULAR, wht);
    glMaterialf(GL_FRONT, GL_SHININESS, 60.);
    glBegin(GL_QUADS);
    glVertex3d(-x,-y,0.0);
    glVertex3d(-x,+y,0.0);
    glVertex3d(+x,+y,0.0);
    glVertex3d(+x,-y,0.0);
    glNormal3d(0.,0.,1.);
    glEnd();
    
}

/*
 This code is redundant because of the following function, but I'm keeping it in case its useful for future reference for drawing cylinders.
 
 void Gdisplay::drawCylinder(std::vector <double> A, std::vector <double> B, double h, double rBase, double rEnd, std::vector <double> col, int res){
 
 
 glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); // fill
 
 GLfloat colr[] = {col[0], col[1], col[2], 1.f};
 GLfloat wht[] = {1.f, 1.f, 1.f, 1.f};
 glMaterialfv(GL_FRONT, GL_DIFFUSE, colr);
 glMaterialfv(GL_FRONT, GL_SPECULAR, wht);
 glMaterialf(GL_FRONT, GL_SHININESS, 60.);
 
 double wA = rBase*acos(-1.0)/((double)res-1);
 double wB = rEnd*acos(-1.0)/((double)res-1);
 
 double xA, xB, yA, yB;
 
 arma::mat Ap (3, res);
 arma::mat An (3, res);
 arma::mat Bp (3, res);
 arma::mat Bn (3, res);
 
 for (int i=0; i<res; i++){
 
 double a = 2.0*acos(-1.0)*(float)i/(float)res;
 
 xA = cos(a);
 yA = sin(a);
 
 glNormal3d(xA,yA,0.);
 
 xB = xA*rEnd;
 yB = yA*rEnd;
 xA *= rBase;
 yA *= rBase;
 
 double aT = a+acos(-1.)*0.5;
 double c = cos(aT);
 double s = sin(aT);
 
 Ap(0,i) = xA + wA * c;
 Ap(1,i) = yA + wA * s;
 Ap(2,i) = 0.;
 Bp(0,i) = xB + wB * c;
 Bp(1,i) = yB + wB * s;
 Bp(2,i) = h;
 Bn(0,i) = xB - wB * c;
 Bn(1,i) = yB - wB * s;
 Bn(2,i) = h;
 An(0,i) = xA - wA * c;
 An(1,i) = yA - wA * s;
 An(2,i) = 0.;
 }
 
 for (int i=0; i<res; i++){
 glBegin(GL_POLYGON);
 glVertex3d(Ap(0,i)+A[0],Ap(1,i)+A[1],Ap(2,i)+A[2]);
 glVertex3d(Bp(0,i)+A[0],Bp(1,i)+A[1],Bp(2,i)+A[2]);
 glVertex3d(Bn(0,i)+A[0],Bn(1,i)+A[1],Bn(2,i)+A[2]);
 glVertex3d(An(0,i)+A[0],An(1,i)+A[1],An(2,i)+A[2]);
 glEnd();
 }
 }
 */
/*
 void Gdisplay::drawCylinder(float x1, float y1, float z1, float x2, float y2, float z2, float radA, float radB, int subdivisions, std::vector <double> col){
 
 GLUquadricObj *quadric=gluNewQuadric();
 gluQuadricNormals(quadric, GLU_SMOOTH);
 
 glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); // fill
 
 GLfloat colr[] = {col[0], col[1], col[2], 1.f};
 GLfloat wht[] = {1.f, 1.f, 1.f, 1.f};
 glMaterialfv(GL_FRONT, GL_DIFFUSE, colr);
 glMaterialfv(GL_FRONT, GL_SPECULAR, wht);
 glMaterialf(GL_FRONT, GL_SHININESS, 60.);
 
 
 float vx = x2-x1;
 float vy = y2-y1;
 float vz = z2-z1;
 float v = sqrt(vx*vx + vy*vy + vz*vz);
 float ax;
 
 if (fabs(vz) < 1.0e-3) {
 ax = 57.2957795*acos( vx/v ); // rotation angle in x-y plane
 if ( vy <= 0.0 )
 ax = -ax;
 }
 else {
 ax = 57.2957795*acos( vz/v ); // rotation angle
 if ( vz <= 0.0 )
 ax = -ax;
 }
 float rx = -vy*vz;
 float ry = vx*vz;
 glPushMatrix();
 
 //draw the cylinder body
 glTranslatef(x1,y1,z1);
 if (fabs(vz) < 1.0e-3){
 glRotated(90.0, 0, 1, 0.0); // Rotate & align with x axis
 glRotated(ax, -1.0, 0.0, 0.0); // Rotate to point 2 in x-y plane
 }
 else {
 glRotated(ax, rx, ry, 0.0); // Rotate about rotation vector
 }
 gluQuadricOrientation(quadric,GLU_OUTSIDE);
 gluCylinder(quadric, radA, radB, v, subdivisions, 1);
 
 //draw the first cap
 gluQuadricOrientation(quadric,GLU_INSIDE);
 gluDisk( quadric, 0.0, radA, subdivisions, 1);
 glTranslatef( 0,0,v );
 
 //draw the second cap
 gluQuadricOrientation(quadric,GLU_OUTSIDE);
 gluDisk( quadric, 0.0, radB, subdivisions, 1);
 glPopMatrix();
 }
 */
/*
 void Gdisplay::drawCylinder(float x1, float y1, float z1, float x2, float y2, float z2, float radA, float radB, int subdivisions, std::vector <double> col){
 
 GLUquadricObj *quadric=gluNewQuadric();
 gluQuadricNormals(quadric, GLU_SMOOTH);
 
 glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); // fill
 
 GLfloat colr[] = {col[0], col[1], col[2], 1.f};
 GLfloat wht[] = {1.f, 1.f, 1.f, 1.f};
 glMaterialfv(GL_FRONT, GL_DIFFUSE, colr);
 glMaterialfv(GL_FRONT, GL_SPECULAR, wht);
 glMaterialf(GL_FRONT, GL_SHININESS, 60.);
 
 double px = x1-x2;
 double py = y1-y2;
 double pz = z1-z2;
 double len = sqrt(px*px + py*py + pz*pz);
 
 glPushMatrix();
 
 if(len>0.&&radB>0.&&radA>0.){
 glTranslatef(x2,y2,z2);
 glRotatef((180./M_PI)*acos(pz/len),-py,+px,0.);
 gluQuadricOrientation(quadric,GLU_OUTSIDE);
 gluCylinder(quadric, radB, radA, len, subdivisions, 1);
 }
 
 //draw the first cap
 gluQuadricOrientation(quadric,GLU_INSIDE);
 gluDisk( quadric, 0.0, radB, subdivisions, 1);
 glTranslatef(0,0,len);
 
 //draw the second cap
 gluQuadricOrientation(quadric,GLU_OUTSIDE);
 gluDisk( quadric, 0.0, radA, subdivisions, 1);
 
 glPopMatrix();
 
 }
 */

// THIS ONE WORKS!
void Gdisplay::drawCylinder(float x1, float y1, float z1, float x2, float y2, float z2, float radA, float radB, int subdivisions, std::vector <double> col){
    
    GLUquadricObj *quadric=gluNewQuadric();
    gluQuadricNormals(quadric, GLU_SMOOTH);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); // fill
    GLfloat colr[] = {col[0], col[1], col[2], 1.f};
    GLfloat wht[] = {1.f, 1.f, 1.f, 1.f};
    glMaterialfv(GL_FRONT, GL_DIFFUSE, colr);
    glMaterialfv(GL_FRONT, GL_SPECULAR, wht);
    glMaterialf(GL_FRONT, GL_SHININESS, 60.);
    
    double px = x1-x2;
    double py = y1-y2;
    double pz = z1-z2;
    double len = sqrt(px*px+py*py+pz*pz);
    
    glPushMatrix();
    
    if (len>0.&&radB>0.&&radA>0.){
        
        glTranslatef(x2,y2,z2);
        glRotatef((180./M_PI)*acos(pz/len),-py+1e-6,+px,0.);
        
        gluQuadricOrientation(quadric,GLU_OUTSIDE);
        gluCylinder(quadric, radB, radA, len, subdivisions, 1);
        
        gluQuadricOrientation(quadric,GLU_INSIDE);
        gluDisk( quadric, 0.0, radB, subdivisions, 1);
        glTranslatef( 0,0,len );
        
        gluQuadricOrientation(quadric,GLU_OUTSIDE);
        gluDisk( quadric, 0.0, radA, subdivisions, 1);
        
    }
    glPopMatrix();
    
}

/*
 void Gdisplay::drawCylinder(float x1, float y1, float z1, float x2, float y2, float z2, float radA, float radB, int subdivisions, std::vector <double> col){
 
 GLUquadricObj *quadric=gluNewQuadric();
 gluQuadricNormals(quadric, GLU_SMOOTH);
 glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); // fill
 GLfloat colr[] = {col[0], col[1], col[2], 1.f};
 GLfloat wht[] = {1.f, 1.f, 1.f, 1.f};
 glMaterialfv(GL_FRONT, GL_DIFFUSE, colr);
 glMaterialfv(GL_FRONT, GL_SPECULAR, wht);
 glMaterialf(GL_FRONT, GL_SHININESS, 60.);
 
 double px = x1-x2;
 double py = y1-y2;
 double pz = z1-z2;
 double len = sqrt(px*px+py*py+pz*pz);
 
 glPushMatrix();
 
 if (len>0.&&radB>0.&&radA>0.){
 
 // obviously this will come out!
 int iN = subdivisions;
 int jN = subdivisions;
 std::vector<std::vector<double> > H(subdivisions);
 for(int i=0;i<iN;i++){
 H[i].resize(jN);
 for(int j=0; j<jN;j++){
 if (i<iN/2){
 H[i][j] = 0.;
 } else{
 H[i][j] = radA*sin(3.*(double)j/(double)jN);
 }
 }
 }
 
 glTranslatef(x2,y2,z2);
 glRotatef((180./M_PI)*acos(pz/len),-py+1e-6,+px,0.);
 
 double sub = 2.*M_PI/(double)subdivisions;
 double ap = 0.;
 double cosap = 1.;
 double sinap = 0.;
 double div = len/(double)jN;
 for (int i=1;i<=iN;i++){
 int imod = i%iN;
 double a = (double)imod*sub;
 double cosa = cos(a);
 double sina = sin(a);
 double jdip = 0.;
 double jdi;
 for (int j=1;j<jN;j++){
 jdi = j*div;
 glBegin(GL_QUADS);
 glVertex3d(H[imod][j-1]*cosa,H[imod][j-1]*sina,jdip);
 glVertex3d(H[imod][j]*cosa,H[imod][j]*sina,jdi);
 glVertex3d(H[imod][j]*cosap,H[imod][j]*sinap,jdi);
 glVertex3d(H[imod][j-1]*cosap,H[imod][j-1]*sinap,jdip);
 glNormal3d(cosa,sina,0.); // close enough for now
 jdip=jdi;
 glEnd();
 }
 
 ap = a;
 cosap = cosa;
 sinap = sina;
 }
 
 //glTranslatef(x2,y2,z2);
 //glRotatef((180./M_PI)*acos(pz/len),-py+1e-6,+px,0.);
 
 //gluQuadricOrientation(quadric,GLU_OUTSIDE);
 //gluCylinder(quadric, radB, radA, len, subdivisions, 1);
 
 //gluQuadricOrientation(quadric,GLU_INSIDE);
 //gluDisk( quadric, 0.0, radB, subdivisions, 1);
 //glTranslatef( 0,0,len );
 
 //gluQuadricOrientation(quadric,GLU_OUTSIDE);
 //gluDisk( quadric, 0.0, radA, subdivisions, 1);
 
 }
 glPopMatrix();
 
 }
 */

void Gdisplay::drawMesh(std::vector< std::vector< std::vector <double> > > X, std::vector<std::vector<std::vector<double> > > C){
    
    //GLUquadricObj *quadric=gluNewQuadric();
    //gluQuadricNormals(quadric, GLU_SMOOTH);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); // fill
    //GLfloat colr[] = {col[0], col[1], col[2], 1.f};
    GLfloat wht[] = {1.f, 1.f, 1.f, 1.f};
    //glMaterialfv(GL_FRONT, GL_DIFFUSE, colr);
    glMaterialfv(GL_FRONT, GL_SPECULAR, wht);
    glMaterialf(GL_FRONT, GL_SHININESS, 60.);
    
    
    // obviously this will come out!
    int iN = X.size();
    int jN = X[0].size();
    
    for (int i=1;i<iN;i++){
        for (int j=1;j<jN;j++){
            GLfloat colr[]={C[i][j][0],C[i][j][1],C[i][j][2],1.f};
            glMaterialfv(GL_FRONT, GL_DIFFUSE, colr);
            glBegin(GL_QUADS);
            double ax = X[i][j-1][0]-X[i-1][j-1][0];
            double ay = X[i][j-1][1]-X[i-1][j-1][1];
            double az = X[i][j-1][2]-X[i-1][j-1][2];
            double bx = X[i-1][j][0]-X[i-1][j-1][0];
            double by = X[i-1][j][1]-X[i-1][j-1][1];
            double bz = X[i-1][j][2]-X[i-1][j-1][2];
            double nx = ay*bz-az*by;
            double ny = az*bx-ax*bz;
            double nz = ax*by-ay*bx;
            double nl = -1./sqrt(nx*nx+ny*ny+nz*nz);
            glNormal3d(nx*nl,ny*nl,nz*nl);
            glVertex3d(X[i-1][j-1][0],X[i-1][j-1][1],X[i-1][j-1][2]);
            glVertex3d(X[i-1][j][0],X[i-1][j][1],X[i-1][j][2]);
            glVertex3d(X[i][j][0],X[i][j][1],X[i][j][2]);
            glVertex3d(X[i][j-1][0],X[i][j-1][1],X[i][j-1][2]);
            glEnd();
            
        }
        
        GLfloat colr[]={C[i][0][0],C[i][0][1],C[i][0][2],1.f};
        glMaterialfv(GL_FRONT, GL_DIFFUSE, colr);
        glBegin(GL_QUADS);
        double ax = X[i][jN-1][0]-X[i-1][jN-1][0];
        double ay = X[i][jN-1][1]-X[i-1][jN-1][1];
        double az = X[i][jN-1][2]-X[i-1][jN-1][2];
        double bx = X[i-1][0][0]-X[i-1][jN-1][0];
        double by = X[i-1][0][1]-X[i-1][jN-1][1];
        double bz = X[i-1][0][2]-X[i-1][jN-1][2];
        double nx = ay*bz-az*by;
        double ny = az*bx-ax*bz;
        double nz = ax*by-ay*bx;
        double nl = -1./sqrt(nx*nx+ny*ny+nz*nz);
        glNormal3d(nx*nl,ny*nl,nz*nl);
        glVertex3d(X[i-1][jN-1][0],X[i-1][jN-1][1],X[i-1][jN-1][2]);
        glVertex3d(X[i-1][0][0],X[i-1][0][1],X[i-1][0][2]);
        glVertex3d(X[i][0][0],X[i][0][1],X[i][0][2]);
        glVertex3d(X[i][jN-1][0],X[i][jN-1][1],X[i][jN-1][2]);
        glEnd();
        
    }
}

void Gdisplay::drawMesh2(std::vector< std::vector< std::vector <double> > > X, std::vector <double> col){
    
    //GLUquadricObj *quadric=gluNewQuadric();
    //gluQuadricNormals(quadric, GLU_SMOOTH);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); // fill
    GLfloat colr[] = {col[0], col[1], col[2], 1.f};
    GLfloat wht[] = {1.f, 1.f, 1.f, 1.f};
    glMaterialfv(GL_FRONT, GL_DIFFUSE, colr);
    glMaterialfv(GL_FRONT, GL_SPECULAR, wht);
    glMaterialf(GL_FRONT, GL_SHININESS, 60.);
    
    // obviously this will come out!
    int iN = X.size();
    int jN = X[0].size();
    
    for (int i=1;i<iN;i++){
        for (int j=1;j<jN;j++){
            glBegin(GL_QUADS);
            double ax = X[i][j-1][0]-X[i-1][j-1][0];
            double ay = X[i][j-1][1]-X[i-1][j-1][1];
            double az = X[i][j-1][2]-X[i-1][j-1][2];
            double bx = X[i-1][j][0]-X[i-1][j-1][0];
            double by = X[i-1][j][1]-X[i-1][j-1][1];
            double bz = X[i-1][j][2]-X[i-1][j-1][2];
            double nx = ay*bz-az*by;
            double ny = az*bx-ax*bz;
            double nz = ax*by-ay*bx;
            double nl = -1./sqrt(nx*nx+ny*ny+nz*nz);
            glNormal3d(nx*nl,ny*nl,nz*nl);
            glVertex3d(X[i-1][j-1][0],X[i-1][j-1][1],X[i-1][j-1][2]);
            glVertex3d(X[i-1][j][0],X[i-1][j][1],X[i-1][j][2]);
            glVertex3d(X[i][j][0],X[i][j][1],X[i][j][2]);
            glVertex3d(X[i][j-1][0],X[i][j-1][1],X[i][j-1][2]);
            glEnd();
            
        }
        
        glBegin(GL_QUADS);
        double ax = X[i][jN-1][0]-X[i-1][jN-1][0];
        double ay = X[i][jN-1][1]-X[i-1][jN-1][1];
        double az = X[i][jN-1][2]-X[i-1][jN-1][2];
        double bx = X[i-1][0][0]-X[i-1][jN-1][0];
        double by = X[i-1][0][1]-X[i-1][jN-1][1];
        double bz = X[i-1][0][2]-X[i-1][jN-1][2];
        double nx = ay*bz-az*by;
        double ny = az*bx-ax*bz;
        double nz = ax*by-ay*bx;
        double nl = -1./sqrt(nx*nx+ny*ny+nz*nz);
        glNormal3d(nx*nl,ny*nl,nz*nl);
        glVertex3d(X[i-1][jN-1][0],X[i-1][jN-1][1],X[i-1][jN-1][2]);
        glVertex3d(X[i-1][0][0],X[i-1][0][1],X[i-1][0][2]);
        glVertex3d(X[i][0][0],X[i][0][1],X[i][0][2]);
        glVertex3d(X[i][jN-1][0],X[i][jN-1][1],X[i][jN-1][2]);
        glEnd();
        
        
    }
    
    
    
}


/*
 void Gdisplay::drawTorus(std::vector< std::vector< std::vector <double> > > X, std::vector<std::vector<std::vector<double> > > C){
 
 //GLUquadricObj *quadric=gluNewQuadric();
 //gluQuadricNormals(quadric, GLU_SMOOTH);
 glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); // fill
 //GLfloat colr[] = {col[0], col[1], col[2], 1.f};
 GLfloat wht[] = {1.f, 1.f, 1.f, 1.f};
 //glMaterialfv(GL_FRONT, GL_DIFFUSE, colr);
 glMaterialfv(GL_FRONT, GL_SPECULAR, wht);
 glMaterialf(GL_FRONT, GL_SHININESS, 60.);
 
 
 // obviously this will come out!
 int iN = X.size();
 int jN = X[0].size();
 
 for (int i=1;i<=iN;i++){
 int I0 = i%iN;
 int I1 = (i-1)%iN;
 for (int j=1;j<=jN;j++){
 int J0 = j%jN;
 int J1 = (j-1)%jN;
 GLfloat colr[]={C[I0][J0][0],C[I0][J0][1],C[I0][J0][2],1.f};
 glMaterialfv(GL_FRONT, GL_DIFFUSE, colr);
 glBegin(GL_QUADS);
 
 double ax = X[I0][J1][0]-X[I1][J1][0];
 double ay = X[I0][J1][1]-X[I1][J1][1];
 double az = X[I0][J1][2]-X[I1][J1][2];
 double bx = X[I1][J0][0]-X[I1][J1][0];
 double by = X[I1][J0][1]-X[I1][J1][1];
 double bz = X[I1][J0][2]-X[I1][J1][2];
 double nx = ay*bz-az*by;
 double ny = az*bx-ax*bz;
 double nz = ax*by-ay*bx;
 double nl = 1./sqrt(nx*nx+ny*ny+nz*nz);
 
 
 glNormal3d(nx*nl,ny*nl,nz*nl);
 
 glVertex3d(X[I1][J1][0],X[I1][J1][1],X[I1][J1][2]);
 glVertex3d(X[I1][J0][0],X[I1][J0][1],X[I1][J0][2]);
 glVertex3d(X[I0][J0][0],X[I0][J0][1],X[I0][J0][2]);
 glVertex3d(X[I0][J1][0],X[I0][J1][1],X[I0][J1][2]);
 
 glEnd();
 
 }
 }
 }*/
void Gdisplay::drawTorus(std::vector< std::vector< std::vector <double> > > X, std::vector<std::vector<std::vector<double> > > C){
    
    //GLUquadricObj *quadric=gluNewQuadric();
    //gluQuadricNormals(quadric, GLU_SMOOTH);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); // fill
    //GLfloat colr[] = {col[0], col[1], col[2], 1.f};
    GLfloat wht[] = {1.f, 1.f, 1.f, 1.f};
    //glMaterialfv(GL_FRONT, GL_DIFFUSE, colr);
    glMaterialfv(GL_FRONT, GL_SPECULAR, wht);
    glMaterialf(GL_FRONT, GL_SHININESS, 60.);
    
    
    // obviously this will come out!
    int iN = X.size();
    int jN = X[0].size();
    
    for (int i=1;i<=iN;i++){
        int I0 = i%iN;
        int I1 = (i-1)%iN;
        for (int j=1;j<=jN;j++){
            int J0 = j%jN;
            int J1 = (j-1)%jN;
            //GLfloat colr[]={C[I0][J0][0],C[I0][J0][1],C[I0][J0][2],1.f};
            GLfloat colr[]={0.,1.,0.5,1.f};
            glMaterialfv(GL_FRONT, GL_DIFFUSE, colr);
            glBegin(GL_QUADS);
            
            double ax = X[I0][J1][0]-X[I1][J1][0];
            double ay = X[I0][J1][1]-X[I1][J1][1];
            double az = X[I0][J1][2]-X[I1][J1][2];
            double bx = X[I1][J0][0]-X[I1][J1][0];
            double by = X[I1][J0][1]-X[I1][J1][1];
            double bz = X[I1][J0][2]-X[I1][J1][2];
            double nx = ay*bz-az*by;
            double ny = az*bx-ax*bz;
            double nz = ax*by-ay*bx;
            double nl = 1./sqrt(nx*nx+ny*ny+nz*nz);
            nx*=nl;
            ny*=nl;
            nz*=nl;
            
            if((C[I0][J0][0]*nx+C[I0][J0][1]*ny+C[I0][J0][2]*nz)>0.){
                glNormal3d(nx,ny,nz);
            } else {
                glNormal3d(-nx,-ny,-nz);
            }
            
            //glNormal3d(C[I0][J0][0],C[I0][J0][1],C[I0][J0][2]);//nx*nl,ny*nl,nz*nl);
            
            glVertex3d(X[I1][J1][0],X[I1][J1][1],X[I1][J1][2]);
            glVertex3d(X[I1][J0][0],X[I1][J0][1],X[I1][J0][2]);
            glVertex3d(X[I0][J0][0],X[I0][J0][1],X[I0][J0][2]);
            glVertex3d(X[I0][J1][0],X[I0][J1][1],X[I0][J1][2]);
            
            glEnd();
            
        }
    }
}



/*
 void Gdisplay::drawCubeSphere(std::vector< std::vector <double> > X){
 
 glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); // fill
 //GLfloat colr[] = {col[0], col[1], col[2], 1.f};
 GLfloat wht[] = {1.f, 1.f, 1.f, 1.f};
 glMaterialfv(GL_FRONT, GL_SPECULAR, wht);
 glMaterialf(GL_FRONT, GL_SHININESS, 60.);
 
 int N = X.size();
 int n = sqrt(N/6);
 int a, b, c, d;
 double ax, ay, az, bx, by, bz, nx, ny, nz, nl;
 for(int i=0;i<6;i++){
 for(int s=1;s<n;s++){
 for(int t=1;t<n;t++){
 GLfloat colr[]={0.,1.,0.5,1.f};
 glMaterialfv(GL_FRONT, GL_DIFFUSE, colr);
 glBegin(GL_QUADS);
 
 a = i*n*n+(s-1)*n+(t-1);
 b = i*n*n+(s-1)*n+(t);
 c = i*n*n+(s)*n+(t);
 d = i*n*n+(s)*n+(t-1);
 
 ax = X[b][0]-X[a][0];
 ay = X[b][1]-X[a][1];
 az = X[b][2]-X[a][2];
 bx = X[d][0]-X[a][0];
 by = X[d][1]-X[a][1];
 bz = X[d][2]-X[a][2];
 nx = ay*bz-az*by;
 ny = az*bx-ax*bz;
 nz = ax*by-ay*bx;
 nl = -1./sqrt(nx*nx+ny*ny+nz*nz);
 
 glNormal3d(nx*nl,ny*nl,nz*nl);
 
 glVertex3d(X[a][0],X[a][1],X[a][2]);
 glVertex3d(X[b][0],X[b][1],X[b][2]);
 glVertex3d(X[c][0],X[c][1],X[c][2]);
 glVertex3d(X[d][0],X[d][1],X[d][2]);
 
 glEnd();
 }
 }
 }
 
 }
 */
/*
void Gdisplay::drawCubeSphere(std::vector< std::vector <double> > X){
    
    
    GLfloat colrZip[]={1.,1.,1.,1.f};
    GLfloat colrTri[]={1.,1.,1.,1.f};
    
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); // fill
    GLfloat wht[] = {1.f, 1.f, 1.f, 1.f};
    glMaterialfv(GL_FRONT, GL_SPECULAR, wht);
    glMaterialf(GL_FRONT, GL_SHININESS, 60.);
    
    int N = X.size();
    int n = sqrt(N/6);
    int k=(n-1);
    int nk=n*k;
    int n2 = n*n;
    int a, b, c, d;
    double ax, ay, az, bx, by, bz, nx, ny, nz, nl;
    for(int i=0;i<6;i++){
        for(int s=1;s<n;s++){
            for(int t=1;t<n;t++){
                GLfloat colr[]={0.,1.,0.5,1.f};
                glMaterialfv(GL_FRONT, GL_DIFFUSE, colr);
                glBegin(GL_QUADS);
                
                a = i*n2+(s-1)*n+(t-1);
                b = i*n2+(s-1)*n+(t);
                c = i*n2+(s)*n+(t);
                d = i*n2+(s)*n+(t-1);
                
                ax = X[b][0]-X[a][0];
                ay = X[b][1]-X[a][1];
                az = X[b][2]-X[a][2];
                bx = X[d][0]-X[a][0];
                by = X[d][1]-X[a][1];
                bz = X[d][2]-X[a][2];
                nx = ay*bz-az*by;
                ny = az*bx-ax*bz;
                nz = ax*by-ay*bx;
                nl = -1./sqrt(nx*nx+ny*ny+nz*nz);
                
                glNormal3d(nx*nl,ny*nl,nz*nl);
                
                glVertex3d(X[a][0],X[a][1],X[a][2]);
                glVertex3d(X[b][0],X[b][1],X[b][2]);
                glVertex3d(X[c][0],X[c][1],X[c][2]);
                glVertex3d(X[d][0],X[d][1],X[d][2]);
                
                glEnd();
            }
        }
    }
    
    //ZIPPING
    for(int i=1;i<n;i++){
        
        glMaterialfv(GL_FRONT, GL_DIFFUSE, colrZip);
        glBegin(GL_QUADS);
        
        //0A->2C
        a=0*n2+(i-1)*n;
        b=0*n2+i*n;
        c=2*n2+i;
        d=2*n2+(i-1);
        ax = X[b][0]-X[a][0];
        ay = X[b][1]-X[a][1];
        az = X[b][2]-X[a][2];
        bx = X[d][0]-X[a][0];
        by = X[d][1]-X[a][1];
        bz = X[d][2]-X[a][2];
        nx = ay*bz-az*by;
        ny = az*bx-ax*bz;
        nz = ax*by-ay*bx;
        nl = -1./sqrt(nx*nx+ny*ny+nz*nz);
        glNormal3d(nx*nl,ny*nl,nz*nl);
        glVertex3d(X[a][0],X[a][1],X[a][2]);
        glVertex3d(X[b][0],X[b][1],X[b][2]);
        glVertex3d(X[c][0],X[c][1],X[c][2]);
        glVertex3d(X[d][0],X[d][1],X[d][2]);
        
        //0B->3A
        a=0*n2+k+n*(i-1);
        b=0*n2+k+n*i;
        c=3*n2+i*n;
        d=3*n2+(i-1)*n;
        ax = X[b][0]-X[a][0];
        ay = X[b][1]-X[a][1];
        az = X[b][2]-X[a][2];
        bx = X[d][0]-X[a][0];
        by = X[d][1]-X[a][1];
        bz = X[d][2]-X[a][2];
        nx = ay*bz-az*by;
        ny = az*bx-ax*bz;
        nz = ax*by-ay*bx;
        nl = +1./sqrt(nx*nx+ny*ny+nz*nz);
        glNormal3d(nx*nl,ny*nl,nz*nl);
        glVertex3d(X[a][0],X[a][1],X[a][2]);
        glVertex3d(X[b][0],X[b][1],X[b][2]);
        glVertex3d(X[c][0],X[c][1],X[c][2]);
        glVertex3d(X[d][0],X[d][1],X[d][2]);
        
        //0C->4A
        a=0*n2+(i-1);
        b=0*n2+i;
        c=4*n2+i*n;
        d=4*n2+(i-1)*n;
        ax = X[b][0]-X[a][0];
        ay = X[b][1]-X[a][1];
        az = X[b][2]-X[a][2];
        bx = X[d][0]-X[a][0];
        by = X[d][1]-X[a][1];
        bz = X[d][2]-X[a][2];
        nx = ay*bz-az*by;
        ny = az*bx-ax*bz;
        nz = ax*by-ay*bx;
        nl = +1./sqrt(nx*nx+ny*ny+nz*nz);
        glNormal3d(nx*nl,ny*nl,nz*nl);
        glVertex3d(X[a][0],X[a][1],X[a][2]);
        glVertex3d(X[b][0],X[b][1],X[b][2]);
        glVertex3d(X[c][0],X[c][1],X[c][2]);
        glVertex3d(X[d][0],X[d][1],X[d][2]);
        
        //0D->5C
        a=0*n2+nk+(i-1);
        b=0*n2+nk+i;
        c=5*n2+i;
        d=5*n2+(i-1);
        ax = X[b][0]-X[a][0];
        ay = X[b][1]-X[a][1];
        az = X[b][2]-X[a][2];
        bx = X[d][0]-X[a][0];
        by = X[d][1]-X[a][1];
        bz = X[d][2]-X[a][2];
        nx = ay*bz-az*by;
        ny = az*bx-ax*bz;
        nz = ax*by-ay*bx;
        nl = -1./sqrt(nx*nx+ny*ny+nz*nz);
        glNormal3d(nx*nl,ny*nl,nz*nl);
        glVertex3d(X[a][0],X[a][1],X[a][2]);
        glVertex3d(X[b][0],X[b][1],X[b][2]);
        glVertex3d(X[c][0],X[c][1],X[c][2]);
        glVertex3d(X[d][0],X[d][1],X[d][2]);
        
        //1A->4B
        a=1*n2+(i-1)*n;
        b=1*n2+i*n;
        c=4*n2+k+n*i;
        d=4*n2+k+n*(i-1);
        ax = X[b][0]-X[a][0];
        ay = X[b][1]-X[a][1];
        az = X[b][2]-X[a][2];
        bx = X[d][0]-X[a][0];
        by = X[d][1]-X[a][1];
        bz = X[d][2]-X[a][2];
        nx = ay*bz-az*by;
        ny = az*bx-ax*bz;
        nz = ax*by-ay*bx;
        nl = -1./sqrt(nx*nx+ny*ny+nz*nz);
        glNormal3d(nx*nl,ny*nl,nz*nl);
        glVertex3d(X[a][0],X[a][1],X[a][2]);
        glVertex3d(X[b][0],X[b][1],X[b][2]);
        glVertex3d(X[c][0],X[c][1],X[c][2]);
        glVertex3d(X[d][0],X[d][1],X[d][2]);
        
        //1B->5D
        a=1*n2+k+n*(i-1);
        b=1*n2+k+n*i;
        c=5*n2+nk+i;
        d=5*n2+nk+(i-1);
        ax = X[b][0]-X[a][0];
        ay = X[b][1]-X[a][1];
        az = X[b][2]-X[a][2];
        bx = X[d][0]-X[a][0];
        by = X[d][1]-X[a][1];
        bz = X[d][2]-X[a][2];
        nx = ay*bz-az*by;
        ny = az*bx-ax*bz;
        nz = ax*by-ay*bx;
        nl = +1./sqrt(nx*nx+ny*ny+nz*nz);
        glNormal3d(nx*nl,ny*nl,nz*nl);
        glVertex3d(X[a][0],X[a][1],X[a][2]);
        glVertex3d(X[b][0],X[b][1],X[b][2]);
        glVertex3d(X[c][0],X[c][1],X[c][2]);
        glVertex3d(X[d][0],X[d][1],X[d][2]);
        
        //1C->2D
        a=1*n2+(i-1);
        b=1*n2+i;
        c=2*n2+nk+i;
        d=2*n2+nk+(i-1);
        ax = X[b][0]-X[a][0];
        ay = X[b][1]-X[a][1];
        az = X[b][2]-X[a][2];
        bx = X[d][0]-X[a][0];
        by = X[d][1]-X[a][1];
        bz = X[d][2]-X[a][2];
        nx = ay*bz-az*by;
        ny = az*bx-ax*bz;
        nz = ax*by-ay*bx;
        nl = +1./sqrt(nx*nx+ny*ny+nz*nz);
        glNormal3d(nx*nl,ny*nl,nz*nl);
        glVertex3d(X[a][0],X[a][1],X[a][2]);
        glVertex3d(X[b][0],X[b][1],X[b][2]);
        glVertex3d(X[c][0],X[c][1],X[c][2]);
        glVertex3d(X[d][0],X[d][1],X[d][2]);
        
        //1D->3B
        a=1*n2+nk+(i-1);
        b=1*n2+nk+i;
        c=3*n2+k+n*i;
        d=3*n2+k+n*(i-1);
        ax = X[b][0]-X[a][0];
        ay = X[b][1]-X[a][1];
        az = X[b][2]-X[a][2];
        bx = X[d][0]-X[a][0];
        by = X[d][1]-X[a][1];
        bz = X[d][2]-X[a][2];
        nx = ay*bz-az*by;
        ny = az*bx-ax*bz;
        nz = ax*by-ay*bx;
        nl = -1./sqrt(nx*nx+ny*ny+nz*nz);
        glNormal3d(nx*nl,ny*nl,nz*nl);
        glVertex3d(X[a][0],X[a][1],X[a][2]);
        glVertex3d(X[b][0],X[b][1],X[b][2]);
        glVertex3d(X[c][0],X[c][1],X[c][2]);
        glVertex3d(X[d][0],X[d][1],X[d][2]);
        
        //2A->4C
        a=2*n2+(i-1)*n;
        b=2*n2+i*n;
        c=4*n2+i;
        d=4*n2+(i-1);
        ax = X[b][0]-X[a][0];
        ay = X[b][1]-X[a][1];
        az = X[b][2]-X[a][2];
        bx = X[d][0]-X[a][0];
        by = X[d][1]-X[a][1];
        bz = X[d][2]-X[a][2];
        nx = ay*bz-az*by;
        ny = az*bx-ax*bz;
        nz = ax*by-ay*bx;
        nl = -1./sqrt(nx*nx+ny*ny+nz*nz);
        glNormal3d(nx*nl,ny*nl,nz*nl);
        glVertex3d(X[a][0],X[a][1],X[a][2]);
        glVertex3d(X[b][0],X[b][1],X[b][2]);
        glVertex3d(X[c][0],X[c][1],X[c][2]);
        glVertex3d(X[d][0],X[d][1],X[d][2]);
        
        //2B->5A
        a=2*n2+k+n*(i-1);
        b=2*n2+k+n*i;
        c=5*n2+i*n;
        d=5*n2+(i-1)*n;
        ax = X[b][0]-X[a][0];
        ay = X[b][1]-X[a][1];
        az = X[b][2]-X[a][2];
        bx = X[d][0]-X[a][0];
        by = X[d][1]-X[a][1];
        bz = X[d][2]-X[a][2];
        nx = ay*bz-az*by;
        ny = az*bx-ax*bz;
        nz = ax*by-ay*bx;
        nl = +1./sqrt(nx*nx+ny*ny+nz*nz);
        glNormal3d(nx*nl,ny*nl,nz*nl);
        glVertex3d(X[a][0],X[a][1],X[a][2]);
        glVertex3d(X[b][0],X[b][1],X[b][2]);
        glVertex3d(X[c][0],X[c][1],X[c][2]);
        glVertex3d(X[d][0],X[d][1],X[d][2]);
        
        //3C->4D
        a=3*n2+(i-1);
        b=3*n2+i;
        c=4*n2+nk+i;
        d=4*n2+nk+(i-1);
        ax = X[b][0]-X[a][0];
        ay = X[b][1]-X[a][1];
        az = X[b][2]-X[a][2];
        bx = X[d][0]-X[a][0];
        by = X[d][1]-X[a][1];
        bz = X[d][2]-X[a][2];
        nx = ay*bz-az*by;
        ny = az*bx-ax*bz;
        nz = ax*by-ay*bx;
        nl = +1./sqrt(nx*nx+ny*ny+nz*nz);
        glNormal3d(nx*nl,ny*nl,nz*nl);
        glVertex3d(X[a][0],X[a][1],X[a][2]);
        glVertex3d(X[b][0],X[b][1],X[b][2]);
        glVertex3d(X[c][0],X[c][1],X[c][2]);
        glVertex3d(X[d][0],X[d][1],X[d][2]);
        
        //3D->5B
        a=3*n2+nk+(i-1);
        b=3*n2+nk+i;
        c=5*n2+n*i+k;
        d=5*n2+n*(i-1)+k;
        ax = X[b][0]-X[a][0];
        ay = X[b][1]-X[a][1];
        az = X[b][2]-X[a][2];
        bx = X[d][0]-X[a][0];
        by = X[d][1]-X[a][1];
        bz = X[d][2]-X[a][2];
        nx = ay*bz-az*by;
        ny = az*bx-ax*bz;
        nz = ax*by-ay*bx;
        nl = -1./sqrt(nx*nx+ny*ny+nz*nz);
        glNormal3d(nx*nl,ny*nl,nz*nl);
        glVertex3d(X[a][0],X[a][1],X[a][2]);
        glVertex3d(X[b][0],X[b][1],X[b][2]);
        glVertex3d(X[c][0],X[c][1],X[c][2]);
        glVertex3d(X[d][0],X[d][1],X[d][2]);
        glEnd();
    }
    
    // CORNER TRIANGLES
    
    glMaterialfv(GL_FRONT, GL_DIFFUSE, colrTri);
    glBegin(GL_QUADS);
    
    a=0*n2+nk;
    b=2*n2+k;
    d=5*n2;
    ax = X[b][0]-X[a][0];
    ay = X[b][1]-X[a][1];
    az = X[b][2]-X[a][2];
    bx = X[d][0]-X[a][0];
    by = X[d][1]-X[a][1];
    bz = X[d][2]-X[a][2];
    nx = ay*bz-az*by;
    ny = az*bx-ax*bz;
    nz = ax*by-ay*bx;
    nl = 1./sqrt(nx*nx+ny*ny+nz*nz);
    glNormal3d(nx*nl,ny*nl,nz*nl);
    glVertex3d(X[a][0],X[a][1],X[a][2]);
    glVertex3d(X[b][0],X[b][1],X[b][2]);
    glVertex3d(X[b][0],X[b][1],X[b][2]);
    glVertex3d(X[d][0],X[d][1],X[d][2]);
    
    a=0*n2+nk+k;
    b=3*n2+nk;
    d=5*n2+k;
    ax = X[b][0]-X[a][0];
    ay = X[b][1]-X[a][1];
    az = X[b][2]-X[a][2];
    bx = X[d][0]-X[a][0];
    by = X[d][1]-X[a][1];
    bz = X[d][2]-X[a][2];
    nx = ay*bz-az*by;
    ny = az*bx-ax*bz;
    nz = ax*by-ay*bx;
    nl = -1./sqrt(nx*nx+ny*ny+nz*nz);
    glNormal3d(nx*nl,ny*nl,nz*nl);
    glVertex3d(X[a][0],X[a][1],X[a][2]);
    glVertex3d(X[b][0],X[b][1],X[b][2]);
    glVertex3d(X[b][0],X[b][1],X[b][2]);
    glVertex3d(X[d][0],X[d][1],X[d][2]);
    
    
    a=0*n2;
    b=2*n2;
    d=4*n2;
    ax = X[b][0]-X[a][0];
    ay = X[b][1]-X[a][1];
    az = X[b][2]-X[a][2];
    bx = X[d][0]-X[a][0];
    by = X[d][1]-X[a][1];
    bz = X[d][2]-X[a][2];
    nx = ay*bz-az*by;
    ny = az*bx-ax*bz;
    nz = ax*by-ay*bx;
    nl = -1./sqrt(nx*nx+ny*ny+nz*nz);
    glNormal3d(nx*nl,ny*nl,nz*nl);
    glVertex3d(X[a][0],X[a][1],X[a][2]);
    glVertex3d(X[b][0],X[b][1],X[b][2]);
    glVertex3d(X[b][0],X[b][1],X[b][2]);
    glVertex3d(X[d][0],X[d][1],X[d][2]);
    
    a=0*n2+k;
    b=3*n2;
    d=4*n2+nk;
    ax = X[b][0]-X[a][0];
    ay = X[b][1]-X[a][1];
    az = X[b][2]-X[a][2];
    bx = X[d][0]-X[a][0];
    by = X[d][1]-X[a][1];
    bz = X[d][2]-X[a][2];
    nx = ay*bz-az*by;
    ny = az*bx-ax*bz;
    nz = ax*by-ay*bx;
    nl = 1./sqrt(nx*nx+ny*ny+nz*nz);
    glNormal3d(nx*nl,ny*nl,nz*nl);
    glVertex3d(X[a][0],X[a][1],X[a][2]);
    glVertex3d(X[b][0],X[b][1],X[b][2]);
    glVertex3d(X[b][0],X[b][1],X[b][2]);
    glVertex3d(X[d][0],X[d][1],X[d][2]);
    
    
    // BOTTON TRIANGLES
    
    a=1*n2+k;
    b=2*n2+nk+k;
    d=5*n2+nk;
    c=b;
    ax = X[b][0]-X[a][0];
    ay = X[b][1]-X[a][1];
    az = X[b][2]-X[a][2];
    bx = X[d][0]-X[a][0];
    by = X[d][1]-X[a][1];
    bz = X[d][2]-X[a][2];
    nx = ay*bz-az*by;
    ny = az*bx-ax*bz;
    nz = ax*by-ay*bx;
    nl = 1./sqrt(nx*nx+ny*ny+nz*nz);
    glNormal3d(nx*nl,ny*nl,nz*nl);
    glVertex3d(X[a][0],X[a][1],X[a][2]);
    glVertex3d(X[b][0],X[b][1],X[b][2]);
    glVertex3d(X[c][0],X[c][1],X[c][2]);
    glVertex3d(X[d][0],X[d][1],X[d][2]);
    
    
    a=1*n2+nk+k;
    b=3*n2+nk+k;
    d=5*n2+nk+k;
    ax = X[b][0]-X[a][0];
    ay = X[b][1]-X[a][1];
    az = X[b][2]-X[a][2];
    bx = X[d][0]-X[a][0];
    by = X[d][1]-X[a][1];
    bz = X[d][2]-X[a][2];
    nx = ay*bz-az*by;
    ny = az*bx-ax*bz;
    nz = ax*by-ay*bx;
    nl = -1./sqrt(nx*nx+ny*ny+nz*nz);
    glNormal3d(nx*nl,ny*nl,nz*nl);
    glVertex3d(X[a][0],X[a][1],X[a][2]);
    glVertex3d(X[b][0],X[b][1],X[b][2]);
    glVertex3d(X[b][0],X[b][1],X[b][2]);
    glVertex3d(X[d][0],X[d][1],X[d][2]);
    
    
    a=1*n2;
    b=2*n2+nk;
    d=4*n2+k;
    ax = X[b][0]-X[a][0];
    ay = X[b][1]-X[a][1];
    az = X[b][2]-X[a][2];
    bx = X[d][0]-X[a][0];
    by = X[d][1]-X[a][1];
    bz = X[d][2]-X[a][2];
    nx = ay*bz-az*by;
    ny = az*bx-ax*bz;
    nz = ax*by-ay*bx;
    nl = -1./sqrt(nx*nx+ny*ny+nz*nz);
    glNormal3d(nx*nl,ny*nl,nz*nl);
    glVertex3d(X[a][0],X[a][1],X[a][2]);
    glVertex3d(X[b][0],X[b][1],X[b][2]);
    glVertex3d(X[b][0],X[b][1],X[b][2]);
    glVertex3d(X[d][0],X[d][1],X[d][2]);
    
    a=1*n2+nk;
    b=3*n2+k;
    d=4*n2+nk+k;
    
    ax = X[b][0]-X[a][0];
    ay = X[b][1]-X[a][1];
    az = X[b][2]-X[a][2];
    bx = X[d][0]-X[a][0];
    by = X[d][1]-X[a][1];
    bz = X[d][2]-X[a][2];
    nx = ay*bz-az*by;
    ny = az*bx-ax*bz;
    nz = ax*by-ay*bx;
    nl = 1./sqrt(nx*nx+ny*ny+nz*nz);
    glNormal3d(nx*nl,ny*nl,nz*nl);
    glVertex3d(X[a][0],X[a][1],X[a][2]);
    glVertex3d(X[b][0],X[b][1],X[b][2]);
    glVertex3d(X[b][0],X[b][1],X[b][2]);
    glVertex3d(X[d][0],X[d][1],X[d][2]);
    glEnd();
    
    
}
 */


void Gdisplay::drawCubeSphere(std::vector< std::vector <double> > X){
    
    
    GLfloat colrZip[]={1.,1.,1.,1.f};
    GLfloat colrTri[]={1.,1.,1.,1.f};
    
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); // fill
    GLfloat wht[] = {1.f, 1.f, 1.f, 1.f};
    glMaterialfv(GL_FRONT, GL_SPECULAR, wht);
    glMaterialf(GL_FRONT, GL_SHININESS, 60.);
    
    int N = X.size();
    int n = sqrt(N/6);
    int k=(n-1);
    int nk=n*k;
    int n2 = n*n;
    int a, b, c, d;
    double ax, ay, az, bx, by, bz, nx, ny, nz, nl;
    for(int i=0;i<6;i++){
        for(int s=1;s<n;s++){
            for(int t=1;t<n;t++){
                GLfloat colr[]={0.,1.,0.5,1.f};
                glMaterialfv(GL_FRONT, GL_DIFFUSE, colr);
                glBegin(GL_QUADS);
                
                a = i*n2+(s-1)*n+(t-1);
                b = i*n2+(s-1)*n+(t);
                c = i*n2+(s)*n+(t);
                d = i*n2+(s)*n+(t-1);
                
                ax = X[b][0]-X[a][0];
                ay = X[b][1]-X[a][1];
                az = X[b][2]-X[a][2];
                bx = X[d][0]-X[a][0];
                by = X[d][1]-X[a][1];
                bz = X[d][2]-X[a][2];
                nx = ay*bz-az*by;
                ny = az*bx-ax*bz;
                nz = ax*by-ay*bx;
                nl = -1./sqrt(nx*nx+ny*ny+nz*nz);
                
                glNormal3d(nx*nl,ny*nl,nz*nl);
                
                glVertex3d(X[a][0],X[a][1],X[a][2]);
                glVertex3d(X[b][0],X[b][1],X[b][2]);
                glVertex3d(X[c][0],X[c][1],X[c][2]);
                glVertex3d(X[d][0],X[d][1],X[d][2]);
                
                glEnd();
            }
        }
    }
    
    //ZIPPING
    for(int i=1;i<n;i++){
        
        glMaterialfv(GL_FRONT, GL_DIFFUSE, colrZip);
        glBegin(GL_QUADS);
        
        //0A->2C
        a=0*n2+(i-1)*n;
        b=0*n2+i*n;
        c=2*n2+i;
        d=2*n2+(i-1);
        
        ax = X[b][0]-X[a][0];
        ay = X[b][1]-X[a][1];
        az = X[b][2]-X[a][2];
        bx = X[d][0]-X[a][0];
        by = X[d][1]-X[a][1];
        bz = X[d][2]-X[a][2];
        nx = ay*bz-az*by;
        ny = az*bx-ax*bz;
        nz = ax*by-ay*bx;
        nl = -1./sqrt(nx*nx+ny*ny+nz*nz);
        glNormal3d(nx*nl,ny*nl,nz*nl);
        glVertex3d(X[a][0],X[a][1],X[a][2]);
        glVertex3d(X[b][0],X[b][1],X[b][2]);
        glVertex3d(X[c][0],X[c][1],X[c][2]);
        glVertex3d(X[d][0],X[d][1],X[d][2]);
        
        //0B->3A
        a=0*n2+k+n*(i-1);
        b=3*n2+(i-1)*n;
        c=3*n2+i*n;
        d=0*n2+k+n*i;

        ax = X[b][0]-X[a][0];
        ay = X[b][1]-X[a][1];
        az = X[b][2]-X[a][2];
        bx = X[d][0]-X[a][0];
        by = X[d][1]-X[a][1];
        bz = X[d][2]-X[a][2];
        nx = ay*bz-az*by;
        ny = az*bx-ax*bz;
        nz = ax*by-ay*bx;
        nl = -1./sqrt(nx*nx+ny*ny+nz*nz); //CHANGED TO -
        glNormal3d(nx*nl,ny*nl,nz*nl);
        glVertex3d(X[a][0],X[a][1],X[a][2]);
        glVertex3d(X[b][0],X[b][1],X[b][2]);
        glVertex3d(X[c][0],X[c][1],X[c][2]);
        glVertex3d(X[d][0],X[d][1],X[d][2]);
        
        //0C->4A
        a=0*n2+(i-1);
        b=4*n2+(i-1)*n;
        c=4*n2+i*n;
        d=0*n2+i;
        
        ax = X[b][0]-X[a][0];
        ay = X[b][1]-X[a][1];
        az = X[b][2]-X[a][2];
        bx = X[d][0]-X[a][0];
        by = X[d][1]-X[a][1];
        bz = X[d][2]-X[a][2];
        nx = ay*bz-az*by;
        ny = az*bx-ax*bz;
        nz = ax*by-ay*bx;
        nl = -1./sqrt(nx*nx+ny*ny+nz*nz); // CHANGED TO -
        glNormal3d(nx*nl,ny*nl,nz*nl);
        glVertex3d(X[a][0],X[a][1],X[a][2]);
        glVertex3d(X[b][0],X[b][1],X[b][2]);
        glVertex3d(X[c][0],X[c][1],X[c][2]);
        glVertex3d(X[d][0],X[d][1],X[d][2]);
        
        //0D->5C
        a=0*n2+nk+(i-1);
        b=0*n2+nk+i;
        c=5*n2+i;
        d=5*n2+(i-1);
        
        ax = X[b][0]-X[a][0];
        ay = X[b][1]-X[a][1];
        az = X[b][2]-X[a][2];
        bx = X[d][0]-X[a][0];
        by = X[d][1]-X[a][1];
        bz = X[d][2]-X[a][2];
        nx = ay*bz-az*by;
        ny = az*bx-ax*bz;
        nz = ax*by-ay*bx;
        nl = -1./sqrt(nx*nx+ny*ny+nz*nz);
        glNormal3d(nx*nl,ny*nl,nz*nl);
        glVertex3d(X[a][0],X[a][1],X[a][2]);
        glVertex3d(X[b][0],X[b][1],X[b][2]);
        glVertex3d(X[c][0],X[c][1],X[c][2]);
        glVertex3d(X[d][0],X[d][1],X[d][2]);
        
        //1A->4B
        a=1*n2+(i-1)*n;
        b=1*n2+i*n;
        c=4*n2+k+n*i;
        d=4*n2+k+n*(i-1);
        
        ax = X[b][0]-X[a][0];
        ay = X[b][1]-X[a][1];
        az = X[b][2]-X[a][2];
        bx = X[d][0]-X[a][0];
        by = X[d][1]-X[a][1];
        bz = X[d][2]-X[a][2];
        nx = ay*bz-az*by;
        ny = az*bx-ax*bz;
        nz = ax*by-ay*bx;
        nl = -1./sqrt(nx*nx+ny*ny+nz*nz);
        glNormal3d(nx*nl,ny*nl,nz*nl);
        glVertex3d(X[a][0],X[a][1],X[a][2]);
        glVertex3d(X[b][0],X[b][1],X[b][2]);
        glVertex3d(X[c][0],X[c][1],X[c][2]);
        glVertex3d(X[d][0],X[d][1],X[d][2]);
        
        //1B->5D
        a=1*n2+k+n*(i-1);
        b=5*n2+nk+(i-1);
        c=5*n2+nk+i;
        d=1*n2+k+n*i;

        
        ax = X[b][0]-X[a][0];
        ay = X[b][1]-X[a][1];
        az = X[b][2]-X[a][2];
        bx = X[d][0]-X[a][0];
        by = X[d][1]-X[a][1];
        bz = X[d][2]-X[a][2];
        nx = ay*bz-az*by;
        ny = az*bx-ax*bz;
        nz = ax*by-ay*bx;
        nl = -1./sqrt(nx*nx+ny*ny+nz*nz); // CHANGED TO -
        glNormal3d(nx*nl,ny*nl,nz*nl);
        glVertex3d(X[a][0],X[a][1],X[a][2]);
        glVertex3d(X[b][0],X[b][1],X[b][2]);
        glVertex3d(X[c][0],X[c][1],X[c][2]);
        glVertex3d(X[d][0],X[d][1],X[d][2]);
        
        //1C->2D
        a=1*n2+(i-1);
        b=2*n2+nk+(i-1);
        c=2*n2+nk+i;
        d=1*n2+i;
        
        ax = X[b][0]-X[a][0];
        ay = X[b][1]-X[a][1];
        az = X[b][2]-X[a][2];
        bx = X[d][0]-X[a][0];
        by = X[d][1]-X[a][1];
        bz = X[d][2]-X[a][2];
        nx = ay*bz-az*by;
        ny = az*bx-ax*bz;
        nz = ax*by-ay*bx;
        nl = -1./sqrt(nx*nx+ny*ny+nz*nz); // CHANGED TO -
        glNormal3d(nx*nl,ny*nl,nz*nl);
        glVertex3d(X[a][0],X[a][1],X[a][2]);
        glVertex3d(X[b][0],X[b][1],X[b][2]);
        glVertex3d(X[c][0],X[c][1],X[c][2]);
        glVertex3d(X[d][0],X[d][1],X[d][2]);
        
        //1D->3B
        a=1*n2+nk+(i-1);
        b=1*n2+nk+i;
        c=3*n2+k+n*i;
        d=3*n2+k+n*(i-1);
        ax = X[b][0]-X[a][0];
        ay = X[b][1]-X[a][1];
        az = X[b][2]-X[a][2];
        bx = X[d][0]-X[a][0];
        by = X[d][1]-X[a][1];
        bz = X[d][2]-X[a][2];
        nx = ay*bz-az*by;
        ny = az*bx-ax*bz;
        nz = ax*by-ay*bx;
        nl = -1./sqrt(nx*nx+ny*ny+nz*nz);
        glNormal3d(nx*nl,ny*nl,nz*nl);
        glVertex3d(X[a][0],X[a][1],X[a][2]);
        glVertex3d(X[b][0],X[b][1],X[b][2]);
        glVertex3d(X[c][0],X[c][1],X[c][2]);
        glVertex3d(X[d][0],X[d][1],X[d][2]);
        
        //2A->4C
        a=2*n2+(i-1)*n;
        b=2*n2+i*n;
        c=4*n2+i;
        d=4*n2+(i-1);
        ax = X[b][0]-X[a][0];
        ay = X[b][1]-X[a][1];
        az = X[b][2]-X[a][2];
        bx = X[d][0]-X[a][0];
        by = X[d][1]-X[a][1];
        bz = X[d][2]-X[a][2];
        nx = ay*bz-az*by;
        ny = az*bx-ax*bz;
        nz = ax*by-ay*bx;
        nl = -1./sqrt(nx*nx+ny*ny+nz*nz);
        glNormal3d(nx*nl,ny*nl,nz*nl);
        glVertex3d(X[a][0],X[a][1],X[a][2]);
        glVertex3d(X[b][0],X[b][1],X[b][2]);
        glVertex3d(X[c][0],X[c][1],X[c][2]);
        glVertex3d(X[d][0],X[d][1],X[d][2]);
        
        //2B->5A
        a=2*n2+k+n*(i-1);
        b=5*n2+(i-1)*n;
        c=5*n2+i*n;
        d=2*n2+k+n*i;

        ax = X[b][0]-X[a][0];
        ay = X[b][1]-X[a][1];
        az = X[b][2]-X[a][2];
        bx = X[d][0]-X[a][0];
        by = X[d][1]-X[a][1];
        bz = X[d][2]-X[a][2];
        nx = ay*bz-az*by;
        ny = az*bx-ax*bz;
        nz = ax*by-ay*bx;
        nl = -1./sqrt(nx*nx+ny*ny+nz*nz);    // CHANGED TO -
        glNormal3d(nx*nl,ny*nl,nz*nl);
        glVertex3d(X[a][0],X[a][1],X[a][2]);
        glVertex3d(X[b][0],X[b][1],X[b][2]);
        glVertex3d(X[c][0],X[c][1],X[c][2]);
        glVertex3d(X[d][0],X[d][1],X[d][2]);
        
        //3C->4D
        a=3*n2+(i-1);
        b=4*n2+nk+(i-1);
        c=4*n2+nk+i;
        d=3*n2+i;
        
        ax = X[b][0]-X[a][0];
        ay = X[b][1]-X[a][1];
        az = X[b][2]-X[a][2];
        bx = X[d][0]-X[a][0];
        by = X[d][1]-X[a][1];
        bz = X[d][2]-X[a][2];
        nx = ay*bz-az*by;
        ny = az*bx-ax*bz;
        nz = ax*by-ay*bx;
        nl = -1./sqrt(nx*nx+ny*ny+nz*nz);   // CHANGED TO -
        glNormal3d(nx*nl,ny*nl,nz*nl);
        glVertex3d(X[a][0],X[a][1],X[a][2]);
        glVertex3d(X[b][0],X[b][1],X[b][2]);
        glVertex3d(X[c][0],X[c][1],X[c][2]);
        glVertex3d(X[d][0],X[d][1],X[d][2]);
        
        //3D->5B
        a=3*n2+nk+(i-1);
        b=3*n2+nk+i;
        c=5*n2+n*i+k;
        d=5*n2+n*(i-1)+k;
        ax = X[b][0]-X[a][0];
        ay = X[b][1]-X[a][1];
        az = X[b][2]-X[a][2];
        bx = X[d][0]-X[a][0];
        by = X[d][1]-X[a][1];
        bz = X[d][2]-X[a][2];
        nx = ay*bz-az*by;
        ny = az*bx-ax*bz;
        nz = ax*by-ay*bx;
        nl = -1./sqrt(nx*nx+ny*ny+nz*nz);
        glNormal3d(nx*nl,ny*nl,nz*nl);
        glVertex3d(X[a][0],X[a][1],X[a][2]);
        glVertex3d(X[b][0],X[b][1],X[b][2]);
        glVertex3d(X[c][0],X[c][1],X[c][2]);
        glVertex3d(X[d][0],X[d][1],X[d][2]);
        glEnd();
    }
    
    // CORNER TRIANGLES
    
    glMaterialfv(GL_FRONT, GL_DIFFUSE, colrTri);
    glBegin(GL_QUADS);
    
    a=0*n2+nk;
    b=5*n2;
    d=2*n2+k;
    ax = X[b][0]-X[a][0];
    ay = X[b][1]-X[a][1];
    az = X[b][2]-X[a][2];
    bx = X[d][0]-X[a][0];
    by = X[d][1]-X[a][1];
    bz = X[d][2]-X[a][2];
    nx = ay*bz-az*by;
    ny = az*bx-ax*bz;
    nz = ax*by-ay*bx;
    nl = -1./sqrt(nx*nx+ny*ny+nz*nz); // CHANGED TO -
    glNormal3d(nx*nl,ny*nl,nz*nl);
    glVertex3d(X[a][0],X[a][1],X[a][2]);
    glVertex3d(X[b][0],X[b][1],X[b][2]);
    glVertex3d(X[b][0],X[b][1],X[b][2]);
    glVertex3d(X[d][0],X[d][1],X[d][2]);
    
    a=0*n2+nk+k;
    b=3*n2+nk;
    d=5*n2+k;
    ax = X[b][0]-X[a][0];
    ay = X[b][1]-X[a][1];
    az = X[b][2]-X[a][2];
    bx = X[d][0]-X[a][0];
    by = X[d][1]-X[a][1];
    bz = X[d][2]-X[a][2];
    nx = ay*bz-az*by;
    ny = az*bx-ax*bz;
    nz = ax*by-ay*bx;
    nl = -1./sqrt(nx*nx+ny*ny+nz*nz);
    glNormal3d(nx*nl,ny*nl,nz*nl);
    glVertex3d(X[a][0],X[a][1],X[a][2]);
    glVertex3d(X[b][0],X[b][1],X[b][2]);
    glVertex3d(X[b][0],X[b][1],X[b][2]);
    glVertex3d(X[d][0],X[d][1],X[d][2]);
    
    
    a=0*n2;
    b=2*n2;
    d=4*n2;
    ax = X[b][0]-X[a][0];
    ay = X[b][1]-X[a][1];
    az = X[b][2]-X[a][2];
    bx = X[d][0]-X[a][0];
    by = X[d][1]-X[a][1];
    bz = X[d][2]-X[a][2];
    nx = ay*bz-az*by;
    ny = az*bx-ax*bz;
    nz = ax*by-ay*bx;
    nl = -1./sqrt(nx*nx+ny*ny+nz*nz);
    glNormal3d(nx*nl,ny*nl,nz*nl);
    glVertex3d(X[a][0],X[a][1],X[a][2]);
    glVertex3d(X[b][0],X[b][1],X[b][2]);
    glVertex3d(X[b][0],X[b][1],X[b][2]);
    glVertex3d(X[d][0],X[d][1],X[d][2]);
    
    a=0*n2+k;
    b=4*n2+nk;
    d=3*n2;
    ax = X[b][0]-X[a][0];
    ay = X[b][1]-X[a][1];
    az = X[b][2]-X[a][2];
    bx = X[d][0]-X[a][0];
    by = X[d][1]-X[a][1];
    bz = X[d][2]-X[a][2];
    nx = ay*bz-az*by;
    ny = az*bx-ax*bz;
    nz = ax*by-ay*bx;
    nl = -1./sqrt(nx*nx+ny*ny+nz*nz);   // CHANGED TO -
    glNormal3d(nx*nl,ny*nl,nz*nl);
    glVertex3d(X[a][0],X[a][1],X[a][2]);
    glVertex3d(X[b][0],X[b][1],X[b][2]);
    glVertex3d(X[b][0],X[b][1],X[b][2]);
    glVertex3d(X[d][0],X[d][1],X[d][2]);
    
    
    // BOTTON TRIANGLES
    
    a=1*n2+k;
    b=5*n2+nk;
    d=2*n2+nk+k;
    c=b;
    ax = X[b][0]-X[a][0];
    ay = X[b][1]-X[a][1];
    az = X[b][2]-X[a][2];
    bx = X[d][0]-X[a][0];
    by = X[d][1]-X[a][1];
    bz = X[d][2]-X[a][2];
    nx = ay*bz-az*by;
    ny = az*bx-ax*bz;
    nz = ax*by-ay*bx;
    nl = -1./sqrt(nx*nx+ny*ny+nz*nz);  // CHANGED TO -
    glNormal3d(nx*nl,ny*nl,nz*nl);
    glVertex3d(X[a][0],X[a][1],X[a][2]);
    glVertex3d(X[b][0],X[b][1],X[b][2]);
    glVertex3d(X[c][0],X[c][1],X[c][2]);
    glVertex3d(X[d][0],X[d][1],X[d][2]);
    
    
    a=1*n2+nk+k;
    b=3*n2+nk+k;
    d=5*n2+nk+k;
    ax = X[b][0]-X[a][0];
    ay = X[b][1]-X[a][1];
    az = X[b][2]-X[a][2];
    bx = X[d][0]-X[a][0];
    by = X[d][1]-X[a][1];
    bz = X[d][2]-X[a][2];
    nx = ay*bz-az*by;
    ny = az*bx-ax*bz;
    nz = ax*by-ay*bx;
    nl = -1./sqrt(nx*nx+ny*ny+nz*nz);
    glNormal3d(nx*nl,ny*nl,nz*nl);
    glVertex3d(X[a][0],X[a][1],X[a][2]);
    glVertex3d(X[b][0],X[b][1],X[b][2]);
    glVertex3d(X[b][0],X[b][1],X[b][2]);
    glVertex3d(X[d][0],X[d][1],X[d][2]);
    
    
    a=1*n2;
    b=2*n2+nk;
    d=4*n2+k;
    ax = X[b][0]-X[a][0];
    ay = X[b][1]-X[a][1];
    az = X[b][2]-X[a][2];
    bx = X[d][0]-X[a][0];
    by = X[d][1]-X[a][1];
    bz = X[d][2]-X[a][2];
    nx = ay*bz-az*by;
    ny = az*bx-ax*bz;
    nz = ax*by-ay*bx;
    nl = -1./sqrt(nx*nx+ny*ny+nz*nz);
    glNormal3d(nx*nl,ny*nl,nz*nl);
    glVertex3d(X[a][0],X[a][1],X[a][2]);
    glVertex3d(X[b][0],X[b][1],X[b][2]);
    glVertex3d(X[b][0],X[b][1],X[b][2]);
    glVertex3d(X[d][0],X[d][1],X[d][2]);
    
    a=1*n2+nk;
    b=4*n2+nk+k;
    d=3*n2+k;
    
    ax = X[b][0]-X[a][0];
    ay = X[b][1]-X[a][1];
    az = X[b][2]-X[a][2];
    bx = X[d][0]-X[a][0];
    by = X[d][1]-X[a][1];
    bz = X[d][2]-X[a][2];
    nx = ay*bz-az*by;
    ny = az*bx-ax*bz;
    nz = ax*by-ay*bx;
    nl = -1./sqrt(nx*nx+ny*ny+nz*nz);   // CHANGED TO -
    glNormal3d(nx*nl,ny*nl,nz*nl);
    glVertex3d(X[a][0],X[a][1],X[a][2]);
    glVertex3d(X[b][0],X[b][1],X[b][2]);
    glVertex3d(X[b][0],X[b][1],X[b][2]);
    glVertex3d(X[d][0],X[d][1],X[d][2]);
    glEnd();
    
    
}


void Gdisplay::drawSphereFromMesh(std::vector<std::vector<double> > X, std::vector<std::vector<int> > M, std::vector<std::vector<double> > C){
    
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); // fill
    GLfloat wht[] = {1.f, 1.f, 1.f, 1.f};
    glMaterialfv(GL_FRONT, GL_SPECULAR, wht);
    glMaterialf(GL_FRONT, GL_SHININESS, 60.);
    
    int a, b, c, d;
    double ax, ay, az, bx, by, bz, nx, ny, nz, nl;
    
    for(int i=0;i<X.size();i++){
        GLfloat col[]={C[i][0],C[i][1],C[i][2],1.f};
        glMaterialfv(GL_FRONT, GL_DIFFUSE, col);
        glBegin(GL_QUADS);
        a = i;
        b = M[i][0];        //+x
        c = M[M[i][0]][1];  //+xy
        d = M[i][1];        //+y
        ax = X[b][0]-X[a][0];
        ay = X[b][1]-X[a][1];
        az = X[b][2]-X[a][2];
        bx = X[d][0]-X[a][0];
        by = X[d][1]-X[a][1];
        bz = X[d][2]-X[a][2];
        nx = ay*bz-az*by;
        ny = az*bx-ax*bz;
        nz = ax*by-ay*bx;
        nl = 1./sqrt(nx*nx+ny*ny+nz*nz);
        glNormal3d(nx*nl,ny*nl,nz*nl);
        glVertex3d(X[a][0],X[a][1],X[a][2]);
        glVertex3d(X[b][0],X[b][1],X[b][2]);
        glVertex3d(X[c][0],X[c][1],X[c][2]);
        glVertex3d(X[d][0],X[d][1],X[d][2]);
        glEnd();
    }
    for(int i=0;i<X.size();i++){
        GLfloat col[]={C[i][0],C[i][1],C[i][2],1.f};
        glMaterialfv(GL_FRONT, GL_DIFFUSE, col);
        glBegin(GL_QUADS);
        a = i;
        b = M[i][2];        //+x
        c = M[M[i][2]][3];  //+xy
        d = M[i][3];        //+y
        ax = X[b][0]-X[a][0];
        ay = X[b][1]-X[a][1];
        az = X[b][2]-X[a][2];
        bx = X[d][0]-X[a][0];
        by = X[d][1]-X[a][1];
        bz = X[d][2]-X[a][2];
        nx = ay*bz-az*by;
        ny = az*bx-ax*bz;
        nz = ax*by-ay*bx;
        nl = 1./sqrt(nx*nx+ny*ny+nz*nz);
        glNormal3d(nx*nl,ny*nl,nz*nl);
        glVertex3d(X[a][0],X[a][1],X[a][2]);
        glVertex3d(X[b][0],X[b][1],X[b][2]);
        glVertex3d(X[c][0],X[c][1],X[c][2]);
        glVertex3d(X[d][0],X[d][1],X[d][2]);
        glEnd();
    }

    
}

void Gdisplay::drawFlatCube(std::vector< std::vector<int> > C, std::vector< std::vector <double> > Col,double X, double Y, double Z){
    
    int n = sqrt(C.size()/6);
    double dn1 = 1./(double)n;
    double dn2 = 0.5*dn1;
    double scale = 0.5;
    double x, y, xp, yp;
    std::vector<double> xoff(6);
    std::vector<double> yoff(6);
    xoff[0]=  0; yoff[0]= 1.5;
    xoff[1]=  0; yoff[1]= 0.5;
    xoff[2]=  0; yoff[2]=-0.5;
    xoff[3]=  0; yoff[3]= -1.5;
    xoff[4]= -1; yoff[4]= 0.5;
    xoff[5]= +1; yoff[5]= 0.5;
    
    for(int i=0;i<C.size();i++){
        int f = C[i][0];
            x = xoff[f]+((((double)C[i][1]+0.5)*dn1)-0.5)-dn2;
            xp = x+dn1;
            x*=scale;
            xp*=scale;
            //for(int t=0;t<n;t++){
                y = yoff[f]+((((double)C[i][2]+0.5)*dn1)-0.5)-dn2;
                yp = y+dn1;
                y*=scale;
                yp*=scale;
                
                GLfloat colr[]={Col[i][0],Col[i][1],Col[i][2],1.f};
                glMaterialfv(GL_FRONT, GL_DIFFUSE, colr);
                glBegin(GL_QUADS);
                glNormal3d(0,0,1);
                glVertex3d(X+x,Y+y,Z+0);
                glVertex3d(X+x,Y+yp,Z+0);
                glVertex3d(X+xp,Y+yp,Z+0);
                glVertex3d(X+xp,Y+y,Z+0);
                glEnd();
 //               k++;
//            }
 //       }
    }
    
}

/*
void Gdisplay::drawFlatCube(std::vector< std::vector <double> > C){
    
    int n = sqrt(C.size()/6);
    double dn1 = 1./(double)n;
    double dn2 = 0.5*dn1;
    double scale = 0.5;
    double x, y, xp, yp;
    std::vector<double> xoff(6);
    std::vector<double> yoff(6);
    xoff[0]=  0; yoff[0]=-1.5;
    xoff[1]=  0; yoff[1]= 0.5;
    xoff[2]=  0; yoff[2]=-0.5;
    xoff[3]=  0; yoff[3]= 1.5;
    xoff[4]= -1; yoff[4]= 0.5;
    xoff[5]= +1; yoff[5]= 0.5;
    
    int k=0;
    for(int i=0;i<6;i++){
        for(int s=0;s<n;s++){
            x = xoff[i]+((((double)s+0.5)*dn1)-0.5)-dn2;
            xp = x+dn1;
            x*=scale;
            xp*=scale;
            for(int t=0;t<n;t++){
                y = yoff[i]+((((double)t+0.5)*dn1)-0.5)-dn2;
                yp = y+dn1;
                y*=scale;
                yp*=scale;
                
                GLfloat colr[]={C[k][0],C[k][1],C[k][2],1.f};
                glMaterialfv(GL_FRONT, GL_DIFFUSE, colr);
                glBegin(GL_QUADS);
                glNormal3d(0,0,1);
                glVertex3d(x,y,0);
                glVertex3d(x,yp,0);
                glVertex3d(xp,yp,0);
                glVertex3d(xp,y,0);
                glEnd();
                k++;
            }
        }
    }
    
}
 */



void Gdisplay::addQuad(std::vector <double> p1, std::vector <double> p2, std::vector <double> p3, std::vector <double> p4, std::vector <double> N, std::vector <double> C){
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); // fill
    glBegin(GL_POLYGON);
    GLfloat col[] = {C[0], C[1], C[2], 1.f};
    GLfloat wht[] = {1.f, 1.f, 1.f, 1.f};
    glMaterialfv(GL_FRONT, GL_DIFFUSE, col);
    //  glMaterialfv(GL_FRONT, GL_SPECULAR, wht); // uncomment for shiny!
    glMaterialf(GL_FRONT, GL_SHININESS, 60.);
    glNormal3d(N[0],N[1],N[2]);
    glVertex3d(p1[0],p1[1],p1[2]);
    glVertex3d(p2[0],p2[1],p2[2]);
    glVertex3d(p3[0],p3[1],p3[2]);
    glVertex3d(p4[0],p4[1],p4[2]);
    glEnd();
}

void Gdisplay::addCrossHairs(double d, double l, int w){
    
    //x col
    drawLine(-d,-d,-d,	-d+l,-d,-d,		1.,0.,0.,w); //-z
    drawLine(-d,+d,-d,	-d+l,+d,-d,		1.,0.,0.,w);
    drawLine(+d,+d,-d,	+d-l,+d,-d,		1.,0.,0.,w);
    drawLine(+d,-d,-d,	+d-l,-d,-d,		1.,0.,0.,w);
    drawLine(-d,-d,0.,	-d+l,-d,0.,		1.,0.,0.,w); //0z
    drawLine(-d,+d,0.,	-d+l,+d,0.,		1.,0.,0.,w);
    drawLine(+d,+d,0.,	+d-l,+d,0.,		1.,0.,0.,w);
    drawLine(+d,-d,0.,	+d-l,-d,0.,		1.,0.,0.,w);
    drawLine(-d,-d,+d,	-d+l,-d,+d,		1.,0.,0.,w); //+z
    drawLine(-d,+d,+d,	-d+l,+d,+d,		1.,0.,0.,w);
    drawLine(+d,+d,+d,	+d-l,+d,+d,		1.,0.,0.,w);
    drawLine(+d,-d,+d,	+d-l,-d,+d,		1.,0.,0.,w);
    
    //y col
    drawLine(-d,-d,-d,	-d,-d+l,-d,		0.,1.,0.,w); //-z
    drawLine(-d,+d,-d,	-d,+d-l,-d,		0.,1.,0.,w);
    drawLine(+d,+d,-d,	+d,+d-l,-d,		0.,1.,0.,w);
    drawLine(+d,-d,-d,	+d,-d+l,-d,		0.,1.,0.,w);
    drawLine(-d,-d,0.,	-d,-d+l,0.,		0.,1.,0.,w); //0z
    drawLine(-d,+d,0.,	-d,+d-l,0.,		0.,1.,0.,w);
    drawLine(+d,+d,0.,	+d,+d-l,0.,		0.,1.,0.,w);
    drawLine(+d,-d,0.,	+d,-d+l,0.,		0.,1.,0.,w);
    drawLine(-d,-d,+d,	-d,-d+l,+d,		0.,1.,0.,w); //+z
    drawLine(-d,+d,+d,	-d,+d-l,+d,		0.,1.,0.,w);
    drawLine(+d,+d,+d,	+d,+d-l,+d,		0.,1.,0.,w);
    drawLine(+d,-d,+d,	+d,-d+l,+d,		0.,1.,0.,w);
    
    //z col
    drawLine(-d,-d,-d,	-d,-d,-d+l,		0.,0.,1.,w); //-z
    drawLine(-d,+d,-d,	-d,+d,-d+l,		0.,0.,1.,w);
    drawLine(+d,+d,-d,	+d,+d,-d+l,		0.,0.,1.,w);
    drawLine(+d,-d,-d,	+d,-d,-d+l,		0.,0.,1.,w);
    drawLine(-d,-d,0.,	-d,-d,0.+l,		0.,0.,1.,w); //-z
    drawLine(-d,+d,0.,	-d,+d,0.+l,		0.,0.,1.,w);
    drawLine(+d,+d,0.,	+d,+d,0.+l,		0.,0.,1.,w);
    drawLine(+d,-d,0.,	+d,-d,0.+l,		0.,0.,1.,w);
    drawLine(-d,-d,+d,	-d,-d,+d-l,		0.,0.,1.,w); //+z
    drawLine(-d,+d,+d,	-d,+d,+d-l,		0.,0.,1.,w);
    drawLine(+d,+d,+d,	+d,+d,+d-l,		0.,0.,1.,w);
    drawLine(+d,-d,+d,	+d,-d,+d-l,		0.,0.,1.,w);
}

void Gdisplay::saveImage(std::string filename){
    glXMakeCurrent(disp, win, glc);
    GLubyte * bits; //RGB bits
    GLint viewport[4]; //current viewport
    glGetIntegerv(GL_VIEWPORT, viewport);
    int w = viewport[2];
    int h = viewport[3];
    bits = new GLubyte[w*3*h];
    glFinish(); //finish all commands of OpenGL
    glPixelStorei(GL_PACK_ALIGNMENT,1); //or glPixelStorei(GL_PACK_ALIGNMENT,4);
    glPixelStorei(GL_PACK_ROW_LENGTH, 0);
    glPixelStorei(GL_PACK_SKIP_ROWS, 0);
    glPixelStorei(GL_PACK_SKIP_PIXELS, 0);
    glReadPixels(0, 0, w, h, GL_BGR_EXT, GL_UNSIGNED_BYTE, bits);
    IplImage * capImg = cvCreateImage(cvSize(w,h), IPL_DEPTH_8U, 3);
    for(int i=0; i<h; ++i){
        for(int j=0; j < w; ++j){
            capImg->imageData[i*capImg->widthStep + j*3+0] = (unsigned char)(bits[(h-i-1)*3*w + j*3+0]);
            capImg->imageData[i*capImg->widthStep + j*3+1] = (unsigned char)(bits[(h-i-1)*3*w + j*3+1]);
            capImg->imageData[i*capImg->widthStep + j*3+2] = (unsigned char)(bits[(h-i-1)*3*w + j*3+2]);
        }
    }
    std::stringstream outFile;
    outFile << filename;
    cvSaveImage(outFile.str().c_str(),capImg);
    cvReleaseImage(&capImg);
    delete[] bits;
}




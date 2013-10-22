#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <iostream>

//include CMINPACK for Levenberg-Marquardt optimization
#include <cminpack.h>

using namespace std;
//using namespace arma;

class Hydrooptics {
    /* Class contais definitions of all hydrooptial equations for 
     * calculation of Remote sensing reflectance and cost function
     * and Jacobian of cost function
     */ 
    public:

    static const double PI =  3.14159265;

    // Rrsw_w = – 0.00036 + 0.110(bb/a) – 0.0447(bb/a)2
    static const double RW0 =-0.00036;
    static const double RW1 = 0.11;
    static const double RW2 =-0.0447;

    // Kd  = (1/mu0)[a2 + ab(0.473mu0 – 0.218)]1/2 ==
    static const double KD0 =-0.218;
    static const double KD1 = 0.473;

    // b = bb / B_ == */
    static const double B_WAT = 0.5;
    static const double B_CHL = 0.0011;
    static const double B_TSM = 0.08;
    static const double B_DOC = 1;

    // water refraction index
    static const double NW = 1.33333;
    
    int bands;
    
    double * aaw;
    double * bbw;
    double * bw;
    double * aam;
    double * bbm;
    double * bm;

    double * s;
    double * al;

    double h;
    double theta;
    double qf;
    double mu0;
    
    //Constructor:
    //Set number of bands and the HO-model
    Hydrooptics(int inBands, double * inModel);

    //Destructor
    //Clean memory from model
    ~Hydrooptics();
    
    //set hydro-optical model in the object
    int set_model(double * model);
    
    //set hydo-optical conditions: reflectance, solar zenith
    int set_params(double inTheta);

    //set hydo-optical conditions: reflectance, solar zenith
    int set_params(double * inS, double inTheta);

    //caluclate Subsurface remote sensing refectance (Rrsw) from given C
    //for deep waters
    //Other parameters (model, solar zenith) are defined
    //in the object (at initialization)
    //Returns vector for the entrie spectrum
    double rrsw (const double * c, int bn) const;

    //Calcualte cost function: difference between measured and
    //reconstructed R for given C
    //Returns vector for the entrie spectrum
    double rs (const double * c, int bn) const;

    //Calcualte sum square error: sum pf squared cost function for all bands
    //reconstructed R for given C
    double sse (const double * c) const;

    //Calculate Jacobian of the cost function:
    //partial derivative of the cost function by each concentration
    //Returns vector for the entrie spectrum
    double jacobian(const double * c, int bn, int vn) const;

    //Calculate Jacobian of the cost function for deep waters
    //Returns single value at one wavelength
    double j_deep(double s,
                double c0, double c1, double c2,
                double a0, double a1, double a2,
                double bb0, double bb1, double bb2,
                double b0, double b1, double b2,
                double aWAT, double bbWAT, double bWAT) const;
};

int startingCPA(double * parameters, double * startC);

//Interface for the LM-optimization library
int fcn(void *p, int m, int n, const double *x, double *fvec, double *fjac, 
	 int ldfjac, int iflag);


#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <iostream>

//include Armadillo: C++ linear algebra library 
//http://arma.sourceforge.net/
#include <armadillo>

//include CMINPACK for Levenberg-Marquardt optimization
//http://devernay.free.fr/hacks/cminpack/index.html
#include <cminpack.h>
#define real __cminpack_real__

using namespace std;
using namespace arma;

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
    
    //parameters for albedo approximation
    static const double G1 = 555;   //peak center
    static const double G2 = 40;    //peak width
    static const double L0 = -0.06; //line 0

    mat m;
    mat s;
    mat al;
    mat lambda;
    double h;
    double theta;
    double qf;
    double mu0;
    
    //Constructor:
    //Set model, albedo, depth, solar zenith
    Hydrooptics(mat inM, mat inAL, double h, double theta);

    //Set       model,   depth,    solar zenith  wavelengths, 
    Hydrooptics(mat inM, double h, double theta, mat inLambda);

    //create matrix [1 c0 c1 c2]
    mat cOne(mat c) const;

    //albedo approximation
    mat albedo(double al1, double al2) const;

    //caluclate Subsurface remote sensing refectance (Rrsw) from given C
    //for deep waters
    //Other parameters (model, albedo, depth, sola zenith) are defined
    //in the object (at initialization)
    //Returns vector for the entrie spectrum
    mat rrsw (mat c) const;

    //caluclate Subsurface remote sensing refectance (Rrsw) from given C
    //and given ALR
    //for deep waters
    //Other parameters (model, albedo, depth, sola zenith) are defined
    //in the object (at initialization)
    //Returns vector for the entrie spectrum
    mat rrsw (mat c, double al1, double al2) const;

    //Calcualte cost function: difference between measured and
    //reconstructed R for given C
    //Returns vector for the entrie spectrum
    mat rs (mat c) const;    

    //Calcualte cost function: difference between measured and
    //reconstructed R for given C and ALR
    //Returns vector for the entrie spectrum
    mat rs (mat c, double al1, double  al2) const;    

    //Calculate Jacobian of the cost function:
    //partial derivative of the cost function by each concentration
    //Returns vector for the entrie spectrum
    mat jacobian(mat c) const;

    //Calculate Jacobian of the cost function:
    //partial derivative of the cost function by each concentration and albedo
    //Returns vector for the entrie spectrum
    mat jacobian(mat c, double al1, double al2) const;
    
    //Calculate Jacobian of the cost function for deep waters
    //Returns single value at one wavelength
    double j_deep(double s,
                double c0, double c1, double c2,
                double a0, double a1, double a2,
                double bb0, double bb1, double bb2,
                double b0, double b1, double b2,
                double aWAT, double bbWAT, double bWAT) const;
    //Calculate Jacobian of the cost function for shallow waters
    //Returns single value at one wavelength
    double j_shallow(double s,
                double c0, double c1, double c2,
                double a0, double a1, double a2,
                double bb0, double bb1, double bb2,
                double b0, double b1, double b2,
                double aWAT, double bbWAT, double bWAT,
                double al) const;

    //Calculate Jacobian of the cost function with additional variable: albedo
    //Returns single value of derivative on C at one wavelength
    double j_shallow_c0(double s,
                double c0, double c1, double c2,
                double a0, double a1, double a2,
                double bb0, double bb1, double bb2,
                double b0, double b1, double b2,
                double aWAT, double bbWAT, double bWAT,
                double al1, double al2, double ll) const;

    //Calculate Jacobian of the cost function with additional variable: albedo
    //Returns single value of derivative on ALR (albedo ratio) at one wavelength
    double j_shallow_al1(double s,
                double c0, double c1, double c2,
                double a0, double a1, double a2,
                double bb0, double bb1, double bb2,
                double b0, double b1, double b2,
                double aWAT, double bbWAT, double bWAT,
                double al1, double al2, double ll) const;

    //Calculate Jacobian of the cost function with additional variable: albedo
    //Returns single value of derivative on ALR (albedo ratio) at one wavelength
    double j_shallow_al2(double s,
                double c0, double c1, double c2,
                double a0, double a1, double a2,
                double bb0, double bb1, double bb2,
                double b0, double b1, double b2,
                double aWAT, double bbWAT, double bWAT,
                double al1, double al2, double ll) const;

};

//Interface for the LM-optimization library
int fcn(void *p, int m, int n, const real *x, real *fvec, real *fjac, 
	 int ldfjac, int iflag);


//Interface for the LM-optimization library
int fcn_al(void *p, int m, int n, const real *x, real *fvec, real *fjac, 
	 int ldfjac, int iflag);

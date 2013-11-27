#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <iostream>

//include CMINPACK for Levenberg-Marquardt optimization
#include <cminpack.h>

using namespace std;

class Hydrooptics {
    /* Class contais definitions of all hydrooptial equations for 
     * calculation of Remote sensing reflectance and cost function
     * and Jacobian of cost function
     */ 
    public:

    static const double PI =  3.14159265;

    // Rrsw_w = – 0.00036 + 0.110(bb/a) – 0.0447(bb/a)2
    static const double RW0 = -0.00036;
    static const double RW1 =  0.11;
    static const double RW2 = -0.0447;

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
    
    //set hydo-optical conditions: solar zenith
    int set_theta(double inTheta);

    //set hydo-optical conditions: reflectance, solar zenith
    int set_s(double * inS);
    
    //caluclate Subsurface remote sensing refectance (Rrsw) from given C
    //for deep waters
    //Other parameters (model, albedo, depth, sola zenith) are defined
    //in the object (at initialization)
    //Returns vector for the entrie spectrum
    virtual double rrsw (const double * c, int bn) const;

    //Calcualte cost function: difference between measured and
    //reconstructed R for given C
    //Returns vector for the entrie spectrum
    virtual double rs (const double * c, int bn) const;

    //Calcualte sum square error: sum pf squared cost function for all bands
    //reconstructed R for given C
    double sse (const double * c) const;

    //Calculate Jacobian of the cost function:
    //partial derivative of the cost function by each concentration
    //Returns vector for the entrie spectrum
    double jacobian(const double * c, int bn, int vn) const;

    //Calculate Jacobian of the cost function for deep waters
    //Returns single value at one wavelength
    virtual double j(int bn, double s,
                     double c0, double c1, double c2,
                     double a0, double a1, double a2,
                     double bb0, double bb1, double bb2,
                     double b0, double b1, double b2,
                     double aWAT, double bbWAT, double bWAT) const;
    
    int retrieve(int startCN, double * startC, double * xBest);
};


class HydroopticsShallow : public Hydrooptics {
    /* Class contais definitions of all hydrooptial equations for 
     * calculation of Remote sensing reflectance and cost function
     * and Jacobian of cost function
     */ 
    public:

    // Kd  = (1/mu0)[a2 + ab(0.473mu0 – 0.218)]1/2 ==
    static const double KD0 = -0.218;
    static const double KD1 =  0.473;
    
    double * al;
    double h;

    //Constructor:
    //Set number of bands and the HO-model
    HydroopticsShallow(int inBands, double * inModel);

    //Destructor
    //Clean memory from model
    ~HydroopticsShallow();

    //set hydo-optical conditions: albedo
    int set_albedo(double * inAL);

    //caluclate Subsurface remote sensing refectance (Rrsw) from given C
    //for deep waters
    //Other parameters (model, albedo, depth, sola zenith) are defined
    //in the object (at initialization)
    //Returns vector for the entrie spectrum
    virtual double rrsw (const double * c, int bn) const;


    //Calculate Jacobian of the cost function for shallow waters
    //Returns single value at one wavelength
    virtual double j(int bn, double s,
                     double c0, double c1, double c2,
                     double a0, double a1, double a2,
                     double bb0, double bb1, double bb2,
                     double b0, double b1, double b2,
                     double aWAT, double bbWAT, double bWAT,
                     double al) const;
    
};

class HydroopticsAlbedo : public HydroopticsShallow {
    /* Class contais definitions of all hydrooptial equations for 
     * calculation of Remote sensing reflectance and cost function
     * and Jacobian of cost function
     */ 
    public:

    double * lambda;

    double h;

    //parameters for albedo approximation
    static const double G1 = 555;   //peak center
    static const double G2 = 40;    //peak width
    static const double L0 = -0.06; //line 0
       
    //Constructor:
    //Set number of bands and the HO-model
    HydroopticsAlbedo(int inBands, double * inModel, double * inLambda);

    //Destructor
    //Clean memory from model
    ~HydroopticsAlbedo();

    //albedo approximation
    double albedo(double al1, double al2, int bn) const;

    //caluclate Subsurface remote sensing refectance (Rrsw) from given C
    //for deep waters
    //Other parameters (model, albedo, depth, sola zenith) are defined
    //in the object (at initialization)
    //Returns vector for the entrie spectrum
    virtual double rrsw (const double * c, int bn) const;
    
    //Calcualte cost function: difference between measured and
    //reconstructed R for given C
    //Returns vector for the entrie spectrum
    virtual double rs (const double * c, int bn) const;

    //Calculate Jacobian of the cost function with additional variable: albedo
    //Returns single value of derivative on C at one wavelength
    virtual double j(int bn, double s,
                    double c0, double c1, double c2,
                    double a0, double a1, double a2,
                    double bb0, double bb1, double bb2,
                    double b0, double b1, double b2,
                    double aWAT, double bbWAT, double bWAT,
                    double al1, double al2, double ll) const;

};

int startingCPA(double * parameters, double * startC);

int startingCPA_al(double * parameters, double * startC);

//Interface for the LM-optimization library
int fcn(void *p, int m, int n, const double *x, double *fvec, double *fjac, 
	 int ldfjac, int iflag);

int compare (const void * v1, const void * v2);
//Interface for the LM-optimization library
//int fcn_al(void *p, int m, int n, const double *x, double *fvec, double *fjac, 
//	 int ldfjac, int iflag);

%module lm
%{
/* Put header files here or function declarations like below */

extern int get_rrsw(
         double *model, int model_n0, int model_n1,
         double inC[3],
         double *albedo, int albedo_n0,
         double h,
         double theta,
         double *outR, int outR_n0);

extern int get_rrsw_al(
         double *model, int model_n0, int model_n1,
         double inC[3],
         double *lambda, int lambda_n0,
         double h,
         double theta,
         double al1,
         double al2,
         double *outR, int outR_n0);

extern int get_c(double parameters[6],
         double *model, int model_n0, int model_n1,
         double *inR, int inR_rows, int inR_cols,
         double *albedo, int albedo_rows, int albedo_cols,
         double *h, int h_rows,
         double *theta, int theta_rows,
         double *outC, int outC_length);

extern int get_c_al(double parameters[6],
         double *model, int model_n0, int model_n1,
         double *inR, int inR_rows, int inR_cols,
         double *lambda, int lambda_n0,
         double *h, int h_rows,
         double *theta, int theta_rows,
         double *outC, int outC_length);
         
#define SWIG_FILE_WITH_INIT

%}
 
%include "numpy.i"

%init %{
    import_array();
%}

%apply (double IN_ARRAY1[ANY]) {
    (double parameters[6]),
    (double inC[3])
};


%apply (double* IN_ARRAY2, int DIM1, int DIM2) {
    (double *model, int model_n0, int model_n1),
    (double *inR, int inR_rows, int inR_cols),
    (double *albedo, int albedo_rows, int albedo_cols)
};

%apply (double* IN_ARRAY1, int DIM1) {
    (double *albedo, int albedo_n0),
    (double *lambda, int lambda_n0),
    (double *h, int h_rows),
    (double *theta, int theta_rows)
};

%apply (double* ARGOUT_ARRAY1, int DIM1) {
    (double *outR, int outR_n0),
    (double *outC, int outC_length)
};

extern int get_rrsw(
         double *model, int model_n0, int model_n1,
         double inC[3],
         double *albedo, int albedo_n0,
         double h,
         double theta,
         double *outR, int outR_n0);

extern int get_rrsw_al(
         double *model, int model_n0, int model_n1,
         double inC[3],
         double *lambda, int lambda_n0,
         double h,
         double theta,
         double al1,
         double al2,
         double *outR, int outR_n0);

extern int get_c(double parameters[6],
         double *model, int model_n0, int model_n1,
         double *inR, int inR_rows, int inR_cols,
         double *albedo, int albedo_rows, int albedo_cols,
         double *h, int h_rows,
         double *theta, int theta_rows,
         double *outC, int outC_length);

extern int get_c_al(double parameters[6],
         double *model, int model_n0, int model_n1,
         double *inR, int inR_rows, int inR_cols,
         double *lambda, int lambda_n0,
         double *h, int h_rows,
         double *theta, int theta_rows,
         double *outC, int outC_length);


#include "fem_functions.h"

#include <iostream>
#include <cstdlib>
#include <cmath>
using namespace std;

static int nqpoints = 1;

/* Quadrature points
   - (x1,y1,x2,y2,...) */
static double qpoints[3] = {
  0.25,
  0.25,
  0.25};


/* Quadrature weights
   - (v1,v2,...) */
static double qweights[1] = {0.166666666666666666667};


/* Nodal basis function evaluations
    - basis component is fastest varying, then basis function, then quad point */
static double Tet4Basis[4] = {
  0.25,    /* N1 */
  0.25,    /* N2 */
  0.25,    /* N3 */
  0.25};   /* N4 */

/* Nodal basis function derivative evaluations,
    - basis component is fastest varying, then derivative direction, then basis function, then quad point */
static double TetBasisDerivatives[4][3] = {
    {-1.0 ,-1.0,-1.0 },
    { 1.0 , 0.0, 0.0 },  /* dN1/dxi   dN1/deta  dN1/dzeta */
    { 0.0 , 1.0, 0.0 },  /* dN2/dxi   dN2/deta  dN2/dzeta */
    { 0.0 , 0.0, 1.0 }   /* dN3/dxi   dN3/deta  dN3/dzeta */
                         /* dN4/dxi   dN4/deta  dN4/dzeta */
  };

#define X(i)      xyz[i*3]
#define Y(i)      xyz[i*3+1]
#define Z(i)      xyz[i*3+2]
#define DPHI(i,j) dphi[i*3+j]

void compute_tet_functions(double *xyz, double *phi, double *dphi, double *JxW)
{
    
    double J[3][3], invJ[3][3];
    double detJ, invDet;
    
     
    J[0][0] = X(1)-X(0);
    J[0][1] = X(2)-X(0);
    J[0][2] = X(3)-X(0);
        
    J[1][0] = Y(1)-Y(0);
    J[1][1] = Y(2)-Y(0);
    J[1][2] = Y(3)-Y(0);
    
    J[2][0] = Z(1)-Z(0);
    J[2][1] = Z(2)-Z(0);
    J[2][2] = Z(3)-Z(0);
    
    detJ =(J[0][0]*( J[1][1]*J[2][2] - J[1][2]*J[2][1]      ) +
           J[0][1]*( J[1][2]*J[2][0] - J[1][0]*J[2][2]      ) +
           J[0][2]*( J[1][0]*J[2][1] - J[1][1]*J[2][0])     ) ;
  
    
    if(detJ < 0.0) {
        cout << "Negative Jacobian!" << endl;
        exit(1);
    }
        
    invDet = 1.0/(detJ);

    invJ[0][0] = invDet*(J[1][1]*J[2][2] - J[1][2]*J[2][1]);
    invJ[0][1] = invDet*(J[0][2]*J[2][1] - J[0][1]*J[2][2]);
    invJ[0][2] = invDet*(J[0][1]*J[1][2] - J[0][2]*J[1][1]);
    invJ[1][0] = invDet*(J[1][2]*J[2][0] - J[1][0]*J[2][2]);
    invJ[1][1] = invDet*(J[0][0]*J[2][2] - J[0][2]*J[2][0]);
    invJ[1][2] = invDet*(J[0][2]*J[1][0] - J[0][0]*J[1][2]);
    invJ[2][0] = invDet*(J[1][0]*J[2][1] - J[1][1]*J[2][0]);
    invJ[2][1] = invDet*(J[0][1]*J[2][0] - J[0][0]*J[2][1]);
    invJ[2][2] = invDet*(J[0][0]*J[1][1] - J[0][1]*J[1][0]);
    
    
    for(int i=0; i < 4; i++)
    {
        for(int j=0; j < 3; j++) {
            DPHI(i,j) = 0.0;
            for(int k=0; k < 3; k++)
                   DPHI(i,j) += invJ[k][j]*TetBasisDerivatives[i][k];
        }
    }
        
        *JxW = detJ*qweights[0];
        phi = &Tet4Basis[0];
}


#if 0
void compute_tet_shape_function(double *qpoints, double *phi, double *dphi, double *w)
{
    /* Quadrature points
   - (x1,y1,x2,y2,...) */

    qpoints[0]= -0.5;
    qpoints[1]= -0.5;
    qpoints[2]= -0.5;

   /* Quadrature weights
   - (v1,v2,...) */
    *w = 1.33333333333;


/* Nodal basis function evaluations
    - basis component is fastest varying, then basis function, then quad point */
    phi[0] = 0.25;
    phi[1] = 0.25;
    phi[2] = 0.25;
    phi[3] = 0.25;

/* Nodal basis function derivative evaluations,
    - basis component is fastest varying, then derivative direction, then basis function, then quad point */
    dphi[0]  = -0.5;
    dphi[1]  = -0.5;
    dphi[2]  = -0.5;
    dphi[3]  =  0.5;
    dphi[4]  =  0.0;
    dphi[5]  =  0.0;
    dphi[6]  =  0.0;
    dphi[7]  =  0.5;
    dphi[8]  =  0.0;
    dphi[9]  =  0.0;
    dphi[10] =  0.0;
    dphi[11] =  0.0;
    dphi[12] =  0.5;

}


void compute_tet_functions(double* xyz, double *dphi, double *VOL)
{
    
    double x1 = xyz[0];
    double y1 = xyz[1];
    double z1 = xyz[2];
    
    double x2 = xyz[3];
    double y2 = xyz[4];
    double z2 = xyz[5];
    
    
    double x3 = xyz[6];
    double y3 = xyz[7];
    double z3 = xyz[8];
    
    double x4 = xyz[9];
    double y4 = xyz[10];
    double z4 = xyz[11];
    
    double x21 = x2 - x1;
    double x31 = x3 - x1;
    double x32 = x3 - x2;
    double x34 = x3 - x4;
    double x41 = x4 - x1;
    double x42 = x4 - x2;

    double y14 = y1 - y4;
    double y12 = y1 - y2;
    double y23 = y2 - y3;
    double y31 = y3 - y1;
    double y34 = y3 - y4;
    double y42 = y4 - y2;

    double z14 = z1 - z4;
    double z12 = z1 - z2;
    double z23 = z2 - z3;
    double z31 = z3 - z1;
    double z34 = z3 - z4;
    double z42 = z4 - z2;

    double x23 = - x32;
    double x43 = - x34;
    double x12 = - x21;

    double y43 = - y34;
    double y41 = - y14;
    double y21 = - y12;

    double z43 = - z34;
    double z32 = - z23;
    double z41 = - z14;
    double z21 = - z12;
    
    double V6 = x21*(z31*y14-y31*z14) + x31*(z12*y14-y12*z14) + x41*(z12*y31-y12*z31);

    double V6INV = 1.0 / V6;

    *VOL = V6/6.0; 

    dphi[0] = ( z43*y23 - z23*y43 ) * V6INV ;
    dphi[1] = ( z42*x32 - z32*x42 ) * V6INV ;
    dphi[2] = ( y43*x23 - x43*y23 ) * V6INV ;

    dphi[3] = ( y31*z41 - y41*z31 ) * V6INV ;
    dphi[4] = ( x41*z31 - x31*z41 ) * V6INV ;
    dphi[5] = ( y41*x31 - y31*x41 ) * V6INV ;

    dphi[6] = ( y41*z21 - y21*z41 ) * V6INV ;
    dphi[7] = ( x21*z41 - x41*z21 ) * V6INV ;
    dphi[8] = ( y42*x12 - x42*y12 ) * V6INV ;

    dphi[9]  = - ( dphi[0] + dphi[3] + dphi[6] );
    dphi[10] = - ( dphi[1] + dphi[4] + dphi[7] );
    dphi[11] = - ( dphi[2] + dphi[5] + dphi[8] );
    
}




#endif


double compute_tet_volume(double *xyz)
{
    double x1 = xyz[0];
    double y1 = xyz[1];
    double z1 = xyz[2];
    
    double x2 = xyz[3];
    double y2 = xyz[4];
    double z2 = xyz[5];
    
    
    double x3 = xyz[6];
    double y3 = xyz[7];
    double z3 = xyz[8];
    
    double x4 = xyz[9];
    double y4 = xyz[10];
    double z4 = xyz[11];
    
    double x21 = x2 - x1;
    double x31 = x3 - x1;
    double x32 = x3 - x2;
    double x34 = x3 - x4;
    double x41 = x4 - x1;
  

    double y14 = y1 - y4;
    double y12 = y1 - y2;
   
    double y31 = y3 - y1;
    double y34 = y3 - y4;
    double z14 = z1 - z4;
    double z12 = z1 - z2;
    double z31 = z3 - z1;

    return (x21*(z31*y14-y31*z14) + x31*(z12*y14-y12*z14) + x41*(z12*y31-y12*z31))/6.0;
    
}

double compute_element_quality(double *xyz)
{
    double x1 = xyz[0];
    double y1 = xyz[1];
    double z1 = xyz[2];
    
    double x2 = xyz[3];
    double y2 = xyz[4];
    double z2 = xyz[5];
    
    
    double x3 = xyz[6];
    double y3 = xyz[7];
    double z3 = xyz[8];
    
    double x4 = xyz[9];
    double y4 = xyz[10];
    double z4 = xyz[11];
    
    double x21 = x2 - x1;
    double x31 = x3 - x1;
    double x32 = x3 - x2;
    double x34 = x3 - x4;
    double x41 = x4 - x1;
    double x42 = x4 - x2;
  

    double y14 = y1 - y4;
    double y12 = y1 - y2;
    double y23 = y2 - y3;
    double y24 = y2 - y4;   
    double y31 = y3 - y1;
    double y34 = y3 - y4;
    
    double z14 = z1 - z4;
    double z12 = z1 - z2;
    double z31 = z3 - z1;
    double z23 = z2 - z3;
    double z34 = z3 - z4;
    double z24 = z2 - z4;

    double vol = (x21*(z31*y14-y31*z14) + x31*(z12*y14-y12*z14) + x41*(z12*y31-y12*z31))/6.0;
    
    double e12 = sqrt(x21*x21+y12*y12+z12*z12);
    double e23 = sqrt(x32*x32+y23*y23+z23*z23);
    double e31 = sqrt(x31*x31+y31*y31+z31*z31);
    double e14 = sqrt(x41*x41+y14*y14+z14*z14);
    double e24 = sqrt(x42*x42+y24*y24+z24*z24);
    double e34 = sqrt(x34*x34+y34*y34+z34*z34);
    
    double et    = e12*e12+e23*e23+e31*e31+e14*e14+e24*e24+e34*e34;
    double cubic = 1.7320508075;
    
    return (72.0*cubic*vol)/(pow(et,1.5));
    
    
    
}
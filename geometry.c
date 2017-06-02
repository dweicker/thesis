//
//  geometry.c
//  
//
//  Created by David Weicker on 7/04/17.
//
//

#include "geometry.h"


/* Quadrilateral mapping for [0;1] */
void quad_mapping_01(double *X,double *Y,double *xsi, double *eta, int n, double *x, double *y){
    int i;
    for(i=0;i<n;i++){
        x[i] = X[0]*(1-xsi[i])*(1-eta[i]) + X[1]*xsi[i]*(1-eta[i]) + X[2]*(1-xsi[i])*eta[i] + X[3]*xsi[i]*eta[i];
        y[i] = Y[0]*(1-xsi[i])*(1-eta[i]) + Y[1]*xsi[i]*(1-eta[i]) + Y[2]*(1-xsi[i])*eta[i] + Y[3]*xsi[i]*eta[i];
    }
}

/* Quadrilateral mapping for [-1;1] */
void quad_mapping_11(double *X,double *Y,double *xsi, double *eta, int n, double *x, double *y){
    int i;
    for(i=0;i<n;i++){
        x[i] = X[0]*(1-xsi[i])*(1-eta[i])/4 + X[1]*(1+xsi[i])*(1-eta[i])/4 + X[2]*(1-xsi[i])*(1+eta[i])/4 + X[3]*(1+xsi[i])*(1+eta[i])/4;
        y[i] = Y[0]*(1-xsi[i])*(1-eta[i])/4 + Y[1]*(1+xsi[i])*(1-eta[i])/4 + Y[2]*(1-xsi[i])*(1+eta[i])/4 + Y[3]*(1+xsi[i])*(1+eta[i])/4;
    }
}

/** Compute the jacobian (from mapping 11) evaluated at a given point
 *
 * \param [in] X,Y          Physical coordinates of the corners of the quadrant
 * \param [in] xsi,eta      Point where we want to evaluate
 * \return                  The value of the jacobian
 */
double jacobian(double *X, double *Y, double xsi, double eta){
    double dxdxsi = (X[1]-X[0])*(1-eta)/4 + (X[3]-X[2])*(1+eta)/4;
    double dxdeta = (X[2]-X[0])*(1-xsi)/4 + (X[3]-X[1])*(1+xsi)/4;
    double dydxsi = (Y[1]-Y[0])*(1-eta)/4 + (Y[3]-Y[2])*(1+eta)/4;
    double dydeta = (Y[2]-Y[0])*(1-xsi)/4 + (Y[3]-Y[1])*(1+xsi)/4;
    return dxdxsi*dydeta - dxdeta*dydxsi;
}


/** Build the derivation matrix for the Lagrange polynomials
 *
 * \param [in] xsi      The 1D GLL points
 * \param [out] H       The derivation matrix
 * \param [in] degree   The degree of the interpolation
 */
void derivation_matrix(double *xsi, double *H, int degree){
    double product = 1.0;
    double sum = 0.0;
    int i,k,j,m;
    for(i=0;i<=degree;i++){
        for(k=0;k<=degree;k++){
            for(j=0;j<=degree;j++){
                if(j!=i){
                    for(m=0;m<=degree;m++){
                        if(m!=j && m!=i)
                            product *= (xsi[k]-xsi[m])/(xsi[i]-xsi[m]);
                    }
                    sum += product/(xsi[i]-xsi[j]);
                    product = 1.0;
                }
            }
            H[i*(degree+1)+k] = sum;
            sum = 0.0;
        }
    }
}

/** Compute the derivatives of the inverse mapping
 *
 * \param [in] X,Y          The corners of the actual element
 * \param [in] xsi,eta      Where we want to evaluate in the reference element
 * \param [out] factor      The geometric factors we want to have (must be allocated to length 3)
 */
void geom_factor(double *X, double *Y, double xsi, double eta, double *factor){
    double dxdxsi = (X[1]-X[0])*(1-eta)/4 + (X[3]-X[2])*(1+eta)/4;
    double dxdeta = (X[2]-X[0])*(1-xsi)/4 + (X[3]-X[1])*(1+xsi)/4;
    double dydxsi = (Y[1]-Y[0])*(1-eta)/4 + (Y[3]-Y[2])*(1+eta)/4;
    double dydeta = (Y[2]-Y[0])*(1-xsi)/4 + (Y[3]-Y[1])*(1+xsi)/4;
    double jac = dxdxsi*dydeta - dxdeta*dydxsi;
    factor[0] = (dydeta*dydeta+dxdeta*dxdeta)/jac;
    factor[1] = (-dydeta*dydxsi-dxdeta*dxdxsi)/jac;
    factor[2] = (dydxsi*dydxsi+dxdxsi*dxdxsi)/jac;
}

/** Compute the generalized weights for the sem scheme
 *
 * \param [in] X,Y                  The corners of the actual element
 * \param [in] gll_points           Gauss-Lobatto-Legendre points
 * \param [in] weights              GLL 1D integration weights
 * \param [in] degree               The degree of the interpolation
 * \param [out] Wee,Wen,Wnn     The different generalized weights
 */
void gen_weights(double *X, double *Y, double *gll_points, double *weights, int degree, double *Wee, double *Wen, double *Wnn){
    int p,q;
    int N = degree+1;
    double factor[3];
    for(q=0;q<=degree;q++){
        for(p=0;p<=degree;p++){
            geom_factor(X,Y,gll_points[p],gll_points[q],factor);
            Wee[q*N+p] = weights[p]*weights[q]*factor[0];
            Wen[q*N+p] = weights[p]*weights[q]*factor[1];
            Wnn[q*N+p] = weights[p]*weights[q]*factor[2];
        }
    }
}

/** Value of the 1D i-th Lagrangian interpolant evaluated at xsi
 *
 * \param[in] gll_points    Gauss-Lobatto-Legendre points
 * \param[in] degree        The degree of the interpolation
 * \param[in] i             The index of the shape function
 * \param[in] xsi           The point where we want to evaluate phi
 * \return                  The value of the 1D i-th Lagrangian interpolant at xsi
 */
double phi(double *gll_points, int degree, int i, double xsi){
    double value = 1.0;
    for(int j=0;j<=degree;j++){
        if(j!=i){
            value *= (xsi-gll_points[j])/(gll_points[i]-gll_points[j]);
        }
    }
    return value;
}

/** Compute the general projection for hanging nodes
 *
 * \param [in] xsi          The gll points in the reference element [-1;1]
 * \param [in] degree       The degree of the interpolation
 * \param [out] proj        The reference matrix for the projection (size : 2*(degree+1) X (degree+1) )
 */
void general_projection(double *xsi, int degree, double *proj){
    int N = degree+1;
    int i,j;
    //we fill the top half of the matrix
    for(j=0;j<N;j++){
        for(i=0;i<N;i++){
            proj[j*N+i] = phi(xsi,degree,i,xsi[j]/2-0.5);
        }
    }
    //we fill the bottom half of the matrix
    for(j=0;j<N;j++){
        for(i=0;i<N;i++){
            proj[(j+N)*N+i] = phi(xsi,degree,i,xsi[j]/2+0.5);
        }
    }
}

/** Compute the transformation matrix for element where there are hanging nodes
 *  The transformation matrix links the nodes in the element with the actual nodes (= identity if no hanging node)
 *
 * \param [in] hanging_corner       The matrix from quad_decode
 * \param [in] degree               The degree of the interpolation
 * \param [in] vnodes               The total number of nodes in an element
 * \param [in] gen_proj             The general projection matrix (1D projection)
 * \param [out] transform           The transformation matrix (size : vnodes X vnodes) - must be calloc !
 * \param [out] interior            Small array defining the interior = [line_beg line_end col_beg col_end]
 */
void transformation_matrix(int *hanging_corner, int degree, int vnodes, double *gen_proj, double *transform, int *interior){
    int i,k;
    int N = degree+1;
    int part;
    int line_beg = 0;
    int col_beg = 0;
    int line_end = N;
    int col_end = N;
    if(hanging_corner[0] == 1 || hanging_corner[1] == 0){
        //the bottom face is hanging
        part = (hanging_corner[0] == 1)? N : 0; //first or second part?
        line_beg = 1;
        for(i=0;i<=degree;i++){
            for(k=0;k<=degree;k++){
                transform[i*vnodes + k] = gen_proj[(i+part)*N + k];
            }
        }
    }
    if(hanging_corner[1] == 3 || hanging_corner[3] == 1){
        //the right face is hanging
        part = (hanging_corner[1] == 3)? N : 0; //first or second part?
        col_end = degree;
        for(i=0;i<=degree;i++){
            for(k=0;k<=degree;k++){
                transform[(i*N+degree)*vnodes + k*N+degree] = gen_proj[(i+part)*N + k];
            }
        }
    }
    if(hanging_corner[0] == 2 || hanging_corner[2] == 0){
        //the left face is hanging
        part = (hanging_corner[0] == 2)? N : 0; //first or second part?
        col_beg = 1;
        for(i=0;i<=degree;i++){
            for(k=0;k<=degree;k++){
                transform[i*N*vnodes + k*N] = gen_proj[(i+part)*N + k];
            }
        }
    }
    if(hanging_corner[2] == 3 || hanging_corner[3] == 2){
        //the top face is hanging
        part = (hanging_corner[2] == 3)? N : 0; //first or second part?
        line_end = degree;
        for(i=0;i<=degree;i++){
            for(k=0;k<=degree;k++){
                transform[(N*degree+i)*vnodes + N*degree+k] = gen_proj[(i+part)*N + k];
            }
        }
    }
    //we finally fill the "interior" with identity
    for(i=line_beg;i<line_end;i++){
        for(k=col_beg;k<col_end;k++){
            transform[(i*N+k)*vnodes + i*N+k] = 1.0;
        }
    }
    //we finish by filling the array interior
    interior[0] = line_beg;
    interior[1] = line_end;
    interior[2] = col_beg;
    interior[3] = col_end;
}



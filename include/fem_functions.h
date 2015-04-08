/* 
 * File:   fem_functions.h
 * Author: camata
 *
 * Created on 18 de Agosto de 2014, 14:45
 */

#ifndef FEM_FUNCTIONS_H
#define	FEM_FUNCTIONS_H


/*
 *   Linear tetrahedra 
 */
void compute_tet_functions(double* xyz, double *phi, double *dphi, double *JxW);


double compute_tet_volume(double *xyz);

double compute_element_quality(double *xyz);

#endif	/* FEM_FUNCTIONS_H */


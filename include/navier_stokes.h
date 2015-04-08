/* 
 * File:   navier_stokes.h
 * Author: camata
 *
 * Created on 18 de Setembro de 2014, 08:39
 */

#ifdef DEVELOP

#ifndef NAVIER_STOKES_H
#define	NAVIER_STOKES_H

#include "equations_system.h"

class NavierStokes: EquationsSystem 
{
    public:
        NavierStokes(Mesh &M);
        solve(double t_step, double tmax);
        
    private:
        void assembly_system();
        double    Reynolds;
        int       n_dof;
};


#endif	/* NAVIER_STOKES_H */

#endif


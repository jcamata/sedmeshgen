/* 
 * File:   fsi.h
 * Author: camata
 *
 * Created on 18 de Agosto de 2014, 14:34
 */

#ifndef FSI_H
#define	FSI_H

#include "equations_system.h"


#define DISPX 0
#define DISPY 1
#define DISPZ 2

class FluidStructure : public EquationsSystem 
{
    
public:
    FluidStructure(Mesh &M);
    void set_disp_from_mapping(const char *fname, NodeTag wall, int dof);
    void solve();
    bool check_volume();
    //void compute_quality();
    //void write_vtk(const char* basename);
    //void write_vtk_legacy(const char* basename);
    
protected:
    void assembly_system();
    void compute_tau();
    void update_mesh();
    
    
    std::vector<int> bc_d;
    
    
    std::vector<double> taue;
    //std::vector<double> quality;
    double vmin;
    double vmax;
};


#endif	/* FSI_H */


/* 
 * File:   equations_system.h
 * Author: camata
 *
 * Created on 18 de Agosto de 2014, 14:28
 */

#ifndef EQUATIONS_SYSTEM_H
#define	EQUATIONS_SYSTEM_H


#include "mesh.h"
#include "petsc.h"


typedef struct {
    NodeTag boundary_wall;
    int     dof;
    double  value;
    std::vector<int> bnd_eq;
} DirichletSurfaceValue;

typedef struct {
    int    node;
    int dof;
    double value;
} DirichletNodalValue;

class EquationsSystem
{
public:
    EquationsSystem(Mesh &M, int n_dof);
    ~EquationsSystem();
    void SetUp();
    void set_number_of_dofs(int ndof);
    void add_dirichlet_surface_value(NodeTag wall, int ndof, double value);
    void add_dirichlet_nodal_value(int node, int dof, double value);
    void apply_dirichlet_condition();
    //void setup_dirichlet_condition();
    int get_equation_number( int local_node_id, int dof_id );
    int get_local_equation_number( int local_node_id, int dof_id );
    //virtual void solve() = 0;
    
   
protected:
    void ScatterSolution();
    
    
    int        n_local_equation;
    int        n_global_equation;
    int        n_dof;
    
    Mat        A;
    Vec        b;
    Vec        x;
    Vec        xcurrent;
    Vec        xold;
    VecScatter scatter;
    KSP        ksp;
    std::vector<int>     dofmap;
    std::vector<double>  dirichlet_values;
    
    std::vector<DirichletSurfaceValue> dirichlet_surface_values;
    std::vector<DirichletNodalValue>   dirichlet_nodal_values;
    
    Mesh *_mesh;
    
};



#endif	/* EQUATIONS_SYSTEM_H */



#ifdef DEVELOP

#include "navier_stokes.h"


#define NDOF 4;
#define UDOF 0;
#define VDOF 1;
#define WDOF 2;
#define PDOF 3;

#define Kuu(i,j) Ke[i*NDOF+UDOF][j*NDOF+UDOF]
#define Kuv(i,j) Ke[i*NDOF+UDOF][j*NDOF+VDOF]
#define Kuw(i,j) Ke[i*NDOF+UDOF][j*NDOF+WDOF]
#define Kup(i,j) Ke[i*NDOF+UDOF][j*NDOF+PDOF]
#define Kvu(i,j) Ke[i*NDOF+VDOF][j*NDOF+UDOF]
#define Kvv(i,j) Ke[i*NDOF+VDOF][j*NDOF+VDOF]
#define Kvw(i,j) Ke[i*NDOF+VDOF][j*NDOF+WDOF]
#define Kvp(i,j) Ke[i*NDOF+VDOF][j*NDOF+PDOF]
#define Kwu(i,j) Ke[i*NDOF+WDOF][j*NDOF+UDOF]
#define Kwv(i,j) Ke[i*NDOF+WDOF][j*NDOF+VDOF]
#define Kww(i,j) Ke[i*NDOF+WDOF][j*NDOF+WDOF]
#define Kwp(i,j) Ke[i*NDOF+WDOF][j*NDOF+PDOF]
#define Kpu(i,j) Ke[i*NDOF+PDOF][j*NDOF+UDOF]
#define Kpv(i,j) Ke[i*NDOF+PDOF][j*NDOF+VDOF]
#define Kpw(i,j) Ke[i*NDOF+PDOF][j*NDOF+WDOF]
#define Kpp(i,j) Ke[i*NDOF+PDOF][j*NDOF+PDOF]


NavierStokes::NavierStokes(Mesh &M):EquationsSystem(M,4){
    
}


void NavierStokes::solve(double t_step, double tmax)
{
    
    double time            = 0.0;
    bool converged          = false;
    int non_linear_step     = 0;
    int max_non_linear_step = 10;
    Vec   N;
    
    VecDuplicate(this->x,N);
    
    SetUp();
    
    while(time < tmax ) 
    {
        time += t_step;
        
        // non linear loop:
        do {
        
            
            
        } while (!converged || non_linear_step > max_non_linear_step);
        
        
    }
    
    
    
    
    
    
}

void NavierStokes::assembly_system()
{
    
    int n_nodes   = this->_mesh->get_n_nodes();
    int n_element = this->_mesh->get_n_elements();
    
    PetscScalar *sol_current;
    PetscScalar *sol_old;
    
    PetscScalar Ke[16][16];
    PetscScalar Fe[16];
    int         eqidx[16];
    double      coords[4][3];
    double      dphi[4][3];
    double      phi[4];
    double      JxW=1.0;
    long        connectivity[4];
    //
    
    VecSet(this->b, 0.0);
    VecGetArray(this->xcurrent, &sol_current);
    
    
    for(int iel = 0; iel < n_element; iel++)
    {
        _mesh->get_element_conn(iel, connectivity);
        
        for(int n=0; n < 4; i++)
        {
             _mesh->get_nodal_coordinate(connectivity[n], &coords[n][0]);
             for(int d = 0; d < this->n_dof; d++)
                 eqidx[n*this->n_dof+d] = this->dofmap[connectivity[n]*this->n_dof+d];
        }
        
        compute_tet_functions(coords,phi,dphi,&JxW);
        
        for(int n = 0; n < 4; n++)
        {
            
            
        }
        
        
        
        
    }
    
    

}

#endif 


#include <assert.h>
#include <vector>
#include <set>

#include "equations_system.h"

EquationsSystem::EquationsSystem(Mesh &M, int n_dof)
{
     this->n_dof       = n_dof;           
     this->_mesh = &M;
     
}

void EquationsSystem::set_number_of_dofs(int ndof)
{
    n_dof = ndof;
}

void EquationsSystem::SetUp()
{
    //IS from;
    //IS to;
    //ISLocalToGlobalMapping ltog;
    
    printf("Setting up equations system...\n");
    int n_nodes = _mesh->get_n_nodes();
    //printf(" % d nodes, %d dofs", n_nodes, this->n_dof);
    this->dofmap.resize(n_nodes*this->n_dof);
    
    if(this->dirichlet_nodal_values.size() == 0)
    {
        printf("No Dirichlet boundary condition was setted up!\n");
        printf("You need to call some add_dirichlet subroutine before set up the equations system!\n");
    }
    
    // set al equations with 0
    //std::fill(this->dofmap.begin(), this->dofmap.end(), 0);
    for(int i = 0; i < this->dirichlet_nodal_values.size(); i++)
    {
        int node_id = this->dirichlet_nodal_values[i].node;
        int dof_id   = this->dirichlet_nodal_values[i].dof;
        
        int   idx = node_id*this->n_dof+dof_id; 
        //printf("%d %d\n", node_id, dof_id);
        if( idx > dofmap.size() ) 
            printf("idx maior\n");
        this->dofmap[idx]      = -1;
    }
    
    //printf("Aqui 1\n");
    
    int n_local_equation = 0;
    for(int i = 0; i < this->dofmap.size(); i++)
    {
        if( this->dofmap[i] == -1 ) continue;
        
        this->dofmap[i] = n_local_equation;
        n_local_equation++;
    }
    
    //printf("Aqui 2\n");
    
  
  
  /*  
    std::vector<PetscInt> gindex;
    
    
    int *fidx   = (int*)malloc(sizeof(int)*n_nodes*n_dof);
    int *tidx   = (int*)malloc(sizeof(int)*n_nodes*n_dof);
    
    // determina numero de equaçoes locais
    n_local_equation = 0;
    for(int i=0; i < n_nodes; i++)
    {
          if(!_mesh->check_nodal_tag(i,REMOTE) ) 
          { 
              n_local_equation++;
              for(int j = 0; j < n_dof; j++) {
                  gindex.push_back(get_equation_number(i,j));
              }
          }
          for(int j = 0; j < n_dof; j++) {
              fidx[i*n_dof+j]  = (_mesh->get_nodal_index(i))*n_dof+j;
              tidx[i*n_dof+j] = i*n_dof+j;
          } 
    }
   * 
   * */
    
    //n_local_equation *=n_dof;
    
    printf("Number of local equations: %d \n", n_local_equation);
    
   // assert(n_local_equation == gindex.size());
    
    /*
    ISLocalToGlobalMappingCreate(PETSC_COMM_WORLD, n_local_equation, &gindex[0], PETSC_COPY_VALUES, &ltog);
    
    ISCreateGeneral(PETSC_COMM_WORLD,n_nodes*n_dof, fidx,PETSC_COPY_VALUES, &from);
    ISCreateGeneral(PETSC_COMM_WORLD,n_nodes*n_dof, tidx,PETSC_COPY_VALUES, &to);
    
    free(fidx);
    free(tidx);
     * */
    
    // Creating parallel vectors
    VecCreateSeq(PETSC_COMM_SELF,n_local_equation, &b);
    VecSetFromOptions(b);
    VecDuplicate(b,&x);
    VecSetFromOptions(x);
    VecSetOption(b,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);
    
    VecCreateSeq(PETSC_COMM_SELF,n_nodes*n_dof, &xcurrent);
    VecDuplicate(xcurrent,&xold);
    //VecScatterCreate(x,from,xcurrent,to,&scatter);
    
    // Prepara a prealocacao da matriz
    std::vector< std::set<int> > diag(n_local_equation);
    
    // Creating parallel matrix
    for(int i=0; i < _mesh->get_n_elements(); i++)
    {
        long conn[4];
        _mesh->get_element_conn(i,conn);
        for(int j=0; j < 4; j++)
        {
            for(int idof = 0; idof < this->n_dof; idof++ )
            {
                int ieq = dofmap[conn[j]*this->n_dof+idof];
                if( ieq == -1) continue;
                 for(int k=0; k < 4; k++)
                {
                    for(int jdof = 0; jdof < this->n_dof; jdof++ )
                    {
                        int jeq = dofmap[conn[k]*this->n_dof+jdof];
                        if( ieq == -1) continue;
                        if( ieq == jeq) continue;
                        diag[ieq].insert(jeq);
                    }
                }
            }
        }
    }
    
    
    int dnz = 0;
    for(int i=0; i < diag.size(); i++) {
        int tmp = diag[i].size();
        if(tmp > dnz) dnz = tmp;
    }
    dnz = dnz + 1;
    
    
    MatCreate(PETSC_COMM_WORLD,&A);
    MatSetType(A, MATSEQAIJ);
    MatSetSizes(A,n_local_equation,n_local_equation, PETSC_DETERMINE, PETSC_DETERMINE );
    MatSetFromOptions(A);
    MatSeqAIJSetPreallocation(A,dnz,PETSC_NULL);
    //MatSetUp(A);
    
    //MatSetLocalToGlobalMapping(A,ltog,ltog);
    
    //ISDestroy(&from);
    //ISDestroy(&to);   
}


EquationsSystem::~EquationsSystem()
{
    if(A != PETSC_NULL)        MatDestroy(&A);
    if(b != PETSC_NULL)        VecDestroy(&b);
    if(x != PETSC_NULL)        VecDestroy(&x);
    if(xcurrent != PETSC_NULL) VecDestroy(&xcurrent);
    if(xold != PETSC_NULL)     VecDestroy(&xold);
    //if(scatter != PETSC_NULL ) VecScatterDestroy(&scatter);
}


 void EquationsSystem::add_dirichlet_surface_value(NodeTag wall, int dof, double value)
{

    
    std::vector<int> nodeset = this->_mesh->getNodesOnBoundary(wall);
    
    for(int i = 0; i < nodeset.size(); i++) {
        DirichletNodalValue dnv;
        dnv.node  = nodeset[i];
        dnv.dof   = dof;
        dnv.value = value;
        this->dirichlet_nodal_values.push_back(dnv);
    }
   
}

void EquationsSystem::add_dirichlet_nodal_value(int node, int dof, double value)
{
    DirichletNodalValue info;
    
    //TODO: Check if wall is valid.
    info.node          = node;
    info.dof           = dof;
    info.value         = value;
    
    dirichlet_nodal_values.push_back(info);
    
}


void EquationsSystem::apply_dirichlet_condition()
{

    
    PetscScalar *xx;
    VecGetArray(this->xcurrent, &xx);
    for(int i = 0; i < this->dirichlet_nodal_values.size(); i++)
    {
        int node = this->dirichlet_nodal_values[i].node;
        int dof  = this->dirichlet_nodal_values[i].dof;
        
        xx[node*this->n_dof+dof] = this->dirichlet_nodal_values[i].value;
    }
    
    VecRestoreArray(this->xcurrent, &xx);
    
}


int EquationsSystem::get_equation_number(int local_node_id, int dof_id)
{
    return (_mesh->get_nodal_index(local_node_id)*n_dof + dof_id);
}

int EquationsSystem::get_local_equation_number(int local_node_id, int dof_id)
{
    return local_node_id*n_dof +dof_id;
}


void EquationsSystem::ScatterSolution()
{
    PetscScalar *to, *from;
    VecGetArray(this->x       , &from);
    VecGetArray(this->xcurrent, &to);
    
    int n_nodes = this->_mesh->get_n_nodes();
    
    for(int node =0; node < n_nodes; node++)
    {
        for(int dof = 0; dof < this->n_dof; dof++)
        {
            int eq = this->dofmap[node*this->n_dof+dof];
            if (eq == -1) continue;
            
            to[node*this->n_dof+dof] = from[eq];
        }
    }
    
    VecRestoreArray(this->xcurrent, &to);
    VecRestoreArray(this->x,        &from);
   
}


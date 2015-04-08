
#include <vector>
#include <mpi.h>
#include <fstream>
#include <string.h>

#include "fsi.h"
#include "fem_functions.h"
#include "mesh.h"

#include <gts.h>
#include <glib.h>

#define DMAX(a,b) a>=b?a:b
#define DMIN(a,b) a<=b?a:b

/*
int openVTKFile(FILE* vtufile, const char* baseName, MPI_Comm comm);
int closeVTKFile(const char* baseName, MPI_Comm comm);
int writeMeshVTKFormatASCII (FILE* vtufile, Mesh *mesh, MPI_Comm comm);
*/

FluidStructure::FluidStructure(Mesh &M):EquationsSystem(M,1) {
}


void FluidStructure::compute_tau()
{
    double        xyz[4][3];
    long          conn[4];
    double        vmax_local = 1.0e-15, vmin_local=1.0e15;
    
    int           n_element = _mesh->get_n_elements();
    taue.resize(n_element);
    
    for(int iel=0; iel < n_element; iel++)
    {
        _mesh->get_element_conn(iel, conn);
        for(int i=0; i < 4; i++)
        {
            _mesh->get_nodal_coordinate(conn[i], &xyz[i][0]);
           
        }
        
        taue[iel] = compute_tet_volume(&xyz[0][0]);
        vmax_local = DMAX(vmax_local, taue[iel]);
        vmin_local = DMIN(vmin_local, taue[iel]);
        
        
    }
    vmax = vmax_local;
    vmin = vmin_local;
    
    //MPI_Allreduce(&vmax_local, &vmax, 1, MPI_DOUBLE, MPI_MAX, PETSC_COMM_WORLD);
    //MPI_Allreduce(&vmax_local, &vmin, 1, MPI_DOUBLE, MPI_MIN, PETSC_COMM_WORLD);
    
    
    double val = vmin/vmax;
    double one_val = 1.0/val;
    
    for(int iel=0; iel < n_element; iel++)
    {
        double t   = taue[iel]*one_val;
        taue[iel]  = (1.0 - val)/t; 
    }
    
}





bool FluidStructure::check_volume()
{
    double        xyz[4][3];
    long          conn[4];
    
    
    int           n_element = _mesh->get_n_elements();
    
    
    for(int iel=0; iel < n_element; iel++)
    {
        _mesh->get_element_conn(iel, conn);
        
        for(int i=0; i < 4; i++)
        {
            _mesh->get_nodal_coordinate(conn[i], &xyz[i][0]);
           
        }
        
        double vol = compute_tet_volume(&xyz[0][0]);
       
        if(vol < 0.0 ) 
        {
            printf("Negative volme found!!\n");
            return false;
        }
    }
    
    return true;
    
}





void FluidStructure::assembly_system()
{
    
    int n_nodes   = _mesh->get_n_nodes();
    int n_element = _mesh->get_n_elements();
    
    PetscScalar *xx;
    PetscScalar Ke[4][4];
    PetscScalar Fe[4];
    int         eqidx[4];
    double      xyz[4][3];
    double      dphi[4][3];
    double      phi[4];
    double      JxW=1.0;
    long        conn[4];
    //
    
    VecSet(b, 0.0);
    
    VecGetArray(this->xcurrent, &xx);
    
    for(int iel=0; iel < n_element; iel++)
    {
        //if(iel % 100000 == 0) printf("assembling matrix %d\n", iel);
        
        
        
        _mesh->get_element_conn(iel, conn);
        
        for(int i=0; i < 4; i++)
        {
           
            _mesh->get_nodal_coordinate(conn[i], &xyz[i][0]);
             eqidx[i] = this->dofmap[conn[i]*this->n_dof];
             for(int j=0; j < 4; j++)
                 Ke[i][j] = 0.0;
        }
        
        
        compute_tet_functions(&xyz[0][0],&phi[0], &dphi[0][0],&JxW);
        
        double diffusivity = 1.0 + taue[iel];
        
        
        for(int i=0; i < 4; i++) {
            for(int j = 0; j < 4; j++)
            {
                
                Ke[i][j] += JxW*diffusivity*( dphi[i][0]*dphi[j][0]  +
                                              dphi[i][1]*dphi[j][1]  + 
                                              dphi[i][2]*dphi[j][2] );
            }
        }
        
        double d0 = xx[conn[0]*this->n_dof];
        double d1 = xx[conn[1]*this->n_dof];
        double d2 = xx[conn[2]*this->n_dof];
        double d3 = xx[conn[3]*this->n_dof];
        
        
        Fe[0] = -(Ke[0][0]*d0 + Ke[0][1]*d1 + Ke[0][2]*d2+Ke[0][3]*d3); 
        Fe[1] = -(Ke[1][0]*d0 + Ke[1][1]*d1 + Ke[1][2]*d2+Ke[1][3]*d3); 
        Fe[2] = -(Ke[2][0]*d0 + Ke[2][1]*d1 + Ke[2][2]*d2+Ke[2][3]*d3); 
        Fe[3] = -(Ke[3][0]*d0 + Ke[3][1]*d1 + Ke[3][2]*d2+Ke[3][3]*d3); 
       
            
         MatSetValues(A,4,eqidx,4,eqidx,&Ke[0][0], ADD_VALUES); 
         VecSetValues(b,4,eqidx, &Fe[0], ADD_VALUES);
         
        
     }   
    
    VecRestoreArray(this->xcurrent, &xx);
    
   
    MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
    VecAssemblyBegin(b);
    
    MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
    VecAssemblyEnd(b);
}

void delete_data(void *data)
{
    if(data) free(data);
}

/*
void readPoints(const char * fname, double scale, kdtree *kd)
{
    std::fstream in;
    double buf[3], x,y,z;
    
    in.open(fname);
    
    while(in.good() )
    {
        if( !(in >> x >> y >> z )) continue;
        buf[0] = x*scale;
        buf[1] = y*scale;
        buf[2] = z*scale;
        double *data = (double*)malloc(sizeof(double));
        *data = buf[2];
        kd_insert(kd,buf, (void*) data);
        
    }
    in.close();
}
 * */


gdouble  distance(GtsPoint *p, gpointer bounded)
{
    GtsTriangle *t = (GtsTriangle*) bounded;
    return gts_point_triangle_distance(p, t);
}


void FluidStructure::set_disp_from_mapping(const char* fname, NodeTag wall, int dof)
{
    double c[3];

    FILE      *gts_file;
    GtsSurface *s;
    GtsPoint   *p;
    GtsFile    *fp;
    GNode      *t;
    
    printf("Setting displacement at bottom wall..\n");

    gts_file = fopen(fname, "r");
    if(!gts_file) 
    {
        fputs ("file does not find!\n", stderr);
	return ;
    }
    fp = gts_file_new(gts_file);
    s  = gts_surface_new(gts_surface_class (),
		         gts_face_class (),
		         gts_edge_class (),
		         gts_vertex_class ());
    if (gts_surface_read (s, fp)) {
        fputs ("file on standard input is not a valid GTS file\n", stderr);
        fprintf (stderr, "stdin:%d:%d: %s\n", fp->line, fp->pos, fp->error);
        return; /* failure */
    }
    
    t = gts_bb_tree_surface(s);
    p = gts_point_new(gts_point_class(), 0.0 , 0.0 ,0.0);
    
    
    std::vector<int> &wall_nodes = _mesh->getNodesOnBoundary(wall);
    for(int i=0; i < wall_nodes.size(); i++)
    {
   
        int node_id  = wall_nodes[i];
        //int global_id = _mesh->get_nodal_index(local_id);
        _mesh->get_nodal_coordinate(node_id, c);

	double factor  = 1.0;
        gts_point_set(p, c[0]*factor, c[1]*factor, c[2]*factor);
        double d = gts_bb_tree_point_distance(t,p,distance,NULL);
        d = 4.6 - (d/factor);
        add_dirichlet_nodal_value(node_id, dof, d);  

              
    }
    
    
    gts_bb_tree_destroy(t, true); 
    printf("%d disp nodes were setted\n", wall_nodes.size()); 
    printf("Setting done. \n");

   
}




void FluidStructure::solve()
{
    printf("Starting Mesh Moviment scheme...\n");
    
    SetUp();
    
    printf("Computing tau^e parameter\n");
    compute_tau();
    
    VecSet(this->xcurrent, 0.0);
    
    printf("Applying dirichlet values...\n");
    apply_dirichlet_condition();
    
    printf("Assembling system...\n");
    assembly_system();
    
    
    MatSetOption(this->A,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);
     
    printf("Solving system\n");
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    
    KSPSetOperators(ksp,A,A);
    KSPSetFromOptions(ksp);
    KSPSolve(ksp,b,x);
     
    ScatterSolution();
    
    KSPDestroy(&ksp);
    
    update_mesh();
     
    
}

void FluidStructure::update_mesh()
{
    double      xyz[3];
    PetscScalar *disp;
    int n_nodes = _mesh->get_n_nodes();
    VecGetArray(xcurrent, &disp);
    for(int i = 0; i < n_nodes; i++)
    {
        _mesh->get_nodal_coordinate(i, xyz);
        
        //xyz[0] += disp[i*n_dof+0];
        //xyz[1] += disp[i*n_dof+1];
        xyz[2] += disp[i*n_dof+0];
        
        _mesh->set_nodal_coordinate(i,xyz[0], xyz[1], xyz[2]);
    }
    
    VecRestoreArray(xcurrent, &disp);
    
    
}


#if 0

void FluidStructure::write_vtk(const char* basename)
{
    FILE* vtufile;
    int mpirank, nprocs;
    char                  vtufilename[100];
    unsigned short       *uint8_data;
    long conn[4];
    int il, jl, sk;
    
    MPI_Comm_size(PETSC_COMM_WORLD, &nprocs);
    MPI_Comm_rank(PETSC_COMM_WORLD, &mpirank);
    
    /* Have each proc write to its own file */
    sprintf (vtufilename, "%s_%04d_%04d.vtu", basename, nprocs, mpirank);
    vtufile = fopen (vtufilename, "w");
    if (vtufile == NULL) {
      printf("Could not open %s for output!\n", vtufilename);
      //return -1;
    }
  
    fprintf (vtufile, "<?xml version=\"1.0\"?>\n");
    fprintf (vtufile, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"");
    fprintf (vtufile, " byte_order=\"LittleEndian\">\n");
    fprintf (vtufile, "  <UnstructuredGrid>\n");
    
    //writeMeshVTKFormatASCII (vtufile, _mesh , PETSC_COMM_WORLD);
    {
        const unsigned long Ncells         = _mesh->get_n_elements();
        const unsigned long Ntotal         = _mesh->get_n_nodes();
        fprintf (vtufile,"    <Piece NumberOfPoints=\"%lld\" NumberOfCells=\"%lld\">\n", (long long) Ntotal, (long long) Ncells);
        fprintf (vtufile, "      <Points>\n");
         /* write point position data */
        fprintf (vtufile, "        <DataArray type=\"%s\" Name=\"Coordinates\" NumberOfComponents=\"3\" format=\"%s\">\n","Float32", "ascii");

        for (unsigned int i =0; i< Ntotal; i++)
        {
            unsigned int d;
            float fx, fy, fz;
            double c[3];
         
            _mesh->get_nodal_coordinate(i,c);
            fx = (float) c[0];
            fy = (float) c[1];
            fz = (float) c[2];
            fprintf(vtufile,"%16.8e %16.8e %16.8e\n",fx,fy,fz);
        }
        fprintf (vtufile, "        </DataArray>\n");
        fprintf (vtufile, "      </Points>\n");
        fprintf (vtufile, "      <Cells>\n");
        /* write connectivity data */
        fprintf (vtufile, "        <DataArray type=\"%s\" Name=\"connectivity\" format=\"%s\">\n", "Int32", "ascii");
        for (il = 0; il < Ncells; ++il)
        {
            _mesh->get_element_conn(il, conn);
            fprintf (vtufile, "          ");
            for(jl=0; jl < 4; jl++)
            fprintf (vtufile, "%lld ",(long long) conn[jl]);
            fprintf(vtufile, "\n");

        }
        fprintf (vtufile, "        </DataArray>\n");
        /* write offset data */
        fprintf (vtufile, "        <DataArray type=\"%s\" Name=\"offsets\" format=\"%s\">\n", "Int32", "ascii");
        fprintf (vtufile, "         ");
        for (il = 1, sk = 1; il <= Ncells; ++il, ++sk) {
            fprintf (vtufile, " %lld", (long long) (4 * il));
            if (!(sk % 4) && il != Ncells)
                fprintf (vtufile, "\n         ");
        }
        fprintf (vtufile, "\n");
        fprintf (vtufile, "        </DataArray>\n");

        /* write type data */
        fprintf (vtufile, "        <DataArray type=\"UInt8\" Name=\"types\" format=\"%s\">\n", "ascii");
        fprintf (vtufile, "         ");
        for (il = 0, sk = 1; il < Ncells; ++il, ++sk) {
            fprintf (vtufile, " 10");
            if (!(sk % 20) && il != (Ncells - 1))
                fprintf (vtufile, "\n         ");
        }
        fprintf (vtufile, "\n");
        fprintf (vtufile, "        </DataArray>\n");
        fprintf (vtufile, "      </Cells>\n");
       
        if( this->quality.size() == Ncells) {
            
            fprintf (vtufile, "      <CellData Scalars=\"Quality\">\n");
            /* write point position data */
            fprintf (vtufile, "        <DataArray type=\"%s\" Name=\"Quality\" format=\"%s\">\n","Float32", "ascii");

            for(int i =0, sk=1; i < Ncells; ++i, ++sk)
            {
                fprintf (vtufile, "%16.8e ", this->quality[i]);
                if (!(sk % 20) && i != (Ncells - 1))
                fprintf (vtufile, "\n         ");
            }
            fprintf (vtufile, "\n");
            fprintf (vtufile, "        </DataArray>\n");
            fprintf (vtufile, "    </CellData>\n");
        }
        fprintf (vtufile, "    </Piece>\n");
        fprintf (vtufile, "  </UnstructuredGrid>\n");
        fprintf (vtufile, "</VTKFile>\n");
        fclose(vtufile);
    
    }
    
    //closeVTKFile(basename,PETSC_COMM_WORLD);
    
}

void FluidStructure::write_vtk_legacy(const char* basename)
{
    int nnodes = _mesh->get_n_nodes();
    int nel    = _mesh->get_n_elements();
    
    FILE* fvtk = fopen(basename, "w");
    fprintf(fvtk,"# vtk DataFile Version 2.0\n");
    fprintf(fvtk,"ASCII string describing the data goes here.\n");
    fprintf(fvtk,"ASCII\n");
    fprintf(fvtk,"DATASET UNSTRUCTURED_GRID\n");
    fprintf(fvtk,"POINTS %d double\n", nnodes);
    for(int i=0; i < nnodes; i++)
    {
        double c[3];
        _mesh->get_nodal_coordinate(i,c);
        fprintf(fvtk,"%16.8e %16.8e %16.8e\n",c[0],c[1],c[2]);
    }
    
    fprintf(fvtk,"\nCELLS %d %d\n", nel, nel*5);
    for(int i =0; i < nel; i++)
    {
        int conn[4];
        _mesh->get_element_conn(i, conn);
        fprintf(fvtk,"%d %d %d %d %d\n",4, conn[0], conn[1], conn[2], conn[3]);
    }
    fprintf(fvtk,"\nCELLS_TYPE %d\n", nel);
    for(int i =0; i < nel; i++)
        fprintf(fvtk,"10\n");
    
    fprintf(fvtk,"\nPOINT_DATA %d\n", nnodes);
    fprintf(fvtk,"SCALARS scalars double 1\n");
    fprintf(fvtk,"LOOKUP_TABLE default\n");
    
    double *scalar;
    VecGetArray(xcurrent, &scalar);
    for(int i=0; i < nnodes; i++)
        fprintf(fvtk,"%16.8e\n", scalar[i]);
    
    VecRestoreArray(xcurrent, &scalar);
     
}
#endif

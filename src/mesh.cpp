/*
 * mesh.cpp
 *
 *  Created on: Mar 5, 2013
 *      Author: camata
 */

#include "petsc.h"


#include "mesh.h"
#include <cassert>
#include <cstdio>
#include <iostream>
using namespace std;


double compute_element_quality(double *xyz);

/*
int openVTKFile(FILE* vtufile, const char* baseName, MPI_Comm comm);
int closeVTKFile(const char* baseName, MPI_Comm comm);
int writeMeshVTKFormatASCII (FILE* vtufile, Mesh *mesh, MPI_Comm comm);
 */

void BoundaryNodes::clean()
{
    for(int i = 0; i < surfaces.size(); i++)
        surfaces[i].clear();
}

Mesh::Mesh()
{
	this->_dim = 3;
	this->_n_elements = 0;
	this->_n_nodes = 0;
}

void    Mesh::set_n_nodes(long n_nodes)
{
	_n_nodes = n_nodes;
	_xyz.resize(_n_nodes*3);
	_node_tag.resize(_n_nodes);
        _node_idx.resize(_n_nodes);
        std::fill(_node_tag.begin(),_node_tag.end(),0);

}

void   Mesh::set_nodal_coordinate(int ith, double x, double y, double z)
{
	assert(ith < _n_nodes);
	_xyz[ith*3  ]  = x;
	_xyz[ith*3+1]  = y;
	_xyz[ith*3+2]  = z;
        _node_idx[ith] = ith;

}


void   Mesh::set_nodal_tag(int ith, NodeTag p)
{
        ith--;
	assert(ith>=0 && ith < _n_nodes);
	_node_tag[ith] = _node_tag[ith] | p;
}

void Mesh::set_nodal_index(int ith, unsigned int index)
{
	this->_node_idx[ith] = index;
}

int Mesh::get_nodal_index(int ith)
{
	return this->_node_idx[ith];
}

void Mesh::get_nodal_coordinate(int ith, double *c)
{
	assert(ith < _n_nodes);
	c[0] = _xyz[ith*3   ];
	c[1] = _xyz[ith*3+1];
	c[2] = _xyz[ith*3+2];
}

bool    Mesh::check_nodal_tag(int ith, NodeTag p)
{
	assert(ith < _n_nodes);
	return ( (_node_tag[ith] & p) == p);
}


unsigned int&  Mesh::nodal_tags(int ith)
{
	return this->_node_tag[ith];
}

unsigned int&  Mesh::element_tags(int ith)
{
	return this->_element_tag[ith];
}

void    Mesh::set_n_element(long n_element)
{
	_n_elements = n_element;
        _elements.resize(_n_elements*4);
        _element_tag.resize(_n_elements);
        std::fill(_element_tag.begin(), _element_tag.end(), 0);

}


void    Mesh::set_element_conn(int ith, int node_0, int node_1, int node_2, int node_3)
{
	assert(ith < _n_elements);
	_elements[ith*4  ]=node_0;
	_elements[ith*4+1]=node_1;
	_elements[ith*4+2]=node_2;
	_elements[ith*4+3]=node_3;

}

void    Mesh::set_element_tag(int ith, ElementTag p)
{
	assert(ith < _n_elements);
	_element_tag[ith] = _element_tag[ith] | p;

}

bool    Mesh::check_element_tag(int ith, ElementTag p)
{
	assert(ith < _n_elements);
	return ( (_element_tag[ith] & p) == p);

}


void    Mesh::get_element_conn(int ith, long *c)
{
    	assert(ith < _n_elements);
	c[0] = _elements[ith*4  ];
	c[1] = _elements[ith*4+1];
	c[2] = _elements[ith*4+2];
	c[3] = _elements[ith*4+3];
}

unsigned long Mesh::get_n_nodes()
{
    return this->_n_nodes;
}

unsigned long Mesh::get_n_elements()
{
    return this->_n_elements;
}

unsigned long* Mesh::get_conn_pointer()
{
	return &this->_elements[0];
}


double* Mesh::get_coords_pointer()
{
    return &_xyz[0];
}

void Mesh::write(const char *fname, IOMode wm)
{
	if (wm == IO_ASCII)
		write_ascii(fname);
	//

}

int Mesh::extract_bnd_nodes()
{
    _bnd_nodes.clean();
    for(int i=0; i < _n_nodes; i++)
    {
        if(check_nodal_tag(i,REMOTE)) continue;
        if(check_nodal_tag(i,NOSLIP)) _bnd_nodes.surfaces[0].push_back(i);
        if(check_nodal_tag(i,SLIPX))    _bnd_nodes.surfaces[1].push_back(i);
        if(check_nodal_tag(i,SLIPY))   _bnd_nodes.surfaces[2].push_back(i);
        if(check_nodal_tag(i,SLIPZ))  _bnd_nodes.surfaces[3].push_back(i);
        if(check_nodal_tag(i,PNULL))   _bnd_nodes.surfaces[4].push_back(i);
        if(check_nodal_tag(i,MOVING))  _bnd_nodes.surfaces[5].push_back(i);
        if(check_nodal_tag(i,SOLID))  _bnd_nodes.surfaces[6].push_back(i);
        if(check_nodal_tag(i,OUTFLOWB))  _bnd_nodes.surfaces[7].push_back(i);
        if(check_nodal_tag(i,OUTFLOWR))  _bnd_nodes.surfaces[8].push_back(i);
        if(check_nodal_tag(i,OUTFLOWL))  _bnd_nodes.surfaces[9].push_back(i);
        if(check_nodal_tag(i,INLET1))  _bnd_nodes.surfaces[10].push_back(i);
        if(check_nodal_tag(i,INLET2))  _bnd_nodes.surfaces[11].push_back(i);
        
        //if(check_nodal_tag(i,MOVING)) _bnd_nodes.surfaces[7].push_back(i);
    }
}


std::vector<int> & Mesh::getNodesOnBoundary(NodeTag wall)
{
        if(wall == NOSLIP )  return _bnd_nodes.surfaces[0];
        if(wall == SLIPX    )  return _bnd_nodes.surfaces[1];
        if(wall == SLIPY   )  return _bnd_nodes.surfaces[2];
        if(wall == SLIPZ  )  return _bnd_nodes.surfaces[3];
        if(wall == PNULL   )  return _bnd_nodes.surfaces[4];
        if(wall == MOVING  )  return _bnd_nodes.surfaces[5];
        if(wall == SOLID  )  return _bnd_nodes.surfaces[6];
        if(wall == OUTFLOWB  )  return _bnd_nodes.surfaces[7];
        if(wall == OUTFLOWR  )  return _bnd_nodes.surfaces[8];
        if(wall == OUTFLOWL  )  return _bnd_nodes.surfaces[9];
        if(wall == INLET1  )  return _bnd_nodes.surfaces[10];
        if(wall == INLET2  )  return _bnd_nodes.surfaces[11];
        
        //if(wall == MOVING  )  return _bnd_nodes.surfaces[7];
        
}


int Mesh::write_ascii(const char * rname)
{
	FILE * fout;
	char fname[80];

        int size;
        int rank;
#if MPI     
        MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
        MPI_Comm_size(PETSC_COMM_WORLD,&size);
#else
        size = 1;
        rank = 0;
#endif        
 
	sprintf(fname,"%s.%d.%d.msh",rname,size, rank);
	fout = fopen(fname, "w");
	if(!fout)
	{
		std::cerr<<"Error: I can not open file: "<<  fname << std::endl;
		return 1;
	}

	fprintf(fout,"$MeshFormat\n");
	fprintf(fout,"%1.1f %d %d %d\n",1.0, IO_ASCII, size, rank);
	fprintf(fout,"$EndMeshFormat\n");
	fprintf(fout,"$Nodes\n");
	fprintf(fout,"%d\n", this->_n_nodes);
	for(int i=0; i < this->_n_nodes; i++)
		fprintf(fout,"%10d %10d %10.10E %10.10E %10.10E\n", this->_node_idx[i], this->_node_tag[i],
				           this->_xyz[i*3],this->_xyz[i*3+1], this->_xyz[i*3+2]);
	fprintf(fout, "$EndNodes\n$");
	fprintf(fout, "$Elements\n");
	fprintf(fout,"%d\n", this->_n_elements);
	for(int i=0; i < this->_n_elements; i++) {
		fprintf(fout,"%10d ", this->_element_tag[i]);
		for(int j=0; j < 4; j++)
			fprintf(fout, "%10d ", this->_elements[i*4+j]);
		fprintf(fout, "\n");
	}
	fprintf(fout,"$EndElements\n");
	fclose(fout);
	return 0;
}



void Mesh::compute_quality()
{
    double        xyz[4][3];
    long          conn[4];
    int qh=0, qmh=0, qmw=0, qw =0; ;
  
    
    
    int           n_element = get_n_elements();
    this->quality.resize(n_element);
    
    
    for(int iel=0; iel < n_element; iel++)
    {
        get_element_conn(iel, conn);
        
        for(int i=0; i < 4; i++)
        {
            get_nodal_coordinate(conn[i], &xyz[i][0]);
           
        }
        
        this->quality[iel] =  compute_element_quality(&xyz[0][0]);
        if(this->quality[iel] > 0.9) { qh++; continue; }
        if(this->quality[iel] > 0.5) { qmh++; continue; }
        if(this->quality[iel] > 0.1) { qmw++; continue; }
        qw++;
       
    }
    
    printf("Element quality Stats: \n");
    printf(" Q_0 > 0.9       : %d elements\n", qh);
    printf(" 0.5 < Q_0 <= 0.9: %d elements\n", qmh);
    printf(" 0.1 < Q_0 <= 0.5: %d elements\n", qmw);
    printf("       Q_0 <= 0.1: %d elements\n", qw);
   
}

void Mesh::write_vtk(const char* basename)
{
    FILE* vtufile;
    int mpirank, nprocs;
    char                  vtufilename[100];
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
        const unsigned long Ncells         = get_n_elements();
        const unsigned long Ntotal         = get_n_nodes();
        fprintf (vtufile,"    <Piece NumberOfPoints=\"%lld\" NumberOfCells=\"%lld\">\n", (long long) Ntotal, (long long) Ncells);
        fprintf (vtufile, "      <Points>\n");
         /* write point position data */
        fprintf (vtufile, "        <DataArray type=\"%s\" Name=\"Coordinates\" NumberOfComponents=\"3\" format=\"%s\">\n","Float32", "ascii");

        for (unsigned int i =0; i< Ntotal; i++)
        {
            unsigned int d;
            float fx, fy, fz;
            double c[3];
         
            get_nodal_coordinate(i,c);
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
            get_element_conn(il, conn);
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
    
}

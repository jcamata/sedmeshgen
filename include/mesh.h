 /*
 * mesh.h
 *
 *  Created on: Mar 5, 2013
 *      Author: camata
 */

#ifndef MESH_H_
#define MESH_H_

#include "define.h"
#include <vector>

/*
typedef enum {
    UNSET_N = 1 << 0,
    BOTTOM  = 1 << 1,
    TOP     = 1 << 2,
    LEFT    = 1 << 3,
    RIGHT   = 1 << 4,
    BACK    = 1 << 5,
    FRONT   = 1 << 6,
    MOVING  = 1 << 7,
    REMOTE  = 1 << 8,
    SOLID   = 1 << 9
} NodeTag;
*/


typedef enum {
    UNSET_N = 1 << 0,
    MOVING  = 1 << 1,  // for moving boundary
    PNULL   = 1 << 2,  // for null pressure
    NOSLIP  = 1 << 3, 
    SLIPX   = 1 << 4,
    SLIPY   = 1 << 5,
    SLIPZ   = 1 << 6,
    SOLID   = 1 << 7,
    REMOTE  = 1 << 8,
    OUTFLOWR = 1 << 9,
    OUTFLOWL = 1 << 10,
    OUTFLOWB = 1 << 11,
    INLET1   = 1 << 12,
    INLET2   = 1 << 13          
} NodeTag;



typedef enum {
    UNSET_E      = 1 << 0,
    FLUID_HEAVY  = 1 << 1,
    FLUID_LIGHT  = 1 << 2,
    FLUID        = 1 << 3,
    SEDIMENT     = 1 << 4
} ElementTag;


typedef enum {
	IO_ASCII  = 1 << 1,
	IO_BINARY = 1 << 2
} IOMode ;


class BoundaryNodes 
{
public:
    BoundaryNodes() {surfaces.resize(12);};
    ~BoundaryNodes() {};
    void clean();
    std::vector< std::vector<int> > surfaces;
    //std::vector<int> top;
    //std::vector<int> bottom;
    //std::vector<int> left;
    //std::vector<int> right;
    //std::vector<int> front;
    //std::vector<int> back;
    //std::vector<int> solid;
};


class Mesh
{
  public:
        Mesh();
        // routines to handle nodes;
	void    set_n_nodes(long n_nodes);
	void    set_nodal_coordinate(int ith, double x, double y, double z);
        void    set_nodal_tag(int ith, NodeTag p);
        unsigned int&  nodal_tags(int ith);
        void    set_nodal_index(int ith, unsigned int idx);
	void    get_nodal_coordinate(int ith, double* c);
        int     get_nodal_index(int ith);
	bool    check_nodal_tag(int ith, NodeTag p);
        unsigned long get_n_nodes();

        void    set_n_element(long n_element);
	void    set_element_conn(int ith, int node_1, int node_2, int node_3, int node_4);
	void    set_element_tag(int ith, ElementTag p);
	bool    check_element_tag(int ith, ElementTag p);
	unsigned int&  element_tags(int ith);
	void    get_element_conn(int ith, long*);
        unsigned long  get_n_elements();
        unsigned long* get_conn_pointer();
        double*  get_coords_pointer();
        void compute_quality();
        void write(const char* fname,  IOMode wm);     
        void write_vtk(const char *fname);
        
        void optimize();
        
        std::vector<int>& getNodesOnBoundary(NodeTag wall);
        int extract_bnd_nodes(); 
       

  private:

    
    int write_ascii(const char* fname);
    int write_binary(const char* binary);


    int                        _dim;         // dimension
    unsigned long               _n_nodes;     // local number of nodes
    unsigned long               _n_elements;  // local number of elements
    std::vector<unsigned int>  _node_tag;    // node tags
    std::vector<unsigned int>  _node_idx;    // global node id
    std::vector<unsigned int>  _element_tag; // element tags
    std::vector<double>        _xyz;         // nodal coordinates
    std::vector<unsigned long> _elements;    // connectivity
    BoundaryNodes              _bnd_nodes;
    std::vector<double>       quality;

};

#endif /* MESH_H_ */

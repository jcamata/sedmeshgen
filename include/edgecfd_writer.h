/* 
 * File:   edgecfd_writer.h
 * Author: camata
 *
 * Created on May 14, 2013, 1:43 PM
 */

#ifndef EDGECFD_WRITER_H
#define	EDGECFD_WRITER_H

#include <vector>
using namespace std;

class Mesh;

typedef struct 
{
    int node;
    int  dof;
    int value;
} InputData;


class EdgeCFDWriter
{
public:
    EdgeCFDWriter(Mesh &M);
    void write(const char* rootname);
   
private:
    
    void write_mesh(const char *fname);
    void write_input(const char* fname);
    Mesh *mesh;
    vector<InputData>   InputDataValues;
};



#endif	/* EDGECFD_WRITER_H */


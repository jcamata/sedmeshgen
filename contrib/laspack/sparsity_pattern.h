/* 
 * File:   sparsity_pattern.h
 * Author: camata
 *
 * Created on 15 de Agosto de 2014, 10:54
 */

#ifndef SPARSITY_PATTERN_H
#define	SPARSITY_PATTERN_H

#include "mesh.h"

#include <vector>
#include <set>
using namespace std;

class SparsityPattern {
public:
    SparsityPattern(Mesh &mesh);
private:
    vector< set<int> > rows;
    int nrows;
};

#endif	/* SPARSITY_PATTERN_H */


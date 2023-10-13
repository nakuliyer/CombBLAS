/****************************************************************/
/* Parallel Combinatorial BLAS Library (for Graph Computations) */
/* version 1.5 -------------------------------------------------*/
/* date: 10/09/2015 ---------------------------------------------*/
/* authors: Ariful Azad, Aydin Buluc, Adam Lugowski ------------*/
/****************************************************************/
/*
 Copyright (c) 2010-2015, The Regents of the University of California
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
 */

#ifndef _COMM_GRID_1D_H_
#define _COMM_GRID_1D_H_

#include <cstdint>
#include <iostream>
#include <cmath>
#include <cassert>
#include <mpi.h>
#include <sstream>
#include <string>
#include <fstream>
#include <stdint.h>
#include "MPIType.h"

namespace combblas {

class CommGrid1D
{
public:
	CommGrid1D(MPI_Comm world)
	{
		MPI_Comm_dup(world, &commWorld);
		MPI_Comm_rank(commWorld, &myrank);
		MPI_Comm_size(commWorld,&worldsize);
	}

	~CommGrid1D(){MPI_Comm_free(&commWorld);}
	CommGrid1D (const CommGrid1D & rhs)
	:worldsize(rhs.worldsize),myrank(rhs.myrank) // copy constructor
	{MPI_Comm_dup(rhs.commWorld, &commWorld);}
	
	CommGrid1D & operator=(const CommGrid1D & rhs)	// assignment operator
	{
		if(this != &rhs)		
		{
			MPI_Comm_free(&commWorld);
			myrank = rhs.myrank;
			worldsize = rhs.worldsize;
			MPI_Comm_dup(rhs.commWorld, &commWorld);
		}
		return *this;
	}
	bool operator== (const CommGrid1D & rhs) const {return *this==rhs;}
	bool operator!= (const CommGrid1D & rhs) const {return (! (*this == rhs));}
	int GetRank(){return myrank;}
	int GetSize(){return worldsize;}
	MPI_Comm GetWorld() const { return commWorld; }
private:
	MPI_Comm commWorld;
	int myrank;
	int worldsize;
	
	
	template <class IT, class NT, class DER>
	friend class SpParMat;
	template <class IT, class NT, class DER>
	friend class SpParMat1D;

	template <class IT, class NT>
	friend class FullyDistSpVec;
};

}

#endif

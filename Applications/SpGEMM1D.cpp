/****************************************************************/
/* Parallel Combinatorial BLAS Library (for Graph Computations) */
/* version 1.6 -------------------------------------------------*/
/* date: 6/15/2017 ---------------------------------------------*/
/* authors: Ariful Azad, Aydin Buluc  --------------------------*/
/****************************************************************/
/*
 Copyright (c) 2010-2017, The Regents of the University of California

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

#include <cstdint>
#include <memory>
#include <mpi.h>
#include <sys/time.h>
#include <iostream>
#include <functional>
#include <algorithm>
#include <vector>
#include <sstream>
#include "CombBLAS/CombBLAS.h"
#include "CombBLAS/SpParMat1D.h"
#include "CombBLAS/ParFriends.h"

using namespace std;
using namespace combblas;

#define EPS 0.0001

#ifdef _OPENMP
int cblas_splits = omp_get_max_threads();
#else
int cblas_splits = 1;
#endif


// Simple helper class for declarations: Just the numerical type is templated
// The index type and the sequential matrix type stays the same for the whole code
// In this case, they are "int" and "SpDCCols"
template <class NT>
class PSpMat
{
public:
    typedef SpDCCols < int64_t, NT > DCCols;
    typedef SpParMat < int64_t, NT, DCCols > MPI_DCCols;
};

int main(int argc, char* argv[])
{
    int nprocs, myrank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

    if(argc < 2){
        if(myrank == 0)
        {
            cout << "Usage: ./<Binary> <MatrixA> " << endl;
        }
        MPI_Finalize();
        return -1;
    }
    else {
        double vm_usage, resident_set;
        string Aname(argv[1]);
        if(myrank == 0){
            fprintf(stderr, "Data: %s\n", argv[1]);
        }
        shared_ptr<CommGrid> fullWorld;
        fullWorld.reset( new CommGrid(MPI_COMM_WORLD, 0, 0) );

        double t0, t1;

        SpParMat<int64_t, double, SpDCCols < int64_t, double >> M(fullWorld);
        
        // Read labelled triple files
        t0 = MPI_Wtime();
        M.ReadGeneralizedTuples(Aname, maximum<double>());
        t1 = MPI_Wtime();
        if(myrank == 0) fprintf(stderr, "Time taken to read file: %lf\n", t1-t0);
        typedef PlusTimesSRing<double, double> PTFF;
    
        // Run 2D multiplication to compare against
        // SpParMat<int64_t, double, SpDCCols < int64_t, double >> A2D(M);
        // SpParMat<int64_t, double, SpDCCols < int64_t, double >> B2D(M);
        // SpParMat<int64_t, double, SpDCCols < int64_t, double >> C2D = 
        // Mult_AnXBn_Synch<PTFF, double, SpDCCols<int64_t, double>, int64_t, double, double, SpDCCols<int64_t, double>, SpDCCols<int64_t, double> >
        // (A2D, B2D);
        
        // do diagonal multi and SB
        SpParMat1D<int64_t, double, SpDCCols<int64_t, double>> A1D_col(M,128,SpParMat1DTYPE::COLWISE); 
        // do (S^TB^T)^T
        SpParMat1D<int64_t, double, SpDCCols<int64_t, double>> A1D_transe(M,128,SpParMat1DTYPE::COLWISE); 
        // do S^2
        SpParMat<int64_t, double, SpDCCols<int64_t, double>> Aoffdiag(M);
        Aoffdiag = Mult_AnXBn_Synch<PTFF, double, SpDCCols<int64_t, double>, int64_t, double, double, SpDCCols<int64_t, double>, SpDCCols<int64_t, double>>(Aoffdiag, Aoffdiag, true, true);

        
        if(myrank == 0) fprintf(stderr, "2D Multiplication done \n");
    }
    MPI_Finalize();
    return 0;
}

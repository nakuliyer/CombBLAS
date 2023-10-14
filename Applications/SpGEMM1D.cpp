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
#include <stdint.h>
#include <sys/time.h>
#include <iostream>
#include <functional>
#include <algorithm>
#include <vector>
#include <sstream>
#include "CombBLAS/CombBLAS.h"
#include "CombBLAS/SpDCCols.h"
#include "CombBLAS/SpDefs.h"
#include "CombBLAS/SpMat.h"
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

typedef SpParMat1D<int64_t, double, SpDCCols<int64_t, double>> Sp1D;
typedef SpParMat<int64_t, double, SpDCCols < int64_t, double >> Sp2D;
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

        SpParMat<int64_t, double, SpDCCols < int64_t, double >> Readingmatrix(fullWorld);
        
        // Read labelled triple files
        t0 = MPI_Wtime();
        // Readingmatrix.ReadGeneralizedTuples(Aname, maximum<double>());
        Readingmatrix.ParallelReadMM(Aname + ".mtx",false,maximum<double>());
        t1 = MPI_Wtime();
        if(myrank == 0) fprintf(stderr, "Time taken to read file: %lf\n", t1-t0);
        cout << Readingmatrix.getnrow() << endl;
        typedef PlusTimesSRing<double, double> PTFF;
        // Run 2D multiplication to compare against
        SpParMat<int64_t, double, SpDCCols < int64_t, double >> A2D(Readingmatrix);
        SpParMat<int64_t, double, SpDCCols < int64_t, double >> B2D(Readingmatrix);
        int64_t totallength = A2D.getnrow();
        const int TMP = (totallength + nprocs-1)/nprocs;
        A2D.KeepDiagBlock(TMP);
        B2D.KeepDiagBlock(TMP);
        SpParMat1D<int64_t, double, SpDCCols<int64_t, double>> A1D(A2D,SpParMat1DTYPE::COLWISE);
        SpParMat1D<int64_t, double, SpDCCols<int64_t, double>> B1D(A2D,SpParMat1DTYPE::COLWISE);
        
        Sp1D C1D = Mult_AnXBn_Diag<PTFF, double, SpDCCols<int64_t, double>, int64_t, double, double, SpDCCols<int64_t, double>, SpDCCols<int64_t, double> >
        (A1D,B1D);
        // // // cout << C1D.getnrow() << ", " << C1D.getncol() << endl;
        Sp2D C2DFrom1D(C1D);
        C2DFrom1D.ParallelWriteMM(Aname+"_From1D.mtx",false);
        cout<< "final in cpp " << C2DFrom1D.getnrow() << ", " << C2DFrom1D.getncol() << endl;
        
        Sp2D C2D = 
        Mult_AnXBn_Synch<PTFF, double, SpDCCols<int64_t, double>, int64_t, double, double, SpDCCols<int64_t, double>, SpDCCols<int64_t, double> >
        (A2D, B2D);
        cout << "in c check " << (C2DFrom1D == C2D) << endl;
        // C2D.ParallelWriteMM(Aname+"_C.mtx",false);
        // eq = (C2D == C2DFrom1D);
        // cout << "diagmulti" << eq << endl;
        // Sp1D C1D = Mult_AnXBn_1D<PTFF, double, SpDCCols<int64_t, double>, int64_t, double, double, SpDCCols<int64_t, double>, SpDCCols<int64_t, double> >
        // (A1D,B1D);
        // SpParMat<int64_t, double, SpDCCols < int64_t, double >> B2D(A2D);
        // t0 = MPI_Wtime();
        // SpParMat<int64_t, double, SpDCCols < int64_t, double >> C2D = 
        // Mult_AnXBn_Synch<PTFF, double, SpDCCols<int64_t, double>, int64_t, double, double, SpDCCols<int64_t, double>, SpDCCols<int64_t, double> >
        // (A2D, B2D);
        // t1 = MPI_Wtime();
        // if(myrank == 0) fprintf(stderr, "Time 1D: %lf\n", t1-t0);
        // Sp2D C2DFrom1D = Sp2D(C1D);
        // eq = (C2D == C2DFrom1D);
        // cout << "mult 1D correct " << eq << endl;
#ifndef NDDEBUG
#endif
        // auto nnz2d = A2D.getnnz();
        // do diagonal and SB
        // Sp1D A1D(A2D,SpParMat1DTYPE::COLWISE); 
        // Sp1D B1D(A2D,SpParMat1DTYPE::COLWISE);
        // Sp1D C1D = Mult_AnXBn_1D<PTFF, double, SpDCCols<int64_t, double>, int64_t, double, double, SpDCCols<int64_t, double>, SpDCCols<int64_t, double> >
        // (A1D, B1D);
        // Sp2D A2D_wodiag(M);
        // A2D_wodiag.RemoveDiagBlock(A1D_col.getblocksize());
        // Sp2D B2D_wodiag(A2D_wodiag);
        // Sp2D C2D_wodiag = Mult_AnXBn_Synch<PTFF, double, SpDCCols<int64_t, double>, int64_t, double, double, SpDCCols<int64_t, double>, SpDCCols<int64_t, double> >
        // (A2D_wodiag, B2D_wodiag);
        // Sp1D A1D_wodiag(C2D_wodiag,SpParMat1DTYPE::COLWISE);
        // A1D_col += A1D_wodiag;
        // A1D_col.allclose(A2D_wodiag);
    }
    MPI_Finalize();
    return 0;
}

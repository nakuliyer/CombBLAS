/****************************************************************/
/* Parallel Combinatorial BLAS Library (for Graph Computations) */
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
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
 */



#include "SpParMat1D.h"
#include "ParFriends.h"
#include "Operations.h"
#include "FileHeader.h"
#include <cassert>
#include <cstdint>
#include <memory>
extern "C" {
#include "mmio.h"
}
#include <sys/types.h>
#include <sys/stat.h>

#include <mpi.h>
#include <fstream>
#include <algorithm>
#include <set>
#include <stdexcept>
#include <string>
#include "CombBLAS/CombBLAS.h"
#include <unistd.h>
#include <memory.h>

namespace combblas
{

    template <class IT, class NT, class DER>
    SpParMat1D<IT, NT, DER>::~SpParMat1D(){
        // No need to delete layermat because it is a smart pointer
        //delete layermat;
    }

    template <class IT, class NT, class DER>
    SpParMat1D< IT,NT,DER >::SpParMat1D (const SpParMat < IT,NT,DER > & A2D, int blocksize, SparMat1DTYPE type)
    :blocksize(blocksize),type(type),total_length(A2D.getnrow())
    {
        typedef typename DER::LocalIT LIT;
        int myrank;
        MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
        commGrid1D = make_shared<CommGrid1D>(MPI_COMM_WORLD);
        auto commGrid2D = A2D.getcommgrid();
        int nprocs = commGrid2D->GetSize();
        IT nrows = A2D.getnrow();
        IT ncols = A2D.getncol();
        assert(nrows == ncols); // 1D only works for sqare matrix
        int pr2d = commGrid2D->GetGridRows();
        int pc2d = commGrid2D->GetGridCols();
        int rowrank2d = commGrid2D->GetRankInProcRow();
        int colrank2d = commGrid2D->GetRankInProcCol();
        IT m_perproc2d = nrows / pr2d;
        IT n_perproc2d = ncols / pc2d;
        DER* spSeq = A2D.seqptr(); // local submatrix
        IT localRowStart2d = colrank2d * m_perproc2d; // first row in this process
        IT localColStart2d = rowrank2d * n_perproc2d; // first col in this process

        LIT lrow1d, lcol1d;
        std::vector<IT> tsendcnt(nprocs,0);
        for(typename DER::SpColIter colit = spSeq->begcol(); colit != spSeq->endcol(); ++colit)
        {
            IT gcol = colit.colid() + localColStart2d;
            for(typename DER::SpColIter::NzIter nzit = spSeq->begnz(colit); nzit != spSeq->endnz(colit); ++nzit)
            {
                IT grow = nzit.rowid() + localRowStart2d;
                int owner = Owner(total_length, grow, gcol, lrow1d, lcol1d, type);
                tsendcnt[owner]++;
            }
        }

        std::vector< std::vector< std::tuple<LIT,LIT, NT> > > sendTuples (nprocs);
        for(typename DER::SpColIter colit = spSeq->begcol(); colit != spSeq->endcol(); ++colit)
        {
            IT gcol = colit.colid() + localColStart2d;
            for(typename DER::SpColIter::NzIter nzit = spSeq->begnz(colit); nzit != spSeq->endnz(colit); ++nzit)
            {
                IT grow = nzit.rowid() + localRowStart2d;
                NT val = nzit.value();
                int owner = Owner(total_length, grow, gcol, lrow1d, lcol1d, type);
                sendTuples[owner].push_back(std::make_tuple(lrow1d, lcol1d, val));
            }
        }

        // LIT datasize;
        // std::tuple<LIT,LIT,NT>* recvTuples = ExchangeData(sendTuples, commGrid2D->GetWorld(), datasize);

        // maxblocksizeperrank = (total_length + nprocs - 1) / nprocs;
        // if (myrank == nprocs-1) maxblocksizeperrank = total_length - maxblocksizeperrank * (nprocs - 1);

    }
    
    /*
     *  Only calculates owner in terms of non-special distribution
     * */
    template <class IT, class NT,class DER>
    template <typename LIT>
    int SpParMat1D<IT,NT,DER>::Owner(IT total_length, IT grow, IT gcol, LIT & lrow, LIT & lcol, SparMat1DTYPE type) const {
        // IT ni = grow / blocksize;
        // IT nj = gcol / blocksize;
        // if (ni != nj) return -1; // the element is not in diag block.
        IT curidx;
        if (type == SparMat1DTYPE::COLWISE)
        {
            curidx = grow;
        }else if(type == SparMat1DTYPE::ROWWISE)
        {
            curidx = gcol;
        }
        IT rowsperrank = (total_length+commGrid1D->worldsize-1) / commGrid1D->worldsize;
        IT ownerrank = curidx / rowsperrank;
        lrow = grow - rowsperrank * commGrid1D->myrank;
        lcol = gcol - rowsperrank * commGrid1D->myrank;
        return ownerrank;
    }

}

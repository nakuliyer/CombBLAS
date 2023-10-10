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
#include "CombBLAS/SpMat.h"
#include "CombBLAS/mtSpGEMM.h"
#include "ParFriends.h"
#include "Operations.h"
#include "FileHeader.h"
#include <cassert>
#include <cstdint>
#include <memory>
#include <tuple>
#include <vector>
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
    using std::cout;
    using std::endl;
    template <class IT, class NT>
    std::tuple<IT,IT,NT>* ExchangeData1D(std::vector<std::vector<std::tuple<IT,IT,NT>>> & tempTuples, MPI_Comm World, IT& datasize)
    {
        /* Create/allocate variables for vector assignment */
        MPI_Datatype MPI_tuple;
        MPI_Type_contiguous(sizeof(std::tuple<IT,IT,NT>), MPI_CHAR, &MPI_tuple);
        MPI_Type_commit(&MPI_tuple);

        int nprocs;
        MPI_Comm_size(World, &nprocs);

        int * sendcnt = new int[nprocs];
        int * recvcnt = new int[nprocs];
        int * sdispls = new int[nprocs]();
        int * rdispls = new int[nprocs]();

        // Set the newly found vector entries
        IT totsend = 0;
        for(IT i=0; i<nprocs; ++i)
        {
            sendcnt[i] = tempTuples[i].size();
            totsend += tempTuples[i].size();
        }

        MPI_Alltoall(sendcnt, 1, MPI_INT, recvcnt, 1, MPI_INT, World);

        std::partial_sum(sendcnt, sendcnt+nprocs-1, sdispls+1);
        std::partial_sum(recvcnt, recvcnt+nprocs-1, rdispls+1);
        IT totrecv = std::accumulate(recvcnt,recvcnt+nprocs, static_cast<IT>(0));

        std::vector< std::tuple<IT,IT,NT> > sendTuples(totsend);
        for(int i=0; i<nprocs; ++i)
        {
            copy(tempTuples[i].begin(), tempTuples[i].end(), sendTuples.data()+sdispls[i]);
            std::vector< std::tuple<IT,IT,NT> >().swap(tempTuples[i]);    // clear memory
        }

        std::tuple<IT,IT,NT>* recvTuples = new std::tuple<IT,IT,NT>[totrecv];
        //std::vector< std::tuple<IT,IT,NT> > recvTuples(totrecv);
        MPI_Alltoallv(sendTuples.data(), sendcnt, sdispls, MPI_tuple, recvTuples, recvcnt, rdispls, MPI_tuple, World);
        DeleteAll(sendcnt, recvcnt, sdispls, rdispls); // free all memory
        MPI_Type_free(&MPI_tuple);
        datasize = totrecv;
        return recvTuples;
    }

    template <class IT, class NT, class DER>
    SpParMat1D<IT, NT, DER>::~SpParMat1D(){
        // No need to delete layermat because it is a smart pointer
        //delete layermat;
    }

    template <class IT, class NT, class DER>
    SpParMat1D< IT,NT,DER >::SpParMat1D (const SpParMat < IT,NT,DER > & A2D, SpParMat1DTYPE mattype)
    :mattype(mattype)
    {
        int myrank;
        MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
        grid1d = std::make_shared<CommGrid1D>(MPI_COMM_WORLD);
        int nprocs = grid1d->worldsize;
        totallength = A2D.getnrow();
        colinmyrank = (totallength + nprocs-1)/nprocs;
        colprefix = colinmyrank * myrank;
        if(myrank == nprocs-1) colinmyrank = totallength - colinmyrank * (nprocs-1);
        blocksize = colinmyrank; // currently blocksize is number of local columns
        auto commGrid2D = A2D.getcommgrid();
        // int nprocs = commGrid2D->GetSize();
        IT nrows = A2D.getnrow();
        IT ncols = A2D.getncol();
        int pr2d = commGrid2D->GetGridRows();
        int pc2d = commGrid2D->GetGridCols();
        int rowrank2d = commGrid2D->GetRankInProcRow();
        int colrank2d = commGrid2D->GetRankInProcCol();
        IT m_perproc2d = nrows / pr2d;
        IT n_perproc2d = ncols / pc2d;
        DER* spSeq = A2D.seqptr(); // local submatrix
        IT localRowStart2d = colrank2d * m_perproc2d; // first row in this process
        IT localColStart2d = rowrank2d * n_perproc2d; // first col in this process
        
        IT lrow1d, lcol1d;
        std::vector<IT> tsendcnt(nprocs,0);
        // for(typename DER::SpColIter colit = spSeq->begcol(); colit != spSeq->endcol(); ++colit)
        // {
        //     IT gcol = colit.colid() + localColStart2d;
        //     for(typename DER::SpColIter::NzIter nzit = spSeq->begnz(colit); nzit != spSeq->endnz(colit); ++nzit)
        //     {
        //         IT grow = nzit.rowid() + localRowStart2d;
        //         int owner = Owner(grow, gcol, lrow1d, lcol1d, mattype);
        //         tsendcnt[owner]++;
        //     }
        // }

        std::vector< std::vector< std::tuple<IT,IT, NT> > > sendTuples (nprocs);
        for(typename DER::SpColIter colit = spSeq->begcol(); colit != spSeq->endcol(); ++colit)
        {
            IT gcol = colit.colid() + localColStart2d;
            for(typename DER::SpColIter::NzIter nzit = spSeq->begnz(colit); nzit != spSeq->endnz(colit); ++nzit)
            {
                IT grow = nzit.rowid() + localRowStart2d;
                NT val = nzit.value();
                int owner = Owner(grow, gcol, lrow1d, lcol1d, mattype);
                sendTuples[owner].push_back(std::make_tuple(lrow1d, lcol1d, val));
            }
        }
        for(int i=0; i<sendTuples.size(); i++) tsendcnt[i] = sendTuples[i].size();
        IT datasize;
        std::tuple<IT,IT,NT>* recvTuples = ExchangeData1D(sendTuples, commGrid2D->GetWorld(), datasize);
        MPI_Barrier(MPI_COMM_WORLD);
        cout << "rank " << myrank << ", " << nrows << ", " << colinmyrank << endl;
        SpTuples<IT, NT>spTuples(datasize, nrows, colinmyrank, recvTuples);
        if(this->spSeq) delete this->spSeq;
        this->spSeq = new DER(spTuples, false);
        cout << "in creation, " << this->spSeq->getnrow() << "," << this->spSeq->getncol() << endl;
    }
    
    /*
     *  Only calculates owner in terms of non-special distribution
     * */
    template <class IT, class NT,class DER>
    int SpParMat1D<IT,NT,DER>::Owner(IT grow, IT gcol, IT & lrow, IT & lcol, SpParMat1DTYPE mattype) const {
        int ownerrank;
        lrow = grow;
        int nprocs = grid1d->worldsize;
        int tmp = (totallength + nprocs - 1) / nprocs;
        ownerrank = gcol / tmp;
        lcol = gcol - ownerrank * tmp;
        assert(mattype == SpParMat1DTYPE::COLWISE); // currently only support colwise split
        return ownerrank;
    }




    template <class IT, class NT,class DER>
    SpParMat1D<IT,NT,DER> SpParMat1D<IT,NT,DER>::Mult_AnXAn_1D()
    {
        int myrank = grid1d->myrank;
        IT diagrowsmax = blocksize * (grid1d->myrank+1);
        // get diag block 
        int nprocs = grid1d->worldsize;
        const int tmp = (totallength + nprocs - 1) / nprocs;

        std::vector< std::tuple<IT,IT, NT>> diagtuple;
        for(typename DER::SpColIter colit = spSeq->begcol(); colit != spSeq->endcol(); ++colit)
        {
            IT gcol = colit.colid() + colprefix;
            for(typename DER::SpColIter::NzIter nzit = spSeq->begnz(colit); nzit != spSeq->endnz(colit); ++nzit)
            {
                IT grow = nzit.rowid();
                if(grow / tmp != gcol / tmp) continue; // by pass offdiag part
                diagtuple.push_back(std::make_tuple(grow - tmp*grid1d->myrank, colit.colid(), nzit.value()));
                if(grow > diagrowsmax)break;
            }
        }
        SpTuples<IT, NT> sptuple(diagtuple.size(),blocksize,blocksize,diagtuple.data());
        sptuple.tuples_deleted = true;

        DER *spmat = new DER(sptuple,false);
        typedef PlusTimesSRing<double, double> PTFF;
        cout << spSeq->getncol() << ", " << spmat->getnrow() << endl;
        auto retSpTuples = LocalSpGEMM<PTFF, double>(*this->spSeq, *spmat, false, false); 
        DER *retSpDER = new DER(*retSpTuples,false);
        // transpose offdiag
        std::vector<IT> sendcnt(nprocs,0);
        std::vector< std::vector< std::tuple<IT,IT, NT> > > sendtuples (nprocs);
        for(typename DER::SpColIter colit = spSeq->begcol(); colit != spSeq->endcol(); ++colit)
        {
            IT gcol = colit.colid() + colprefix;
            for(typename DER::SpColIter::NzIter nzit = spSeq->begnz(colit); nzit != spSeq->endnz(colit); ++nzit)
            {
                IT grow = nzit.rowid();
                if(grow / tmp == gcol / tmp) continue;
                int transposeowner = grow / tmp;
                sendtuples[transposeowner].push_back(std::make_tuple(gcol, grow - transposeowner * tmp, nzit.value()));
            }
        }
        for(int i=0;i<sendtuples.size();i++)sendcnt[i] = sendtuples[i].size();
        IT datasize;
        std::tuple<IT,IT,NT>* recvTuples = ExchangeData1D(sendtuples, MPI_COMM_WORLD, datasize);
        MPI_Barrier(MPI_COMM_WORLD);
        SpTuples<IT, NT> offspTuples(datasize, totallength, colinmyrank, recvTuples);
        DER *offdiagtrans = new DER(offspTuples,false);
        
        spmat->Transpose();
        auto retSpTuples2 = LocalSpGEMM<PTFF, double>(*offdiagtrans, *spmat, false, false);
        DER *retSpDER2 = new DER(*retSpTuples2,false);
        // transpose 
        sendcnt.clear();
        sendtuples.clear();
        for(typename DER::SpColIter colit = spSeq->begcol(); colit != spSeq->endcol(); ++colit)
        {
            IT gcol = colit.colid() + colprefix;
            for(typename DER::SpColIter::NzIter nzit = spSeq->begnz(colit); nzit != spSeq->endnz(colit); ++nzit)
            {
                IT grow = nzit.rowid();
                int transposeowner = grow / tmp;
                sendtuples[transposeowner].push_back(std::make_tuple(gcol, grow - transposeowner * tmp, nzit.value()));
            }
        }
        for(int i=0;i<sendtuples.size();i++)sendcnt[i] = sendtuples[i].size();
        std::tuple<IT,IT,NT>* recvTuples2 = ExchangeData1D(sendtuples, MPI_COMM_WORLD, datasize);
        MPI_Barrier(MPI_COMM_WORLD);
        SpTuples<IT, NT> offspTuples2(datasize, totallength, colinmyrank, recvTuples2);
        offspTuples2.tuples_deleted = false;
        DER * dtrans2 = new DER(offspTuples2, false);
        *retSpDER2 += *dtrans2;
        spSeq = retSpDER2;
        return *this;
    }


    template <class IT, class NT,class DER>
    SpParMat1D< IT,NT,DER > & SpParMat1D< IT,NT,DER >::operator+=(const SpParMat1D< IT,NT,DER > & rhs)
    {
        *this->spSeq += *(rhs.spSeq);
        return *this;
    }
}




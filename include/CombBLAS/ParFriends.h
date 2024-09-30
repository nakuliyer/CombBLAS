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


#ifndef _PAR_FRIENDS_H_
#define _PAR_FRIENDS_H_

#include "mpi.h"
#include <iostream>
#include <chrono>
#include <cstdarg>
#include "SpParMat.h"	
#include "SpParHelper.h"
#include "MPIType.h"
#include "Friends.h"
#include "OptBuf.h"
#include "mtSpGEMM.h"
#include "MultiwayMerge.h"
#include "DnParMat.h"

// #include "mkl.h"

#include <type_traits>

namespace combblas {

template <class IT, class NT, class DER>
class SpParMat;

/*************************************************************************************************/
/**************************** FRIEND FUNCTIONS FOR PARALLEL CLASSES ******************************/
/*************************************************************************************************/


/**
 ** Concatenate all the FullyDistVec<IT,NT> objects into a single one
 **/
template <typename IT, typename NT>
FullyDistVec<IT,NT> Concatenate ( std::vector< FullyDistVec<IT,NT> > & vecs)
{
	if(vecs.size() < 1)
	{
		SpParHelper::Print("Warning: Nothing to concatenate, returning empty ");
		return FullyDistVec<IT,NT>();
	}
	else if (vecs.size() < 2)
	{
		return vecs[1];
	
	}
	else 
	{
		typename std::vector< FullyDistVec<IT,NT> >::iterator it = vecs.begin();
		std::shared_ptr<CommGrid> commGridPtr = it->getcommgrid();
		MPI_Comm World = commGridPtr->GetWorld();
		
		IT nglen = it->TotalLength();	// new global length
		IT cumloclen = it->MyLocLength();	// existing cumulative local lengths 
		++it;
		for(; it != vecs.end(); ++it)
		{
			if(*(commGridPtr) != *(it->getcommgrid()))
			{
				SpParHelper::Print("Grids are not comparable for FullyDistVec<IT,NT>::EWiseApply\n");
				MPI_Abort(MPI_COMM_WORLD, GRIDMISMATCH);
			}
			nglen += it->TotalLength();
			cumloclen += it->MyLocLength();
		}
		FullyDistVec<IT,NT> ConCat (commGridPtr, nglen, NT());	
		int nprocs = commGridPtr->GetSize();
		
		std::vector< std::vector< NT > > data(nprocs);
		std::vector< std::vector< IT > > inds(nprocs);
		IT gloffset = 0;
		for(it = vecs.begin(); it != vecs.end(); ++it)
		{
			IT loclen = it->LocArrSize();
			for(IT i=0; i < loclen; ++i)
			{
				IT locind;
				IT loffset = it->LengthUntil();
				int owner = ConCat.Owner(gloffset+loffset+i, locind);	
				data[owner].push_back(it->arr[i]);
				inds[owner].push_back(locind);
			}
			gloffset += it->TotalLength();
		}
		
		int * sendcnt = new int[nprocs];
		int * sdispls = new int[nprocs];
		for(int i=0; i<nprocs; ++i)
			sendcnt[i] = (int) data[i].size();
		
		int * rdispls = new int[nprocs];
		int * recvcnt = new int[nprocs];
		MPI_Alltoall(sendcnt, 1, MPI_INT, recvcnt, 1, MPI_INT, World);  // share the request counts
		sdispls[0] = 0;
		rdispls[0] = 0;
		for(int i=0; i<nprocs-1; ++i)
		{
			sdispls[i+1] = sdispls[i] + sendcnt[i];
			rdispls[i+1] = rdispls[i] + recvcnt[i];
		}
		IT totrecv = std::accumulate(recvcnt,recvcnt+nprocs,static_cast<IT>(0));
		NT * senddatabuf = new NT[cumloclen];
		for(int i=0; i<nprocs; ++i)
		{
      std::copy(data[i].begin(), data[i].end(), senddatabuf+sdispls[i]);
			std::vector<NT>().swap(data[i]);	// delete data vectors
		}
		NT * recvdatabuf = new NT[totrecv];
		MPI_Alltoallv(senddatabuf, sendcnt, sdispls, MPIType<NT>(), recvdatabuf, recvcnt, rdispls, MPIType<NT>(), World);  // send data
		delete [] senddatabuf;
		
		IT * sendindsbuf = new IT[cumloclen];
		for(int i=0; i<nprocs; ++i)
		{
      std::copy(inds[i].begin(), inds[i].end(), sendindsbuf+sdispls[i]);
			std::vector<IT>().swap(inds[i]);	// delete inds vectors
		}
		IT * recvindsbuf = new IT[totrecv];
		MPI_Alltoallv(sendindsbuf, sendcnt, sdispls, MPIType<IT>(), recvindsbuf, recvcnt, rdispls, MPIType<IT>(), World);  // send new inds
		DeleteAll(sendindsbuf, sendcnt, sdispls);

		for(int i=0; i<nprocs; ++i)
		{
			for(int j = rdispls[i]; j < rdispls[i] + recvcnt[i]; ++j)			
			{
				ConCat.arr[recvindsbuf[j]] = recvdatabuf[j];
			}
		}
		DeleteAll(recvindsbuf, recvcnt, rdispls);
		return ConCat;
	}
}

template <typename MATRIXA, typename MATRIXB>
bool CheckSpGEMMCompliance(const MATRIXA & A, const MATRIXB & B)
{
	if(A.getncol() != B.getnrow())
	{
		std::ostringstream outs;
		outs << "Can not multiply, dimensions does not match"<< std::endl;
		outs << A.getncol() << " != " << B.getnrow() << std::endl;
		SpParHelper::Print(outs.str());
		MPI_Abort(MPI_COMM_WORLD, DIMMISMATCH);
		return false;
	}	
	if((void*) &A == (void*) &B)
	{
		std::ostringstream outs;
		outs << "Can not multiply, inputs alias (make a temporary copy of one of them first)"<< std::endl;
		SpParHelper::Print(outs.str());
		MPI_Abort(MPI_COMM_WORLD, MATRIXALIAS);
		return false;
	}	
	return true;
}	


// Combined logic for prune, recovery, and select
template <typename IT, typename NT, typename DER>
void MCLPruneRecoverySelect(SpParMat<IT,NT,DER> & A, NT hardThreshold, IT selectNum, IT recoverNum, NT recoverPct, int kselectVersion)
{
    
#ifdef TIMING
    double t0, t1;
#endif
    // Prune and create a new pruned matrix
    SpParMat<IT,NT,DER> PrunedA = A.Prune(std::bind2nd(std::less_equal<NT>(), hardThreshold), false);
    // column-wise statistics of the pruned matrix
    FullyDistVec<IT,NT> colSums = PrunedA.Reduce(Column, std::plus<NT>(), 0.0);
    FullyDistVec<IT,NT> nnzPerColumn = PrunedA.Reduce(Column, std::plus<NT>(), 0.0, [](NT val){return 1.0;});
    FullyDistVec<IT,NT> pruneCols(A.getcommgrid(), A.getncol(), hardThreshold);
    PrunedA.FreeMemory();
    
    
    // Check if we need recovery
    // columns with nnz < recoverNum (r)
    FullyDistSpVec<IT,NT> recoverCols(nnzPerColumn, std::bind2nd(std::less<NT>(), recoverNum));
    recoverCols = recoverPct;
    // columns with nnz < r AND sum < recoverPct (pct)
    recoverCols = EWiseApply<NT>(recoverCols, colSums,
                                 [](NT spval, NT dval){return spval;},
                                 [](NT spval, NT dval){return dval < spval;},
                                 false, NT());
    
    IT nrecover = recoverCols.getnnz();
    if(nrecover > 0)
    {
#ifdef TIMING
        t0=MPI_Wtime();
#endif
        A.Kselect(recoverCols, recoverNum, kselectVersion);

#ifdef TIMING
        t1=MPI_Wtime();
        mcl_kselecttime += (t1-t0);
#endif

        pruneCols.Set(recoverCols);

#ifdef COMBBLAS_DEBUG
        std::ostringstream outs;
        outs << "Number of columns needing recovery: " << nrecover << std::endl;
        SpParHelper::Print(outs.str());
#endif
        
    }
    
    
    if(selectNum>0)
    {
        // remaining columns will be up for selection
        FullyDistSpVec<IT,NT> selectCols = EWiseApply<NT>(recoverCols, colSums,
                                                          [](NT spval, NT dval){return spval;},
                                                          [](NT spval, NT dval){return spval==-1;},
                                                          true, static_cast<NT>(-1));
        
        selectCols = selectNum;
        selectCols = EWiseApply<NT>(selectCols, nnzPerColumn,
                                    [](NT spval, NT dval){return spval;},
                                    [](NT spval, NT dval){return dval > spval;},
                                    false, NT());
        IT nselect = selectCols.getnnz();
        
        if(nselect > 0 )
        {
#ifdef TIMING
            t0=MPI_Wtime();
#endif
            A.Kselect(selectCols, selectNum, kselectVersion); // PrunedA would also work
#ifdef TIMING
            t1=MPI_Wtime();
            mcl_kselecttime += (t1-t0);
#endif
        
            pruneCols.Set(selectCols);
#ifdef COMBBLAS_DEBUG
            std::ostringstream outs;
            outs << "Number of columns needing selection: " << nselect << std::endl;
            SpParHelper::Print(outs.str());
#endif
#ifdef TIMING
            t0=MPI_Wtime();
#endif
            SpParMat<IT,NT,DER> selectedA = A.PruneColumn(pruneCols, std::less<NT>(), false);
#ifdef TIMING
            t1=MPI_Wtime();
            mcl_prunecolumntime += (t1-t0);
#endif
            if(recoverNum>0 ) // recovery can be attempted after selection
            {

                FullyDistVec<IT,NT> nnzPerColumn1 = selectedA.Reduce(Column, std::plus<NT>(), 0.0, [](NT val){return 1.0;});
                FullyDistVec<IT,NT> colSums1 = selectedA.Reduce(Column, std::plus<NT>(), 0.0);
                selectedA.FreeMemory();
  
                // slected columns with nnz < recoverNum (r)
                selectCols = recoverNum;
                selectCols = EWiseApply<NT>(selectCols, nnzPerColumn1,
                                            [](NT spval, NT dval){return spval;},
                                            [](NT spval, NT dval){return dval < spval;},
                                            false, NT());
                
                // selected columns with sum < recoverPct (pct)
                selectCols = recoverPct;
                selectCols = EWiseApply<NT>(selectCols, colSums1,
                                            [](NT spval, NT dval){return spval;},
                                            [](NT spval, NT dval){return dval < spval;},
                                            false, NT());
                
                IT n_recovery_after_select = selectCols.getnnz();
                if(n_recovery_after_select>0)
                {
                    // mclExpandVector2 does it on the original vector
                    // mclExpandVector1 does it one pruned vector
#ifdef TIMING
                    t0=MPI_Wtime();
#endif
                    A.Kselect(selectCols, recoverNum, kselectVersion); // Kselect on PrunedA might give different result
#ifdef TIMING
                    t1=MPI_Wtime();
                    mcl_kselecttime += (t1-t0);
#endif
                    pruneCols.Set(selectCols);
                    
#ifdef COMBBLAS_DEBUG
                    std::ostringstream outs1;
                    outs1 << "Number of columns needing recovery after selection: " << nselect << std::endl;
                    SpParHelper::Print(outs1.str());
#endif
                }
                
            }
        }
    }
    

    // final prune
#ifdef TIMING
    t0=MPI_Wtime();
#endif
    A.PruneColumn(pruneCols, std::less<NT>(), true);
#ifdef TIMING
    t1=MPI_Wtime();
    mcl_prunecolumntime += (t1-t0);
#endif
    // Add loops for empty columns
    if(recoverNum<=0 ) // if recoverNum>0, recovery would have added nonzeros in empty columns
    {
        FullyDistVec<IT,NT> nnzPerColumnA = A.Reduce(Column, std::plus<NT>(), 0.0, [](NT val){return 1.0;});
        FullyDistSpVec<IT,NT> emptyColumns(nnzPerColumnA, std::bind2nd(std::equal_to<NT>(), 0.0));
        emptyColumns = 1.00;
        //Ariful: We need a selective AddLoops function with a sparse vector
        //A.AddLoops(emptyColumns);
    }
}




/**
 * Broadcasts A multiple times (#phases) in order to save storage in the output
 * Only uses 1/phases of C memory if the threshold/max limits are proper
 */
template <typename SR, typename NUO, typename UDERO, typename IU, typename NU1, typename NU2, typename UDERA, typename UDERB>
SpParMat<IU,NUO,UDERO> MemEfficientSpGEMM (SpParMat<IU,NU1,UDERA> & A, SpParMat<IU,NU2,UDERB> & B,
                                           int phases, NUO hardThreshold, IU selectNum, IU recoverNum, NUO recoverPct, int kselectVersion, int64_t perProcessMemory)
{
    typedef typename UDERA::LocalIT LIA;
    typedef typename UDERB::LocalIT LIB;
    typedef typename UDERO::LocalIT LIC;
    
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    if(A.getncol() != B.getnrow())
    {
        std::ostringstream outs;
        outs << "Can not multiply, dimensions does not match"<< std::endl;
        outs << A.getncol() << " != " << B.getnrow() << std::endl;
        SpParHelper::Print(outs.str());
        MPI_Abort(MPI_COMM_WORLD, DIMMISMATCH);
        return SpParMat< IU,NUO,UDERO >();
    }
    if(phases <1 || phases >= A.getncol())
    {
        SpParHelper::Print("MemEfficientSpGEMM: The value of phases is too small or large. Resetting to 1.\n");
        phases = 1;
    }
    
    int stages, dummy; 	// last two parameters of ProductGrid are ignored for Synch multiplication
    std::shared_ptr<CommGrid> GridC = ProductGrid((A.commGrid).get(), (B.commGrid).get(), stages, dummy, dummy);
    
    
    if(perProcessMemory>0) // estimate the number of phases permitted by memory
    {
        int p;
        MPI_Comm World = GridC->GetWorld();
        MPI_Comm_size(World,&p);
        
        int64_t perNNZMem_in = sizeof(IU)*2 + sizeof(NU1);
        int64_t perNNZMem_out = sizeof(IU)*2 + sizeof(NUO);
        
        // max nnz(A) in a porcess
        int64_t lannz = A.getlocalnnz();
        int64_t gannz;
        MPI_Allreduce(&lannz, &gannz, 1, MPIType<int64_t>(), MPI_MAX, World);
        int64_t inputMem = gannz * perNNZMem_in * 4; // for four copies (two for SUMMA)
        
        // max nnz(A^2) stored by summa in a porcess
        int64_t asquareNNZ = EstPerProcessNnzSUMMA(A,B);
		int64_t asquareMem = asquareNNZ * perNNZMem_out * 2; // an extra copy in multiway merge and in selection/recovery step
        
        
        // estimate kselect memory
        int64_t d = ceil( (asquareNNZ * sqrt(p))/ B.getlocalcols() ); // average nnz per column in A^2 (it is an overestimate because asquareNNZ is estimated based on unmerged matrices)
        // this is equivalent to (asquareNNZ * p) / B.getcol()
        int64_t k = std::min(int64_t(std::max(selectNum, recoverNum)), d );
        int64_t kselectmem = B.getlocalcols() * k * 8 * 3;
        
        // estimate output memory
        int64_t outputNNZ = (B.getlocalcols() * k)/sqrt(p);
        int64_t outputMem = outputNNZ * perNNZMem_in * 2;
        
        //inputMem + outputMem + asquareMem/phases + kselectmem/phases < memory
        int64_t remainingMem = perProcessMemory*1000000000 - inputMem - outputMem;
        if(remainingMem > 0)
        {
            phases = 1 + (asquareMem+kselectmem) / remainingMem;
        }
        
        
        if(myrank==0)
        {
            if(remainingMem < 0)
            {
                std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n Warning: input and output memory requirement is greater than per-process avaiable memory. Keeping phase to the value supplied at the command line. The program may go out of memory and crash! \n !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
            }
#ifdef SHOW_MEMORY_USAGE
            int64_t maxMemory = kselectmem/phases + inputMem + outputMem + asquareMem / phases;
            if(maxMemory>1000000000)
            std::cout << "phases: " << phases << ": per process memory: " << perProcessMemory << " GB asquareMem: " << asquareMem/1000000000.00 << " GB" << " inputMem: " << inputMem/1000000000.00 << " GB" << " outputMem: " << outputMem/1000000000.00 << " GB" << " kselectmem: " << kselectmem/1000000000.00 << " GB" << std::endl;
            else
            std::cout << "phases: " << phases << ": per process memory: " << perProcessMemory << " GB asquareMem: " << asquareMem/1000000.00 << " MB" << " inputMem: " << inputMem/1000000.00 << " MB" << " outputMem: " << outputMem/1000000.00 << " MB" << " kselectmem: " << kselectmem/1000000.00 << " MB" << std::endl;
#endif
            
        }
    }
    
    LIA C_m = A.spSeq->getnrow();
    LIB C_n = B.spSeq->getncol();
    
    std::vector< UDERB > PiecesOfB;
    UDERB CopyB = *(B.spSeq); // we allow alias matrices as input because of this local copy
    
    CopyB.ColSplit(phases, PiecesOfB); // CopyB's memory is destroyed at this point
    MPI_Barrier(GridC->GetWorld());

    
    LIA ** ARecvSizes = SpHelper::allocate2D<LIA>(UDERA::esscount, stages);
    LIB ** BRecvSizes = SpHelper::allocate2D<LIB>(UDERB::esscount, stages);
    
    static_assert(std::is_same<LIA, LIB>::value, "local index types for both input matrices should be the same");
    static_assert(std::is_same<LIA, LIC>::value, "local index types for input and output matrices should be the same");
    
    
    SpParHelper::GetSetSizes( *(A.spSeq), ARecvSizes, (A.commGrid)->GetRowWorld());
    
    // Remotely fetched matrices are stored as pointers
    UDERA * ARecv;
    UDERB * BRecv;
    
    std::vector< UDERO > toconcatenate;
    
    int Aself = (A.commGrid)->GetRankInProcRow();
    int Bself = (B.commGrid)->GetRankInProcCol();

    for(int p = 0; p< phases; ++p)
    {
        SpParHelper::GetSetSizes( PiecesOfB[p], BRecvSizes, (B.commGrid)->GetColWorld());
        std::vector< SpTuples<LIC,NUO>  *> tomerge;
        for(int i = 0; i < stages; ++i)
        {
            std::vector<LIA> ess;
            if(i == Aself)  ARecv = A.spSeq;	// shallow-copy
            else
            {
                ess.resize(UDERA::esscount);
                for(int j=0; j< UDERA::esscount; ++j)
                    ess[j] = ARecvSizes[j][i];		// essentials of the ith matrix in this row
                ARecv = new UDERA();				// first, create the object
            }
            
#ifdef TIMING
            double t0=MPI_Wtime();
#endif
            SpParHelper::BCastMatrix(GridC->GetRowWorld(), *ARecv, ess, i);	// then, receive its elements
#ifdef TIMING
            double t1=MPI_Wtime();
            mcl_Abcasttime += (t1-t0);
#endif
            ess.clear();

            if(i == Bself)  BRecv = &(PiecesOfB[p]);	// shallow-copy
            else
            {
                ess.resize(UDERB::esscount);
                for(int j=0; j< UDERB::esscount; ++j)
                    ess[j] = BRecvSizes[j][i];
                BRecv = new UDERB();
            }
#ifdef TIMING
            double t2=MPI_Wtime();
#endif
            SpParHelper::BCastMatrix(GridC->GetColWorld(), *BRecv, ess, i);	// then, receive its elements
#ifdef TIMING
            double t3=MPI_Wtime();
            mcl_Bbcasttime += (t3-t2);
#endif
            
            
#ifdef TIMING
            double t4=MPI_Wtime();
#endif
            //SpTuples<LIC,NUO> * C_cont = LocalSpGEMM<SR, NUO>(*ARecv, *BRecv,i != Aself, i != Bself);
            SpTuples<LIC,NUO> * C_cont = LocalHybridSpGEMM<SR, NUO>(*ARecv, *BRecv,i != Aself, i != Bself);

#ifdef TIMING
            double t5=MPI_Wtime();
            mcl_localspgemmtime += (t5-t4);
#endif
            
            if(!C_cont->isZero())
                tomerge.push_back(C_cont);
            else
                delete C_cont;
            
        }   // all stages executed
        
#ifdef SHOW_MEMORY_USAGE
        int64_t gcnnz_unmerged, lcnnz_unmerged = 0;
         for(size_t i = 0; i < tomerge.size(); ++i)
         {
              lcnnz_unmerged += tomerge[i]->getnnz();
         }
        MPI_Allreduce(&lcnnz_unmerged, &gcnnz_unmerged, 1, MPIType<int64_t>(), MPI_MAX, MPI_COMM_WORLD);
        int64_t summa_memory = gcnnz_unmerged*20;//(gannz*2 + phase_nnz + gcnnz_unmerged + gannz + gannz/phases) * 20; // last two for broadcasts
        
        if(myrank==0)
        {
            if(summa_memory>1000000000)
                std::cout << p+1 << ". unmerged: " << summa_memory/1000000000.00 << "GB " ;
            else
                std::cout << p+1 << ". unmerged: " << summa_memory/1000000.00 << " MB " ;
            
        }
#endif

#ifdef TIMING
        double t6=MPI_Wtime();
#endif
        //UDERO OnePieceOfC(MergeAll<SR>(tomerge, C_m, PiecesOfB[p].getncol(),true), false);
        // TODO: MultiwayMerge can directly return UDERO inorder to avoid the extra copy
        SpTuples<LIC,NUO> * OnePieceOfC_tuples = MultiwayMerge<SR>(tomerge, C_m, PiecesOfB[p].getncol(),true);
        
#ifdef SHOW_MEMORY_USAGE
        int64_t gcnnz_merged, lcnnz_merged ;
        lcnnz_merged = OnePieceOfC_tuples->getnnz();
        MPI_Allreduce(&lcnnz_merged, &gcnnz_merged, 1, MPIType<int64_t>(), MPI_MAX, MPI_COMM_WORLD);
       
        // TODO: we can remove gcnnz_merged memory here because we don't need to concatenate anymore
        int64_t merge_memory = gcnnz_merged*2*20;//(gannz*2 + phase_nnz + gcnnz_unmerged + gcnnz_merged*2) * 20;
        
        if(myrank==0)
        {
            if(merge_memory>1000000000)
                std::cout << " merged: " << merge_memory/1000000000.00 << "GB " ;
            else
                std::cout << " merged: " << merge_memory/1000000.00 << " MB " ;
            
        }
#endif
        
        
#ifdef TIMING
        double t7=MPI_Wtime();
        mcl_multiwaymergetime += (t7-t6);
#endif
        UDERO * OnePieceOfC = new UDERO(* OnePieceOfC_tuples, false);
        delete OnePieceOfC_tuples;
        
        SpParMat<IU,NUO,UDERO> OnePieceOfC_mat(OnePieceOfC, GridC);
        MCLPruneRecoverySelect(OnePieceOfC_mat, hardThreshold, selectNum, recoverNum, recoverPct, kselectVersion);
        
#ifdef SHOW_MEMORY_USAGE
        int64_t gcnnz_pruned, lcnnz_pruned ;
        lcnnz_pruned = OnePieceOfC_mat.getlocalnnz();
        MPI_Allreduce(&lcnnz_pruned, &gcnnz_pruned, 1, MPIType<int64_t>(), MPI_MAX, MPI_COMM_WORLD);
        
        
        // TODO: we can remove gcnnz_merged memory here because we don't need to concatenate anymore
        int64_t prune_memory = gcnnz_pruned*2*20;//(gannz*2 + phase_nnz + gcnnz_pruned*2) * 20 + kselectmem; // 3 extra copies of OnePieceOfC_mat, we can make it one extra copy!
        //phase_nnz += gcnnz_pruned;
        
        if(myrank==0)
        {
            if(prune_memory>1000000000)
                std::cout << "Prune: " << prune_memory/1000000000.00 << "GB " << std::endl ;
            else
                std::cout << "Prune: " << prune_memory/1000000.00 << " MB " << std::endl ;
            
        }
#endif
        
        // ABAB: Change this to accept pointers to objects
        toconcatenate.push_back(OnePieceOfC_mat.seq());
    }
    
    
    UDERO * C = new UDERO(0,C_m, C_n,0);
    C->ColConcatenate(toconcatenate); // ABAB: Change this to accept a vector of pointers to pointers to DER objects

    
    SpHelper::deallocate2D(ARecvSizes, UDERA::esscount);
    SpHelper::deallocate2D(BRecvSizes, UDERA::esscount);
    return SpParMat<IU,NUO,UDERO> (C, GridC);
}



/**
 * Parallel C = A*B routine that uses a double buffered broadcasting scheme 
 * @pre { Input matrices, A and B, should not alias }
 * Most memory efficient version available. Total stages: 2*sqrt(p)
 * Memory requirement during first sqrt(p) stages: <= (3/2)*(nnz(A)+nnz(B))+(1/2)*nnz(C)
 * Memory requirement during second sqrt(p) stages: <= nnz(A)+nnz(B)+nnz(C)
 * Final memory requirement: nnz(C) if clearA and clearB are true 
 **/  
template <typename SR, typename NUO, typename UDERO, typename IU, typename NU1, typename NU2, typename UDERA, typename UDERB> 
SpParMat<IU,NUO,UDERO> Mult_AnXBn_DoubleBuff
		(SpParMat<IU,NU1,UDERA> & A, SpParMat<IU,NU2,UDERB> & B, bool clearA = false, bool clearB = false )

{
	if(!CheckSpGEMMCompliance(A,B) )
	{
		return SpParMat< IU,NUO,UDERO >();
	}
	typedef typename UDERA::LocalIT LIA;
    	typedef typename UDERB::LocalIT LIB;
	typedef typename UDERO::LocalIT LIC;

	static_assert(std::is_same<LIA, LIB>::value, "local index types for both input matrices should be the same");
    	static_assert(std::is_same<LIA, LIC>::value, "local index types for input and output matrices should be the same");

	int stages, dummy; 	// last two parameters of ProductGrid are ignored for Synch multiplication
	std::shared_ptr<CommGrid> GridC = ProductGrid((A.commGrid).get(), (B.commGrid).get(), stages, dummy, dummy);
	LIA C_m = A.spSeq->getnrow();
	LIB C_n = B.spSeq->getncol();
    
	UDERA * A1seq = new UDERA();
	UDERA * A2seq = new UDERA(); 
	UDERB * B1seq = new UDERB();
	UDERB * B2seq = new UDERB();
	(A.spSeq)->Split( *A1seq, *A2seq); 
	const_cast< UDERB* >(B.spSeq)->Transpose();
	(B.spSeq)->Split( *B1seq, *B2seq);
    
    	// Transpose back for the column-by-column algorithm
    	const_cast< UDERB* >(B1seq)->Transpose();
    	const_cast< UDERB* >(B2seq)->Transpose();
    
	LIA ** ARecvSizes = SpHelper::allocate2D<LIA>(UDERA::esscount, stages);
	LIB ** BRecvSizes = SpHelper::allocate2D<LIB>(UDERB::esscount, stages);

	SpParHelper::GetSetSizes( *A1seq, ARecvSizes, (A.commGrid)->GetRowWorld());
	SpParHelper::GetSetSizes( *B1seq, BRecvSizes, (B.commGrid)->GetColWorld());

	// Remotely fetched matrices are stored as pointers
	UDERA * ARecv; 
	UDERB * BRecv;
	std::vector< SpTuples<LIC,NUO>  *> tomerge;

	int Aself = (A.commGrid)->GetRankInProcRow();
	int Bself = (B.commGrid)->GetRankInProcCol();	

	for(int i = 0; i < stages; ++i) 
	{
		std::vector<LIA> ess;	
		if(i == Aself)
		{	
			ARecv = A1seq;	// shallow-copy 
		}
		else
		{
			ess.resize(UDERA::esscount);
			for(int j=0; j< UDERA::esscount; ++j)	
			{
				ess[j] = ARecvSizes[j][i];		// essentials of the ith matrix in this row	
			}
			ARecv = new UDERA();				// first, create the object
		}
		SpParHelper::BCastMatrix(GridC->GetRowWorld(), *ARecv, ess, i);	// then, receive its elements	
		ess.clear();	
		if(i == Bself)
		{
			BRecv = B1seq;	// shallow-copy
		}
		else
		{
			ess.resize(UDERB::esscount);		
			for(int j=0; j< UDERB::esscount; ++j)	
			{
				ess[j] = BRecvSizes[j][i];	
			}	
			BRecv = new UDERB();
		}
		SpParHelper::BCastMatrix(GridC->GetColWorld(), *BRecv, ess, i);	// then, receive its elements
		
		// before activating this remove transposing B1seq
        	/*
		SpTuples<LIC,NUO> * C_cont = MultiplyReturnTuples<SR, NUO>
						(*ARecv, *BRecv, // parameters themselves
						false, true,	// transpose information (B is transposed)
						i != Aself, 	// 'delete A' condition
						i != Bself);	// 'delete B' condition
        
        	*/
        
        	SpTuples<LIC,NUO> * C_cont = LocalHybridSpGEMM<SR, NUO>
                        (*ARecv, *BRecv, // parameters themselves
                        i != Aself,    // 'delete A' condition
                        i != Bself);   // 'delete B' condition
        
        
        
		
		if(!C_cont->isZero())
			tomerge.push_back(C_cont);
		else
			delete C_cont;
	}
	if(clearA) delete A1seq;
	if(clearB) delete B1seq;
	
	// Set the new dimensions
	SpParHelper::GetSetSizes( *A2seq, ARecvSizes, (A.commGrid)->GetRowWorld());
	SpParHelper::GetSetSizes( *B2seq, BRecvSizes, (B.commGrid)->GetColWorld());

	// Start the second round
	for(int i = 0; i < stages; ++i) 
	{
		std::vector<LIA> ess;	
		if(i == Aself)
		{	
			ARecv = A2seq;	// shallow-copy 
		}
		else
		{
			ess.resize(UDERA::esscount);
			for(int j=0; j< UDERA::esscount; ++j)	
			{
				ess[j] = ARecvSizes[j][i];		// essentials of the ith matrix in this row	
			}
			ARecv = new UDERA();				// first, create the object
		}

		SpParHelper::BCastMatrix(GridC->GetRowWorld(), *ARecv, ess, i);	// then, receive its elements	
		ess.clear();	
		
		if(i == Bself)
		{
			BRecv = B2seq;	// shallow-copy
		}
		else
		{
			ess.resize(UDERB::esscount);		
			for(int j=0; j< UDERB::esscount; ++j)	
			{
				ess[j] = BRecvSizes[j][i];	
			}	
			BRecv = new UDERB();
		}
		SpParHelper::BCastMatrix(GridC->GetColWorld(), *BRecv, ess, i);	// then, receive its elements

        	// before activating this remove transposing B2seq
        	/*
		SpTuples<LIC,NUO> * C_cont = MultiplyReturnTuples<SR, NUO>
						(*ARecv, *BRecv, // parameters themselves
						false, true,	// transpose information (B is transposed)
						i != Aself, 	// 'delete A' condition
						i != Bself);	// 'delete B' condition
		
        
        	*/
        
        	SpTuples<LIC,NUO> * C_cont = LocalHybridSpGEMM<SR, NUO>
                	(*ARecv, *BRecv, // parameters themselves
                 	i != Aself,    // 'delete A' condition
                 	i != Bself);   // 'delete B' condition
        
		if(!C_cont->isZero())
			tomerge.push_back(C_cont);
		else
			delete C_cont;
	}
	SpHelper::deallocate2D(ARecvSizes, UDERA::esscount);
	SpHelper::deallocate2D(BRecvSizes, UDERB::esscount);
	if(clearA) 
	{
		delete A2seq;
		delete A.spSeq;
		A.spSeq = NULL;
	}
	else
	{
		(A.spSeq)->Merge(*A1seq, *A2seq);
		delete A1seq;
		delete A2seq;
	}
	if(clearB) 
	{
		delete B2seq;
		delete B.spSeq;
		B.spSeq = NULL;	
	}
	else
	{
		B1seq->Transpose();
		B2seq->Transpose();
		(B.spSeq)->Merge(*B1seq, *B2seq);	
		delete B1seq;
		delete B2seq;
		const_cast< UDERB* >(B.spSeq)->Transpose();	// transpose back to original
	}
			
	UDERO * C = new UDERO(MergeAll<SR>(tomerge, C_m, C_n,true), false);
	return SpParMat<IU,NUO,UDERO> (C, GridC);		// return the result object
}


/**
 * Parallel A = B*C routine that uses only MPI-1 features
 * Relies on simple blocking broadcast
 * @pre { Input matrices, A and B, should not alias }
 **/  
template <typename SR, typename NUO, typename UDERO, typename IU, typename NU1, typename NU2, typename UDERA, typename UDERB> 
SpParMat<IU, NUO, UDERO> Mult_AnXBn_Synch 
		(SpParMat<IU,NU1,UDERA> & A, SpParMat<IU,NU2,UDERB> & B, bool clearA = false, bool clearB = false )

{
	if(!CheckSpGEMMCompliance(A,B) )
	{
		return SpParMat< IU,NUO,UDERO >();
	}
	int stages, dummy; 	// last two parameters of ProductGrid are ignored for Synch multiplication
	std::shared_ptr<CommGrid> GridC = ProductGrid((A.commGrid).get(), (B.commGrid).get(), stages, dummy, dummy);		
	IU C_m = A.spSeq->getnrow();
	IU C_n = B.spSeq->getncol();
	
	//const_cast< UDERB* >(B.spSeq)->Transpose(); // do not transpose for colum-by-column multiplication

	IU ** ARecvSizes = SpHelper::allocate2D<IU>(UDERA::esscount, stages);
	IU ** BRecvSizes = SpHelper::allocate2D<IU>(UDERB::esscount, stages);
	
	SpParHelper::GetSetSizes( *(A.spSeq), ARecvSizes, (A.commGrid)->GetRowWorld());
	SpParHelper::GetSetSizes( *(B.spSeq), BRecvSizes, (B.commGrid)->GetColWorld());

	// Remotely fetched matrices are stored as pointers
	UDERA * ARecv; 
	UDERB * BRecv;
	std::vector< SpTuples<IU,NUO>  *> tomerge;

	int Aself = (A.commGrid)->GetRankInProcRow();
	int Bself = (B.commGrid)->GetRankInProcCol();	
	
	for(int i = 0; i < stages; ++i) 
	{
		std::vector<IU> ess;	
		if(i == Aself)
		{	
			ARecv = A.spSeq;	// shallow-copy 
		}
		else
		{
			ess.resize(UDERA::esscount);
			for(int j=0; j< UDERA::esscount; ++j)	
			{
				ess[j] = ARecvSizes[j][i];		// essentials of the ith matrix in this row	
			}
			ARecv = new UDERA();				// first, create the object
		}

		SpParHelper::BCastMatrix(GridC->GetRowWorld(), *ARecv, ess, i);	// then, receive its elements	
		ess.clear();	
		
		if(i == Bself)
		{
			BRecv = B.spSeq;	// shallow-copy
		}
		else
		{
			ess.resize(UDERB::esscount);		
			for(int j=0; j< UDERB::esscount; ++j)	
			{
				ess[j] = BRecvSizes[j][i];	
			}	
			BRecv = new UDERB();
		}
		
		SpParHelper::BCastMatrix(GridC->GetColWorld(), *BRecv, ess, i);	// then, receive its elements

		
		 // before activating this transpose B first
		/*SpTuples<IU,NUO> * C_cont = MultiplyReturnTuples<SR, NUO>
						(*ARecv, *BRecv, // parameters themselves
						false, true,	// transpose information (B is transposed)
						i != Aself, 	// 'delete A' condition
						i != Bself);	// 'delete B' condition
						
         */
		/*
		SpTuples<IU,NUO> * C_cont = LocalSpGEMM<SR, NUO>
						(*ARecv, *BRecv, // parameters themselves
						i != Aself, 	// 'delete A' condition
						i != Bself);	// 'delete B' condition
		
        */
        SpTuples<IU,NUO> * C_cont = LocalHybridSpGEMM<SR, NUO>
                    (*ARecv, *BRecv, // parameters themselves
                     i != Aself,    // 'delete A' condition
                     i != Bself);   // 'delete B' condition
    
		if(!C_cont->isZero()) 
			tomerge.push_back(C_cont);

		#ifdef COMBBLAS_DEBUG
   		std::ostringstream outs;
		outs << i << "th SUMMA iteration"<< std::endl;
		SpParHelper::Print(outs.str());
		#endif
	}
	if(clearA && A.spSeq != NULL) 
	{	
		delete A.spSeq;
		A.spSeq = NULL;
	}	
	if(clearB && B.spSeq != NULL) 
	{
		delete B.spSeq;
		B.spSeq = NULL;
	}

	SpHelper::deallocate2D(ARecvSizes, UDERA::esscount);
	SpHelper::deallocate2D(BRecvSizes, UDERB::esscount);
		
	//UDERO * C = new UDERO(MergeAll<SR>(tomerge, C_m, C_n,true), false);
	// First get the result in SpTuples, then convert to UDER
	// the last parameter to MergeAll deletes tomerge arrays
	
	SpTuples<IU,NUO> * C_tuples = MultiwayMerge<SR>(tomerge, C_m, C_n,true);
	UDERO * C = new UDERO(*C_tuples, false);
    delete C_tuples;

	//if(!clearB)
	//	const_cast< UDERB* >(B.spSeq)->Transpose();	// transpose back to original

	return SpParMat<IU,NUO,UDERO> (C, GridC);		// return the result object
}
    

    
/**
  * Estimate the maximum nnz needed to store in a process from all stages of SUMMA before reduction
  * @pre { Input matrices, A and B, should not alias }
  **/
template <typename IU, typename NU1, typename NU2, typename UDERA, typename UDERB>
int64_t EstPerProcessNnzSUMMA(SpParMat<IU,NU1,UDERA> & A, SpParMat<IU,NU2,UDERB> & B)  
{
    	typedef typename UDERA::LocalIT LIA;
    	typedef typename UDERB::LocalIT LIB;
	static_assert(std::is_same<LIA, LIB>::value, "local index types for both input matrices should be the same");


        int64_t nnzC_SUMMA = 0;
        
        if(A.getncol() != B.getnrow())
        {
            std::ostringstream outs;
            outs << "Can not multiply, dimensions does not match"<< std::endl;
            outs << A.getncol() << " != " << B.getnrow() << std::endl;
            SpParHelper::Print(outs.str());
            MPI_Abort(MPI_COMM_WORLD, DIMMISMATCH);
            return nnzC_SUMMA;
        }
       
        int stages, dummy;     // last two parameters of ProductGrid are ignored for Synch multiplication
        std::shared_ptr<CommGrid> GridC = ProductGrid((A.commGrid).get(), (B.commGrid).get(), stages, dummy, dummy);
  
        MPI_Barrier(GridC->GetWorld());
        
        LIA ** ARecvSizes = SpHelper::allocate2D<LIA>(UDERA::esscount, stages);
        LIB ** BRecvSizes = SpHelper::allocate2D<LIB>(UDERB::esscount, stages);
        SpParHelper::GetSetSizes( *(A.spSeq), ARecvSizes, (A.commGrid)->GetRowWorld());
        SpParHelper::GetSetSizes( *(B.spSeq), BRecvSizes, (B.commGrid)->GetColWorld());
        
        // Remotely fetched matrices are stored as pointers
        UDERA * ARecv;
        UDERB * BRecv;

        int Aself = (A.commGrid)->GetRankInProcRow();
        int Bself = (B.commGrid)->GetRankInProcCol();
        
        
        for(int i = 0; i < stages; ++i)
        {
            std::vector<LIA> ess;
            if(i == Aself)
            {
                ARecv = A.spSeq;    // shallow-copy
            }
            else
            {
                ess.resize(UDERA::esscount);
                for(int j=0; j< UDERA::esscount; ++j)
                {
                    ess[j] = ARecvSizes[j][i];        // essentials of the ith matrix in this row
                }
                ARecv = new UDERA();                // first, create the object
            }
            
            SpParHelper::BCastMatrix(GridC->GetRowWorld(), *ARecv, ess, i);    // then, receive its elements
            ess.clear();
            
            if(i == Bself)
            {
                BRecv = B.spSeq;    // shallow-copy
            }
            else
            {
                ess.resize(UDERB::esscount);
                for(int j=0; j< UDERB::esscount; ++j)
                {
                    ess[j] = BRecvSizes[j][i];
                }
                BRecv = new UDERB();
            }
            
            SpParHelper::BCastMatrix(GridC->GetColWorld(), *BRecv, ess, i);    // then, receive its elements
            
	    // no need to keep entries of colnnzC in larger precision 
	    // because colnnzC is of length nzc and estimates nnzs per column
			// @OGUZ-EDIT Using hash spgemm for estimation
            //LIB * colnnzC = estimateNNZ(*ARecv, *BRecv);
			LIB* flopC = estimateFLOP(*ARecv, *BRecv);
			LIB* colnnzC = estimateNNZ_Hash(*ARecv, *BRecv, flopC);
            if (flopC) delete [] flopC;

            LIB nzc = BRecv->GetDCSC()->nzc;
            int64_t nnzC_stage = 0;
#ifdef THREADED
#pragma omp parallel for reduction (+:nnzC_stage)
#endif
            for (LIB k=0; k<nzc; k++)
            {
                nnzC_stage = nnzC_stage + colnnzC[k];
            }
            nnzC_SUMMA += nnzC_stage;
            if(colnnzC) delete [] colnnzC;

			// sampling-based estimation (comment the estimation above, and
			// comment out below to use)			
			// int64_t nnzC_stage = estimateNNZ_sampling(*ARecv, *BRecv);
			// nnzC_SUMMA += nnzC_stage;
            
            // delete received data
            if(i != Aself)
                delete ARecv;
            if(i != Bself)
                delete BRecv;
        }
        
        SpHelper::deallocate2D(ARecvSizes, UDERA::esscount);
        SpHelper::deallocate2D(BRecvSizes, UDERB::esscount);
        
        int64_t nnzC_SUMMA_max = 0;
        MPI_Allreduce(&nnzC_SUMMA, &nnzC_SUMMA_max, 1, MPIType<int64_t>(), MPI_MAX, GridC->GetWorld());
        
        return nnzC_SUMMA_max;
}



/**
  * Estimate the maximum nnz needed to store in a process from all stages of SUMMA before reduction
  * @pre { Input matrices, A and B, should not alias }
  **/
template <typename IU, typename NU1, typename NU2, typename UDERA, typename UDERB>
int64_t EstPerProcessNnzSUMMAx(SpParMat<IU,NU1,UDERA> & A, SpParMat<IU,NU2,UDERB> & B,
							  int nrounds,
							  std::vector<std::pair<int64_t, double> > &stage_stats,
							  int iter)
{
		int myrank;
    	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
		
    	typedef typename UDERA::LocalIT LIA;
    	typedef typename UDERB::LocalIT LIB;
		static_assert(std::is_same<LIA, LIB>::value, "local index types for both input matrices should be the same");

        int64_t nnzC_SUMMA = 0;
		int64_t nnzC_SUMMA_sampled = 0;
        
        if(A.getncol() != B.getnrow())
        {
            std::ostringstream outs;
            outs << "Can not multiply, dimensions does not match"<< std::endl;
            outs << A.getncol() << " != " << B.getnrow() << std::endl;
            SpParHelper::Print(outs.str());
            MPI_Abort(MPI_COMM_WORLD, DIMMISMATCH);
            return nnzC_SUMMA;
        }
       
        int stages, dummy;     // last two parameters of ProductGrid are ignored for Synch multiplication
        std::shared_ptr<CommGrid> GridC = ProductGrid((A.commGrid).get(), (B.commGrid).get(), stages, dummy, dummy);
  
        MPI_Barrier(GridC->GetWorld());
        
        LIA ** ARecvSizes = SpHelper::allocate2D<LIA>(UDERA::esscount, stages);
        LIB ** BRecvSizes = SpHelper::allocate2D<LIB>(UDERB::esscount, stages);
        SpParHelper::GetSetSizes( *(A.spSeq), ARecvSizes, (A.commGrid)->GetRowWorld());
        SpParHelper::GetSetSizes( *(B.spSeq), BRecvSizes, (B.commGrid)->GetColWorld());
        
        // Remotely fetched matrices are stored as pointers
        UDERA * ARecv;
        UDERB * BRecv;

        int Aself = (A.commGrid)->GetRankInProcRow();
        int Bself = (B.commGrid)->GetRankInProcCol();
        
        
        for(int i = 0; i < stages; ++i)
        {
			double t_tmp = MPI_Wtime();
			
            std::vector<LIA> ess;
            if(i == Aself)
            {
                ARecv = A.spSeq;    // shallow-copy
            }
            else
            {
                ess.resize(UDERA::esscount);
                for(int j=0; j< UDERA::esscount; ++j)
                {
                    ess[j] = ARecvSizes[j][i];        // essentials of the ith matrix in this row
                }
                ARecv = new UDERA();                // first, create the object
            }
            
            SpParHelper::BCastMatrix(GridC->GetRowWorld(), *ARecv, ess, i);    // then, receive its elements
            ess.clear();
            
            if(i == Bself)
            {
                BRecv = B.spSeq;    // shallow-copy
            }
            else
            {
                ess.resize(UDERB::esscount);
                for(int j=0; j< UDERB::esscount; ++j)
                {
                    ess[j] = BRecvSizes[j][i];
                }
                BRecv = new UDERB();
            }
            
            SpParHelper::BCastMatrix(GridC->GetColWorld(), *BRecv, ess, i);    // then, receive its elements

			// *t_est_comm += MPI_Wtime() - t_tmp;
			
            
	    	// no need to keep entries of colnnzC in larger precision 
	    	// because colnnzC is of length nzc and estimates nnzs per column
			// @OGUZ-EDIT Using hash spgemm for estimation
            // LIB * colnnzC = estimateNNZ(*ARecv, *BRecv);
			LIB*	flopC		= estimateFLOP(*ARecv, *BRecv);
			LIB		nzc			= BRecv->GetDCSC()->nzc;
			int64_t flopC_stage = 0;
			int64_t nnzC_stage	= 0;
			if (iter >= 10)				// !!!
			{
				t_tmp = MPI_Wtime();
				LIB* colnnzC = estimateNNZ_Hash(*ARecv, *BRecv, flopC);
				
				#ifdef THREADED
				#pragma omp parallel for reduction (+:nnzC_stage, flopC_stage)
				#endif
            	for (LIB k=0; k<nzc; k++)
            	{
                	nnzC_stage	= nnzC_stage + colnnzC[k];
					flopC_stage = flopC_stage + flopC[k];
            	}
				// *t_est_comp += MPI_Wtime() - t_tmp;

				if(colnnzC)
					delete [] colnnzC;
			}
			else
			{
				t_tmp = MPI_Wtime();
				#ifdef THREADED
				#pragma omp parallel for reduction (+:flopC_stage)
				#endif
				for (LIB k=0; k<nzc; k++)
				{
					flopC_stage = flopC_stage + flopC[k];
				}
				nnzC_stage = estimateNNZ_sampling(*ARecv, *BRecv, nrounds);
				// if (myrank == 0)
				// 	std::cout << "nnzC_stage " << nnzC_stage << "\n" << std::flush;
				// *t_est_comp += MPI_Wtime() - t_tmp;
			}

			stage_stats[i] =
				std::make_pair(flopC_stage,
							   static_cast<double>(flopC_stage) /
							   static_cast<double>(nnzC_stage));

            nnzC_SUMMA += nnzC_stage;
            
            // delete received data
            if(i != Aself)
                delete ARecv;
            if(i != Bself)
                delete BRecv;
        }
        
        SpHelper::deallocate2D(ARecvSizes, UDERA::esscount);
        SpHelper::deallocate2D(BRecvSizes, UDERB::esscount);
        
        int64_t nnzC_SUMMA_max = 0;
        MPI_Allreduce(&nnzC_SUMMA, &nnzC_SUMMA_max, 1, MPIType<int64_t>(), MPI_MAX, GridC->GetWorld());
		int64_t nnzC_SUMMA_tot = 0;
        MPI_Allreduce(&nnzC_SUMMA, &nnzC_SUMMA_tot, 1, MPIType<int64_t>(), MPI_SUM, GridC->GetWorld());
		/* if (myrank == 0) */
		/* 	std::cout << "est max " << nnzC_SUMMA_max */
		/* 			  << " est tot " << nnzC_SUMMA_tot << std::endl; */
		return nnzC_SUMMA_max;
}
    
    
template <typename MATRIX, typename VECTOR>
void CheckSpMVCompliance(const MATRIX & A, const VECTOR & x)
{
	if(A.getncol() != x.TotalLength())
	{
		std::ostringstream outs;
		outs << "Can not multiply, dimensions does not match"<< std::endl;
		outs << A.getncol() << " != " << x.TotalLength() << std::endl;
		SpParHelper::Print(outs.str());
		MPI_Abort(MPI_COMM_WORLD, DIMMISMATCH);
	}
	if(! ( *(A.getcommgrid()) == *(x.getcommgrid())) ) 		
	{
		std::cout << "Grids are not comparable for SpMV" << std::endl; 
		MPI_Abort(MPI_COMM_WORLD, GRIDMISMATCH);
	}
}			


template <typename SR, typename IU, typename NUM, typename UDER> 
FullyDistSpVec<IU,typename promote_trait<NUM,IU>::T_promote>  SpMV 
	(const SpParMat<IU,NUM,UDER> & A, const FullyDistSpVec<IU,IU> & x, bool indexisvalue, OptBuf<int32_t, typename promote_trait<NUM,IU>::T_promote > & optbuf);

template <typename SR, typename IU, typename NUM, typename UDER> 
FullyDistSpVec<IU,typename promote_trait<NUM,IU>::T_promote>  SpMV 
	(const SpParMat<IU,NUM,UDER> & A, const FullyDistSpVec<IU,IU> & x, bool indexisvalue)
{
	typedef typename promote_trait<NUM,IU>::T_promote T_promote;
	OptBuf<int32_t, T_promote > optbuf = OptBuf<int32_t, T_promote >();
	return SpMV<SR>(A, x, indexisvalue, optbuf);
}

/**
 * Step 1 of the sparse SpMV algorithm 
 * @param[in,out]   trxlocnz, lenuntil,trxinds,trxnums  { set or allocated }
 * @param[in] 	indexisvalue	
 **/
template<typename IU, typename NV>
void TransposeVector(MPI_Comm & World, const FullyDistSpVec<IU,NV> & x, int32_t & trxlocnz, IU & lenuntil, int32_t * & trxinds, NV * & trxnums, bool indexisvalue)
{
	int32_t xlocnz = (int32_t) x.getlocnnz();	
	int32_t roffst = (int32_t) x.RowLenUntil();	// since trxinds is int32_t
	int32_t roffset;
	IU luntil = x.LengthUntil();
	int diagneigh = x.commGrid->GetComplementRank();

	MPI_Status status;
	MPI_Sendrecv(&roffst, 1, MPIType<int32_t>(), diagneigh, TROST, &roffset, 1, MPIType<int32_t>(), diagneigh, TROST, World, &status);
	MPI_Sendrecv(&xlocnz, 1, MPIType<int32_t>(), diagneigh, TRNNZ, &trxlocnz, 1, MPIType<int32_t>(), diagneigh, TRNNZ, World, &status);
	MPI_Sendrecv(&luntil, 1, MPIType<IU>(), diagneigh, TRLUT, &lenuntil, 1, MPIType<IU>(), diagneigh, TRLUT, World, &status);
	
	// ABAB: Important observation is that local indices (given by x.ind) is 32-bit addressible
	// Copy them to 32 bit integers and transfer that to save 50% of off-node bandwidth
	trxinds = new int32_t[trxlocnz];
	int32_t * temp_xind = new int32_t[xlocnz];
#ifdef THREADED
#pragma omp parallel for
#endif
	for(int i=0; i< xlocnz; ++i)
        temp_xind[i] = (int32_t) x.ind[i];
	MPI_Sendrecv(temp_xind, xlocnz, MPIType<int32_t>(), diagneigh, TRI, trxinds, trxlocnz, MPIType<int32_t>(), diagneigh, TRI, World, &status);
	delete [] temp_xind;
	if(!indexisvalue)
	{
		trxnums = new NV[trxlocnz];
		MPI_Sendrecv(const_cast<NV*>(SpHelper::p2a(x.num)), xlocnz, MPIType<NV>(), diagneigh, TRX, trxnums, trxlocnz, MPIType<NV>(), diagneigh, TRX, World, &status);
	}
    
  std::transform(trxinds, trxinds+trxlocnz, trxinds, std::bind2nd(std::plus<int32_t>(), roffset)); // fullydist indexing (p pieces) -> matrix indexing (sqrt(p) pieces)
}


/**
 * Step 2 of the sparse SpMV algorithm 
 * @param[in,out]   trxinds, trxnums { deallocated }
 * @param[in,out]   indacc, numacc { allocated }
 * @param[in,out]	accnz { set }
 * @param[in] 		trxlocnz, lenuntil, indexisvalue
 **/
template<typename IU, typename NV>
void AllGatherVector(MPI_Comm & ColWorld, int trxlocnz, IU lenuntil, int32_t * & trxinds, NV * & trxnums, 
					 int32_t * & indacc, NV * & numacc, int & accnz, bool indexisvalue)
{
    int colneighs, colrank;
	MPI_Comm_size(ColWorld, &colneighs);
	MPI_Comm_rank(ColWorld, &colrank);
	int * colnz = new int[colneighs];
	colnz[colrank] = trxlocnz;
	MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, colnz, 1, MPI_INT, ColWorld);
	int * dpls = new int[colneighs]();	// displacements (zero initialized pid) 
	std::partial_sum(colnz, colnz+colneighs-1, dpls+1);
	accnz = std::accumulate(colnz, colnz+colneighs, 0);
	indacc = new int32_t[accnz];
	numacc = new NV[accnz];
	
	// ABAB: Future issues here, colnz is of type int (MPI limitation)
	// What if the aggregate vector size along the processor row/column is not 32-bit addressible?
	// This will happen when n/sqrt(p) > 2^31
	// Currently we can solve a small problem (scale 32) with 4096 processor
	// For a medium problem (scale 35), we'll need 32K processors which gives sqrt(p) ~ 180
	// 2^35 / 180 ~ 2^29 / 3 which is not an issue !
	
#ifdef TIMING
	double t0=MPI_Wtime();
#endif
	MPI_Allgatherv(trxinds, trxlocnz, MPIType<int32_t>(), indacc, colnz, dpls, MPIType<int32_t>(), ColWorld);
	
	delete [] trxinds;
	if(indexisvalue)
	{
		IU lenuntilcol;
		if(colrank == 0)  lenuntilcol = lenuntil;
		MPI_Bcast(&lenuntilcol, 1, MPIType<IU>(), 0, ColWorld);
		for(int i=0; i< accnz; ++i)	// fill numerical values from indices
		{
			numacc[i] = indacc[i] + lenuntilcol;
		}
	}
	else
	{
		MPI_Allgatherv(trxnums, trxlocnz, MPIType<NV>(), numacc, colnz, dpls, MPIType<NV>(), ColWorld);
		delete [] trxnums;
	}	
#ifdef TIMING
	double t1=MPI_Wtime();
	cblas_allgathertime += (t1-t0);
#endif
	DeleteAll(colnz,dpls);
}	



/**
  * Step 3 of the sparse SpMV algorithm, with the semiring 
  * @param[in,out] optbuf {scratch space for all-to-all (fold) communication}
  * @param[in,out] indacc, numacc {index and values of the input vector, deleted upon exit}
  * @param[in,out] sendindbuf, sendnumbuf {index and values of the output vector, created}
 **/
template<typename SR, typename IVT, typename OVT, typename IU, typename NUM, typename UDER>
void LocalSpMV(const SpParMat<IU,NUM,UDER> & A, int rowneighs, OptBuf<int32_t, OVT > & optbuf, int32_t * & indacc, IVT * & numacc, 
			   int32_t * & sendindbuf, OVT * & sendnumbuf, int * & sdispls, int * sendcnt, int accnz, bool indexisvalue, PreAllocatedSPA<OVT> & SPA)
{
    if(optbuf.totmax > 0)	// graph500 optimization enabled
	{ 
		if(A.spSeq->getnsplit() > 0)
		{
			// optbuf.{inds/nums/dspls} and sendcnt are all pre-allocated and only filled by dcsc_gespmv_threaded
			generic_gespmv_threaded_setbuffers<SR> (*(A.spSeq), indacc, numacc, accnz, optbuf.inds, optbuf.nums, sendcnt, optbuf.dspls, rowneighs);
		}
		else
		{
			generic_gespmv<SR> (*(A.spSeq), indacc, numacc, accnz, optbuf.inds, optbuf.nums, sendcnt, optbuf.dspls, rowneighs, indexisvalue);
		}
		DeleteAll(indacc,numacc);
	}
	else
	{
		if(A.spSeq->getnsplit() > 0)
		{
			// sendindbuf/sendnumbuf/sdispls are all allocated and filled by dcsc_gespmv_threaded
			int totalsent = generic_gespmv_threaded<SR> (*(A.spSeq), indacc, numacc, accnz, sendindbuf, sendnumbuf, sdispls, rowneighs, SPA);
			
			DeleteAll(indacc, numacc);
			for(int i=0; i<rowneighs-1; ++i)
				sendcnt[i] = sdispls[i+1] - sdispls[i];
			sendcnt[rowneighs-1] = totalsent - sdispls[rowneighs-1];
		}
		else
		{
            // default SpMSpV
            std::vector< int32_t > indy;
            std::vector< OVT >  numy;
            generic_gespmv<SR>(*(A.spSeq), indacc, numacc, accnz, indy, numy, SPA);	
            
            DeleteAll(indacc, numacc);
            
            int32_t bufsize = indy.size();	// as compact as possible
            sendindbuf = new int32_t[bufsize];
            sendnumbuf = new OVT[bufsize];
            int32_t perproc = A.getlocalrows() / rowneighs;
            
            int k = 0;	// index to buffer
            for(int i=0; i<rowneighs; ++i)
            {
                int32_t end_this = (i==rowneighs-1) ? A.getlocalrows(): (i+1)*perproc;
                while(k < bufsize && indy[k] < end_this)
                {
                    sendindbuf[k] = indy[k] - i*perproc;
                    sendnumbuf[k] = numy[k];
                    ++sendcnt[i];
                    ++k;
                }
            }
            sdispls = new int[rowneighs]();
            std::partial_sum(sendcnt, sendcnt+rowneighs-1, sdispls+1);
            
//#endif

		}
	}

}



// non threaded
template <typename SR, typename IU, typename OVT>
void MergeContributions(int*  listSizes, std::vector<int32_t *> & indsvec, std::vector<OVT *> & numsvec, std::vector<IU>& mergedind, std::vector<OVT>& mergednum)
{
    
    int nlists = indsvec.size();
    // this condition is checked in the caller SpMV function.
    // I am still putting it here for completeness
    if(nlists == 1)
    {
        // simply copy data
        int veclen = listSizes[0];
        mergedind.resize(veclen);
        mergednum.resize(veclen);
        for(int i=0; i<veclen; i++)
        {
            mergedind[i] = indsvec[0][i];
            mergednum[i] = numsvec[0][i];
        }
        return;
    }

    int32_t hsize = 0;
    int32_t inf = std::numeric_limits<int32_t>::min();
    int32_t sup = std::numeric_limits<int32_t>::max();
    KNHeap< int32_t, int32_t > sHeap(sup, inf);
    int * processed = new int[nlists]();
    for(int i=0; i<nlists; ++i)
    {
        if(listSizes[i] > 0)
        {
            // key, list_id
            sHeap.insert(indsvec[i][0], i);
            ++hsize;
        }
    }
    int32_t key, locv;
    if(hsize > 0)
    {
        sHeap.deleteMin(&key, &locv);
        mergedind.push_back( static_cast<IU>(key));
        mergednum.push_back(numsvec[locv][0]);	// nothing is processed yet
        
        if( (++(processed[locv])) < listSizes[locv] )
            sHeap.insert(indsvec[locv][processed[locv]], locv);
        else
            --hsize;
    }
    while(hsize > 0)
    {
        sHeap.deleteMin(&key, &locv);
        if(mergedind.back() == static_cast<IU>(key))
        {
            mergednum.back() = SR::add(mergednum.back(), numsvec[locv][processed[locv]]);
            // ABAB: Benchmark actually allows us to be non-deterministic in terms of parent selection
            // We can just skip this addition operator (if it's a max/min select)
        }
        else
        {
            mergedind.push_back(static_cast<IU>(key));
            mergednum.push_back(numsvec[locv][processed[locv]]);
        }
        
        if( (++(processed[locv])) < listSizes[locv] )
            sHeap.insert(indsvec[locv][processed[locv]], locv);
        else
            --hsize;
    }
    DeleteAll(processed);
}



template <typename SR, typename IU, typename OVT>
void MergeContributions_threaded(int * & listSizes, std::vector<int32_t *> & indsvec, std::vector<OVT *> & numsvec, std::vector<IU> & mergedind, std::vector<OVT> & mergednum, IU maxindex)
{
    
    int nlists = indsvec.size();
    // this condition is checked in the caller SpMV function.
    // I am still putting it here for completeness
    if(nlists == 1)
    {
        // simply copy data
        int veclen = listSizes[0];
        mergedind.resize(veclen);
        mergednum.resize(veclen);
        
#ifdef THREADED
#pragma omp parallel for
#endif
        for(int i=0; i<veclen; i++)
        {
            mergedind[i] = indsvec[0][i];
            mergednum[i] = numsvec[0][i];
        }
        return;
    }
    
    int nthreads=1;
#ifdef THREADED
#pragma omp parallel
    {
        nthreads = omp_get_num_threads();
    }
#endif
    int nsplits = 4*nthreads; // oversplit for load balance
    nsplits = std::min(nsplits, (int)maxindex);
    std::vector< std::vector<int32_t> > splitters(nlists);
    for(int k=0; k< nlists; k++)
    {
        splitters[k].resize(nsplits+1);
        splitters[k][0] = static_cast<int32_t>(0);
#pragma omp parallel for
        for(int i=1; i< nsplits; i++)
        {
            IU cur_idx = i * (maxindex/nsplits);
            auto it = std::lower_bound (indsvec[k], indsvec[k] + listSizes[k], cur_idx);
            splitters[k][i] = (int32_t) (it - indsvec[k]);
        }
        splitters[k][nsplits] = listSizes[k];
    }
    
    // ------ perform merge in parallel ------
    std::vector<std::vector<IU>> indsBuf(nsplits);
    std::vector<std::vector<OVT>> numsBuf(nsplits);
    //TODO: allocate these vectors here before calling MergeContributions
#pragma omp parallel for schedule(dynamic)
    for(int i=0; i< nsplits; i++)
    {
        std::vector<int32_t *> tIndsVec(nlists);
        std::vector<OVT *> tNumsVec(nlists);
        std::vector<int> tLengths(nlists);
        for(int j=0; j< nlists; ++j)
        {
            tIndsVec[j] = indsvec[j] + splitters[j][i];
            tNumsVec[j] = numsvec[j] + splitters[j][i];
            tLengths[j]= splitters[j][i+1] - splitters[j][i];
            
        }
        MergeContributions<SR>(tLengths.data(), tIndsVec, tNumsVec, indsBuf[i], numsBuf[i]);
    }

    // ------ concatenate merged tuples processed by threads ------
    std::vector<IU> tdisp(nsplits+1);
    tdisp[0] = 0;
    for(int i=0; i<nsplits; ++i)
    {
        tdisp[i+1] = tdisp[i] + indsBuf[i].size();
    }
    
    mergedind.resize(tdisp[nsplits]);
    mergednum.resize(tdisp[nsplits]);
    
    
#pragma omp parallel for schedule(dynamic)
    for(int i=0; i< nsplits; i++)
    {
        std::copy(indsBuf[i].data() , indsBuf[i].data() + indsBuf[i].size(), mergedind.data() + tdisp[i]);
        std::copy(numsBuf[i].data() , numsBuf[i].data() + numsBuf[i].size(), mergednum.data() + tdisp[i]);
    }
}


/** 
  * This version is the most flexible sparse matrix X sparse vector [Used in KDT]
  * It accepts different types for the matrix (NUM), the input vector (IVT) and the output vector (OVT)
  * without relying on automatic type promotion
  * Input (x) and output (y) vectors can be ALIASED because y is not written until the algorithm is done with x.
  */
template <typename SR, typename IVT, typename OVT, typename IU, typename NUM, typename UDER>
void SpMV (const SpParMat<IU,NUM,UDER> & A, const FullyDistSpVec<IU,IVT> & x, FullyDistSpVec<IU,OVT> & y, 
			bool indexisvalue, OptBuf<int32_t, OVT > & optbuf, PreAllocatedSPA<OVT> & SPA)
{
	CheckSpMVCompliance(A,x);
	optbuf.MarkEmpty();
    y.glen = A.getnrow(); // in case it is not set already
	
	MPI_Comm World = x.commGrid->GetWorld();
	MPI_Comm ColWorld = x.commGrid->GetColWorld();
	MPI_Comm RowWorld = x.commGrid->GetRowWorld();
	
	int accnz;
	int32_t trxlocnz;
	IU lenuntil;
	int32_t *trxinds, *indacc;
	IVT *trxnums, *numacc;
	
#ifdef TIMING
    double t0=MPI_Wtime();
#endif
    
	TransposeVector(World, x, trxlocnz, lenuntil, trxinds, trxnums, indexisvalue);
    
#ifdef TIMING
    double t1=MPI_Wtime();
    cblas_transvectime += (t1-t0);
#endif
    
    if(x.commGrid->GetGridRows() > 1)
    {
        AllGatherVector(ColWorld, trxlocnz, lenuntil, trxinds, trxnums, indacc, numacc, accnz, indexisvalue);   // trxindS/trxnums deallocated, indacc/numacc allocated, accnz set
    }
    else
    {
        accnz = trxlocnz;
        indacc = trxinds;   // aliasing ptr
        numacc = trxnums;   // aliasing ptr
    }
	
	int rowneighs;
	MPI_Comm_size(RowWorld, &rowneighs);
	int * sendcnt = new int[rowneighs]();	
	int32_t * sendindbuf;	
	OVT * sendnumbuf;
	int * sdispls;
    
#ifdef TIMING
    double t2=MPI_Wtime();
#endif
    
	LocalSpMV<SR>(A, rowneighs, optbuf, indacc, numacc, sendindbuf, sendnumbuf, sdispls, sendcnt, accnz, indexisvalue, SPA);	// indacc/numacc deallocated, sendindbuf/sendnumbuf/sdispls allocated

#ifdef TIMING
    double t3=MPI_Wtime();
    cblas_localspmvtime += (t3-t2);
#endif
	

    if(x.commGrid->GetGridCols() == 1)
    {
        y.ind.resize(sendcnt[0]);
        y.num.resize(sendcnt[0]);


		if(optbuf.totmax > 0 )	// graph500 optimization enabled
		{
#ifdef THREADED
#pragma omp parallel for
#endif
			for(int i=0; i<sendcnt[0]; i++)
			{
				y.ind[i] = optbuf.inds[i];
				y.num[i] = optbuf.nums[i];
			}
		}
		else
		{
#ifdef THREADED
#pragma omp parallel for
#endif
			for(int i=0; i<sendcnt[0]; i++)
			{
				y.ind[i] = sendindbuf[i];
				y.num[i] = sendnumbuf[i];
			}
			DeleteAll(sendindbuf, sendnumbuf,sdispls);
		}
		delete [] sendcnt;
        return;
    }
	int * rdispls = new int[rowneighs];
	int * recvcnt = new int[rowneighs];
	MPI_Alltoall(sendcnt, 1, MPI_INT, recvcnt, 1, MPI_INT, RowWorld);       // share the request counts
	
	// receive displacements are exact whereas send displacements have slack
	rdispls[0] = 0;
	for(int i=0; i<rowneighs-1; ++i)
	{
		rdispls[i+1] = rdispls[i] + recvcnt[i];
	}
	
	int totrecv = std::accumulate(recvcnt,recvcnt+rowneighs,0);	
	int32_t * recvindbuf = new int32_t[totrecv];
	OVT * recvnumbuf = new OVT[totrecv];
	
#ifdef TIMING
	double t4=MPI_Wtime();
#endif
	if(optbuf.totmax > 0 )	// graph500 optimization enabled
	{
		MPI_Alltoallv(optbuf.inds, sendcnt, optbuf.dspls, MPIType<int32_t>(), recvindbuf, recvcnt, rdispls, MPIType<int32_t>(), RowWorld);
		MPI_Alltoallv(optbuf.nums, sendcnt, optbuf.dspls, MPIType<OVT>(), recvnumbuf, recvcnt, rdispls, MPIType<OVT>(), RowWorld);
		delete [] sendcnt;
	}
	else
    {
		MPI_Alltoallv(sendindbuf, sendcnt, sdispls, MPIType<int32_t>(), recvindbuf, recvcnt, rdispls, MPIType<int32_t>(), RowWorld);
		MPI_Alltoallv(sendnumbuf, sendcnt, sdispls, MPIType<OVT>(), recvnumbuf, recvcnt, rdispls, MPIType<OVT>(), RowWorld);
		DeleteAll(sendindbuf, sendnumbuf, sendcnt, sdispls);
	}
#ifdef TIMING
	double t5=MPI_Wtime();
	cblas_alltoalltime += (t5-t4);
#endif
	
#ifdef TIMING
    double t6=MPI_Wtime();
#endif
    //MergeContributions<SR>(y,recvcnt, rdispls, recvindbuf, recvnumbuf, rowneighs);
    // free memory of y, in case it was aliased
    std::vector<IU>().swap(y.ind);
    std::vector<OVT>().swap(y.num);
    
    std::vector<int32_t *> indsvec(rowneighs);
    std::vector<OVT *> numsvec(rowneighs);
    
#ifdef THREADED
#pragma omp parallel for
#endif
    for(int i=0; i<rowneighs; i++)
    {
        indsvec[i] = recvindbuf+rdispls[i];
        numsvec[i] = recvnumbuf+rdispls[i];
    }
#ifdef THREADED
    MergeContributions_threaded<SR>(recvcnt, indsvec, numsvec, y.ind, y.num, y.MyLocLength());
#else
    MergeContributions<SR>(recvcnt, indsvec, numsvec, y.ind, y.num);
#endif
    
    DeleteAll(recvcnt, rdispls,recvindbuf, recvnumbuf);
#ifdef TIMING
    double t7=MPI_Wtime();
    cblas_mergeconttime += (t7-t6);
#endif
    
}


template <typename SR, typename IVT, typename OVT, typename IU, typename NUM, typename UDER>
void SpMV (const SpParMat<IU,NUM,UDER> & A, const FullyDistSpVec<IU,IVT> & x, FullyDistSpVec<IU,OVT> & y, bool indexisvalue, PreAllocatedSPA<OVT> & SPA)
{
	OptBuf< int32_t, OVT > optbuf = OptBuf< int32_t,OVT >(); 
	SpMV<SR>(A, x, y, indexisvalue, optbuf, SPA);
}

template <typename SR, typename IVT, typename OVT, typename IU, typename NUM, typename UDER>
void SpMV (const SpParMat<IU,NUM,UDER> & A, const FullyDistSpVec<IU,IVT> & x, FullyDistSpVec<IU,OVT> & y, bool indexisvalue)
{
    OptBuf< int32_t, OVT > optbuf = OptBuf< int32_t,OVT >();
    PreAllocatedSPA<OVT> SPA;
    SpMV<SR>(A, x, y, indexisvalue, optbuf, SPA);
}

template <typename SR, typename IVT, typename OVT, typename IU, typename NUM, typename UDER>
void SpMV (const SpParMat<IU,NUM,UDER> & A, const FullyDistSpVec<IU,IVT> & x, FullyDistSpVec<IU,OVT> & y, bool indexisvalue, OptBuf<int32_t, OVT > & optbuf)
{
	PreAllocatedSPA<OVT> SPA;
	SpMV<SR>(A, x, y, indexisvalue, optbuf, SPA);
}


/**
 * Automatic type promotion is ONLY done here, all the callee functions (in Friends.h and below) are initialized with the promoted type
 * If indexisvalues = true, then we do not need to transfer values for x (happens for BFS iterations with boolean matrices and integer rhs vectors)
 **/
template <typename SR, typename IU, typename NUM, typename UDER>
FullyDistSpVec<IU,typename promote_trait<NUM,IU>::T_promote>  SpMV 
(const SpParMat<IU,NUM,UDER> & A, const FullyDistSpVec<IU,IU> & x, bool indexisvalue, OptBuf<int32_t, typename promote_trait<NUM,IU>::T_promote > & optbuf)
{		
	typedef typename promote_trait<NUM,IU>::T_promote T_promote;
	FullyDistSpVec<IU, T_promote> y ( x.getcommgrid(), A.getnrow());	// identity doesn't matter for sparse vectors
	SpMV<SR>(A, x, y, indexisvalue, optbuf);
	return y;
}

/**
 * Parallel dense SpMV
 **/ 
template <typename SR, typename IU, typename NUM, typename NUV, typename UDER> 
FullyDistVec<IU,typename promote_trait<NUM,NUV>::T_promote>  SpMV 
	(const SpParMat<IU,NUM,UDER> & A, const FullyDistVec<IU,NUV> & x )
{
	typedef typename promote_trait<NUM,NUV>::T_promote T_promote;
	CheckSpMVCompliance(A, x);

	MPI_Comm World = x.commGrid->GetWorld();
	MPI_Comm ColWorld = x.commGrid->GetColWorld();
	MPI_Comm RowWorld = x.commGrid->GetRowWorld();

	int xsize = (int) x.LocArrSize();
	int trxsize = 0;

	int diagneigh = x.commGrid->GetComplementRank();
	MPI_Status status;
	MPI_Sendrecv(&xsize, 1, MPI_INT, diagneigh, TRX, &trxsize, 1, MPI_INT, diagneigh, TRX, World, &status);
	
	NUV * trxnums = new NUV[trxsize];
	MPI_Sendrecv(const_cast<NUV*>(SpHelper::p2a(x.arr)), xsize, MPIType<NUV>(), diagneigh, TRX, trxnums, trxsize, MPIType<NUV>(), diagneigh, TRX, World, &status);

        int colneighs, colrank;
	MPI_Comm_size(ColWorld, &colneighs);
	MPI_Comm_rank(ColWorld, &colrank);
	int * colsize = new int[colneighs];
	colsize[colrank] = trxsize;
	MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, colsize, 1, MPI_INT, ColWorld);
	int * dpls = new int[colneighs]();	// displacements (zero initialized pid) 
	std::partial_sum(colsize, colsize+colneighs-1, dpls+1);
	int accsize = std::accumulate(colsize, colsize+colneighs, 0);
	NUV * numacc = new NUV[accsize];

	MPI_Allgatherv(trxnums, trxsize, MPIType<NUV>(), numacc, colsize, dpls, MPIType<NUV>(), ColWorld);
	delete [] trxnums;

	// serial SpMV with dense vector
	T_promote id = SR::id();
	IU ysize = A.getlocalrows();
	T_promote * localy = new T_promote[ysize];
	std::fill_n(localy, ysize, id);		

#ifdef THREADED
	dcsc_gespmv_threaded<SR>(*(A.spSeq), numacc, localy);
#else
	dcsc_gespmv<SR>(*(A.spSeq), numacc, localy);	
#endif
	

	DeleteAll(numacc,colsize, dpls);

	// FullyDistVec<IT,NT>(shared_ptr<CommGrid> grid, IT globallen, NT initval, NT id)
	FullyDistVec<IU, T_promote> y ( x.commGrid, A.getnrow(), id);
	
	int rowneighs;
	MPI_Comm_size(RowWorld, &rowneighs);

	IU begptr, endptr;
	for(int i=0; i< rowneighs; ++i)
	{
		begptr = y.RowLenUntil(i);
		if(i == rowneighs-1)
		{
			endptr = ysize;
		}
		else
		{
			endptr = y.RowLenUntil(i+1);
		}
		MPI_Reduce(localy+begptr, SpHelper::p2a(y.arr), endptr-begptr, MPIType<T_promote>(), SR::mpi_op(), i, RowWorld);
	}
	delete [] localy;
	return y;
}

	
/**
 * \TODO: Old version that is no longer considered optimal
 * Kept for legacy purposes
 * To be removed when other functionals are fully tested.
 **/ 
template <typename SR, typename IU, typename NUM, typename NUV, typename UDER> 
FullyDistSpVec<IU,typename promote_trait<NUM,NUV>::T_promote>  SpMV 
	(const SpParMat<IU,NUM,UDER> & A, const FullyDistSpVec<IU,NUV> & x)
{
	typedef typename promote_trait<NUM,NUV>::T_promote T_promote;
	CheckSpMVCompliance(A, x);

	MPI_Comm World = x.commGrid->GetWorld();
	MPI_Comm ColWorld = x.commGrid->GetColWorld();
	MPI_Comm RowWorld = x.commGrid->GetRowWorld();

	int xlocnz = (int) x.getlocnnz();
	int trxlocnz = 0;
	int roffst = x.RowLenUntil();
	int offset;

	int diagneigh = x.commGrid->GetComplementRank();
	MPI_Status status;
	MPI_Sendrecv(&xlocnz, 1, MPI_INT, diagneigh, TRX, &trxlocnz, 1, MPI_INT, diagneigh, TRX, World, &status);
	MPI_Sendrecv(&roffst, 1, MPI_INT, diagneigh, TROST, &offset, 1, MPI_INT, diagneigh, TROST, World, &status);
	
	IU * trxinds = new IU[trxlocnz];
	NUV * trxnums = new NUV[trxlocnz];
	MPI_Sendrecv(const_cast<IU*>(SpHelper::p2a(x.ind)), xlocnz, MPIType<IU>(), diagneigh, TRX, trxinds, trxlocnz, MPIType<IU>(), diagneigh, TRX, World, &status);
	MPI_Sendrecv(const_cast<NUV*>(SpHelper::p2a(x.num)), xlocnz, MPIType<NUV>(), diagneigh, TRX, trxnums, trxlocnz, MPIType<NUV>(), diagneigh, TRX, World, &status);
  std::transform(trxinds, trxinds+trxlocnz, trxinds, std::bind2nd(std::plus<IU>(), offset)); // fullydist indexing (n pieces) -> matrix indexing (sqrt(p) pieces)

        int colneighs, colrank;
	MPI_Comm_size(ColWorld, &colneighs);
	MPI_Comm_rank(ColWorld, &colrank);
	int * colnz = new int[colneighs];
	colnz[colrank] = trxlocnz;
	MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, colnz, 1, MPI_INT, ColWorld);
	int * dpls = new int[colneighs]();	// displacements (zero initialized pid) 
	std::partial_sum(colnz, colnz+colneighs-1, dpls+1);
	int accnz = std::accumulate(colnz, colnz+colneighs, 0);
	IU * indacc = new IU[accnz];
	NUV * numacc = new NUV[accnz];

	// ABAB: Future issues here, colnz is of type int (MPI limitation)
	// What if the aggregate vector size along the processor row/column is not 32-bit addressible?
	MPI_Allgatherv(trxinds, trxlocnz, MPIType<IU>(), indacc, colnz, dpls, MPIType<IU>(), ColWorld);
	MPI_Allgatherv(trxnums, trxlocnz, MPIType<NUV>(), numacc, colnz, dpls, MPIType<NUV>(), ColWorld);
	DeleteAll(trxinds, trxnums);

	// serial SpMV with sparse vector
	std::vector< int32_t > indy;
	std::vector< T_promote >  numy;
	
        int32_t * tmpindacc = new int32_t[accnz];
        for(int i=0; i< accnz; ++i) tmpindacc[i] = indacc[i];
	delete [] indacc;

	dcsc_gespmv<SR>(*(A.spSeq), tmpindacc, numacc, accnz, indy, numy);	// actual multiplication

	DeleteAll(tmpindacc, numacc);
	DeleteAll(colnz, dpls);

	FullyDistSpVec<IU, T_promote> y ( x.commGrid, A.getnrow());	// identity doesn't matter for sparse vectors
	IU yintlen = y.MyRowLength();

	int rowneighs;
	MPI_Comm_size(RowWorld,&rowneighs);
	std::vector< std::vector<IU> > sendind(rowneighs);
	std::vector< std::vector<T_promote> > sendnum(rowneighs);
	typename std::vector<int32_t>::size_type outnz = indy.size();
	for(typename std::vector<IU>::size_type i=0; i< outnz; ++i)
	{
		IU locind;
		int rown = y.OwnerWithinRow(yintlen, static_cast<IU>(indy[i]), locind);
		sendind[rown].push_back(locind);
		sendnum[rown].push_back(numy[i]);
	}

	IU * sendindbuf = new IU[outnz];
	T_promote * sendnumbuf = new T_promote[outnz];
	int * sendcnt = new int[rowneighs];
	int * sdispls = new int[rowneighs];
	for(int i=0; i<rowneighs; ++i)
		sendcnt[i] = sendind[i].size();

	int * rdispls = new int[rowneighs];
	int * recvcnt = new int[rowneighs];
	MPI_Alltoall(sendcnt, 1, MPI_INT, recvcnt, 1, MPI_INT, RowWorld);       // share the request counts

	sdispls[0] = 0;
	rdispls[0] = 0;
	for(int i=0; i<rowneighs-1; ++i)
	{
		sdispls[i+1] = sdispls[i] + sendcnt[i];
		rdispls[i+1] = rdispls[i] + recvcnt[i];
	}
	int totrecv = std::accumulate(recvcnt,recvcnt+rowneighs,0);
	IU * recvindbuf = new IU[totrecv];
	T_promote * recvnumbuf = new T_promote[totrecv];

	for(int i=0; i<rowneighs; ++i)
	{
    std::copy(sendind[i].begin(), sendind[i].end(), sendindbuf+sdispls[i]);
		std::vector<IU>().swap(sendind[i]);
	}
	for(int i=0; i<rowneighs; ++i)
	{
    std::copy(sendnum[i].begin(), sendnum[i].end(), sendnumbuf+sdispls[i]);
		std::vector<T_promote>().swap(sendnum[i]);
	}
	MPI_Alltoallv(sendindbuf, sendcnt, sdispls, MPIType<IU>(), recvindbuf, recvcnt, rdispls, MPIType<IU>(), RowWorld);
	MPI_Alltoallv(sendnumbuf, sendcnt, sdispls, MPIType<T_promote>(), recvnumbuf, recvcnt, rdispls, MPIType<T_promote>(), RowWorld);
	
	DeleteAll(sendindbuf, sendnumbuf);
	DeleteAll(sendcnt, recvcnt, sdispls, rdispls);
		
	// define a SPA-like data structure
	IU ysize = y.MyLocLength();
	T_promote * localy = new T_promote[ysize];
	bool * isthere = new bool[ysize];
	std::vector<IU> nzinds;	// nonzero indices		
  std::fill_n(isthere, ysize, false);
	
	for(int i=0; i< totrecv; ++i)
	{
		if(!isthere[recvindbuf[i]])
		{
			localy[recvindbuf[i]] = recvnumbuf[i];	// initial assignment
			nzinds.push_back(recvindbuf[i]);
			isthere[recvindbuf[i]] = true;
		} 
		else
		{
			localy[recvindbuf[i]] = SR::add(localy[recvindbuf[i]], recvnumbuf[i]);	
		}
	}
	DeleteAll(isthere, recvindbuf, recvnumbuf);
	sort(nzinds.begin(), nzinds.end());
	int nnzy = nzinds.size();
	y.ind.resize(nnzy);
	y.num.resize(nnzy);
	for(int i=0; i< nnzy; ++i)
	{
		y.ind[i] = nzinds[i];
		y.num[i] = localy[nzinds[i]]; 	
	}
	delete [] localy;
	return y;
}
	

template <typename IU, typename NU1, typename NU2, typename UDERA, typename UDERB> 
SpParMat<IU,typename promote_trait<NU1,NU2>::T_promote,typename promote_trait<UDERA,UDERB>::T_promote> EWiseMult 
	(const SpParMat<IU,NU1,UDERA> & A, const SpParMat<IU,NU2,UDERB> & B , bool exclude)
{
	typedef typename promote_trait<NU1,NU2>::T_promote N_promote;
	typedef typename promote_trait<UDERA,UDERB>::T_promote DER_promote;

	if(*(A.commGrid) == *(B.commGrid))	
	{
		DER_promote * result = new DER_promote( EWiseMult(*(A.spSeq),*(B.spSeq),exclude) );
		return SpParMat<IU, N_promote, DER_promote> (result, A.commGrid);
	}
	else
	{
		std::cout << "Grids are not comparable elementwise multiplication" << std::endl; 
		MPI_Abort(MPI_COMM_WORLD, GRIDMISMATCH);
		return SpParMat< IU,N_promote,DER_promote >();
	}
}
	
template <typename RETT, typename RETDER, typename IU, typename NU1, typename NU2, typename UDERA, typename UDERB, typename _BinaryOperation> 
SpParMat<IU,RETT,RETDER> EWiseApply 
	(const SpParMat<IU,NU1,UDERA> & A, const SpParMat<IU,NU2,UDERB> & B, _BinaryOperation __binary_op, bool notB, const NU2& defaultBVal)
{
	if(*(A.commGrid) == *(B.commGrid))	
	{
		RETDER * result = new RETDER( EWiseApply<RETT>(*(A.spSeq),*(B.spSeq), __binary_op, notB, defaultBVal) );
		return SpParMat<IU, RETT, RETDER> (result, A.commGrid);
	}
	else
	{
		std::cout << "Grids are not comparable elementwise apply" << std::endl; 
		MPI_Abort(MPI_COMM_WORLD, GRIDMISMATCH);
		return SpParMat< IU,RETT,RETDER >();
	}
}

template <typename RETT, typename RETDER, typename IU, typename NU1, typename NU2, typename UDERA, typename UDERB, typename _BinaryOperation, typename _BinaryPredicate> 
SpParMat<IU,RETT,RETDER> EWiseApply
	(const SpParMat<IU,NU1,UDERA> & A, const SpParMat<IU,NU2,UDERB> & B, _BinaryOperation __binary_op, _BinaryPredicate do_op, bool allowANulls, bool allowBNulls, const NU1& ANullVal, const NU2& BNullVal, const bool allowIntersect, const bool useExtendedBinOp)
{
	if(*(A.commGrid) == *(B.commGrid))	
	{
		RETDER * result = new RETDER( EWiseApply<RETT>(*(A.spSeq),*(B.spSeq), __binary_op, do_op, allowANulls, allowBNulls, ANullVal, BNullVal, allowIntersect) );
		return SpParMat<IU, RETT, RETDER> (result, A.commGrid);
	}
	else
	{
		std::cout << "Grids are not comparable elementwise apply" << std::endl; 
		MPI_Abort(MPI_COMM_WORLD, GRIDMISMATCH);
		return SpParMat< IU,RETT,RETDER >();
	}
}

// plain adapter
template <typename RETT, typename RETDER, typename IU, typename NU1, typename NU2, typename UDERA, typename UDERB, typename _BinaryOperation, typename _BinaryPredicate> 
SpParMat<IU,RETT,RETDER>
EWiseApply (const SpParMat<IU,NU1,UDERA> & A, const SpParMat<IU,NU2,UDERB> & B, _BinaryOperation __binary_op, _BinaryPredicate do_op, bool allowANulls, bool allowBNulls, const NU1& ANullVal, const NU2& BNullVal, const bool allowIntersect = true)
{
	return EWiseApply<RETT, RETDER>(A, B,
				EWiseExtToPlainAdapter<RETT, NU1, NU2, _BinaryOperation>(__binary_op),
				EWiseExtToPlainAdapter<bool, NU1, NU2, _BinaryPredicate>(do_op),
				allowANulls, allowBNulls, ANullVal, BNullVal, allowIntersect, true);
}
// end adapter

/**
 * if exclude is true, then we prune all entries W[i] != zero from V 
 * if exclude is false, then we perform a proper elementwise multiplication
**/
template <typename IU, typename NU1, typename NU2>
FullyDistSpVec<IU,typename promote_trait<NU1,NU2>::T_promote> EWiseMult 
	(const FullyDistSpVec<IU,NU1> & V, const FullyDistVec<IU,NU2> & W , bool exclude, NU2 zero)
{
	typedef typename promote_trait<NU1,NU2>::T_promote T_promote;

	if(*(V.commGrid) == *(W.commGrid))	
	{
		FullyDistSpVec< IU, T_promote> Product(V.commGrid);
		if(V.glen != W.glen)
		{
			std::cerr << "Vector dimensions don't match for EWiseMult\n";
			MPI_Abort(MPI_COMM_WORLD, DIMMISMATCH);
		}
		else
		{
			Product.glen = V.glen;
			IU size= V.getlocnnz();
			if(exclude)
			{
				#if defined(_OPENMP) && defined(CBLAS_EXPERIMENTAL)	// not faster than serial
				int actual_splits = cblas_splits * 1;	// 1 is the parallel slackness
        std::vector <IU> tlosizes (actual_splits, 0);
        std::vector < std::vector<IU> > tlinds(actual_splits);
        std::vector < std::vector<T_promote> > tlnums(actual_splits);
				IU tlsize = size / actual_splits;
				#pragma omp parallel for //schedule(dynamic, 1)
				for(IU t = 0; t < actual_splits; ++t)
				{
					IU tlbegin = t*tlsize;
					IU tlend = (t==actual_splits-1)? size : (t+1)*tlsize;
					for(IU i=tlbegin; i<tlend; ++i)
					{
						if(W.arr[V.ind[i]] == zero) 	// keep only those
						{
							tlinds[t].push_back(V.ind[i]);
							tlnums[t].push_back(V.num[i]);
							tlosizes[t]++;
						}
					}
				}
        std::vector<IU> prefix_sum(actual_splits+1,0);
        std::partial_sum(tlosizes.begin(), tlosizes.end(), prefix_sum.begin()+1); 
				Product.ind.resize(prefix_sum[actual_splits]);
				Product.num.resize(prefix_sum[actual_splits]);
			
				#pragma omp parallel for //schedule(dynamic, 1)
				for(IU t=0; t< actual_splits; ++t)
				{
          std::copy(tlinds[t].begin(), tlinds[t].end(), Product.ind.begin()+prefix_sum[t]);
          std::copy(tlnums[t].begin(), tlnums[t].end(), Product.num.begin()+prefix_sum[t]);
				}
				#else
				for(IU i=0; i<size; ++i)
				{
					if(W.arr[V.ind[i]] == zero)     // keep only those
					{
                        	       		Product.ind.push_back(V.ind[i]);
                                		Product.num.push_back(V.num[i]);
                                      	}	
				}
				#endif
			}
			else
			{
				for(IU i=0; i<size; ++i)
				{
					if(W.arr[V.ind[i]] != zero) 	// keep only those
					{
						Product.ind.push_back(V.ind[i]);
						Product.num.push_back(V.num[i] * W.arr[V.ind[i]]);
					}
				}
			}
		}
		return Product;
	}
	else
	{
		std::cout << "Grids are not comparable elementwise multiplication" << std::endl; 
		MPI_Abort(MPI_COMM_WORLD, GRIDMISMATCH);
		return FullyDistSpVec< IU,T_promote>();
	}
}


/**
 Threaded EWiseApply. Only called internally from EWiseApply.
**/
template <typename RET, typename IU, typename NU1, typename NU2, typename _BinaryOperation, typename _BinaryPredicate>
FullyDistSpVec<IU,RET> EWiseApply_threaded
	(const FullyDistSpVec<IU,NU1> & V, const FullyDistVec<IU,NU2> & W , _BinaryOperation _binary_op, _BinaryPredicate _doOp, bool allowVNulls, NU1 Vzero, const bool useExtendedBinOp)
{
	typedef RET T_promote; //typedef typename promote_trait<NU1,NU2>::T_promote T_promote;
	if(*(V.commGrid) == *(W.commGrid))	
	{
		FullyDistSpVec< IU, T_promote> Product(V.commGrid);
		if(V.TotalLength() != W.TotalLength())
		{
			std::ostringstream outs;
			outs << "Vector dimensions don't match (" << V.TotalLength() << " vs " << W.TotalLength() << ") for EWiseApply (short version)\n";
			SpParHelper::Print(outs.str());
			MPI_Abort(MPI_COMM_WORLD, DIMMISMATCH);
		}
		else
		{
            int nthreads=1;
#ifdef _OPENMP
#pragma omp parallel
            {
                nthreads = omp_get_num_threads();
            }
#endif

			Product.glen = V.glen;
			IU size= W.LocArrSize();
			IU spsize = V.getlocnnz();
            
            // temporary result vectors per thread
            std::vector<std::vector<IU>> tProductInd(nthreads);
            std::vector<std::vector<T_promote>> tProductVal(nthreads);
            IU perthread; //chunk of tProductInd or tProductVal allocated to each thread
            if (allowVNulls)
                perthread = size/nthreads;
            else
                perthread = spsize/nthreads;
            
#ifdef _OPENMP
#pragma omp parallel
#endif
            {
                int curthread = 0;
#ifdef _OPENMP
                curthread = omp_get_thread_num();
#endif
                IU tStartIdx = perthread * curthread;
                IU tNextIdx = perthread * (curthread+1);
                
                if (allowVNulls)
                {
                    if(curthread == nthreads-1) tNextIdx = size;
                    
                    // get sparse part for the current thread
                    auto it = std::lower_bound (V.ind.begin(), V.ind.end(), tStartIdx);
                    IU tSpIdx = (IU) std::distance(V.ind.begin(), it);
                    
                    // iterate over the dense vector
                    for(IU tIdx=tStartIdx; tIdx < tNextIdx; ++tIdx)
                    {
                        if(tSpIdx < spsize && V.ind[tSpIdx] < tNextIdx && V.ind[tSpIdx] == tIdx)
                        {
                            if (_doOp(V.num[tSpIdx], W.arr[tIdx], false, false))
                            {
                                tProductInd[curthread].push_back(tIdx);
                                tProductVal[curthread].push_back (_binary_op(V.num[tSpIdx], W.arr[tIdx], false, false));
                            }
                            tSpIdx++;
                        }
                        else
                        {
                            if (_doOp(Vzero, W.arr[tIdx], true, false))
                            {
                                tProductInd[curthread].push_back(tIdx);
                                tProductVal[curthread].push_back (_binary_op(Vzero, W.arr[tIdx], true, false));
                            }
                        }
                    }
                }
                else // iterate over the sparse vector
                {
                    if(curthread == nthreads-1) tNextIdx = spsize;
                    for(IU tSpIdx=tStartIdx; tSpIdx < tNextIdx; ++tSpIdx)
                    {
                        if (_doOp(V.num[tSpIdx], W.arr[V.ind[tSpIdx]], false, false))
                        {
                            
                            tProductInd[curthread].push_back( V.ind[tSpIdx]);
                            tProductVal[curthread].push_back (_binary_op(V.num[tSpIdx], W.arr[V.ind[tSpIdx]], false, false));
                        }
                    }
                }
            }
            
            std::vector<IU> tdisp(nthreads+1);
            tdisp[0] = 0;
            for(int i=0; i<nthreads; ++i)
            {
                tdisp[i+1] = tdisp[i] + tProductInd[i].size();
            }
            
            // copy results from temporary vectors
            Product.ind.resize(tdisp[nthreads]);
            Product.num.resize(tdisp[nthreads]);
            
#ifdef _OPENMP
#pragma omp parallel
#endif
            {
                int curthread = 0;
#ifdef _OPENMP
                curthread = omp_get_thread_num();
#endif
                std::copy(tProductInd[curthread].begin(), tProductInd[curthread].end(), Product.ind.data() + tdisp[curthread]);
                std::copy(tProductVal[curthread].begin() , tProductVal[curthread].end(), Product.num.data() + tdisp[curthread]);
            }
		}
		return Product;
	}
	else
	{
		std::cout << "Grids are not comparable for EWiseApply" << std::endl;
		MPI_Abort(MPI_COMM_WORLD, GRIDMISMATCH);
		return FullyDistSpVec< IU,T_promote>();
	}
}



/**
 * Performs an arbitrary binary operation _binary_op on the corresponding elements of two vectors with the result stored in a return vector ret.
 * The binary operatiation is only performed if the binary predicate _doOp returns true for those elements. Otherwise the binary operation is not
 * performed and ret does not contain an element at that position.
 * More formally the operation is defined as:
 * if (_doOp(V[i], W[i]))
 *    ret[i] = _binary_op(V[i], W[i])
 * else
 *    // ret[i] is not set
 * Hence _doOp can be used to implement a filter on either of the vectors.
 *
 * The above is only defined if both V[i] and W[i] exist (i.e. an intersection). To allow a union operation (ex. when V[i] doesn't exist but W[i] does)
 * the allowVNulls flag is set to true and the Vzero argument is used as the missing V[i] value.
 *
 * The type of each element of ret must not necessarily be related to the types of V or W, so the return type must be explicitly specified as a template parameter:
 * FullyDistSpVec<int, double> r = EWiseApply<double>(V, W, plus, retTrue, false, 0)
 **/
template <typename RET, typename IU, typename NU1, typename NU2, typename _BinaryOperation, typename _BinaryPredicate>
FullyDistSpVec<IU,RET> EWiseApply
(const FullyDistSpVec<IU,NU1> & V, const FullyDistVec<IU,NU2> & W , _BinaryOperation _binary_op, _BinaryPredicate _doOp, bool allowVNulls, NU1 Vzero, const bool useExtendedBinOp)
{
    
#ifdef _OPENMP
    return EWiseApply_threaded<RET>(V, W, _binary_op, _doOp, allowVNulls, Vzero, useExtendedBinOp);
    
#else
    typedef RET T_promote; //typedef typename promote_trait<NU1,NU2>::T_promote T_promote;
    if(*(V.commGrid) == *(W.commGrid))
    {
        FullyDistSpVec< IU, T_promote> Product(V.commGrid);
        //FullyDistVec< IU, NU1> DV (V); // Ariful: I am not sure why it was there??
        if(V.TotalLength() != W.TotalLength())
        {
          std::ostringstream outs;
            outs << "Vector dimensions don't match (" << V.TotalLength() << " vs " << W.TotalLength() << ") for EWiseApply (short version)\n";
            SpParHelper::Print(outs.str());
            MPI_Abort(MPI_COMM_WORLD, DIMMISMATCH);
        }
        else
        {
            Product.glen = V.glen;
            IU size= W.LocArrSize();
            IU spsize = V.getlocnnz();
            IU sp_iter = 0;
            if (allowVNulls)
            {
                // iterate over the dense vector
                for(IU i=0; i<size; ++i)
                {
                    if(sp_iter < spsize && V.ind[sp_iter] == i)
                    {
                        if (_doOp(V.num[sp_iter], W.arr[i], false, false))
                        {
                            Product.ind.push_back(i);
                            Product.num.push_back(_binary_op(V.num[sp_iter], W.arr[i], false, false));
                        }
                        sp_iter++;
                    }
                    else
                    {
                        if (_doOp(Vzero, W.arr[i], true, false))
                        {
                            Product.ind.push_back(i);
                            Product.num.push_back(_binary_op(Vzero, W.arr[i], true, false));
                        }
                    }
                }
            }
            else
            {
                // iterate over the sparse vector
                for(sp_iter = 0; sp_iter < spsize; ++sp_iter)
                {
                    if (_doOp(V.num[sp_iter], W.arr[V.ind[sp_iter]], false, false))
                    {
                        Product.ind.push_back(V.ind[sp_iter]);
                        Product.num.push_back(_binary_op(V.num[sp_iter], W.arr[V.ind[sp_iter]], false, false));
                    }
                }
                
            }
        }
        return Product;
    }
    else
    {
      std::cout << "Grids are not comparable for EWiseApply" << std::endl; 
        MPI_Abort(MPI_COMM_WORLD, GRIDMISMATCH);
        return FullyDistSpVec< IU,T_promote>();
    }
#endif
}



/**
 * Performs an arbitrary binary operation _binary_op on the corresponding elements of two vectors with the result stored in a return vector ret. 
 * The binary operatiation is only performed if the binary predicate _doOp returns true for those elements. Otherwise the binary operation is not 
 * performed and ret does not contain an element at that position.
 * More formally the operation is defined as:
 * if (_doOp(V[i], W[i]))
 *    ret[i] = _binary_op(V[i], W[i])
 * else
 *    // ret[i] is not set
 * Hence _doOp can be used to implement a filter on either of the vectors.
 *
 * The above is only defined if both V[i] and W[i] exist (i.e. an intersection). To allow a union operation (ex. when V[i] doesn't exist but W[i] does) 
 * the allowVNulls flag is set to true and the Vzero argument is used as the missing V[i] value.
 * !allowVNulls && !allowWNulls => intersection
 * !allowVNulls &&  allowWNulls => operate on all elements of V 
 *  allowVNulls && !allowWNulls => operate on all elements of W
 *  allowVNulls &&  allowWNulls => union
 *
 * The type of each element of ret must not necessarily be related to the types of V or W, so the return type must be explicitly specified as a template parameter:
 * FullyDistSpVec<int, double> r = EWiseApply<double>(V, W, plus, ...)
 * For intersection, Vzero and Wzero are irrelevant
 * ABAB: \todo: Should allowIntersect be "false" for all SetDifference uses?
**/
template <typename RET, typename IU, typename NU1, typename NU2, typename _BinaryOperation, typename _BinaryPredicate>
FullyDistSpVec<IU,RET> EWiseApply 
	(const FullyDistSpVec<IU,NU1> & V, const FullyDistSpVec<IU,NU2> & W , _BinaryOperation _binary_op, _BinaryPredicate _doOp, bool allowVNulls, bool allowWNulls, NU1 Vzero, NU2 Wzero, const bool allowIntersect, const bool useExtendedBinOp)
{

	typedef RET T_promote; // typename promote_trait<NU1,NU2>::T_promote T_promote;
	if(*(V.commGrid) == *(W.commGrid))	
	{
		FullyDistSpVec< IU, T_promote> Product(V.commGrid);
		if(V.glen != W.glen)
		{
			std::ostringstream outs;
			outs << "Vector dimensions don't match (" << V.glen << " vs " << W.glen << ") for EWiseApply (full version)\n";
			SpParHelper::Print(outs.str());
			MPI_Abort(MPI_COMM_WORLD, DIMMISMATCH);
		}
		else
		{
			Product.glen = V.glen;
			typename std::vector< IU  >::const_iterator indV = V.ind.begin();
			typename std::vector< NU1 >::const_iterator numV = V.num.begin();
			typename std::vector< IU  >::const_iterator indW = W.ind.begin();
			typename std::vector< NU2 >::const_iterator numW = W.num.begin();
			
			while (indV < V.ind.end() && indW < W.ind.end())
			{
				if (*indV == *indW)
				{
					// overlap
					if (allowIntersect)
					{
						if (_doOp(*numV, *numW, false, false))
						{
							Product.ind.push_back(*indV);
							Product.num.push_back(_binary_op(*numV, *numW, false, false));
						}
					}
					indV++; numV++;
					indW++; numW++;
				}
				else if (*indV < *indW)
				{
					// V has value but W does not
					if (allowWNulls)
					{
						if (_doOp(*numV, Wzero, false, true))
						{
							Product.ind.push_back(*indV);
							Product.num.push_back(_binary_op(*numV, Wzero, false, true));
						}
					}
					indV++; numV++;
				}
				else //(*indV > *indW)
				{
					// W has value but V does not
					if (allowVNulls)
					{
						if (_doOp(Vzero, *numW, true, false))
						{
							Product.ind.push_back(*indW);
							Product.num.push_back(_binary_op(Vzero, *numW, true, false));
						}
					}
					indW++; numW++;
				}
			}
			// clean up
			while (allowWNulls && indV < V.ind.end())
			{
				if (_doOp(*numV, Wzero, false, true))
				{
					Product.ind.push_back(*indV);
					Product.num.push_back(_binary_op(*numV, Wzero, false, true));
				}
				indV++; numV++;
			}
			while (allowVNulls && indW < W.ind.end())
			{
				if (_doOp(Vzero, *numW, true, false))
				{
					Product.ind.push_back(*indW);
					Product.num.push_back(_binary_op(Vzero, *numW, true, false));
				}
				indW++; numW++;
			}
		}
		return Product;
	}
	else
	{
		std::cout << "Grids are not comparable for EWiseApply" << std::endl; 
		MPI_Abort(MPI_COMM_WORLD, GRIDMISMATCH);
		return FullyDistSpVec< IU,T_promote>();
	}
}

// plain callback versions
template <typename RET, typename IU, typename NU1, typename NU2, typename _BinaryOperation, typename _BinaryPredicate>
FullyDistSpVec<IU,RET> EWiseApply 
	(const FullyDistSpVec<IU,NU1> & V, const FullyDistVec<IU,NU2> & W , _BinaryOperation _binary_op, _BinaryPredicate _doOp, bool allowVNulls, NU1 Vzero)
{


	return EWiseApply<RET>(V, W,
					EWiseExtToPlainAdapter<RET, NU1, NU2, _BinaryOperation>(_binary_op),
					EWiseExtToPlainAdapter<bool, NU1, NU2, _BinaryPredicate>(_doOp),
					allowVNulls, Vzero, true);
}



template <typename RET, typename IU, typename NU1, typename NU2, typename _BinaryOperation, typename _BinaryPredicate>
FullyDistSpVec<IU,RET> EWiseApply 
	(const FullyDistSpVec<IU,NU1> & V, const FullyDistSpVec<IU,NU2> & W , _BinaryOperation _binary_op, _BinaryPredicate _doOp, bool allowVNulls, bool allowWNulls, NU1 Vzero, NU2 Wzero, const bool allowIntersect = true)
{
	return EWiseApply<RET>(V, W,
					EWiseExtToPlainAdapter<RET, NU1, NU2, _BinaryOperation>(_binary_op),
					EWiseExtToPlainAdapter<bool, NU1, NU2, _BinaryPredicate>(_doOp),
					allowVNulls, allowWNulls, Vzero, Wzero, allowIntersect, true);
}


////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// SpMM /////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
inline void
print_spmm_stats (spmm_stats &stats, int nruns)
{
	stats.t_g_h2d_memcpy	   /= nruns;
	stats.t_g_spmm_buf		   /= nruns;
	stats.t_g_spmm			   /= nruns;
	stats.t_g_d2h_memcpy	   /= nruns;
	stats.t_g_mtx_create	   /= nruns;
	stats.t_comp			   /= nruns;
	stats.t_tot				   /= nruns;
	stats.t_sA_comm_pre		   /= nruns;
	stats.t_sA_comm_post	   /= nruns;
	stats.t_sC_comm_bcastA	   /= nruns;
	stats.t_sC_comm_bcastX	   /= nruns;
	stats.t_sA_2D_comm_bcastX  /= nruns;
	stats.t_sA_2D_comm_reduceY /= nruns;

	float	t_g_h2d_memcpy[3];
	float	t_g_spmm_buf[3];
	float	t_g_spmm[3];
	float	t_g_d2h_memcpy[3];
	float	t_g_mtx_create[3];
	double	t_comp[3];
	double	t_tot[3];
	double	t_sA_comm_pre[3];
	double	t_sA_comm_post[3];	
	double	t_sC_comm_bcastA[3];
	double	t_sC_comm_bcastX[3];
	double	t_sA_2D_comm_bcastX[3];
	double	t_sA_2D_comm_reduceY[3];	

	std::string tmp;
	int np, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	MPI_Reduce(&stats.t_g_h2d_memcpy, &t_g_h2d_memcpy[0], 1,
			   MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&stats.t_g_h2d_memcpy, &t_g_h2d_memcpy[1], 1,
			   MPI_FLOAT, MPI_MIN, 0, MPI_COMM_WORLD);
	MPI_Reduce(&stats.t_g_h2d_memcpy, &t_g_h2d_memcpy[2], 1,
			   MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
	tmp += "[G] h2d-copy       " +
		std::to_string(t_g_h2d_memcpy[1]/1e3) + " " +
		std::to_string(t_g_h2d_memcpy[0]/1e3/np) + " " +
		std::to_string(t_g_h2d_memcpy[2]/1e3) + "\n";

	MPI_Reduce(&stats.t_g_spmm, &t_g_spmm[0], 1,
			   MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&stats.t_g_spmm, &t_g_spmm[1], 1,
			   MPI_FLOAT, MPI_MIN, 0, MPI_COMM_WORLD);
	MPI_Reduce(&stats.t_g_spmm, &t_g_spmm[2], 1,
			   MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
	tmp += "[G] spmm           " +
		std::to_string(t_g_spmm[1]/1e3) + " " +
		std::to_string(t_g_spmm[0]/1e3/np) + " " +
		std::to_string(t_g_spmm[2]/1e3) + "\n";
	
	MPI_Reduce(&stats.t_g_d2h_memcpy, &t_g_d2h_memcpy[0], 1,
			   MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&stats.t_g_d2h_memcpy, &t_g_d2h_memcpy[1], 1,
			   MPI_FLOAT, MPI_MIN, 0, MPI_COMM_WORLD);
	MPI_Reduce(&stats.t_g_d2h_memcpy, &t_g_d2h_memcpy[2], 1,
			   MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
	tmp += "[G] d2h-copy       " +
		std::to_string(t_g_d2h_memcpy[1]/1e3) + " " +
		std::to_string(t_g_d2h_memcpy[0]/1e3/np) + " " +
		std::to_string(t_g_d2h_memcpy[2]/1e3) + "\n";

	MPI_Reduce(&stats.t_g_mtx_create, &t_g_mtx_create[0], 1,
			   MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&stats.t_g_mtx_create, &t_g_mtx_create[1], 1,
			   MPI_FLOAT, MPI_MIN, 0, MPI_COMM_WORLD);
	MPI_Reduce(&stats.t_g_mtx_create, &t_g_mtx_create[2], 1,
			   MPI_FLOAT, MPI_MAX, 0, MPI_COMM_WORLD);
	tmp += "[G] mtx-create     " +
		std::to_string(t_g_mtx_create[1]/1e3) + " " +
		std::to_string(t_g_mtx_create[0]/1e3/np) + " " +
		std::to_string(t_g_mtx_create[2]/1e3) + "\n";

	MPI_Reduce(&stats.t_sA_comm_pre, &t_sA_comm_pre[0], 1,
			   MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&stats.t_sA_comm_pre, &t_sA_comm_pre[1], 1,
			   MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
	MPI_Reduce(&stats.t_sA_comm_pre, &t_sA_comm_pre[2], 1,
			   MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	tmp += "sA-comm-pre        " +
		std::to_string(t_sA_comm_pre[1]) + " " +
		std::to_string(t_sA_comm_pre[0]/np) + " " +
		std::to_string(t_sA_comm_pre[2]) + "\n";

	MPI_Reduce(&stats.t_sA_comm_post, &t_sA_comm_post[0], 1,
			   MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&stats.t_sA_comm_post, &t_sA_comm_post[1], 1,
			   MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
	MPI_Reduce(&stats.t_sA_comm_post, &t_sA_comm_post[2], 1,
			   MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	tmp += "sA-comm-post       " +
		std::to_string(t_sA_comm_post[1]) + " " +
		std::to_string(t_sA_comm_post[0]/np) + " " +
		std::to_string(t_sA_comm_post[2]) + "\n";
	
	MPI_Reduce(&stats.t_sC_comm_bcastA, &t_sC_comm_bcastA[0], 1,
			   MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&stats.t_sC_comm_bcastA, &t_sC_comm_bcastA[1], 1,
			   MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
	MPI_Reduce(&stats.t_sC_comm_bcastA, &t_sC_comm_bcastA[2], 1,
			   MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	tmp += "sC-comm-bcastS     " +
		std::to_string(t_sC_comm_bcastA[1]) + " " +
		std::to_string(t_sC_comm_bcastA[0]/np) + " " +
		std::to_string(t_sC_comm_bcastA[2]) + "\n";

	MPI_Reduce(&stats.t_sC_comm_bcastX, &t_sC_comm_bcastX[0], 1,
			   MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&stats.t_sC_comm_bcastX, &t_sC_comm_bcastX[1], 1,
			   MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
	MPI_Reduce(&stats.t_sC_comm_bcastX, &t_sC_comm_bcastX[2], 1,
			   MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	tmp += "sC-comm-bcastD     " +
		std::to_string(t_sC_comm_bcastX[1]) + " " +
		std::to_string(t_sC_comm_bcastX[0]/np) + " " +
		std::to_string(t_sC_comm_bcastX[2]) + "\n";

	MPI_Reduce(&stats.t_sA_2D_comm_bcastX, &t_sA_2D_comm_bcastX[0], 1,
			   MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&stats.t_sA_2D_comm_bcastX, &t_sA_2D_comm_bcastX[1], 1,
			   MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
	MPI_Reduce(&stats.t_sA_2D_comm_bcastX, &t_sA_2D_comm_bcastX[2], 1,
			   MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	tmp += "sA-2D-comm-bcastX  " +
		std::to_string(t_sA_2D_comm_bcastX[1]) + " " +
		std::to_string(t_sA_2D_comm_bcastX[0]/np) + " " +
		std::to_string(t_sA_2D_comm_bcastX[2]) + "\n";

	MPI_Reduce(&stats.t_sA_2D_comm_reduceY, &t_sA_2D_comm_reduceY[0], 1,
			   MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&stats.t_sA_2D_comm_reduceY, &t_sA_2D_comm_reduceY[1], 1,
			   MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
	MPI_Reduce(&stats.t_sA_2D_comm_reduceY, &t_sA_2D_comm_reduceY[2], 1,
			   MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	tmp += "sA-2D-comm-reduceY " +
		std::to_string(t_sA_2D_comm_reduceY[1]) + " " +
		std::to_string(t_sA_2D_comm_reduceY[0]/np) + " " +
		std::to_string(t_sA_2D_comm_reduceY[2]) + "\n";

	MPI_Reduce(&stats.t_comp, &t_comp[0], 1,
			   MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&stats.t_comp, &t_comp[1], 1,
			   MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
	MPI_Reduce(&stats.t_comp, &t_comp[2], 1,
			   MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	tmp += "comp               " +
		std::to_string(t_comp[1]) + " " +
		std::to_string(t_comp[0]/np) + " " +
		std::to_string(t_comp[2]) + "\n";

	MPI_Reduce(&stats.t_tot, &t_tot[0], 1,
			   MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&stats.t_tot, &t_tot[1], 1,
			   MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
	MPI_Reduce(&stats.t_tot, &t_tot[2], 1,
			   MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	tmp += "tot                " +
		std::to_string(t_tot[1]) + " " +
		std::to_string(t_tot[0]/np) + " " +
		std::to_string(t_tot[2]) + "\n";

	
	if (rank == 0)
	{
		std::cout << "\n";
		std::cout << std::string(80, '*') << "\n";
		std::cout << tmp;
		std::cout << std::string(80, '*') << "\n";
		std::cout << std::endl;
	}	
}


	
template <typename SR,
		  typename IU,
		  typename NUM,
		  typename NUV,
		  typename UDER> 
FullyDistDMat<IU, typename promote_trait<NUM, NUV>::T_promote> 
SpMM_sA
(
    const SpParMat<IU, NUM, UDER>	&A,
	const FullyDistDMat<IU, NUV>	&X,
	spmm_stats						&stats
)
{
	auto t_tot_beg = std::chrono::high_resolution_clock::now();
	
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	typedef typename promote_trait<NUM, NUV>::T_promote T_promote;


	auto t_beg = std::chrono::high_resolution_clock::now();	

	MPI_Comm	world  = X.commGrid->GetWorld();
	MPI_Comm	cworld = X.commGrid->GetColWorld();
	MPI_Comm	rworld = X.commGrid->GetRowWorld();

	int			xsize	  = (int) X.LocArrSize();
	int			trxsize	  = 0;
	int			diagneigh = X.commGrid->GetComplementRank();
	MPI_Status	status;
	MPI_Sendrecv(&xsize, 1, MPI_INT, diagneigh, TRX,
				 &trxsize, 1, MPI_INT, diagneigh, TRX,
				 world, &status);

	
	NUV *trxnums = new NUV[trxsize];
	MPI_Sendrecv(const_cast<NUV *>(SpHelper::p2a(X.arr_)),
				 xsize, MPIType<NUV>(),
				 diagneigh, TRX,
				 trxnums, trxsize, MPIType<NUV>(),
				 diagneigh, TRX, world, &status);

	int cneighs, crank;
	MPI_Comm_size(cworld, &cneighs);
	MPI_Comm_rank(cworld, &crank);
	int *csize = new int[cneighs];
	csize[crank] = trxsize;
	MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, csize, 1, MPI_INT, cworld);
	
	int *dpls = new int[cneighs]();
	std::partial_sum(csize, csize+cneighs-1, dpls+1);
	int accsize = std::accumulate(csize, csize+cneighs, 0);
	NUV *numacc = new NUV[accsize];
	
	MPI_Allgatherv(trxnums, trxsize, MPIType<NUV>(),
				   numacc, csize, dpls, MPIType<NUV>(), cworld);
	delete [] trxnums;

	auto t_end = std::chrono::high_resolution_clock::now();
	stats.t_sA_comm_pre += static_cast<std::chrono::duration<double> >
		(t_end-t_beg).count();
	
	// local SpMM
	T_promote	 id	   = SR::id();
	IU			 ysize = (int64_t)A.getlocalrows() * X.getncol();
	T_promote	*y_loc = new T_promote[ysize];	
		
	std::fill_n(y_loc, ysize, id);


	t_beg = std::chrono::high_resolution_clock::now();
	
	csc_gespmm_cusparse<SR>(*(A.spSeq), numacc, y_loc, (int)X.getncol(),
							static_cast<NUM>(0.0), stats);
	
	t_end = std::chrono::high_resolution_clock::now();
	stats.t_comp += static_cast<std::chrono::duration<double> >
		(t_end-t_beg).count();


	DeleteAll(numacc, csize, dpls);
		
	t_beg = std::chrono::high_resolution_clock::now();
	

	// reduce output - MPI_Reduce_scatter version
	FullyDistDMat<IU, T_promote> Y(X.commGrid, A.getnrow(), X.getncol(), id);
	int rneighs;
	MPI_Comm_size(rworld, &rneighs);
	
	int *recvcounts = new int[rneighs];
	for (int i = 0; i < rneighs; ++i)
	{
		IU tmp;
		if (i == rneighs-1)
			tmp = ysize;
		else
			tmp = Y.RowLenUntil(i+1) * Y.getncol();
		
		recvcounts[i] = tmp - (Y.RowLenUntil(i) * Y.getncol());
	}
	
	MPI_Reduce_scatter(y_loc, SpHelper::p2a(Y.arr_), recvcounts,
					   MPIType<T_promote>(), SR::mpi_op(), rworld);
	
	t_end = std::chrono::high_resolution_clock::now();
	stats.t_sA_comm_post += static_cast<std::chrono::duration<double> >
		(t_end-t_beg).count();


	auto t_tot_end = std::chrono::high_resolution_clock::now();
	stats.t_tot += static_cast<std::chrono::duration<double> >
		(t_tot_end-t_tot_beg).count();

	delete [] y_loc;
	delete [] recvcounts;
	
	return Y;
}



template <typename SR,
		  typename IU,
		  typename NUM,
		  typename NUV,
		  typename UDER> 
DnParMat<IU, typename promote_trait<NUM, NUV>::T_promote> 
SpMM_sA_2D
(
    const SpParMat<IU, NUM, UDER>	&A,
	const DnParMat<IU, NUV>			&X,
	spmm_stats						&stats
)
{
	auto t_tot_beg = std::chrono::high_resolution_clock::now();

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	typedef typename promote_trait<NUM, NUV>::T_promote NT_Y;
	typedef typename UDER::LocalIT LIA;

	int nstages, dummy;
	std::shared_ptr<CommGrid> cgr_Y =
		ProductGrid(A.getcommgrid().get(), (X.cgr).get(), nstages,
					dummy, dummy);
	LIA		Y_m = A.getnrow();
	IU		Y_n = X.getgncol();
	NT_Y	id	= SR::id();
	DnParMat<IU, NT_Y> Y(cgr_Y, Y_m, Y_n, id);

	MPI_Comm	world  = X.cgr->GetWorld();
	MPI_Comm	cworld = X.cgr->GetColWorld();
	MPI_Comm	rworld = X.cgr->GetRowWorld();

	// exchange local matrices along the diagonal
	int	xsize = (int)X.getnels();
	IU	trx_nrows, trx_ncols, trxsize;
	X.get_local_length(X.cgr->GetRankInProcRow(), X.cgr->GetRankInProcCol(),
					   trx_nrows, trx_ncols, trxsize);
		
	int			diagneigh = X.cgr->GetComplementRank();
	MPI_Status	status;
	NUV *trxnums = new NUV[trxsize];
	MPI_Sendrecv(const_cast<NUV *>(SpHelper::p2a(X.arr)),
				 xsize, MPIType<NUV>(),
				 diagneigh, TRX,
				 trxnums, trxsize, MPIType<NUV>(),
				 diagneigh, TRX, world, &status);

	IU *XRecvSizes = new IU[nstages];
	IU	max_xrecv_size = 0;
	for (int s = 0; s < nstages; ++s)
	{
		IU tmp1, tmp2;
		X.get_local_length(X.cgr->GetRankInProcRow(), s,
						   tmp1, tmp2, XRecvSizes[s]);
		max_xrecv_size = max(max_xrecv_size, XRecvSizes[s]);
	}

	IU *YSendSizes = new IU[nstages];
	IU	max_ysend_size = 0;
	for (int s = 0; s < nstages; ++s)
	{
		IU tmp1, tmp2;
		Y.get_local_length(Y.cgr->GetRankInProcCol(), s,
						   tmp1, tmp2, YSendSizes[s]);
		max_ysend_size = max(max_ysend_size, YSendSizes[s]);
	}
	

	NUV	*XRecv = new NUV[max_xrecv_size];
	NUV	*YSend = new NUV[max_ysend_size];
	for (int s = 0; s < nstages; ++s)
	{
		auto t_beg = std::chrono::high_resolution_clock::now();
		
		// broadcast X
		NUV *XRecv_cur = XRecv;
		if (s == X.cgr->GetRankInProcCol())
			XRecv_cur = trxnums;

		MPI_Bcast(XRecv_cur, XRecvSizes[s], MPIType<NUV>(),
				  s, X.cgr->GetColWorld());

		auto t_end = std::chrono::high_resolution_clock::now();
		
		stats.t_sA_2D_comm_bcastX += static_cast<std::chrono::duration<double> >
			(t_end-t_beg).count();
		
		t_beg = std::chrono::high_resolution_clock::now();

		// init Y		
		NUV *YSend_cur = YSend;
		
		std::fill_n(YSend_cur, YSendSizes[s], id);

		csc_gespmm_cusparse<SR>(*(A.spSeq), XRecv_cur, YSend_cur,
								(int)X.getncol(), static_cast<NUM>(0.0), stats);

		t_end = std::chrono::high_resolution_clock::now();
		
		stats.t_comp += static_cast<std::chrono::duration<double> >
			(t_end-t_beg).count();
		
		t_beg = std::chrono::high_resolution_clock::now();

		// reduce
		MPI_Reduce(YSend_cur, SpHelper::p2a(Y.arr), YSendSizes[s],
				   MPIType<NT_Y>(), SR::mpi_op(), s, Y.cgr->GetRowWorld());

		t_end = std::chrono::high_resolution_clock::now();
		stats.t_sA_2D_comm_reduceY +=
			static_cast<std::chrono::duration<double> >(t_end-t_beg).count();
	}

	delete [] XRecv;
	delete [] YSend;

	auto t_tot_end = std::chrono::high_resolution_clock::now();
	stats.t_tot += static_cast<std::chrono::duration<double> >
		(t_tot_end-t_tot_beg).count();
	
	return Y;
}



template <typename SR,
		  typename IU,
		  typename NUM,
		  typename NUV,
		  typename UDER> 
DnParMat<IU, typename promote_trait<NUM, NUV>::T_promote> 
SpMM_sC
(
    const SpParMat<IU, NUM, UDER>	&A,
	const DnParMat<IU, NUV>			&X,
	spmm_stats						&stats
)
{
	auto t_tot_beg = std::chrono::high_resolution_clock::now();
	
	
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	std::string dbgfile = std::string("spmm-debug-P") + std::to_string(rank);
	std::ofstream ofs;
	// ofs.open(dbgfile.c_str());
	// ofs << "in SpMM" << std::endl;
	
	typedef typename promote_trait<NUM, NUV>::T_promote NT_Y;
	typedef typename UDER::LocalIT LIA;

	int nstages, dummy;
	std::shared_ptr<CommGrid> cgr_Y =
		ProductGrid((A.commGrid).get(), (X.cgr).get(), nstages, dummy, dummy);
	LIA		Y_m = A.getnrow();
	IU		Y_n = X.getgncol();
	NT_Y	id	= SR::id();
	DnParMat<IU, NT_Y> Y(cgr_Y, Y_m, Y_n, id);
	int Aself		   = (A.commGrid)->GetRankInProcRow();
    int Xself		   = X.cgr->GetRankInProcCol();
	IU	max_xrecv_size = 0;

	// ofs << "matrix Y sizes: " << Y_m << " " << Y_n << "\n";	

	LIA **ARecvSizes = SpHelper::allocate2D<LIA>(UDER::esscount, nstages);
	IU	 *XRecvSizes = new IU[nstages];
	
	XRecvSizes[Xself] = X.getnels();

	SpParHelper::GetSetSizes(*(A.spSeq), ARecvSizes,
							 (A.commGrid)->GetRowWorld());
	MPI_Allgather(MPI_IN_PLACE, 1, MPIType<IU>(),
				  XRecvSizes, 1, MPIType<IU>(), X.cgr->GetColWorld());

	for (int i = 0; i < nstages; ++i)
		max_xrecv_size = max(max_xrecv_size, XRecvSizes[i]);

	// ofs << "XRecvSizes:\n";
	// for (int i = 0; i < nstages; ++i)
	// 	ofs << XRecvSizes[i] << " ";
	// ofs << "\n";
	// ofs << "max X recv size: " << max_xrecv_size << "\n";

	UDER	*ARecv = NULL;
	NUV		*XRecv = new NUV[max_xrecv_size];

	for (int s = 0; s < nstages; ++s)
	{
		auto t_beg = std::chrono::high_resolution_clock::now();
		
		
		// broadcast A
		std::vector<LIA> ess;
		if (s == Aself)
			ARecv = A.spSeq;
		else
		{
			ess.resize(UDER::esscount);
			for (int j = 0; j < UDER::esscount; ++j)
				ess[j] = ARecvSizes[j][s];
			ARecv = new UDER();
		}

		SpParHelper::BCastMatrix(cgr_Y->GetRowWorld(), *ARecv, ess, s);


		auto t_end = std::chrono::high_resolution_clock::now();
		stats.t_sC_comm_bcastA += static_cast<std::chrono::duration<double> >
			(t_end-t_beg).count();
		
		t_beg = std::chrono::high_resolution_clock::now();
		

		// broadcast X
		NUV *XRecv_cur = XRecv;
		if (s == Xself)
			XRecv_cur = (NUV *)X.arr.data();

		MPI_Bcast(XRecv_cur, XRecvSizes[s], MPIType<NUV>(),
				  s, cgr_Y->GetColWorld());


		t_end = std::chrono::high_resolution_clock::now();
		stats.t_sC_comm_bcastX += static_cast<std::chrono::duration<double> >
			(t_end-t_beg).count();

		t_beg = std::chrono::high_resolution_clock::now();
		

		// ofs << "stage " << s << " ARecv " << ARecv->getnrow() << " "
		// 	<< ARecv->getncol() << " " << ARecv->getnnz() << "\n";
		// ofs << "stage " << s << " XRecv size " << XRecvSizes[s] << "\n";
		// ofs << "stage " << s << " Y " << Y.m_ << " " << Y.n_ << std::endl;
		// ofs << std::flush;

		// assert (XRecv_cur != NULL && Y.arr_.data() != NULL);
		csc_gespmm_cusparse<SR>(*ARecv, XRecv_cur, Y.arr.data(),
								(int)X.getncol(), static_cast<NUM>(1.0), stats);

		t_end = std::chrono::high_resolution_clock::now();
		stats.t_comp += static_cast<std::chrono::duration<double> >
			(t_end-t_beg).count();

		// Y.PrintToFile("SpMM-Y");
	}


	delete [] XRecv;

	auto t_tot_end = std::chrono::high_resolution_clock::now();
	stats.t_tot += static_cast<std::chrono::duration<double> >
		(t_tot_end-t_tot_beg).count();
	
	// ofs.close();
	
	return Y;
}

	

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// sampling-based nnz estimation via SpMV
// @OGUZ-NOTE This is not based on SUMMA, do not use. Estimates the number of
// nonzeros in the final output matrix.


#define NROUNDS 5
typedef std::array<float, NROUNDS> samparr_t;

template <typename NZT>
struct promote_trait<NZT, samparr_t>
{
	typedef samparr_t T_promote;
};



class SamplesSaveHandler
{
public:
	template<typename c, typename t, typename V>
	void save(std::basic_ostream<c, t> &os,
			  std::array<V, NROUNDS> &sample_vec,
			  int64_t index)
	{
		for (auto it = sample_vec.begin(); it != sample_vec.end(); ++it)
			os << *it << " ";
	}
};



template<typename NZT>
struct SelectMinxSR
{
	static samparr_t id()
	{
		samparr_t arr;
		for (auto it = arr.begin(); it != arr.end(); ++it)
			*it = std::numeric_limits<float>::max();
		return arr;
	}


	static bool returnedSAID()
	{
		return false;
	}


	static samparr_t
	add (const samparr_t &arg1, const samparr_t &arg2)
	{
		samparr_t out;
		for (int i = 0; i < NROUNDS; ++i)
			out[i] = std::min(arg1[i], arg2[i]);
		return out;
	}


	static samparr_t
	multiply (const NZT arg1, const samparr_t &arg2)
	{
		return arg2;
	}


	static void axpy (const NZT a, const samparr_t &x, samparr_t &y)
	{
		y = add(y, multiply(a, x));
	}


	static MPI_Op mpi_op()
	{
		static MPI_Op mpiop;
		static bool exists = false;
		if (exists)
			return mpiop;
		else
		{
			MPI_Op_create(MPI_func, true, &mpiop);
			exists = true;
			return mpiop;
		}
	}


	static void
	MPI_func(void *invec, void *inoutvec, int *len, MPI_Datatype *datatype)
	{
		samparr_t   *in    = static_cast<samparr_t *>(invec);
		samparr_t   *inout = static_cast<samparr_t *>(inoutvec);
		for (int i = 0; i < *len; ++i)
			inout[i] = add(inout[i], in[i]);
	}
};



template <typename IU, typename NU1, typename NU2,
		  typename UDERA, typename UDERB>
int64_t
EstPerProcessNnzSpMV(
    SpParMat<IU, NU1, UDERA> &A, SpParMat<IU, NU2, UDERB> &B
	)  
{
	int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	float lambda = 1.0f;

	int nthds = 1;
	#ifdef THREADED
	#pragma omp parallel
	#endif
	{
		nthds = omp_get_num_threads();
	}

	if (myrank == 0)
		std::cout << "taking transposes." << std::endl;
	
	A.Transpose();
	B.Transpose();

	if (myrank == 0)
		std::cout << "setting initial samples." << std::endl;
	
	samparr_t sa;
	FullyDistVec<IU, samparr_t> samples_init(A.getcommgrid(), A.getncol(), sa);

	#ifdef THREADED
	#pragma omp parallel
	#endif
	{
		std::default_random_engine gen;
		std::exponential_distribution<float> exp_dist(lambda);

		#ifdef THREADED
		#pragma omp parallel for
		#endif
		for (IU i = 0; i < samples_init.LocArrSize(); ++i)
		{
			samparr_t tmp;
			for (auto it = tmp.begin(); it != tmp.end(); ++it)
				*it = exp_dist(gen);
			samples_init.SetLocalElement(i, tmp);
		}
	}

	// std::string fname("samples_init");
	// samples_init.ParallelWrite(fname, 1, SamplesSaveHandler(), true);
	
	if (myrank == 0)
		std::cout << "computing mid samples." << std::endl;

	FullyDistVec<IU, samparr_t> samples_mid =
		SpMV<SelectMinxSR<NU1> > (A, samples_init);

	// fname = "samples_mid";
	// samples_mid.ParallelWrite(fname, 1, SamplesSaveHandler(), true);

	if (myrank == 0)
		std::cout << "computing final samples." << std::endl;

	FullyDistVec<IU, samparr_t> samples_final =
		SpMV<SelectMinxSR<NU2> > (B, samples_mid);

	// fname = "samples_final";
	// samples_final.ParallelWrite(fname, 1, SamplesSaveHandler(), true);
	
	if (myrank == 0)
		std::cout << "computing nnz estimation." << std::endl;
	
	float nnzest = 0.0f;

	std::cout << myrank << "samples_final loc size: "
			  << samples_final.LocArrSize() << std::endl;

	const samparr_t *lsamples = samples_final.GetLocArr();
	
	#ifdef THREADED
	#pragma omp parallel for reduction (+:nnzest)
	#endif
	for (IU i = 0; i < samples_final.LocArrSize(); ++i)
	{
		float tmp = 0.0f;
		for (auto it = lsamples[i].begin(); it != lsamples[i].end(); ++it)
			tmp += *it;
		nnzest += static_cast<float>(NROUNDS - 1) / tmp;
	}

	if (myrank == 0)
		std::cout << "taking transposes again." << std::endl;

	int64_t nnzC_est = nnzest;
	int64_t nnzC_tot = 0;
	MPI_Allreduce(&nnzC_est, &nnzC_tot, 1, MPIType<int64_t>(), MPI_SUM,
				  (B.commGrid)->GetWorld());
	
	if (myrank == 0)
		std::cout << "sampling-based spmv est tot: " << nnzC_tot << std::endl;

	// revert back
	A.Transpose();
	B.Transpose();

	return nnzC_tot;
	
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


}


#endif


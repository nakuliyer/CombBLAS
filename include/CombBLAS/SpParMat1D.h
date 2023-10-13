#ifndef _SP_PAR_MAT_1D_H_
#define _SP_PAR_MAT_1D_H_

#include <iostream>
#include <fstream>
#include <cmath>
#include <memory>
#include <mpi.h>
#include <tuple>
#include <vector>
#include <iterator>

#include "SpMat.h"
#include "SpTuples.h"
#include "SpDCCols.h"
#include "CommGrid.h"
#include "CommGrid1D.h"

#include "MPIType.h"
#include "LocArr.h"
#include "SpDefs.h"
#include "Deleter.h"
#include "SpHelper.h"
#include "SpParHelper.h"
#include "FullyDistVec.h"
#include "Friends.h"
#include "Operations.h"
#include "DistEdgeList.h"
#include "mtSpGEMM.h"
#include "MultiwayMerge.h"
#include "CombBLAS.h"

namespace combblas
{
template <class IT, class NT>
std::tuple<IT,IT,NT>* ExchangeDataGeneral(std::vector<std::vector<std::tuple<IT,IT,NT>>> & tempTuples, MPI_Comm World, IT& datasize);

    enum SpParMat1DTYPE{
        ROWWISE = 1,
        COLWISE = 2,
    };
    // 1D only works for square matrix for now
    template <class IT, class NT, class DER>
    class SpParMat1D{
    public:
        typedef typename DER::LocalIT LocalIT;
        typedef typename DER::LocalNT LocalNT;
        typedef IT GlobalIT;
        typedef NT GlobalNT;
        typedef std::vector<std::tuple<IT,IT,NT>> vtuple;
        typedef std::vector<std::vector<std::tuple<IT,IT,NT>>> vvtuple;
        // Constructors
        SpParMat1D (const SpParMat < IT,NT,DER > & A2D, SpParMat1DTYPE type);
        SpParMat1D (std::shared_ptr<CommGrid1D> commgrid, SpParMat1DTYPE type=SpParMat1DTYPE::COLWISE);
        ~SpParMat1D ();
        int Owner(IT grow, IT gcol, IT & lrow, IT & lcol,SpParMat1DTYPE type) const;
        SpParMat1D< IT,NT,DER > & operator+=(const SpParMat1D< IT,NT,DER > & rhs);
        void Prune();
        IT getblocksize(){return blocksize;}
        bool allclose(const SpParMat1D< IT,NT,DER > & rhs);
        SpParMat<IT, NT, DER> ConvertTo2D();
        template <typename SR, typename NUO, typename UDERO, typename IU, typename NU1, typename NU2, typename UDERA, typename UDERB> 
        friend SpParMat1D<IU, NUO, UDERO> Mult_AnXBn_1D(SpParMat1D<IU,NU1,UDERA> & A, SpParMat1D<IU,NU2,UDERB> & B, bool clearA, bool clearB );
#ifndef NDDEBUG
        template <typename SR, typename NUO, typename UDERO, typename IU, typename NU1, typename NU2, typename UDERA, typename UDERB> 
        friend SpParMat1D<IU, NUO, UDERO> Mult_AnXBn_Diag(SpParMat1D<IU,NU1,UDERA> & A, SpParMat1D<IU,NU2,UDERB> & B, bool clearA, bool clearB );
#endif
        IT getnrow() const;
        IT getncol() const;
        friend class SpParMat<IT, NT, DER>;
    private:
        DER * spSeq;
        std::shared_ptr<CommGrid1D> grid1d;
        SpParMat1DTYPE mattype;
        int blocksize;
        IT totallength;
        IT colinmyrank;
        IT colprefix;
    };
}

#include "SpParMat1D.cpp"

#endif


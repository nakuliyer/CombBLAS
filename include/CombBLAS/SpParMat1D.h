#ifndef _SP_PAR_MAT_1D_H_
#define _SP_PAR_MAT_1D_H_

#include <iostream>
#include <fstream>
#include <cmath>
#include <mpi.h>
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

    enum SparMat1DTYPE{
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
        
        // Constructors
        SpParMat1D (const SpParMat < IT,NT,DER > & A2D, int blocksize, SparMat1DTYPE type);
        ~SpParMat1D ();
        
        template <typename LIT>
        int Owner(IT total_length, IT grow, IT gcol, LIT & lrow, LIT & lcol,SparMat1DTYPE type) const;

        template <typename SR, typename NUO, typename UDERO, typename IU, typename NU1, typename NU2, typename UDER1, typename UDER2>
        friend SpParMat3D<IU,NUO,UDERO> Mult_AnXBn_1D(SpParMat3D<IU,NU1,UDER1> & A, SpParMat3D<IU,NU2,UDER2> & B);
        
    private:
        std::shared_ptr<CommGrid1D> commGrid1D;
        SparMat1DTYPE type;
        int blocksize;
        int total_length;
        int maxblocksizeperrank; 
    };
}

#include "SpParMat1D.cpp"

#endif


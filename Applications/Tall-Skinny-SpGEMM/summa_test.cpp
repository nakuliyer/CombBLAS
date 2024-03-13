#include <mpi.h>
#include <sys/time.h>
#include <iostream>
#include <functional>
#include <algorithm>
#include <vector>
#include <sstream>
#include "CombBLAS/CombBLAS.h"
#include "CombBLAS/CommGrid3D.h"
#include "CombBLAS/SpParMat3D.h"
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

  if(argc < 7){
    if(myrank == 0)
    {
      cout << "Usage: ./<Binary> -A <MatrixA> -B <MatrixB> -C <MatrixC> -permute < sym | nonsym | no >" << endl;
    }
    MPI_Finalize();
    return -1;
  }
  else {
    string Aname;
    string Bname;
    string Cname = "";
    string perm;
    for (int i = 1; i < argc; i++)
    {
        if (strcmp(argv[i],"-A")==0){
            Aname = string(argv[i+1]);
            if(myrank == 0) printf("A Matrix filename: %s\n", Aname.c_str());
        }
        if (strcmp(argv[i],"-B")==0){
            Bname = string(argv[i+1]);
            if(myrank == 0) printf("B Matrix filename: %s\n", Bname.c_str());
        }
        if (strcmp(argv[i],"-C")==0){
            Cname = string(argv[i+1]);
            if(myrank == 0) printf("C Matrix filename: %s\n", Cname.c_str());
        }
        if (strcmp(argv[i],"-permute")==0){
            perm = string(argv[i+1]);
            if(myrank == 0) printf("Random permutation: %s\n", perm.c_str());
        }
    }
    shared_ptr<CommGrid> fullWorld;
    fullWorld.reset( new CommGrid(MPI_COMM_WORLD, 0, 0) );

    double t0, t1, t2;

    SpParMat<int64_t, double, SpDCCols < int64_t, double >> A2D(fullWorld);
    SpParMat<int64_t, double, SpDCCols < int64_t, double >> B2D(fullWorld);

    // Read matrix market
    t0 = MPI_Wtime();
    A2D.ParallelReadMM(Aname,true, maximum<double>());
    t1 = MPI_Wtime();
    if(myrank == 0) fprintf(stdout, "Time taken to read A: %lf\n", t1-t0);
    A2D.PrintInfo();
    t0 = MPI_Wtime();
    B2D.ParallelReadMM(Bname,true, maximum<double>());
    t1 = MPI_Wtime();
    if(myrank == 0) fprintf(stdout, "Time taken to read B: %lf\n", t1-t0);
    A2D.PrintInfo();

	double load_imbalance_A = A2D.LoadImbalance();
	double load_imbalance_B = B2D.LoadImbalance();
    if(myrank == 0) fprintf(stdout, "Load imbalance of A: %lf\n", load_imbalance_A);
    if(myrank == 0) fprintf(stdout, "Load imbalance of B: %lf\n", load_imbalance_B);
    
    if(perm != "no"){
        FullyDistVec<int64_t, int64_t> p( fullWorld );
        FullyDistVec<int64_t, int64_t> q( fullWorld );
        FullyDistVec<int64_t, int64_t> r( fullWorld );
        p.iota(A2D.getnrow(), 0);
        q.iota(A2D.getncol(), 0);
        r.iota(B2D.getncol(), 0);
        p.RandPerm(123);
        q.RandPerm(456);
        if(perm == "unsym") r.RandPerm(789);

        t0 = MPI_Wtime();
        if(perm == "unsym") (A2D)(p,q,true);
        else (A2D)(p,p,true);
        t1 = MPI_Wtime();
        if(myrank == 0) fprintf(stdout, "Time taken to permuate A: %lf\n", t1-t0);

        t0 = MPI_Wtime();
        (B2D)(p,r,true);// in-place permute to save memory
        t1 = MPI_Wtime();
        if(myrank == 0) fprintf(stdout, "Time taken to permuate B: %lf\n", t1-t0);

        load_imbalance_A = A2D.LoadImbalance();
        load_imbalance_B = B2D.LoadImbalance();
        if(myrank == 0) fprintf(stdout, "After permutation - Load imbalance of A: %lf\n", load_imbalance_A);
        if(myrank == 0) fprintf(stdout, "After permutation - Load imbalance of B: %lf\n", load_imbalance_B);
    }

    typedef PlusTimesSRing<double, double> PTFF;

    //Run 2D multiplication to compare against
    t0 = MPI_Wtime();
    SpParMat<int64_t, double, SpDCCols < int64_t, double >> C2D =
        Mult_AnXBn_Synch<PTFF, double, SpDCCols<int64_t, double>, int64_t, double, double, SpDCCols<int64_t, double>, SpDCCols<int64_t, double> >
        (A2D, B2D);
    t1 = MPI_Wtime();
    if(myrank == 0) fprintf(stdout, "Time taken for Mult_AnXBn_Synch: %lf\n", t1-t0);
    if(Cname != "") C2D.ParallelWriteMM(Cname, 1);
    if(myrank == 0) fprintf(stdout, "2D Multiplication done \n");
  }
  MPI_Finalize();
  return 0;
}

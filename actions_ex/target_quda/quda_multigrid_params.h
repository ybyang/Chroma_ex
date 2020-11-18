#ifndef __QUDA_MULTIGRID_PARAMS_v2_H__
#define __QUDA_MULTIGRID_PARAMS_v2_H__

#include "chromabase.h"
#include <string>
#include "actions/ferm/invert/quda_solvers/enum_quda_io.h"


namespace Chroma 
{

  struct MULTIGRIDSolverParams_2 {
    
    Real tol;
    int  maxIterations;
    QudaSolverType smootherType;
    bool verbosity;
    QudaPrecisionType prec;
    QudaReconsType reconstruct;
    QudaSchwarzMethod schwarzType;
    multi1d<int> nvec;
    int mg_levels; 
    bool generate_nullspace;
    bool generate_all_levels;
    bool check_multigrid_setup;
    multi1d<bool> setup_on_gpu;
    multi1d<int> nu_pre;
    multi1d<int> nu_post;
    multi1d< multi1d<int> > blocking;
    int outer_gcr_nkrylov;
    int precond_gcr_nkrylov;
    std::string cycle_type;
    Real relaxationOmegaMG;
    Real relaxationOmegaOuter;
    multi1d<int> maxIterSubspaceCreate;
    multi1d<Real> rsdTargetSubspaceCreate;
    multi1d<int> maxIterSubspaceRefresh;
    std::string nullspace_file;
    
    MULTIGRIDSolverParams_2(XMLReader& xml, const std::string& path);
    MULTIGRIDSolverParams_2() {
      tol = .000001;
      relaxationOmegaMG =Real(1.0);
      relaxationOmegaOuter = Real(1.0);
      maxIterations = 10;
      smootherType = MR;
      verbosity = false;
      prec = DEFAULT;
      reconstruct = RECONS_NONE;
      schwarzType = ADDITIVE_SCHWARZ;

      mg_levels = 2;
      setup_on_gpu.resize(mg_levels-1);
      // Default is to set up on gpu
      for(int i=0; i < mg_levels - 1; ++i) { 
        setup_on_gpu[i] = true;
      }
      blocking.resize(mg_levels-1);
      nvec.resize(mg_levels-1);
      nu_pre.resize(mg_levels-1);
      nu_post.resize(mg_levels-1);
      maxIterSubspaceCreate.resize(mg_levels-1);
      maxIterSubspaceRefresh.resize(mg_levels-1);
      rsdTargetSubspaceCreate.resize(mg_levels-1);

      generate_nullspace = true;
      for(int l = 0; l < mg_levels - 1; l++) 
      {
	blocking[l].resize(4);
        blocking[l][0] = blocking[l][1] = blocking[l][2] = blocking[l][3] = 4;
	nu_pre[l] = 2;
	nu_post[l] = 2;
	nvec[l] = 16;

	// Default params: 
	maxIterSubspaceCreate[l] = 500;
	rsdTargetSubspaceCreate[l] = 5.0e-6;
	maxIterSubspaceRefresh[l] = maxIterSubspaceCreate[l];
 
      }
      outer_gcr_nkrylov = 12;
      precond_gcr_nkrylov = 12;

      generate_all_levels = true;
      check_multigrid_setup= true;
      cycle_type = "MG_VCYCLE";

    };

  };
  void read(XMLReader& xml, const std::string& path, MULTIGRIDSolverParams_2& p);
 
  void write(XMLWriter& xml, const std::string& path, 
	     const MULTIGRIDSolverParams_2& param);

}
#endif

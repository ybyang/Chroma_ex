
#include "chromabase.h"
#include <string>
#include "chroma_config.h"
#ifdef BUILD_QUDA
#include "actions_ex/target_quda/quda_multigrid_params.h"

using namespace QDP;

namespace Chroma {
  MULTIGRIDSolverParams_2::MULTIGRIDSolverParams_2(XMLReader& xml, 
					     const std::string& path)
  {
    XMLReader paramtop(xml, path);
    read(paramtop, "Residual", tol);
    read(paramtop, "MaxIterations", maxIterations);
    read(paramtop, "SmootherType", smootherType);
    read(paramtop, "Verbosity", verbosity);
    read(paramtop, "Precision", prec);
    read(paramtop, "Reconstruct", reconstruct);
    
    read(paramtop, "Blocking", blocking);
    mg_levels = blocking.size()+1;
    nvec.resize(mg_levels-1);
    nu_pre.resize(mg_levels-1);
    nu_post.resize(mg_levels-1);
    maxIterSubspaceCreate.resize(mg_levels-1);
    maxIterSubspaceRefresh.resize(mg_levels-1);
    rsdTargetSubspaceCreate.resize(mg_levels-1);

    read(paramtop, "NullVectors", nvec);
    read(paramtop, "Pre-SmootherApplications", nu_pre);
    read(paramtop, "Post-SmootherApplications", nu_post);
    if (nvec.size() != mg_levels-1 ) {
 
      QDPIO::cout<<"Warning. There are "<< blocking.size() 
		 << " blockings but only " << nvec.size() << " sets of NullVectors" << std::endl;
      QDP_abort(1);
    }

    if (nu_pre.size() != mg_levels-1 ) {
 
      QDPIO::cout<<"Error. There are "<< (mg_levels-1)  
		 << " blockings but only " << nu_pre.size() << " sets pre-smoothing iterations" << std::endl;
      QDP_abort(1);
    }


    if( paramtop.count("./RsdTargetSubspaceCreate") == 1 ) { 
      read(paramtop, "RsdTargetSubspaceCreate", rsdTargetSubspaceCreate);
    }
    else {
      // Default values.
      for(int l=0; l < mg_levels-1; ++l) { 
	rsdTargetSubspaceCreate[l] = 5.0e-6;
      }
    }

    if( paramtop.count("./MaxIterSubspaceCreate") == 1) {
      read(paramtop, "./MaxIterSubspaceCreate", maxIterSubspaceCreate);
    }
    else {
      for(int l=0; l < mg_levels-1; ++l) {
	maxIterSubspaceCreate[l] = 500;
      }
    }

    if( paramtop.count("./MaxIterSubspaceRefresh") == 1) { 
      read(paramtop, "./MaxIterSubspaceRefresh", maxIterSubspaceRefresh);
    }
    else {
      for(int l=0; l < mg_levels-1; ++l) { 
	maxIterSubspaceRefresh[l] = maxIterSubspaceCreate[l];
      }
    }


    if( paramtop.count("./OuterGCRNKrylov") == 1 ) {
	read(paramtop, "OuterGCRNKrylov", outer_gcr_nkrylov);
    }
    else { 
	outer_gcr_nkrylov = 12;
    }

    if( paramtop.count("./PrecondGCRNKrylov") == 1 ) { 
	read(paramtop, "PrecondGCRNKrylov", precond_gcr_nkrylov);
    }
    else {
        precond_gcr_nkrylov = 12;
    }

    if (nu_post.size() != mg_levels-1 ) {
 
      QDPIO::cout<<"Warning. There are "<< blocking.size() 
		 << " blockings but only " << nu_post.size() << " sets post-smoothing iterations " << std::endl;
      QDP_abort(1);
    }

    read(paramtop, "GenerateNullspace", generate_nullspace);

    if( paramtop.count("./CheckMultigridSetup") == 1 ) {
      read(paramtop, "CheckMultigridSetup", check_multigrid_setup);
    }
    else {
      check_multigrid_setup = true; // Default: Checks for suspace and coarss op ar one!
    }

    read(paramtop, "GenerateAllLevels", generate_all_levels);

    //FIXME: Should be an enum
    read(paramtop, "CycleType", cycle_type);
  

    // Read optional params otherwise use defaults
    if( paramtop.count("SchwarzType") > 0 ) { 
      read(paramtop, "SchwarzType", schwarzType);
    }

    // Read optional relaxation param
    if( paramtop.count("RelaxationOmegaMG") > 0 ) {
      read(paramtop, "RelaxationOmegaMG", relaxationOmegaMG);
    }
    else { 
	relaxationOmegaMG = Real(1.0);
    }

   // Read optional relaxation param
    if( paramtop.count("RelaxationOmegaOuter") > 0 ) {
      read(paramtop, "RelaxationOmegaOuter", relaxationOmegaOuter);
    }
    else {
        relaxationOmegaOuter = Real(1.0);
    }

    if(paramtop.count("SetupOnGPU") == 1) {
      read(paramtop, "SetupOnGPU", setup_on_gpu);
      if ( setup_on_gpu.size() != mg_levels - 1 ) { 

        // if size != n_levels - 1 then it is an error 
        // unless size is 1 in which case we broadcast it.
        //
	if ( setup_on_gpu.size() == 1 ) { 

          // Only one value was entered broadcast it
          bool value = setup_on_gpu[0];
          setup_on_gpu.resize(mg_levels-1);
          for(int l=0; l < mg_levels-1; ++l) { 
            setup_on_gpu[l] = value;
          } 
        }
        else { 
	   // Size mismatch and size is not 1 error. 
	   QDPIO::cerr << "setup_on_gpu has size = " << setup_on_gpu.size() << 
              " but it should be either 1 (broadcast to all levels) or " << mg_levels -1 
	      << " (mg_levels -1) \n";
           QDP_abort(1);
        } // setup_on_gpu_size == 1
      } // setup_on_gpu_size != mg_levels - 1
    } // paramtop.count("SetupOnGPU") == 1
    else {

      // No specification in XML. Default
      // behaviour is that all levels set up on GPU
      setup_on_gpu.resize(mg_levels-1);
      for(int l=0; l < mg_levels-1; ++l) {
        setup_on_gpu[l] = true;
      }
    }
    
  };

  void read(XMLReader& xml, const std::string& path, 
	    MULTIGRIDSolverParams_2& p)
  {
    MULTIGRIDSolverParams_2 tmp(xml, path);
    p = tmp;
  }

  
  
  void write(XMLWriter& xml, const std::string& path, 
	     const MULTIGRIDSolverParams_2& p) {
    push(xml, path);
    write(xml, "Residual", p.tol);
    write(xml, "MaxIterations", p.maxIterations);
    write(xml, "SmootherType", p.smootherType);
    write(xml, "RelaxationOmegaMG", p.relaxationOmegaMG);
    write(xml, "RelaxationOmegaOuter", p.relaxationOmegaOuter);
    write(xml, "Verbosity", p.verbosity);
    write(xml, "Precision", p.prec);
    write(xml, "Reconstruct", p.reconstruct);
    write(xml, "SchwarzType", p.schwarzType);
    write(xml, "NullVectors", p.nvec);
    write(xml, "MultiGridLevels", p.mg_levels);
    write(xml, "GenerateNullSpace", p.generate_nullspace);
    write(xml, "GenerateAllLevels", p.generate_all_levels);
    write(xml, "CheckMultigridSetup", p.check_multigrid_setup);
    write(xml, "CycleType", p.cycle_type);
    write(xml, "Pre-SmootherApplications", p.nu_pre);
    write(xml, "Post-SmootherApplications", p.nu_post);
    /* FIXME: This should go in the general solver interface, and work for all GCR solvers, not just GCR inner params */
    write(xml, "OuterGCRNKrylov", p.outer_gcr_nkrylov);
    write(xml, "PrecondGCRNKrylov", p.precond_gcr_nkrylov);
    write(xml, "Blocking", p.blocking);
    write(xml, "MaxIterSubspaceCreate", p.maxIterSubspaceCreate);
    write(xml, "MaxIterSubspaceRefresh", p.maxIterSubspaceRefresh);
    write(xml, "RsdTargetSubspaceCreate", p.rsdTargetSubspaceCreate);
    write(xml, "SetupOnGPU", p.setup_on_gpu);
    pop(xml);

  }
}
#endif

// -*- C++ -*-
// $Id: gaus_quark_smearing.h,v 3.2 2008-11-04 18:43:57 edwards Exp $
/*! \file
 *  \brief Gaussian smearing of color std::vector and propagator
 */

#ifndef __mom_quark_smearing_h__
#define __mom_quark_smearing_h__

#include "chroma_config.h"

#ifdef BUILD_QUDA
#include <quda.h>
#include "actions/ferm/invert/quda_solvers/enum_quda_io.h"
#endif

#include "meas/smear/quark_smearing.h"
#include "actions/boson/operator/klein_gord.h"

namespace Chroma {
//! Name and registration
/*! @ingroup smear */
namespace MomQuarkSmearingEnv {
bool registerAll();

//! Return the name
std::string getName();

//! Params for Gauss quark smearing
/*! @ingroup smear */
struct Params {
  Params() {
  }
  Params(XMLReader& in, const std::string& path);
  void writeXML(XMLWriter& in, const std::string& path) const;

  Real wvf_param; /*!< Smearing width */
  int wvfIntPar; /*!< Number of smearing hits */
  int no_smear_dir; /*!< No smearing in this direction */
   int smear_dir; /*!< Only apply smearing in this direction */
  int smear_tr_size; /*!< the size of operator in the transverse direction */
#ifdef BUILD_QUDA  
  QudaReconsType cudaReconstruct;
  QudaPrecisionType cudaPrecision;
  bool tuneDslashP, verboseP, checkP;
#endif  
  bool qudaSmearingP;
  

  // !!! THIS CODE NOW ONLY WORKS FOR MOMENTUM-INJECTED SMEARING !!! //
  multi1d<Real> mom; /*<! quasi-momentum */
};

//! Gaussian quark smearing
/*! @ingroup smear
 *
 * Gaussian quark smearing object
 */
 
 
template<typename T>
class QuarkSmear: public QuarkSmearing<T> {
public:

  //! Full constructor
  QuarkSmear(const Params& p) :
        params(p) {
#ifdef BUILD_QUDA
           QudaPrecision_s cpu_prec;
           QudaPrecision_s gpu_prec;

           int s = sizeof( WordType<Real>::Type_t );
           if (s == 4) {
             cpu_prec = QUDA_SINGLE_PRECISION;
           }
           else {
             cpu_prec = QUDA_DOUBLE_PRECISION;
           }
      
           switch( params.cudaPrecision ) {
           case HALF:
             gpu_prec = QUDA_HALF_PRECISION;
             break;
           case SINGLE:
             gpu_prec = QUDA_SINGLE_PRECISION;
             break;
           case DOUBLE:
             gpu_prec = QUDA_DOUBLE_PRECISION;
             break;
           default:
             gpu_prec = cpu_prec;
             break;
           }
      
           q_gauge_param = newQudaGaugeParam();
           
           const multi1d<int>& latdims = Layout::subgridLattSize();
           q_gauge_param.X[0] = latdims[0];
           q_gauge_param.X[1] = latdims[1];
           q_gauge_param.X[2] = latdims[2];
           q_gauge_param.X[3] = latdims[3];
           q_gauge_param.type = QUDA_SU3_LINKS;
           q_gauge_param.gauge_order = QUDA_QDP_GAUGE_ORDER; 
           q_gauge_param.t_boundary = QUDA_PERIODIC_T;//smearing is independent to BC
           q_gauge_param.cpu_prec = cpu_prec;
           q_gauge_param.cuda_prec = gpu_prec;
                 
           switch( params.cudaReconstruct ) {
           case RECONS_NONE:
             q_gauge_param.reconstruct = QUDA_RECONSTRUCT_NO;
             break;
           case RECONS_8:
             q_gauge_param.reconstruct = QUDA_RECONSTRUCT_8;
             break;
           case RECONS_12:
             q_gauge_param.reconstruct = QUDA_RECONSTRUCT_12;
             break;
           default:
             q_gauge_param.reconstruct = QUDA_RECONSTRUCT_12;
             break;
           };
      
           q_gauge_param.cuda_prec_sloppy = gpu_prec;
           q_gauge_param.reconstruct_sloppy = QUDA_RECONSTRUCT_12;
      
           q_gauge_param.gauge_fix = QUDA_GAUGE_FIXED_NO;
           q_gauge_param.anisotropy = 1.0;//smearing don't care this;
      
        // etup padding
           multi1d<int> face_size(4);
           face_size[0] = latdims[1]*latdims[2]*latdims[3]/2;
           face_size[1] = latdims[0]*latdims[2]*latdims[3]/2;
           face_size[2] = latdims[0]*latdims[1]*latdims[3]/2;
           face_size[3] = latdims[0]*latdims[1]*latdims[2]/2;
      
           int max_face = face_size[0];
           for(int i=1; i <=3; i++) {
             if ( face_size[i] > max_face ) {
               max_face = face_size[i];
             }
           }
      
           q_gauge_param.ga_pad = max_face;
           q_gauge_param.cuda_prec_precondition = gpu_prec;
           q_gauge_param.reconstruct_precondition = QUDA_RECONSTRUCT_12;
      
           quda_inv_param = newQudaInvertParam();
           quda_inv_param.cpu_prec = cpu_prec;
           quda_inv_param.cuda_prec = gpu_prec;
           quda_inv_param.cuda_prec_sloppy = gpu_prec;
           quda_inv_param.preserve_source = QUDA_PRESERVE_SOURCE_NO;
           quda_inv_param.dirac_order = QUDA_DIRAC_ORDER;            
           quda_inv_param.input_location = QUDA_CPU_FIELD_LOCATION;
           quda_inv_param.output_location = QUDA_CPU_FIELD_LOCATION;
           
           quda_inv_param.dslash_type = QUDA_WILSON_DSLASH;
           quda_inv_param.gamma_basis = QUDA_DEGRAND_ROSSI_GAMMA_BASIS;
           quda_inv_param.Ls = 1;
           quda_inv_param.solution_type = QUDA_MAT_SOLUTION ;
           
           if( params.tuneDslashP ) {
             QDPIO::cout << "Enabling Dslash Autotuning" << std::endl;
      
             quda_inv_param.tune = QUDA_TUNE_YES;
           }
           else {
             QDPIO::cout << "Disabling Dslash Autotuning" << std::endl;
      
             quda_inv_param.tune = QUDA_TUNE_NO;
           }
           quda_inv_param.sp_pad = 0;
           quda_inv_param.cl_pad = 0;
           if( params.verboseP ) {   
//             quda_inv_param.verbosity = QUDA_VERBOSE;
             quda_inv_param.verbosity = QUDA_DEBUG_VERBOSE;
           }
           else {
             quda_inv_param.verbosity = QUDA_SUMMARIZE;
           }
#endif        
  }

  void gausSmear_plane(const multi1d<LatticeColorMatrix>& u, T& chi, Real width,
      int ItrGaus, int j_decay, int j_smear, int tr_size) {
    T psi;
  
    Real ftmp = -(width * width) / Real(4 * ItrGaus);
    /* The Klein-Gordon operator is (Lapl + mass_sq), where Lapl = -d^2/dx^2.. */
    /* We want (1 + ftmp * Lapl ) = (Lapl + 1/ftmp)*ftmp */
    Real ftmpi = Real(1) / ftmp;
    Real ftmp2 = Real(2 * Nd) + ftmpi;
  
    if (j_decay < Nd) ftmp2 = Real(2 * Nd - 2) + ftmpi;
    if(j_smear>=0) ftmp2=2+ftmpi;

    multi1d<Real> mom(Nd);
    for(int mu=0;mu<Nd;mu++) mom[mu]=0.0;
  
    for (int n = 0; n < ItrGaus; ++n) {
      psi = chi * ftmp;
      chi = psi * ftmp2;
      for (int mu = 0; mu < Nd; ++mu)
        if ((mu != j_decay&&j_smear<0)||mu==j_smear) {
          chi -= u[mu] * shift(psi, FORWARD, mu)
                 + shift(adj(u[mu]) * psi, BACKWARD, mu);
        }
    }
    
    if(tr_size>1)
    {
       T vp1,vm1;
       for (int mu =0; mu<Nd; mu++)
       if(mu!=j_smear && mu!=j_decay)
       {
           vp1=chi;vm1=chi;
           for (int n = 1; n < tr_size; ++n)
           {
               psi=shift(vp1, FORWARD, mu);vp1=psi;
               psi=shift(vm1, BACKWARD, mu);vm1=psi;
               chi+=vp1;
               chi+=vm1;
           }
       }
    }
      
  }    

  //! Smear the quark
//  void operator()(T& quark, const multi1d<LatticeColorMatrix>& u) const ;
   void operator()(T& quark, const multi1d<LatticeColorMatrix>& u) const{
      StopWatch Timer;
      Timer.reset();
      Timer.start();
      multi1d<LatticeColorMatrix> links_single(Nd);
      const Real twopi = 6.283185307179586476925286;
      if(params.wvfIntPar>0)
      for(int mu=0; mu < Nd; mu++) 
      {
          Real k = params.mom[mu]*twopi/Layout::lattSize()[mu];
          QDPIO::cout << mu << ":" << cos(k) << "+" << -sin(k) << "I;"  << std::endl;
          links_single[mu] = u[mu]*cmplx(cos(k),-sin(k));
      }
        
      LatticeComplex corr=trace(adj(quark)*quark);
      double value=sum(corr).elem().elem().elem().real();
      QDPIO::cout << "Norm2 is (before smearing)" << value << std::endl;
      Timer.stop();
      QDPIO::cout << "time to prepare = "
              << Timer.getTimeInSeconds() 
              << " seconds" << std::endl;

#ifdef BUILD_QUDA
      if(params.qudaSmearingP) {

           void* gauge[4];
           for(int mu=0; mu < Nd; mu++) {
#ifndef BUILD_QUDA_DEVIFACE_GAUGE
               gauge[mu] = (void *)&(links_single[mu].elem(all.start()).elem().elem(0,0).real());
#else
               gauge[mu] = QDPCache::Instance().getDevicePtr( links_single[mu].getId() );
               QDPIO::cout << "MDAGM CUDA gauge[" << mu << "] in = " << gauge[mu] << "\n";
#endif
           }
           QDPIO::cout << "smearing on GPU" << std::endl;

             Timer.reset();
             Timer.start();
           loadGaugeQuda((void *)gauge, &q_gauge_param);
           Timer.stop();
           QDPIO::cout << "time to load gauge = "
              << Timer.getTimeInSeconds() 
              << " seconds" << std::endl;

             Timer.reset();
             Timer.start();
           if(params.checkP)
           {
              T quark_cpu=quark;
              if(params.smear_dir>=0||params.smear_tr_size>0) 
                   QDPIO::cout << "block smearing is not supported at present" << std::endl;
              gausSmear_plane(links_single, quark_cpu, params.wvf_param, params.wvfIntPar, params.no_smear_dir,params.smear_dir, params.smear_tr_size);
              apply(quark);
              quark_cpu=quark_cpu-quark;
              corr=trace(adj(quark_cpu)*quark_cpu);
              value=sum(corr).elem().elem().elem().real();
              corr=trace(adj(quark)*quark);
              Real value2=sum(corr).elem().elem().elem().real();
              
              QDPIO::cout << "Norm of diff/vec is" << sqrt(value/value2) << std::endl;
              fflush(stdout);
           }
           else
              apply(quark);
           Timer.stop();
           QDPIO::cout << "time to smearing on GPU = "
              << Timer.getTimeInSeconds() 
              << " seconds" << std::endl;
      }
      else {
#endif   
      gausSmear_plane(links_single, quark, params.wvf_param, params.wvfIntPar, params.no_smear_dir,params.smear_dir, params.smear_tr_size);
#ifdef BUILD_QUDA
      }
#endif   
      corr=trace(adj(quark)*quark);
      value=sum(corr).elem().elem().elem().real();
      QDPIO::cout << "Norm2 is (after smearing)" << value << std::endl;
      fflush(stdout);
   }
  
private:
  //! Hide partial constructor
  QuarkSmear() {
  }

#ifdef BUILD_QUDA
  QudaInvertParam quda_inv_param;
  QudaGaugeParam q_gauge_param;
  void apply(T& quark);
#endif

private:
  Params params; /*!< smearing params */
};

}  // end namespace

//! Reader
/*! @ingroup smear */
void read(XMLReader& xml, const std::string& path,
    MomQuarkSmearingEnv::Params& param);

//! Writer
/*! @ingroup smear */
void write(XMLWriter& xml, const std::string& path,
    const MomQuarkSmearingEnv::Params& param);

}  // end namespace Chroma

#endif

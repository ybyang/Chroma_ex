// $Id: gaus_quark_smearing.cc,v 3.3 2008-11-04 18:43:57 edwards Exp $
/*! \file
 *  \brief Gaussian smearing of color std::vector
 */

#include "chromabase.h"

#include "meas/smear/quark_smearing_factory.h"
#include "meas_ex/mom_quark_smearing.h"
#include "meas/smear/gaus_smear.h"
#include "util/ferm/transf.h"

#ifdef BUILD_QUDA
#include <comm_quda.h>
#endif

namespace Chroma {

// Read parameters
void read(XMLReader& xml, const std::string& path,
    MomQuarkSmearingEnv::Params& param) {
  MomQuarkSmearingEnv::Params tmp(xml, path);
  param = tmp;
}

//! Parameters for running code
void write(XMLWriter& xml, const std::string& path,
    const MomQuarkSmearingEnv::Params& param) {
  param.writeXML(xml, path);
}

//! Hooks to register the class
namespace MomQuarkSmearingEnv {
namespace {
//! Callback function
QuarkSmearing<LatticePropagator>* createProp(XMLReader& xml_in,
    const std::string& path) {
  return new QuarkSmear<LatticePropagator>(Params(xml_in, path));
}

//! Callback function
QuarkSmearing<LatticeStaggeredPropagator>* createStagProp(XMLReader& xml_in,
    const std::string& path) {
  return new QuarkSmear<LatticeStaggeredPropagator>(Params(xml_in, path));
}

//! Callback function
QuarkSmearing<LatticeFermion>* createFerm(XMLReader& xml_in,
    const std::string& path) {
  return new QuarkSmear<LatticeFermion>(Params(xml_in, path));
}

//! Callback function
QuarkSmearing<LatticeColorVector>* createColorVec(XMLReader& xml_in,
    const std::string& path) {
  return new QuarkSmear<LatticeColorVector>(Params(xml_in, path));
}

//! Local registration flag
bool registered = false;

//! Name to be used
const std::string name = "MOM_GAUSSIAN";
}

//! Return the name
std::string getName() {
  return name;
}

//! Register all the factories
bool registerAll() {
  bool success = true;
  if (!registered) {
    success &= Chroma::ThePropSmearingFactory::Instance().registerObject(name,
      createProp);
    success &= Chroma::TheStagPropSmearingFactory::Instance().registerObject(
      name, createStagProp);
    success &= Chroma::TheFermSmearingFactory::Instance().registerObject(name,
      createFerm);
    success &= Chroma::TheColorVecSmearingFactory::Instance().registerObject(
      name, createColorVec);
    registered = true;
  }
  return success;
}

//! Parameters for running code
Params::Params(XMLReader& xml, const std::string& path) {
  XMLReader paramtop(xml, path);

  read(paramtop, "wvf_param", wvf_param);
  read(paramtop, "wvfIntPar", wvfIntPar);
  read(paramtop, "no_smear_dir", no_smear_dir);
  if( paramtop.count("smear_dir") > 0 ){
      read(paramtop, "smear_dir", smear_dir);
      if( paramtop.count("smear_tr_size") > 0 )
          read(paramtop, "smear_tr_size", smear_tr_size);
      else
          smear_tr_size=1;
  }
  else
  {
      smear_tr_size=0;
      smear_dir=-1;
  }

  if( paramtop.count("qudaSmearingP") > 0 ) {
       read(paramtop, "qudaSmearingP", qudaSmearingP);
#ifdef BUILD_QUDA
            if( paramtop.count("CudaPrecision") > 0 && qudaSmearingP) {
               read(paramtop, "CudaPrecision", cudaPrecision);
             }
             else {
               cudaPrecision = DEFAULT;
             }
             if ( paramtop.count("Verbose") > 0 && qudaSmearingP) {
               read(paramtop, "Verbose", verboseP);
             }
             else {
               verboseP = false;
             }
             if( paramtop.count("AutotuneDslash") > 0 && qudaSmearingP) {
               read(paramtop, "AutotuneDslash", tuneDslashP);
             }
             else {
               tuneDslashP = false;
             }
             if( paramtop.count("CudaReconstruct") > 0 && qudaSmearingP) {
               read(paramtop, "CudaReconstruct", cudaReconstruct);
             }
             else {
               cudaReconstruct = RECONS_NONE;
             }
             if( paramtop.count("checkP") > 0 && qudaSmearingP) {
               read(paramtop, "checkP", checkP);
             }
             else {
               checkP = false;
             }
#endif 
 }
 else {
   qudaSmearingP = false;
 }
   
  read(paramtop, "mom", mom);
  if (mom.size() != Nd) {
    QDPIO::cerr << name << ": wrong size of mom array: expected length=" << Nd
    << std::endl;
    QDP_abort(1);
  }
}

//! Parameters for running code
void Params::writeXML(XMLWriter& xml, const std::string& path) const {
  push(xml, path);

  write(xml, "wvf_kind", MomQuarkSmearingEnv::getName());
  write(xml, "wvf_param", wvf_param);
  write(xml, "wvfIntPar", wvfIntPar);
  write(xml, "no_smear_dir", no_smear_dir);
  write(xml, "mom", mom);

  pop(xml);
}

#ifdef BUILD_QUDA


void apply_sm(LatticeFermion& quark,QudaInvertParam &quda_inv_param, Params &params)
{
           if(params.wvfIntPar==0) return;
           
           LatticeFermion mod_chi,psi_s;
           mod_chi[rb[0]] = zero;
           mod_chi[rb[1]] = quark;
#ifndef BUILD_QUDA_DEVIFACE_SPINOR           
           void* spinorIn =(void *)&(quark.elem(0).elem(0).elem(0).real());
           void* spinorOut =(void *)&(psi_s.elem(0).elem(0).elem(0).real());
#else
           void* spinorIn = GetMemoryPtr( quark.getId() );
           void* spinorOut = GetMemoryPtr( psi_s.getId() );
#endif        
           double wvf_param=params.wvf_param.elem().elem().elem().elem();
           double ftmp=-(wvf_param*wvf_param)/(4.0*params.wvfIntPar);
           double alpha=-ftmp/(1+6*ftmp);
           performWuppertalnStep(spinorOut,spinorIn,&quda_inv_param,params.wvfIntPar,alpha);
           quark=psi_s;
}

void apply_sm(LatticeColorVector& quark,QudaInvertParam &quda_inv_param, Params &params)
{

          LatticeFermion tmp=zero;
          pokeSpin(tmp,quark,0);
          apply_sm(tmp,quda_inv_param,params);
          quark=peekSpin(tmp, 0);
}

template<>
void QuarkSmear<LatticeColorVector>::apply(LatticeColorVector& quark)
{
          apply_sm(quark,quda_inv_param,params);
}

template<>
void QuarkSmear<LatticeFermion>::apply(LatticeFermion& quark)
{
          apply_sm(quark,quda_inv_param,params);
}

template<>
void QuarkSmear<LatticePropagator>::apply(LatticePropagator& quark)
{
      LatticeFermion tmp;
      for(int is=0;is<4;is++)
      for(int ic=0;ic<3;ic++)
      {
           PropToFerm(quark,tmp,ic,is);
           apply_sm(tmp,quda_inv_param,params);
           FermToProp(tmp,quark,ic,is);
      }
}

template<>
void QuarkSmear<LatticeStaggeredPropagator>::apply(LatticeStaggeredPropagator& quark)
{
      QDPIO::cout << "Staggered fermion is not supported here" << std::endl;
}

#endif
}  // end namespace
}  // end namespace Chroma


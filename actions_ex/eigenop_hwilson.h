// -*- C++ -*-
// $Id: hwilson_eigenop.h, v0.9 2020-07-22$
/*! \file
 * \brief The H-Wilson eigensystem.
 *
 * Calculates the eigensystem 
 */

#ifndef __Hwilson_eigenOP_h__
#define __Hwilson_eigenOP_h__

#include "actions_ex/eigenop.h"
#include "actions/ferm/linop/dslash_w.h"
#ifdef BUILD_QUDA
#include <quda.h>
#endif

namespace Chroma 
{ 

  class HwilsonEigenOperator: public EigenOperator<LatticeFermion>
  {
  public:
     
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;
    
    HwilsonEigenOperator(GroupXML_t &fermact,multi1d<LatticeColorMatrix> &u,int neig=0):EigenOperator<LatticeFermion>(fermact,u,neig)
    {
         //read param;
         std::istringstream  is(fermact.xml);
         XMLReader  paramtop(is);
         Real kappa;
         read(paramtop,"kappa",kappa);
         rho = 4-0.5/kappa;

         AnisoParam_t anisoParam;
         D.create(fs,anisoParam);
    }
     
    void operator()(LatticeFermion& chi, const LatticeFermion& psi, enum PlusMinus isign) const
    {
         LatticeFermion tmp,tmp2;
         D(tmp,psi,isign);
    	 chi=Gamma(15)*((4-rho)*psi-0.5*tmp);
    }

    const Real& Rho() const {return rho;}
    
    const FermState<T,P,Q>& getFermState() const {return *fs;}

    int create_eigen(InlineEigenMakerParams &Params)
    {
          
    }
     

  protected:    
    Real rho;
    WilsonDslash D;
  };

}
#endif

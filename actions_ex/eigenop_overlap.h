// -*- C++ -*-
// $Id: hwilson_eigenop.h, v0.9 2020-07-22$
/*! \file
 * \brief The H-Wilson eigensystem.
 *
 * Calculates the eigensystem 
 */

#ifndef __Overlap_eigenOP_h__
#define __Overlap_eigenOP_h__

#include "actions_ex/eigenop_hwilson.h"

namespace Chroma 
{ 

  class OverlapEigenOperator: public EigenOperator<LatticeFermion>
  {
  public:
     
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;
    
    OverlapEigenOperator(GroupXML_t &fermact,multi1d<LatticeColorMatrix> &u,int neig=0):EigenOperator<LatticeFermion>(fermact,u,neig)
    {
        
    
    }
    
    void eps(LatticeFermion& chi, const LatticeFermion& psi, Real prec) const
    {
         T high_src,low;
         
         // deflate psi;
         // low=\sum_i hw.vec_i*sign(hw.val_i)*conj(hw.vec_i); high_src=(1-\sum_i hw.vec_i*conj(hw.vec_i))psi;
         
         // generate high:
         // high=adj(rotate)*sign(Hw->oper_chi)*tmp;
         
         //psi=low+high;
    
    }
    
    void eps5(LatticeFermion& chi, const LatticeFermion& psi, Real prec) const
    {
         T tmp;
         eps(tmp,psi,prec);
         chi=Gamma(15)*tmp; //Gamma(8) in the chiral basis is gamma_5
    }
     
    void operator()(LatticeFermion& chi, const LatticeFermion& psi, enum PlusMinus isign) const 
    {
         T tmp;
         eps5(tmp,psi,1e-12);
         chi=rho*(psi+tmp);
    }
    
    int create_eigen(InlineEigenMakerParams &Params)
    {
          
    }
     
  protected:    
    Real rho;
    HwilsonEigenOperator *Hw;
  };

};

#endif

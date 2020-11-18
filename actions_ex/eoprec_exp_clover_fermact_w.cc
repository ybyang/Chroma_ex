/*! \file
 *  \brief Luscher's Even-odd preconditioned Exponential Clover fermion action
 
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"

#include "actions_ex/eoprec_exp_clover_linop_w.h"
#include "actions_ex/eoprec_exp_clover_fermact_w.h"

//#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermstates/ferm_createstate_reader_w.h"

namespace Chroma
{

  //! Hooks to register the class with the fermact factory
  namespace EvenOddPrecExpCloverFermActEnv
  {
    //! Callback function
    WilsonTypeFermAct<LatticeFermion,
		      multi1d<LatticeColorMatrix>,
		      multi1d<LatticeColorMatrix> >* createFermAct4D(XMLReader& xml_in,
								     const std::string& path)
    {
      return new EvenOddPrecExpCloverFermAct(CreateFermStateEnv::reader(xml_in, path), 
					  CloverFermActParams(xml_in, path));
    }

    //! Callback function
    /*! Differs in return type */
    FermionAction<LatticeFermion,
		  multi1d<LatticeColorMatrix>,
		  multi1d<LatticeColorMatrix> >* createFermAct(XMLReader& xml_in,
							       const std::string& path)
    {
      return createFermAct4D(xml_in, path);
    }

    //! Name to be used
    const std::string name = "EXP_CLOVER";

    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= Chroma::TheFermionActionFactory::Instance().registerObject(name, createFermAct);
	success &= Chroma::TheWilsonTypeFermActFactory::Instance().registerObject(name, createFermAct4D);
	registered = true;
      }
      return success;
    }
  }

  //! Produce a linear operator for this action
  /*!
   * The operator acts on the odd subset
   *
   * \param state	    gauge field     	       (Read)
   */
  EvenOddPrecLogDetLinearOperator<LatticeFermion,
				  multi1d<LatticeColorMatrix>,
				  multi1d<LatticeColorMatrix> >* 
  EvenOddPrecExpCloverFermAct::linOp(Handle< FermState<T,P,Q> > state) const
  {
    return new EvenOddPrecExpCloverLinOp(state,param);
  }

}


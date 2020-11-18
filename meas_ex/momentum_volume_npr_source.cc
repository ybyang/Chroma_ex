// $Id: mom_source_const.cc,v 3.7 2008-11-04 18:43:58 edwards Exp $
/*! \file
 *  \brief Momentum (wall) source construction
 */

#include "chromabase.h"

#include "meas/sources/source_const_factory.h"
#include "meas_ex/momentum_volume_npr_source.h"
#include "util/ft/sftmom.h"
#include "util/ferm/transf.h"

namespace Chroma
{
  // Read parameters
  void read(XMLReader& xml, const std::string& path, MomWallNPRQuarkSourceConstEnv::Params& param)
  {
    MomWallNPRQuarkSourceConstEnv::Params tmp(xml, path);
    param = tmp;
  }

  // Writer
  void write(XMLWriter& xml, const std::string& path, const MomWallNPRQuarkSourceConstEnv::Params& param)
  {
    param.writeXML(xml, path);
  }

  //! Fill a specific color and spin index with 1.0 within a volume
  void boxfil(LatticeFermion& a, int color_index, int spin_index);
/*
  {
    START_CODE();

    if (color_index >= Nc || color_index < 0)
      QDP_error_exit("invalid color index", color_index);

    if (spin_index >= Ns || spin_index < 0)
      QDP_error_exit("invalid spin index", spin_index);

    // Write ONE to all field
    Real one = 1;
    Complex sitecomp = cmplx(one,0);
    ColorVector sitecolor = zero;
    Fermion sitefield = zero;

    pokeSpin(sitefield,
	     pokeColor(sitecolor,sitecomp,color_index),
	     spin_index);

    // Broadcast to all sites
    a = sitefield;  // QDP (not installed version) now supports   construct OLattice = OScalar
      
    END_CODE();
  }
*/

  //! Hooks to register the class
  namespace MomWallNPRQuarkSourceConstEnv
  {
    namespace
    {
      //! Callback function
      QuarkSourceConstruction<LatticePropagator>* createProp(XMLReader& xml_in,
							     const std::string& path)
      {
	return new SourceConst<LatticePropagator>(Params(xml_in, path));
      }

      //! Name to be used
      const std::string name("MOMENTUM_VOLUME_NPR_SOURCE");

      //! Local registration flag
      bool registered = false;
    }

    //! Return the name
    std::string getName() {return name;}

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(name, createProp);
	registered = true;
      }
      return success;
    }


    //! Initialize
    Params::Params()
    {
       j_decay=-1;
    }


    //! Read parameters
    Params::Params(XMLReader& xml, const std::string& path)
    {
      XMLReader paramtop(xml, path);

      int version;
      read(paramtop, "version", version);

      switch (version) 
      {
      case 1:
	break;

      default:
	QDPIO::cerr << __func__ << ": parameter version " << version 
		    << " unsupported." << std::endl;
	QDP_abort(1);
      }

      read(paramtop, "mom", mom);

      if (mom.size() != Nd)
      {
	QDPIO::cerr << name << ": wrong size of mom array: expected length=" << Nd << std::endl;
	QDP_abort(1);
      }
    }


    // Writer
    void Params::writeXML(XMLWriter& xml, const std::string& path) const
    {
      push(xml, path);

      int version = 1;
      write(xml, "version", version);

      write(xml, "mom", mom);

      pop(xml);
    }


    //! Construct the source
    template<>
    LatticePropagator
    SourceConst<LatticePropagator>::operator()(const multi1d<LatticeColorMatrix>& u) const
    {
      QDPIO::cout << "Volume Momentum Source" << std::endl;

      // Initialize the slow Fourier transform phases
        LatticeReal theta=zero;
        for(int id=0;id<4;id++)
            theta+=twopi*Real(params.mom[id])/Real(Layout::lattSize()[id])* QDP::Layout::latticeCoordinate(id);
        LatticeComplex  phase=cmplx(cos(theta),sin(theta));

	QDPIO::cout<<params.mom[0]<<params.mom[1]<<params.mom[2]<<params.mom[3]<<std::endl;

      // Create the quark source
      LatticePropagator quark_source;
      for(int color_source = 0; color_source < Nc; ++color_source)
      {
	for(int spin_source = 0; spin_source < Ns; ++spin_source)
	{
	  // MomWall fill a fermion source. Insert it into the propagator source
	  LatticeFermion chi;
	  boxfil(chi, color_source, spin_source);
	  // Multiply in the time direction phases (not handled in sftmom)
	  chi *= phase;
	  FermToProp(chi, quark_source, color_source, spin_source);
	}
      }

      return quark_source;
    }

  }

}

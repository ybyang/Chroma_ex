// $Id: simple_baryon_seqsrc_w.cc,v 3.7 2009-03-19 17:17:20 mcneile Exp $
/*! \file
 *  \brief Construct baryon sequential sources.
 */

#include "qdp_config.h"

#include "meas/hadron/simple_baryon_seqsrc_w.h"
#include "meas_ex/simple_baryon_seqsrc_w.h"
#include "meas/hadron/seqsource_factory_w.h"
#include "meas/hadron/barspinmat_w.h"
#include "meas/hadron/barhqlq_w.h"
#include "meas/hadron/baryon_spinmat_funcmap_w.h"
#include "util/ft/sftmom.h"
#include "io/xml_group_reader.h"

namespace Chroma 
{

  //! Baryon sequential sources
  /*! \ingroup hadron */
  namespace SimpleBaryonSeqSourceEnv
  { 

    //! Anonymous namespace
    namespace
    {

       void checkArgs(const char* name, const multi1d<LatticePropagator>& quark_propagators)
       {
         if (quark_propagators.size() != 1)
         {
//         QDPIO::cerr << __func__ << ": expect only 1 prop" << std::endl;
//         QDP_abort(1);

           std::ostringstream s;
           s << name << ": expecting 1 prop, instead passed = " << quark_propagators.size() << std::endl;
           throw s.str();
         }
       }
    }

//! Nucleon-Nucleon ISO piece with general projector and Cg5
    LatticePropagator
    BarNucl_ISO_TCg5::operator()(const multi1d<LatticeColorMatrix>& u,
                             const multi1d<ForwardProp_t>& forward_headers,
                             const multi1d<LatticePropagator>& quark_propagators)
    {
      START_CODE();
      if ( Nc != 3 ){    /* Code is specific to Ns=4 and Nc=3. */
        QDPIO::cerr<<" code only works for Nc=3 and Ns=4\n";
        QDP_abort(111) ;
      }
#if QDP_NC == 3

      checkArgs("BarNucl_ISO_TCg5", quark_propagators);

      LatticePropagator src_prop_tmp;
      LatticePropagator q1_tmp;
      LatticePropagator q2_tmp;

      LatticePropagator di_quark;
      LatticeColorMatrix col_mat;

      /* "\bar u O u" insertion in NR proton, ie. 
       * "(u Cg5 d) u" */
      /* Some generic T */

      // Use precomputed Cg5
      q1_tmp = quark_propagators[0] * Cg5;
      q2_tmp = Cg5 * quark_propagators[0];
      di_quark = quarkContract24(q1_tmp, q2_tmp);

      // First term
      src_prop_tmp = T * di_quark;

      // Now the second term
      src_prop_tmp += traceSpin(di_quark) * T;

      // The third term...
      q1_tmp = q2_tmp * Cg5;
      q2_tmp = quark_propagators[0] * T;

      src_prop_tmp -= quarkContract13(q1_tmp, q2_tmp) + transposeSpin(quarkContract12(q2_tmp, q1_tmp));

      /* "\bar d O d" insertion in NR proton, ie. 
       * "(u Cg5 d) u" */
      /* Some generic T */

      // First term
      q2_tmp = quark_propagators[0] * Cg5;
      q1_tmp = T * q2_tmp;

      q2_tmp = Cg5 * quark_propagators[0];
      src_prop_tmp+= quarkContract14(q1_tmp, q2_tmp);

      // Second term
      q1_tmp = q2_tmp * Cg5;
      q2_tmp = quark_propagators[0] * T;

      src_prop_tmp+= transposeSpin(quarkContract12(q2_tmp, q1_tmp));

      END_CODE();

      return projectBaryon(src_prop_tmp,
                           forward_headers);
#else
      LatticePropagator q1_tmp;

      q1_tmp = zero ;
      return q1_tmp ;
#endif
    }


    Complex
    BarNucl_ISO_TCg5::twoPtSink(const multi1d<LatticeColorMatrix>& u,
                            const multi1d<ForwardProp_t>& forward_headers,
                            const multi1d<LatticePropagator>& quark_propagators,
                            int gamma_insertion)
    {
      checkArgs("BarNucl_ISO_TCg5", quark_propagators);
      setTSrce(forward_headers);
      setBC(forward_headers);

      // Constructor the 2pt
      LatticeComplex b_prop = Baryon2PtContractions::sigma2pt(quark_propagators[0],
                                                              quark_propagators[0],
                                                              T, Cg5);

      // Extract the sink at the appropriate momenta
      SftMom sft(0, getTSrce(), getSinkMom(), false, getDecayDir());
      multi2d<DComplex> hsum;
      hsum = sft.sft(b_prop);

      // D-quark rescaling factor is 1 for this type of seqsource
      return Real(1) * hsum[0][getTSink()];
    }

    //! Anonymous namespace
    namespace
    {


      //-------------------- callback functions ---------------------------------------


      //! \bar u O u" insertion in NR proton, ie. "(u C gamma_5 (1/2)(1 + gamma_4)  d) u" 
      /*!
       * \ingroup hadron
       * 
       * \f$C g_5 NR = (1/2)*C gamma_5 * ( 1 + g_4 )\f$
       * 
       * \f$T = (1 + \Sigma_3)*(1 + gamma_4) / 2 
       *      = (1 + Gamma(8) - i G(3) - i G(11)) / 2\f$
       */
      HadronSeqSource<LatticePropagator>* barNuclUFullMixedNR(XMLReader& xml_in,
							  const std::string& path)
      {
        SpinMatrix g_one = 1.0;
        
       // (1 + A_1 + A_2 + A_3)*(1 + V4) / 2 =(1+ g8 - I*g3 - I*g11 - I*g6 - I*g14 + I*g5 + I*g13)
      
	return new BarNuclUTCg5(SimpleBaryonSeqSourceEnv::Params(xml_in, path), 
	      SpinMatrix(0.5*(g_one + Gamma(8)*g_one + timesMinusI(Gamma(3)*g_one + Gamma(6)*g_one + Gamma(11)*g_one + Gamma(14)*g_one)
	             +timesI(Gamma(5)*g_one + Gamma(13)*g_one))),
 	      BaryonSpinMats::Cg5NR());
      }


      //! "\bar d O d" insertion in NR proton, ie. "(u C gamma_5 (1/2)(1 + gamma_4)  d) u"
      /*!
       * \ingroup hadron
       * 
       * C g_5 NR = (1/2)*C gamma_5 * ( 1 + g_4 )
       *
       * T = (1 + \Sigma_3)*(1 + gamma_4) / 2 
       *   = (1 + Gamma(8) - i G(3) - i G(11)) / 2
       */
      HadronSeqSource<LatticePropagator>* barNuclDFullMixedNR(XMLReader& xml_in,
							  const std::string& path)
      {
        SpinMatrix g_one = 1.0;
        
       // (1 + A_1 + A_2 + A_3)*(1 + V4) / 2 =(1+ g8 - I*g3 - I*g11 - I*g6 - I*g14 + I*g5 + I*g13)
	return new BarNuclDTCg5(SimpleBaryonSeqSourceEnv::Params(xml_in, path), 
	      SpinMatrix(0.5*(g_one + Gamma(8)*g_one + timesMinusI(Gamma(3)*g_one + Gamma(6)*g_one + Gamma(11)*g_one + Gamma(14)*g_one)
	             +timesI(Gamma(5)*g_one + Gamma(13)*g_one))),
 	      BaryonSpinMats::Cg5NR());
      }
      
      HadronSeqSource<LatticePropagator>* barNucl_ISO_MixedNR(XMLReader& xml_in,
							  const std::string& path)
      {
        SpinMatrix g_one = 1.0;
        
	return new BarNucl_ISO_TCg5(SimpleBaryonSeqSourceEnv::Params(xml_in, path), 
	             BaryonSpinMats::Tmixed(),BaryonSpinMats::Cg5NR());
      }
      
      HadronSeqSource<LatticePropagator>* barNucl_ISO_UnpolNR(XMLReader& xml_in,
							  const std::string& path)
      {
        SpinMatrix g_one = 1.0;
        
	return new BarNucl_ISO_TCg5(SimpleBaryonSeqSourceEnv::Params(xml_in, path), 
	             BaryonSpinMats::Tunpol(),BaryonSpinMats::Cg5NR());
      }

      HadronSeqSource<LatticePropagator>* barNucl_ISO_Pol_xNR(XMLReader&
							  xml_in, const
							  std::string& path)
      {
        // A_1*(1 + V4) / 2 =-I(g6+g14)/2
        SpinMatrix g_one = 1.0;
        
	return new BarNucl_ISO_TCg5(SimpleBaryonSeqSourceEnv::Params(xml_in, path), 
	             SpinMatrix(0.5*timesMinusI(Gamma(6)*g_one + Gamma(14)*g_one)),BaryonSpinMats::Cg5NR());
      }

      HadronSeqSource<LatticePropagator>* barNuclUPol_xNR(XMLReader&
							  xml_in, const
							  std::string& path)
      {
        SpinMatrix g_one = 1.0;
        
	return new BarNuclUTCg5(SimpleBaryonSeqSourceEnv::Params(xml_in, path), 
	             SpinMatrix(0.5*timesMinusI(Gamma(6)*g_one + Gamma(14)*g_one)),BaryonSpinMats::Cg5NR());
      }

      HadronSeqSource<LatticePropagator>* barNuclDPol_xNR(XMLReader&
							  xml_in, const
							  std::string& path)
      {
        SpinMatrix g_one = 1.0;
        
	return new BarNuclDTCg5(SimpleBaryonSeqSourceEnv::Params(xml_in, path), 
	             SpinMatrix(0.5*timesMinusI(Gamma(6)*g_one + Gamma(14)*g_one)),BaryonSpinMats::Cg5NR());
      }

      HadronSeqSource<LatticePropagator>* barNucl_ISO_Mixed(XMLReader& xml_in,
							  const std::string& path)
      {
        SpinMatrix g_one = 1.0;
        
	return new BarNucl_ISO_TCg5(SimpleBaryonSeqSourceEnv::Params(xml_in, path), 
	             BaryonSpinMats::Tmixed(),BaryonSpinMats::Cg5());
      }
      
      HadronSeqSource<LatticePropagator>* barNucl_ISO_Unpol(XMLReader& xml_in,
							  const std::string& path)
      {
        SpinMatrix g_one = 1.0;
        
	return new BarNucl_ISO_TCg5(SimpleBaryonSeqSourceEnv::Params(xml_in, path), 
	             BaryonSpinMats::Tunpol(),BaryonSpinMats::Cg5());
      }

      HadronSeqSource<LatticePropagator>* barNucl_ISO_Pol_x(XMLReader&
							  xml_in, const
							  std::string& path)
      {
        // A_1*(1 + V4) / 2 =-I(g6+g14)/2
        SpinMatrix g_one = 1.0;
        
	return new BarNucl_ISO_TCg5(SimpleBaryonSeqSourceEnv::Params(xml_in, path), 
	             SpinMatrix(0.5*timesMinusI(Gamma(6)*g_one + Gamma(14)*g_one)),BaryonSpinMats::Cg5());
      }

      HadronSeqSource<LatticePropagator>* barNuclUPol_x(XMLReader&
							  xml_in, const
							  std::string& path)
      {
        SpinMatrix g_one = 1.0;
        
	return new BarNuclUTCg5(SimpleBaryonSeqSourceEnv::Params(xml_in, path), 
	             SpinMatrix(0.5*timesMinusI(Gamma(6)*g_one + Gamma(14)*g_one)),BaryonSpinMats::Cg5());
      }

      HadronSeqSource<LatticePropagator>* barNuclDPol_x(XMLReader&
							  xml_in, const
							  std::string& path)
      {
        SpinMatrix g_one = 1.0;
        
	return new BarNuclDTCg5(SimpleBaryonSeqSourceEnv::Params(xml_in, path), 
	             SpinMatrix(0.5*timesMinusI(Gamma(6)*g_one + Gamma(14)*g_one)),BaryonSpinMats::Cg5());
      }


      
      bool registered_extra = false;

    }// end anonymous namespace

    //! Register all the factories
    bool registerAll_extra() 
    {
      bool success = true; 
      if (! registered_extra)
      {
	//! Register needed stuff

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(std::string("NUCL_U_FULL_MIXED_NONREL"),
										      barNuclUFullMixedNR);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(std::string("NUCL_D_FULL_MIXED_NONREL"),  
										      barNuclDFullMixedNR);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(std::string("NUCL_ISO_MIXED_NONREL"),  
										      barNucl_ISO_MixedNR);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(std::string("NUCL_ISO_UNPOL_NONREL"),  
										      barNucl_ISO_UnpolNR);

        success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(std::string("NUCL_U_POL_x_NONREL"),
                                                                                      barNuclUPol_xNR);

        success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(std::string("NUCL_D_POL_x_NONREL"),
                                                                                      barNuclDPol_xNR);

        success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(std::string("NUCL_ISO_POL_x_NONREL"),
                                                                                      barNucl_ISO_Pol_xNR);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(std::string("NUCL_ISO_MIXED"),  
										      barNucl_ISO_Mixed);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(std::string("NUCL_ISO_UNPOL"),  
										      barNucl_ISO_Unpol);

        success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(std::string("NUCL_U_POL_x"),
                                                                                      barNuclUPol_x);

        success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(std::string("NUCL_D_POL_x"),
                                                                                      barNuclDPol_x);

        success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(std::string("NUCL_ISO_POL_x"),
                                                                                      barNucl_ISO_Pol_x);

	registered_extra = true;
      }
      return success;
    }

  } // namespace BaryonSeqSourceCallMapEnv


}  // end namespace Chroma


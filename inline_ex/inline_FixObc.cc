/*! \file
 *  \brief Inline plaquette
 */

#include "inline_ex/inline_FixObc.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/make_xml_file.h"
#include "meas/glue/mesplq.h"
#include "meas/inline/io/named_objmap.h"
#include "util/ft/sftmom.h"
#include "meas/smear/hyp_smear.h"

#include "meas/inline/io/default_gauge_field.h"

#include "actions/gauge/gaugestates/gauge_createstate_factory.h"
#include "actions/gauge/gaugestates/gauge_createstate_aggregate.h"

namespace Chroma 
{  

  //! OBC input
  void read(XMLReader& xml, const std::string& path, InlineOBCEnv::Params::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "gauge_id", input.gauge_id);
  }

  //! OBC output
  void write(XMLWriter& xml, const std::string& path, const InlineOBCEnv::Params::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "gauge_id", input.gauge_id);

    pop(xml);
  }


  namespace InlineOBCEnv 
  { 
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineMeas(Params(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "FIX_OBC";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= CreateGaugeStateEnv::registerAll();
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	registered = true;
      }
      return success;
    }


    // Params
    Params::Params() 
    { 
      frequency = 0; 
      named_obj.gauge_id = InlineDefaultGaugeField::getId();
      xml_file ="";
    }

    Params::Params(XMLReader& xml_in, const std::string& path) 
    {
      try 
      {
	XMLReader paramtop(xml_in, path);

	if (paramtop.count("Frequency") == 1)
	  read(paramtop, "Frequency", frequency);
	else
	  frequency = 1;

	// Ids
	read(paramtop, "NamedObject", named_obj);

	// Possible alternate XML file pattern
	if (paramtop.count("xml_file") != 0) {
	  read(paramtop, "xml_file", xml_file);
	}
      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << "Caught Exception reading XML: " << e << std::endl;
	QDP_abort(1);
      }
    }


    void 
    InlineMeas::operator()(unsigned long update_no,
			   XMLWriter& xml_out) 
    {
      if( params.xml_file != "" ) 
      {
	std::string xml_file = makeXMLFileName(params.xml_file, update_no);
	push( xml_out, "OBC");
	write(xml_out, "update_no", update_no);
	write(xml_out, "xml_file", xml_file);
	pop(xml_out);

	XMLFileWriter xml(xml_file);
	func(update_no, xml);
      }
      else 
      {
	func(update_no, xml_out);
      }

    }

    void 
    InlineMeas::func(const unsigned long update_no, 
		     XMLWriter& xml_out) 
    {
      START_CODE();
    
      // Test and grab a reference to the gauge field
      multi1d<LatticeColorMatrix> &u=
             TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);

      push(xml_out, "OBC");
      write(xml_out, "update_no", update_no);
      
      Double w_plaq, s_plaq, t_plaq, link; 
      multi2d<Double> plane_plaq;
      SftMom phases_nomom(0, true, 3);
      Set timeslice= phases_nomom.getSet();
      
      LatticeComplex site_act;
      multi1d<DComplex> hsum;
      size_t factor=Layout::lattSize()[0]*Layout::lattSize()[1]*Layout::lattSize()[2]*3;      
      
      MesPlq(u, w_plaq, s_plaq, t_plaq, plane_plaq, link);
      
      if(Layout::primaryNode())
         printf("Before fix: s_plaq=%20.10f, t_plaq=%20.10f, aver=%20.10f\n",toDouble(s_plaq),toDouble(t_plaq),toDouble(w_plaq));
         
      site_act= trace(u[3]*shift(u[2],FORWARD,3)*adj(shift(u[3],FORWARD,2))*adj(u[2]));
      hsum = sumMulti(site_act,timeslice);
      
      if(Layout::primaryNode())
      {
          for(int it=0;it<Layout::lattSize()[3];it++)
          if(it<Layout::lattSize()[3]/10||it>9*Layout::lattSize()[3]/10)
              printf("%4d%20.10f\n",it,toDouble(real(hsum[it]/factor)));
      }
                
      if(toDouble(real(hsum[Layout::lattSize()[3]-1]/factor))<1e-7)
      {
         LatticeColorMatrix u_tmp=where
           (Layout::latticeCoordinate(3)!=Layout::lattSize()[3]-1,
             u[3],LatticeColorMatrix(1.0));
         u[3]=u_tmp;
         if(Layout::primaryNode()) printf("Fix the configuration with OBC for HYP smearing...\n");
      }
      
      MesPlq(u, w_plaq, s_plaq, t_plaq, plane_plaq, link);
      
      if(Layout::primaryNode())
         printf("After fix: s_plaq=%20.10f, t_plaq=%20.10f, aver=%20.10f\n",toDouble(s_plaq),toDouble(t_plaq),toDouble(w_plaq));

      site_act= trace(u[3]*shift(u[2],FORWARD,3)*adj(shift(u[3],FORWARD,2))*adj(u[2]));
      hsum = sumMulti(site_act,timeslice);
      
      if(Layout::primaryNode())
      {
          for(int it=0;it<Layout::lattSize()[3];it++)
          if(it<Layout::lattSize()[3]/10||it>9*Layout::lattSize()[3]/10)
              printf("%4d%20.10f\n",it,toDouble(real(hsum[it]/factor)));
      }
      
      
/*
      multi1d<LatticeColorMatrix> u_hyp(Nd);
      Hyp_Smear(u, u_hyp, 
                      0.75, 0.6, 0.3, 1.0e-5, 100);

      site_act= trace(u_hyp[3]*shift(u_hyp[2],FORWARD,3)*adj(shift(u_hyp[3],FORWARD,2))*adj(u_hyp[2]));
      hsum = sumMulti(site_act,timeslice);
      
      if(Layout::primaryNode())
      {
          for(int it=0;it<Layout::lattSize()[3];it++)
          if(it<Layout::lattSize()[3]/10||it>9*Layout::lattSize()[3]/10)
              printf("%4d%13.5f\n",it,toDouble(real(hsum[it]/factor)));
      }
*/
      pop(xml_out); // pop("OBC");
    
      END_CODE();
    } 

  } // namespace InlineOBCEnv

} // namespace Chroma

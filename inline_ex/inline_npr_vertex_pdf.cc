/*! \file
 * \brief Inline construction of NPR vertices
 *
 * NPR vertices on  props
 */

#include "inline_ex/inline_npr_vertex_pdf.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "util/ft/sftmom.h"
// #include "meas/hadron/BuildingBlocks_w.h"
#include "util/info/proginfo.h"
#include "meas/inline/make_xml_file.h"

#include "meas_ex/npr_vertex_pdf.h"

#include "meas/inline/io/named_objmap.h"

#include "actions/ferm/fermstates/ferm_createstate_factory_w.h"
#include "actions/ferm/fermstates/ferm_createstate_aggregate_w.h"
#include <string>


namespace Chroma 
{ 
  namespace InlineNprVertexPDFEnv 
  { 
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineNprVertexPDF(InlineNprVertexPDFParams(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "NPR_VERTEX_PDF";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= CreateFermStateEnv::registerAll();
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	registered = true;
      }
      return success;
    }
  }

  DPropagator FTpropagator2(const LatticePropagator& prop,
			   const multi1d<int> mom,
			   const multi1d<int> t_src){
    // Initialize the slow Fourier transform phases
    // Makes no sence to do momemtum average here...

    LatticeReal theta=zero;
    for(int id=0;id<4;id++) 
        theta+=twopi*Real(mom[id])/Real(Layout::lattSize()[id])* QDP::Layout::latticeCoordinate(id); 
    LatticeComplex phase=cmplx(cos(theta),sin(theta));
    return sum(prop*phase);    
  }

  //! Param input
  void read(XMLReader& xml, const std::string& path, InlineNprVertexPDFParams::Param_t& input)
  {
    XMLReader paramtop(xml, path);

    int version;
    read(paramtop, "version", version);

    if (paramtop.count("FermState") != 0)
      input.cfs = readXMLGroup(paramtop, "FermState", "Name");
    else
      input.cfs = CreateFermStateEnv::nullXMLGroup();

    switch (version) 
    {
    case 1:
      break;

    default :
      QDPIO::cerr << InlineNprVertexPDFEnv::name << ": input parameter version " 
		  << version << " unsupported." << std::endl;
      QDP_abort(1);
    }
    
    read(paramtop, "links_max", input.links_max);
    if (paramtop.count("links_prep_max")>0) {
       read(paramtop, "links_prep_max", input.links_prep_max);
    }
    else
       input.links_prep_max=0;
    read(paramtop, "file_name", input.file_name);
    read(paramtop, "save_prop", input.save_prop);
    read(paramtop, "cfg_serial", input.cfg_serial);
  }


  //! Param write
  void write(XMLWriter& xml, const std::string& path, const InlineNprVertexPDFParams::Param_t& input)
  {
    push(xml, path);

    int version = 1;
    write(xml, "version", version);
    write(xml, "links_max", input.links_max);
    write(xml, "links_prep_max", input.links_prep_max);
    write(xml, "file_name", input.file_name);    
    write(xml, "save_prop", input.save_prop);
    write(xml, "cfg_serial", input.cfg_serial);  
    xml << input.cfs.xml;

    pop(xml);
  }

  //! Propagator input
  void read(XMLReader& xml, const std::string& path, InlineNprVertexPDFParams::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "gauge_id", input.gauge_id);
    read(inputtop, "prop_id", input.prop_id);
  }

  //! Propagator output
  void write(XMLWriter& xml, const std::string& path, const InlineNprVertexPDFParams::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "gauge_id", input.gauge_id);
    write(xml, "prop_id", input.prop_id);

    pop(xml);
  }


  // Param stuff
  InlineNprVertexPDFParams::InlineNprVertexPDFParams() {frequency = 0;}

  InlineNprVertexPDFParams::InlineNprVertexPDFParams(XMLReader& xml_in, const std::string& path) 
  {
    try 
    {
      XMLReader paramtop(xml_in, path);

      if (paramtop.count("Frequency") == 1)
	read(paramtop, "Frequency", frequency);
      else
	frequency = 1;

      // Read program parameters
      read(paramtop, "Param", param);

      // Read in the output propagator/source configuration info
      read(paramtop, "NamedObject", named_obj);

      // Possible alternate XML file pattern
      if (paramtop.count("xml_file") != 0) 
      {
	read(paramtop, "xml_file", xml_file);
      }
    }
    catch(const std::string& e) 
    {
      QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e << std::endl;
      QDP_abort(1);
    }
  }


  void
  InlineNprVertexPDFParams::write(XMLWriter& xml_out, const std::string& path) 
  {
    push(xml_out, path);
    
    Chroma::write(xml_out, "Param", param); 
    Chroma::write(xml_out, "NamedObject", named_obj);
    QDP::write(xml_out, "xml_file", xml_file);

    pop(xml_out);
  }


  //###################################################################################//
  // Accept All Link Patterns                                                          //
  //###################################################################################//

  void AllLinkPatterns( bool &                          DoThisPattern,
			bool &                          DoFurtherPatterns,
			multi1d< int > & LinkPattern );
/*  {
    DoThisPattern     = true;
    DoFurtherPatterns = true;
    
    return;
  }*/


  // Function call
  void 
  InlineNprVertexPDF::operator()(unsigned long update_no,
				   XMLWriter& xml_out) 
  {
    // If xml file not empty, then use alternate
    if (params.xml_file != "")
    {
      std::string xml_file = makeXMLFileName(params.xml_file, update_no);

      push(xml_out, "NprVertexPDF");
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

  void write_FT_prop(const LatticePropagator& prop,
//                     const multi1d<int> mom,
                     std::string name,int serial)
  {
//    multi1d<int> t_src(mom.size());
//    DPropagator FTprop(FTpropagator2(prop,mom,t_src));
    DPropagator FTprop(sum(prop));
  if(Layout::primaryNode())
  {
      general_data_base io_prop;
      sprintf(io_prop.name,"%s",name.c_str());
	io_prop.add_dimension(dim_conf, 1,&serial);
	io_prop.add_dimension(dim_temporary, 144);
	io_prop.add_dimension(dim_complex, 2);
	io_prop.initialize();
	
        SpinMatrixD U = KYToDRMat();
	FTprop= adj(U) * FTprop * U;
	
      int data_index=0;
      for(int is2=0; is2!=Ns; is2++)
      for(int ic2=0; ic2!=Nc; ic2++)
      for(int is1=0; is1!=Ns; is1++)
      for(int ic1=0; ic1!=Nc; ic1++)
      {
		io_prop.data[data_index]=FTprop.elem().elem(is1,is2).elem(ic1,ic2).real();
		data_index++;
		io_prop.data[data_index]=FTprop.elem().elem(is1,is2).elem(ic1,ic2).imag();
		data_index++;
      }
      io_prop.save();
    }
  
  }
                     

  // Function call
  void 
  InlineNprVertexPDF::func(unsigned long update_no,
			     XMLWriter& XmlOut) 
  {
    START_CODE();

    StopWatch snoop;
    snoop.reset();
    snoop.start();

    push(XmlOut, "NprVertexPDF");
    write(XmlOut, "update_no", update_no);

    QDPIO::cout << " ExampleNprVertexPDF" << std::endl;
    QDPIO::cout << "     volume: " << QDP::Layout::lattSize()[0];
    for (int i=1; i<Nd; ++i) {
      QDPIO::cout << " x " << QDP::Layout::lattSize()[i];
    }
    QDPIO::cout << std::endl;

    //#################################################################################//
    // XML output
    //#################################################################################//

    proginfo(XmlOut);    // Print out basic program info

    push(XmlOut, "Output_version");
    write(XmlOut, "out_version", 2);
    pop(XmlOut);

    //###############################################################################//
    // Read Gauge Field                                                              //
    //###############################################################################//

    QDPIO::cout << "Attempt to initialize the gauge field" << std::endl << std::flush ;

    // Grab the gauge field
    multi1d<LatticeColorMatrix> U;
    XMLBufferWriter gauge_xml;

    try
    {
      U = TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);
      TheNamedObjMap::Instance().get(params.named_obj.gauge_id).getRecordXML(gauge_xml);

      // Set the construct state and modify the fields
      {
	QDPIO::cout << "cfs=XX" << params.param.cfs.xml << "XX" << std::endl;
	std::istringstream  xml_s(params.param.cfs.xml);
	XMLReader  fermtop(xml_s);

	Handle<CreateFermState< LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > 
	  cfs(TheCreateFermStateFactory::Instance().createObject(params.param.cfs.id,
								 fermtop,
								 params.param.cfs.path));

	Handle<FermState< LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > 
	  state((*cfs)(U));

	// Pull the u fields back out from the state since they might have been
	// munged with fermBC's
	U = state->getLinks();
      }
    
    }
    catch( std::bad_cast ) 
    {
      QDPIO::cerr << InlineNprVertexPDFEnv::name << ": caught dynamic cast error" 
		  << std::endl;
      QDP_abort(1);
    }
    catch (const std::string& e) 
    {
      QDPIO::cerr << InlineNprVertexPDFEnv::name << ": std::map call failed: " << e 
		  << std::endl;
      QDP_abort(1);
    }
    catch( ... )
    {
      QDPIO::cerr << InlineNprVertexPDFEnv::name << ": caught generic exception "
		  << std::endl;
      QDP_abort(1);
    }

    // Write out the input
    params.write(XmlOut, "Input");

    // Write out the config info
    write(XmlOut, "Config_info", gauge_xml);

    // Calculate some gauge invariant observables just for info.
    MesPlq(XmlOut, "Observables", U);

    //#################################################################################//
    // Read Forward Propagator                                                         //
    //#################################################################################//

    SftMom phases_nomom( 0, true, Nd-1 );  // used to check props. Fix to Nd-1 direction.

    LatticePropagator F;
    ChromaProp_t prop_header;
    PropSourceConst_t source_header;
    QDPIO::cout << "Attempt to parse forward propagator" << std::endl;
    QDPIO::cout << "parsing forward propagator " << params.named_obj.prop_id << " ... " << std::endl << std::flush;
    

    try
    {
      // Snarf a copy
      F = TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.prop_id);
	
      // Snarf the frwd prop info. This is will throw if the frwd prop id is not there
      XMLReader PropXML, PropRecordXML;
      TheNamedObjMap::Instance().get(params.named_obj.prop_id).getFileXML(PropXML);
      TheNamedObjMap::Instance().get(params.named_obj.prop_id).getRecordXML(PropRecordXML);

      // Try to invert this record XML into a ChromaProp struct
      {
	read(PropRecordXML, "/Propagator/ForwardProp", prop_header);
	read(PropRecordXML, "/Propagator/PropSource", source_header);
      }

      // Sanity check - write out the norm2 of the forward prop in the j_decay direction
      // Use this for any possible verification
      {
	multi1d<Double> PropCheck = 
	  sumMulti( localNorm2( F ), phases_nomom.getSet() );

	QDPIO::cout << "forward propagator check = " << PropCheck[0] << std::endl;

	// Write out the forward propagator header
	push(XmlOut, "ForwardProp");
	write(XmlOut, "PropXML", PropXML);
	write(XmlOut, "PropRecordXML", PropRecordXML);
	write(XmlOut, "PropCheck", PropCheck);
	pop(XmlOut);
      }
    }
    catch( std::bad_cast ) 
    {
      QDPIO::cerr << InlineNprVertexPDFEnv::name << ": caught dynamic cast error" 
		  << std::endl;
      QDP_abort(1);
    }
    catch (const std::string& e) 
    {
      QDPIO::cerr << InlineNprVertexPDFEnv::name << ": forward prop: error message: " << e 
		  << std::endl;
      QDP_abort(1);
    }

    QDPIO::cout << "Forward propagator successfully parsed" << std::endl;


    // Get the momentum from the header

    multi1d<int> mom  ;
    multi1d<int> t_src ;
    int       t_dir = source_header.j_decay ;

    try{
      mom = source_header.getMom() ;
      t_src = source_header.getTSrce() ;
    }
    catch (const std::string& e){
      QDPIO::cerr << InlineNprVertexPDFEnv::name << ": propagator does not have a momentum source or t_src not present: error message: " << e << std::endl;
      QDP_abort(1);
    }

    QDP::StopWatch swatch;
    swatch.reset();
  if(params.param.save_prop==true)    
  {
    //Fourier transform the propagator 
    QDPIO::cout << "Fourier Transforming propagator" << std::endl;

    swatch.start();
    std::string name_prop=params.param.file_name+".ky_prop";
    writeKYUQprop2(F,name_prop);
    swatch.stop();
    QDPIO::cout << "finished write KY propagator"
                << "  time= "
                << swatch.getTimeInSeconds()
                << " secs" << std::endl;
  }

    LatticeReal theta=zero;
    for(int id=0;id<4;id++) 
        theta+=twopi*Real(mom[id])/Real(Layout::lattSize()[id])* QDP::Layout::latticeCoordinate(id); 
    //need the oposit momentum on the sink
    LatticeComplex phase=cmplx(cos(theta),-sin(theta));
    F=F*phase;//project propagator to the given momentum.
    
    std::string name_prop=params.param.file_name+".ftp";
    write_FT_prop(F,name_prop,params.param.cfg_serial);

//    LatticeComplex phase2=cmplx(cos(theta),sin(theta));
//    F=F*phase;

/*
    sprintf(name,"../gs8t8iwa165_052.pdf.wall.%03d%03d%03d%03d.m0.010000",mom[0],mom[1],mom[2],mom[3]);
    readKYUQprop2(F,name);
    std::string ky_prop_name=std::string(name)+".test";
    writeKYUQprop2(F,ky_prop_name);
    
    name_chroma=params.param.file_name+".g";
    write_FT_prop(F,neg_mom,name_chroma);
*/

    //#################################################################################//
    // Construct Building Blocks                                                       //
    //#################################################################################//
    
    XMLBufferWriter file_xml;
    push(file_xml, "NprVertexPDF");
    write(file_xml, "Param", params.param);
    write(file_xml, "ForwardProp", prop_header);
    write(file_xml, "PropSource", source_header);
    write(file_xml, "Config", gauge_xml);
    pop(file_xml);


    XMLBufferWriter prop_xml;
    push(prop_xml,"QuarkPropagator");
    write(prop_xml,"mom",mom);
    write(prop_xml,"origin",t_src);
    write(prop_xml,"t_dir",t_dir);
    pop(prop_xml) ;

    QDPIO::cout << "Calculating building blocks" << std::endl;
    swatch.reset();
    swatch.start();
	std::string str=params.param.file_name+ ".io";

    NprVertexPDF(F, U, params.param.links_max, params.param.links_prep_max, AllLinkPatterns, str.c_str(),params.param.cfg_serial);
    swatch.stop();
      
    QDPIO::cout << "finished calculating NprVertexPDF"
		<< "  time= "
		<< swatch.getTimeInSeconds() 
		<< " secs" << std::endl;

    pop(XmlOut);   // NprVertexPDF

    snoop.stop();
    QDPIO::cout << InlineNprVertexPDFEnv::name << ": total time = "
		<< snoop.getTimeInSeconds() 
		<< " secs" << std::endl;

    QDPIO::cout << InlineNprVertexPDFEnv::name << ": ran successfully" << std::endl;

    END_CODE();
  } 

}

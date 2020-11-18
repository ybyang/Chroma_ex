/*! \file
 * \brief Inline construction of hadron spectrum
 *
 * Spectrum calculations
 */

#include "inline_hadspec_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "io/param_io.h"
#include "io/qprop_io.h"
#include "meas/hadron/mesons2_w.h"
#include "meas/hadron/barhqlq_w.h"
#include "meas/hadron/curcor2_w.h"
#include "meas/inline/make_xml_file.h"
#include "meas/inline/io/named_objmap.h"
#include "meas/smear/no_quark_displacement.h"
#include "meas/hadron/barspinmat_w.h"
#include "io_ex/io_general_class.h"

namespace Chroma 
{ 
  namespace InlineHadSpec2Env 
  { 
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineHadSpec2(InlineHadSpec2Params(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "HADRON_SPECTRUM_v2";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	registered = true;
      }
      return success;
    }
  }



  //! Reader for parameters
  void read(XMLReader& xml, const std::string& path, InlineHadSpec2Params::Param_t& param)
  {
    XMLReader paramtop(xml, path);

    int version;
    read(paramtop, "version", version);

    switch (version) 
    {
    case 1:
      break;

    default:
      QDPIO::cerr << "Input parameter version " << version << " unsupported." << std::endl;
      QDP_abort(1);
    }

    if (paramtop.count("boundary_fac")>0){
       read(paramtop, "boundary_fac", param.boundary_fac);
    }
    else
      param.boundary_fac=1;
    read(paramtop, "hadron_list", param.hadron_list);
    param.prj_type=0;
    read(paramtop, "prj_type", param.prj_type);
    param.translate=true;
    if (paramtop.count("translate")>0) {
      read(paramtop, "translate", param.translate);
    }
    param.save_vec=false;
    if (paramtop.count("save_vec")>0) {
      read(paramtop, "save_vec", param.save_vec);
      if(param.save_vec==true)
          param.translate=false;
    }
    
  
    if (paramtop.count("mom_list")>0) {
      read(paramtop, "mom_list", param.mom_list);
      param.mom2_max = -1;
    }
    else {
      read(paramtop, "mom2_max", param.mom2_max);
    }
    if (paramtop.count("src")>0) {
      read(paramtop, "src", param.src);
    }
    else {
      param.src=multi1d<int>(3);
      param.src[0] = 0;
      param.src[1] = 0;
      param.src[2] = 0;
    }
    
    read(paramtop, "cfg_serial", param.cfg_serial);

    read(paramtop, "avg_equiv_mom", param.avg_equiv_mom);
  }


  //! Writer for parameters
  void write(XMLWriter& xml, const std::string& path, const InlineHadSpec2Params::Param_t& param)
  {
    push(xml, path);

    int version = 1;
    write(xml, "version", version);

    write(xml, "prj_type", param.prj_type);

    if (param.mom2_max >= 0) {
      write(xml, "mom2_max", param.mom2_max);
    }
    else {
      write(xml, "mom_list", param.mom_list);
    }
    write(xml, "avg_equiv_mom", param.avg_equiv_mom);
    
    write(xml, "cfg_serial", param.cfg_serial);

    pop(xml);
  }


  //! Propagator input
  void read(XMLReader& xml, const std::string& path, InlineHadSpec2Params::NamedObject_t::Props_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "first_id", input.first_id);
    read(inputtop, "second_id", input.second_id);
  }

  //! Propagator output
  void write(XMLWriter& xml, const std::string& path, const InlineHadSpec2Params::NamedObject_t::Props_t& input)
  {
    push(xml, path);

    write(xml, "first_id", input.first_id);
    write(xml, "second_id", input.second_id);

    pop(xml);
  }


  //! Propagator input
  void read(XMLReader& xml, const std::string& path, InlineHadSpec2Params::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "gauge_id", input.gauge_id);
    read(inputtop, "sink_pairs", input.sink_pairs);
  }

  //! Propagator output
  void write(XMLWriter& xml, const std::string& path, const InlineHadSpec2Params::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "gauge_id", input.gauge_id);
    write(xml, "sink_pairs", input.sink_pairs);

    pop(xml);
  }


  // Param stuff
  InlineHadSpec2Params::InlineHadSpec2Params()
  { 
    frequency = 0; 
  }

  InlineHadSpec2Params::InlineHadSpec2Params(XMLReader& xml_in, const std::string& path) 
  {
    try 
    {
      XMLReader paramtop(xml_in, path);

      if (paramtop.count("Frequency") == 1)
	read(paramtop, "Frequency", frequency);
      else
	frequency = 1;

      // Parameters for source construction
      read(paramtop, "Param", param);

      // Read in the output propagator/source configuration info
      read(paramtop, "NamedObject", named_obj);

	read(paramtop, "output", output);
    }
    catch(const std::string& e) 
    {
      QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e << std::endl;
      QDP_abort(1);
    }
  }


  void
  InlineHadSpec2Params::write(XMLWriter& xml_out, const std::string& path) 
  {
    push(xml_out, path);
    
    Chroma::write(xml_out, "Param", param);
    Chroma::write(xml_out, "NamedObject", named_obj);
//    QDP::write(xml_out, "xml_file", xml_file);

    pop(xml_out);
  }



  // Anonymous namespace
  namespace 
  {
    //! Useful structure holding sink props
    struct SinkPropContainer_t
    {
      ForwardProp_t prop_header;
      std::string quark_propagator_id;
      Real Mass;
    
      multi1d<int> bc; 
    
      std::string source_type;
      std::string source_disp_type;
      std::string sink_type;
      std::string sink_disp_type;
    };


    //! Useful structure holding sink props
    struct AllSinkProps_t
    {
      SinkPropContainer_t  sink_prop_1;
      SinkPropContainer_t  sink_prop_2;
    };


    //! Read a sink prop
    void readSinkProp(SinkPropContainer_t& s, const std::string& id)
    {
      try
      {
	// Try a cast to see if it succeeds
	const LatticePropagator& foo = 
	  TheNamedObjMap::Instance().getData<LatticePropagator>(id);

	// Snarf the data into a copy
	s.quark_propagator_id = id;
	
	// Snarf the prop info. This is will throw if the prop_id is not there
	XMLReader prop_file_xml, prop_record_xml;
	TheNamedObjMap::Instance().get(id).getFileXML(prop_file_xml);
	TheNamedObjMap::Instance().get(id).getRecordXML(prop_record_xml);
   
	// Try to invert this record XML into a ChromaProp struct
	// Also pull out the id of this source
	{
	  std::string xpath;
	  read(prop_record_xml, "/SinkSmear", s.prop_header);
	  
	  read(prop_record_xml, "/SinkSmear/PropSource/Source/SourceType", s.source_type);
	  xpath = "/SinkSmear/PropSource/Source/Displacement/DisplacementType";
	  if (prop_record_xml.count(xpath) != 0)
	    read(prop_record_xml, xpath, s.source_disp_type);
	  else
	    s.source_disp_type = NoQuarkDisplacementEnv::getName();

	  read(prop_record_xml, "/SinkSmear/PropSink/Sink/SinkType", s.sink_type);
	  xpath = "/SinkSmear/PropSink/Sink/Displacement/DisplacementType";
	  if (prop_record_xml.count(xpath) != 0)
	    read(prop_record_xml, xpath, s.sink_disp_type);
	  else
	    s.sink_disp_type = NoQuarkDisplacementEnv::getName();
	}
      }
      catch( std::bad_cast ) 
      {
	QDPIO::cerr << InlineHadSpec2Env::name << ": caught dynamic cast error" 
		    << std::endl;
	QDP_abort(1);
      }
      catch (const std::string& e) 
      {
	QDPIO::cerr << InlineHadSpec2Env::name << ": error message: " << e 
		    << std::endl;
	QDP_abort(1);
      }


      // Derived from input prop
      // Hunt around to find the mass
      // NOTE: this may be problematic in the future if actions are used with no
      // clear def. of a Mass
      QDPIO::cout << "Try action and mass" << std::endl;
      s.Mass = getMass(s.prop_header.prop_header.fermact);

      // Only baryons care about boundaries
      // Try to find them. If not present, assume dirichlet.
      // This turns off any attempt to time reverse which is the
      // only thing that the BC are affecting.
      s.bc.resize(Nd);
      s.bc = 0;
    
      try
      {
	s.bc = getFermActBoundary(s.prop_header.prop_header.fermact);
      }
      catch (const std::string& e) 
      {
	QDPIO::cerr << InlineHadSpec2Env::name 
		    << ": caught exception. No BC found in these headers. Will assume dirichlet: " << e 
		    << std::endl;
      }

      QDPIO::cout << "FermAct = " << s.prop_header.prop_header.fermact.id << std::endl;
      QDPIO::cout << "Mass = " << s.Mass << std::endl;
    }


    //! Read all sinks
    void readAllSinks(AllSinkProps_t& s, 
		      InlineHadSpec2Params::NamedObject_t::Props_t sink_pair)
    {
      QDPIO::cout << "Attempt to parse forward propagator = " << sink_pair.first_id << std::endl;
      readSinkProp(s.sink_prop_1, sink_pair.first_id);
      QDPIO::cout << "Forward propagator successfully parsed" << std::endl;

      QDPIO::cout << "Attempt to parse forward propagator = " << sink_pair.second_id << std::endl;
      readSinkProp(s.sink_prop_2, sink_pair.second_id);
      QDPIO::cout << "Forward propagator successfully parsed" << std::endl;
    }

  } // namespace anonymous



  // Function call
  void 
  InlineHadSpec2::operator()(unsigned long update_no,
			    XMLWriter& xml_out) 
  {
    // If xml file not empty, then use alternate
    if (params.xml_file != "")
    {
      std::string xml_file = makeXMLFileName(params.xml_file, update_no);

      push(xml_out, "hadspec");
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

/*  
  class QDPFileWriter_iog:public QDPFileWriter
  {
    public:
//       XMLBufferWriter file_xml;
       DML_Checksum checksum;
       uint64_t nbytes;
       int msg_begin, msg_end;
       
//       QDPFileWriter_iog(const std::string& path, QDP_volfmt_t qdp_volfmt,
       QDPFileWriter_iog( XMLBufferWriter & file_xml,const std::string& path, QDP_volfmt_t qdp_volfmt,
                   QDP_serialparallel_t qdp_serpar)
//                   :QDPFileWriter(file_xml,path,qdp_volfmt,qdp_serpar){}
                   :QDPFileWriter(file_xml,path,qdp_volfmt,qdp_serpar,QDPIO_OPEN){}
       
       void write_head(datatype &res_vec)
       {
            char scidac_type[] = QIO_LIMETYPE_BINARY_DATA; 
            LIME_type lime_type=scidac_type;                                            
            int status=QIO_write_field(get(),msg_begin,msg_end,&(QDPOScalarFactoryGet<char>),1,
                                         N_MAX_HEADLENGTH*sizeof(char),sizeof(char),res_vec.type.blank,
                                         &checksum, &nbytes,lime_type);
       }
       
       template<class T>
       void write_data(const OLattice<T>& s1)
       {
            if(get.layout.this_node==get.layout.master_io_node)
         
            char scidac_type[] = QIO_LIMETYPE_BINARY_DATA; 
            LIME_type lime_type=scidac_type;                                            
            int status=QIO_write_field(get(),msg_begin,msg_end,&(QDPOLatticeFactoryGet<T>),1,
                                         sizeof(T),sizeof(typename WordType<T>::Type_t),(void *)s1.getF(),
                                         &checksum, &nbytes,lime_type);
       }
  
       template<class T>
       void write_data(const multi1d< OLattice<T> >& s1)
       {
            char scidac_type[] = QIO_LIMETYPE_BINARY_DATA; 
            LIME_type lime_type=scidac_type;                                            
            int status=QIO_write_field(get(),msg_begin,msg_end,&(QDPOLatticeFactoryGet<T>),s1.size(),
                                         sizeof(T),sizeof(typename WordType<T>::Type_t),(void *)&s1,
                                         &checksum, &nbytes,lime_type);
       }
  
  
  };
*/

  // Real work done here
  void 
  InlineHadSpec2::func(unsigned long update_no,
		      XMLWriter& xml_out) 
  {
    START_CODE();

    StopWatch snoop;
    snoop.reset();
    snoop.start();

    StopWatch Tparts;
    Tparts.reset();
    Tparts.start();

    // Test and grab a reference to the gauge field
    XMLBufferWriter gauge_xml;
    try
    {
      TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);
      TheNamedObjMap::Instance().get(params.named_obj.gauge_id).getRecordXML(gauge_xml);
    }
    catch( std::bad_cast ) 
    {
      QDPIO::cerr << InlineHadSpec2Env::name << ": caught dynamic cast error" 
		  << std::endl;
      QDP_abort(1);
    }
    catch (const std::string& e) 
    {
      QDPIO::cerr << InlineHadSpec2Env::name << ": std::map call failed: " << e 
		  << std::endl;
      QDP_abort(1);
    }
    const multi1d<LatticeColorMatrix>& u = 
      TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);

    push(xml_out, "hadspec");
    write(xml_out, "update_no", update_no);

    QDPIO::cout << " HADSPEC: Spectroscopy for Wilson-like fermions" << std::endl;
    QDPIO::cout << std::endl << "     Gauge group: SU(" << Nc << ")" << std::endl;
    QDPIO::cout << "     volume: " << Layout::lattSize()[0];
    for (int i=1; i<Nd; ++i) {
      QDPIO::cout << " x " << Layout::lattSize()[i];
    }
    QDPIO::cout << std::endl;

    proginfo(xml_out);    // Print out basic program info

    // Write out the input
    params.write(xml_out, "Input");

    // Write out the config info
    write(xml_out, "Config_info", gauge_xml);

    push(xml_out, "Output_version");
    write(xml_out, "out_version", 15);
    pop(xml_out);


    // First calculate some gauge invariant observables just for info.
    MesPlq(xml_out, "Observables", u);

    // Keep an array of all the xml output buffers
    push(xml_out, "Wilson_hadron_measurements");
   
    Tparts.stop();
    QDPIO::cout << InlineHadSpec2Env::name << ": setup time = "
		<< Tparts.getTimeInSeconds() 
		<< " secs" << std::endl;
    Tparts.reset();
    Tparts.start();

    std::vector<int> mom_serial;
    int mom_max=5;
    const Real twopi = 6.283185307179586476925286;
    if(params.param.mom2_max>=0)
    {
        int mom_max=sqrt(params.param.mom2_max)+1;
        int count=0;
        for(int i=-mom_max;i<=mom_max;i++)
        for(int j=-mom_max;j<=mom_max;j++)
        for(int k=-mom_max;k<=mom_max;k++)
        if(i*i+j*j+k*k<=params.param.mom2_max)
           mom_serial.push_back((50+i)*10000+(50+j)*100+50+k);
    }
    else
    for(int ik=0;ik<params.param.mom_list.size();ik++)
       mom_serial.push_back(params.param.mom_list[ik]);
   
    if(mom_serial.size()>1024) mom_serial.resize(1024);
    
    multi1d<LatticeComplex> phase(mom_serial.size());
    
    //to do: set the mom_serial and phase;

    // initialize the output data;    
    std::vector<int> hadron;
    for(int iop=0;iop<params.param.hadron_list.size();iop++)
       hadron.push_back(params.param.hadron_list[iop]);
    general_data_base res(params.output.c_str());
    int offset=0;
    
    bool initialized=false;
    BinaryFileWriter *bin;
    if(params.param.save_vec==true) bin=new BinaryFileWriter(params.output+".vec");

/*    XMLBufferWriter file_xml;
    QDPFileWriter_iog *bin;
    if(params.param.save_vec==true) bin=new QDPFileWriter_iog(file_xml,params.output+".vec",QDPIO_SINGLEFILE,QDPIO_PARALLEL);*/
    
    Tparts.stop();
    QDPIO::cout << InlineHadSpec2Env::name << ": set param time = "
		<< Tparts.getTimeInSeconds() 
		<< " secs" << std::endl;
    Tparts.reset();
    Tparts.start();

    // Now loop over the various fermion pairs
    for(int lpair=0; lpair < params.named_obj.sink_pairs.size(); ++lpair)
    {
      const InlineHadSpec2Params::NamedObject_t::Props_t named_obj = params.named_obj.sink_pairs[lpair];

      push(xml_out, "elem");

      AllSinkProps_t all_sinks;
      readAllSinks(all_sinks, named_obj);

      // Derived from input prop
      multi1d<int> t_srce
                  = all_sinks.sink_prop_1.prop_header.source_header.getTSrce();
      int j_decay = all_sinks.sink_prop_1.prop_header.source_header.j_decay;
      int t0      = all_sinks.sink_prop_1.prop_header.source_header.t_source;
      if(initialized==false)
      {
        for(int ik=0;ik<mom_serial.size();ik++)
        {
           int p[3];
//           p[0]=mom_serial[ik]%100-50;
//           p[1]=(mom_serial[ik]/100)%100-50;
//           p[2]=(mom_serial[ik]/10000)%100-50;
           p[2]=mom_serial[ik]%100-50;
           p[1]=(mom_serial[ik]/100)%100-50;
           p[0]=(mom_serial[ik]/10000)%100-50;
           LatticeReal p_dot_x=0.;
           const Real twopi = 6.283185307179586476925286;
           for(int mu=0,j=0;mu<Nd;++mu)
           if(mu!=j_decay)
           {
                p_dot_x += LatticeReal(Layout::latticeCoordinate(mu) - t_srce[mu]) * twopi * 
                     p[j]/Layout::lattSize()[mu];
                ++j;
           }
           phase[ik]=cmplx(cos(p_dot_x),sin(p_dot_x));
        }

        res.add_dimension(dim_conf, 1, &params.param.cfg_serial);
        res.add_dimension(dim_t2, 1,&t0);
        int isp=params.param.src[2]+1000*(params.param.src[1]+1000*params.param.src[0]);
        res.add_dimension(dim_temporary, 1,&isp);
        res.add_dimension(dim_channel, params.named_obj.sink_pairs.size());
        res.add_dimension(dim_operator,hadron.size(),hadron.data());
        
        general_data_base res_vec;
        res_vec.type=res.type;

        res.add_dimension(dim_momentum,mom_serial.size(),mom_serial.data());
        res.add_dimension(dim_t,Layout::lattSize()[j_decay]);
        res.add_dimension(dim_complex, 2);
        if(Layout::primaryNode()) res.initialize();

        res_vec.add_dimension(dim_t,Layout::lattSize()[3]);
        res_vec.add_dimension(dim_z,Layout::lattSize()[2]);
        res_vec.add_dimension(dim_y,Layout::lattSize()[1]);
        res_vec.add_dimension(dim_x,Layout::lattSize()[0]);
        res_vec.add_dimension(dim_complex, 2);
        
        if(params.param.save_vec==true) bin->writeArray(res_vec.type.blank,1,N_MAX_HEADLENGTH);
//        if(params.param.save_vec==true) bin->write_head(res_vec);
        initialized=true;
      }
    Tparts.stop();
    QDPIO::cout << InlineHadSpec2Env::name << ": setting data time = "
		<< Tparts.getTimeInSeconds() 
		<< " secs" << std::endl;
    Tparts.reset();
    Tparts.start();

      // Sanity checks
      {
	if (all_sinks.sink_prop_2.prop_header.source_header.j_decay != j_decay)
	{
	  QDPIO::cerr << "Error!! j_decay must be the same for all propagators " << std::endl;
	  QDP_abort(1);
	}
	if (all_sinks.sink_prop_2.prop_header.source_header.t_source != 
	    all_sinks.sink_prop_1.prop_header.source_header.t_source)
	{
	  QDPIO::cerr << "Error!! t_source must be the same for all propagators " << std::endl;
	  QDP_abort(1);
	}
	if (all_sinks.sink_prop_1.source_type != all_sinks.sink_prop_2.source_type)
	{
	  QDPIO::cerr << "Error!! source_type must be the same in a pair " << std::endl;
	  QDP_abort(1);
	}
	if (all_sinks.sink_prop_1.sink_type != all_sinks.sink_prop_2.sink_type)
	{
	  QDPIO::cerr << "Error!! source_type must be the same in a pair " << std::endl;
	  QDP_abort(1);
	}
      }

      int bc_spec = 0;
	bc_spec = all_sinks.sink_prop_1.bc[j_decay] ;
	if (all_sinks.sink_prop_2.bc[j_decay] != bc_spec)
	{
	  QDPIO::cerr << "Error!! bc must be the same for all propagators " << std::endl;
	  QDP_abort(1);
	}

      // Initialize the slow Fourier transform phases
//      SftMom phases(params.param.mom2_max, t_srce, params.param.avg_equiv_mom,
//                    j_decay);

      // Keep a copy of the phases with NO momenta
      SftMom phases_nomom(0, true, j_decay);
      Set timeslice= phases_nomom.getSet();

      // Masses
      write(xml_out, "Mass_1", all_sinks.sink_prop_1.Mass);
      write(xml_out, "Mass_2", all_sinks.sink_prop_2.Mass);
      write(xml_out, "t0", t0);

      // Save prop input
      push(xml_out, "Forward_prop_headers");
      write(xml_out, "First_forward_prop", all_sinks.sink_prop_1.prop_header);
      write(xml_out, "Second_forward_prop", all_sinks.sink_prop_2.prop_header);
      pop(xml_out);

      // Sanity check - write out the norm2 of the forward prop in the j_decay direction
      // Use this for any possible verification
/*
      push(xml_out, "Forward_prop_correlator");
      {
	const LatticePropagator& sink_prop_1 = 
	  TheNamedObjMap::Instance().getData<LatticePropagator>(all_sinks.sink_prop_1.quark_propagator_id);
	const LatticePropagator& sink_prop_2 = 
	  TheNamedObjMap::Instance().getData<LatticePropagator>(all_sinks.sink_prop_2.quark_propagator_id);

	write(xml_out, "forward_prop_corr_1", sumMulti(localNorm2(sink_prop_1), phases.getSet()));
	write(xml_out, "forward_prop_corr_2", sumMulti(localNorm2(sink_prop_2), phases.getSet()));
      }
      pop(xml_out);
*/

      push(xml_out, "SourceSinkType");
      {
	QDPIO::cout << "Source_type_1 = " << all_sinks.sink_prop_1.source_type << std::endl;
	QDPIO::cout << "Sink_type_1 = " << all_sinks.sink_prop_1.sink_type << std::endl;
	QDPIO::cout << "Source_type_2 = " << all_sinks.sink_prop_2.source_type << std::endl;
	QDPIO::cout << "Sink_type_2 = " << all_sinks.sink_prop_2.sink_type << std::endl;

	write(xml_out, "source_type_1", all_sinks.sink_prop_1.source_type);
	write(xml_out, "source_disp_type_1", all_sinks.sink_prop_1.source_disp_type);
	write(xml_out, "sink_type_1", all_sinks.sink_prop_1.sink_type);
	write(xml_out, "sink_disp_type_1", all_sinks.sink_prop_1.sink_disp_type);

	write(xml_out, "source_type_2", all_sinks.sink_prop_2.source_type);
	write(xml_out, "source_disp_type_2", all_sinks.sink_prop_2.source_disp_type);
	write(xml_out, "sink_type_2", all_sinks.sink_prop_2.sink_type);
	write(xml_out, "sink_disp_type_2", all_sinks.sink_prop_2.sink_disp_type);
      }
      pop(xml_out);


      // References for use later
      const LatticePropagator& sink_prop_1 = 
	TheNamedObjMap::Instance().getData<LatticePropagator>(all_sinks.sink_prop_1.quark_propagator_id);
      const LatticePropagator& sink_prop_2 = 
	TheNamedObjMap::Instance().getData<LatticePropagator>(all_sinks.sink_prop_2.quark_propagator_id);


      // Construct group name for output
      std::string src_type;
      if (all_sinks.sink_prop_1.source_type == "POINT_SOURCE")
	src_type = "Point";
      else if (all_sinks.sink_prop_1.source_type == "SF_POINT_SOURCE")
	src_type = "Point";
      else if (all_sinks.sink_prop_1.source_type == "SHELL_SOURCE")
	src_type = "Shell";
      else if (all_sinks.sink_prop_1.source_type == "NORM_SHELL_SOURCE")
	src_type = "Shell";
      else if (all_sinks.sink_prop_1.source_type == "SF_SHELL_SOURCE")
	src_type = "Shell";
      else if (all_sinks.sink_prop_1.source_type == "WALL_SOURCE")
	src_type = "Wall";
      else if (all_sinks.sink_prop_1.source_type == "SF_WALL_SOURCE")
	src_type = "Wall";
      else if (all_sinks.sink_prop_1.source_type == "RAND_ZN_WALL_SOURCE")
	src_type = "Wall";
      else if (all_sinks.sink_prop_1.source_type == "MOM_GRID_SOURCE")
	src_type = "grid";
      else
      {
	QDPIO::cerr << "Unsupported source type = " << all_sinks.sink_prop_1.source_type << std::endl;
	QDP_abort(1);
      }

      std::string snk_type;
      if (all_sinks.sink_prop_1.sink_type == "POINT_SINK")
	snk_type = "Point";
      else if (all_sinks.sink_prop_1.sink_type == "SHELL_SINK")
	snk_type = "Shell";
      else if (all_sinks.sink_prop_1.sink_type == "NORM_SHELL_SINK")
	snk_type = "Shell";
      else if (all_sinks.sink_prop_1.sink_type == "WALL_SINK")
	snk_type = "Wall";
      else
      {
	QDPIO::cerr << "Unsupported sink type = " << all_sinks.sink_prop_1.sink_type << std::endl;
	QDP_abort(1);
      }

      std::string source_sink_type = src_type + "_" + snk_type;
      QDPIO::cout << "Source type = " << src_type << std::endl;
      QDPIO::cout << "Sink type = "   << snk_type << std::endl;
      
    Tparts.stop();
    QDPIO::cout << InlineHadSpec2Env::name << ": prepare time = "
		<< Tparts.getTimeInSeconds() 
		<< " secs" << std::endl;
    Tparts.reset();
    Tparts.start();
      
      //do the contraction;
      LatticePropagator anti_sink_prop =  Gamma(15) * sink_prop_2 * Gamma(15);
      SpinMatrix Cg5NR = BaryonSpinMats::Cg5NR();
      if(params.param.prj_type>0) Cg5NR = BaryonSpinMats::Cg5();//patch for the folded projection.
      //sink_prop_1: d type sink_prop_2: u type
      LatticePropagator di_quark = quarkContract13(sink_prop_1 * Cg5NR,
                                                Cg5NR * sink_prop_2);
      if(params.param.prj_type>0) di_quark=0.5*di_quark;
      
      multi1d<SpinMatrix> projs(4);
      SpinMatrix g_one = 1.0, g_mone = -1.0;
      projs[0] = Gamma(0 )*g_one; 
      projs[1] = Gamma(14)*g_one;
      projs[2] = Gamma(13)*g_mone;
      projs[3] = Gamma(11)*g_one;
      LatticeSpinMatrix T_unpol = BaryonSpinMats::Tunpol();
      if(params.param.prj_type>0)
      {
          int nt_r=Layout::lattSize()[3]/params.param.prj_type;
          SpinMatrix prj_p(0.5 * (g_one + (g_one * Gamma(8))));
          SpinMatrix prj_m(0.5 * (g_one - (g_one * Gamma(8))));
          T_unpol =
             where(((Layout::latticeCoordinate(3))-t0+Layout::lattSize()[3])%nt_r<nt_r/2,
                 prj_p,prj_m);
          QDPIO::cout << "prj_type = " << params.param.prj_type
                      << "t=" << t0 << std::endl; 
      }
      
      
    Tparts.stop();
    QDPIO::cout << InlineHadSpec2Env::name << ": checking, di_quark and projection time = "
		<< Tparts.getTimeInSeconds() 
		<< " secs" << std::endl;
    Tparts.reset();
    Tparts.start();

      for(int iop=0;iop<hadron.size();iop++)
      {
           LatticeComplex corr=zero;
           if(hadron[iop]>=100000000)//1e8, do meson
           {
                int ig1=(hadron[iop]/100)%100;
                int ig2=hadron[iop]%100;
                if(ig1<0||ig1>15||ig2<0||ig2>15) continue;
                corr = trace(adj(anti_sink_prop)*Gamma(ig1)*sink_prop_1*Gamma(ig2));
           }
           else // do bayron 
           {                
              int iprj=(hadron[iop]/10000)%10;
              if(hadron[iop]==200505||hadron[iop]==210505||hadron[iop]==220505||hadron[iop]==230505||
                 hadron[iop]==203131||hadron[iop]==213131||hadron[iop]==223131||hadron[iop]==233131) // with Cg=ga1.ga3
                  corr=LatticeComplex(trace(projs[iprj]* T_unpol * traceColor(sink_prop_2 * traceSpin(di_quark)))
                                              + trace(projs[iprj]* T_unpol * traceColor(sink_prop_2 * di_quark)));
              if(hadron[iop]==6)// mixed prj nuecleon
              {
                  SpinMatrix T_mixed = BaryonSpinMats::Tmixed();
                  corr=LatticeComplex(trace(T_mixed * traceColor(sink_prop_2 * traceSpin(di_quark)))
                                              + trace(T_mixed * traceColor(sink_prop_2 * di_quark)));
              }
              if(hadron[iop]==240505)// fully mixed prj nuecleon
              {
                  SpinMatrix T_mixed(0.5*(g_one + Gamma(8)*g_one + timesMinusI(Gamma(3)*g_one + Gamma(6)*g_one + Gamma(11)*g_one + Gamma(14)*g_one)
                                       +timesI(Gamma(5)*g_one + Gamma(13)*g_one)));
                  corr=LatticeComplex(trace(T_mixed * traceColor(sink_prop_2 * traceSpin(di_quark)))
                                              + trace(T_mixed * traceColor(sink_prop_2 * di_quark)));
              }

              if(hadron[iop]==2)// omega and delta
              {
                  SpinMatrix T_mixed = BaryonSpinMats::Tmixed();
                  corr=Baryon2PtContractions::sigmast2pt(sink_prop_1, sink_prop_2,T_mixed, BaryonSpinMats::Cgm());
              }
              if(hadron[iop]==5)// delta ver.2
              {
                  SpinMatrix T_mixed = BaryonSpinMats::Tmixed();
                  corr=Baryon2PtContractions::sigmast2pt(sink_prop_1, sink_prop_2,T_mixed, BaryonSpinMats::Cg4m());
              }
                  
           }
           int nt=Layout::lattSize()[j_decay];
           for(int ik=0;ik<mom_serial.size();ik++)
           {
               multi1d<DComplex> hsum = sumMulti(phase[ik]*corr,timeslice);
               if(Layout::primaryNode())
               for (int t=0; t < nt; ++t) 
               {
                  int t_r=t;
                  if(params.param.translate==true)t_r=(t-t_srce[j_decay]+nt)%nt;
                  int fac=(t<t_srce[j_decay]&&hadron[iop]<100000000)?params.param.boundary_fac:1;
                  res.data[offset*2*nt+2*t_r]=hsum[t].elem().elem().elem().real()*fac;
                  res.data[offset*2*nt+2*t_r+1]=hsum[t].elem().elem().elem().imag()*fac;
               }
               offset++;
           }
    Tparts.stop();
    QDPIO::cout << InlineHadSpec2Env::name << ": calculation time = "
		<< Tparts.getTimeInSeconds() 
		<< " secs" << std::endl;
    Tparts.reset();
    Tparts.start();
           // save corr;
           if(params.param.save_vec==true)
           {
              LatticeComplexD corrD;
              for(int s(0); s < Layout::sitesOnNode(); s++)
              {
                      corrD.elem(s).elem()=corr.elem(s).elem();
                      QDPUtil::byte_swap((void *)&(corrD.elem(s).elem()),sizeof(RealD),2);
              }
              write(*bin,corrD);
//              bin->write_data(corrD);
           }
    Tparts.stop();
    QDPIO::cout << InlineHadSpec2Env::name << ": IO time = "
		<< Tparts.getTimeInSeconds() 
		<< " secs" << std::endl;
    Tparts.reset();
    Tparts.start();
      }
      
      pop(xml_out);  // array element
    }
    if(params.param.save_vec==true) 
    {  
       bin->close();
       delete bin;
    }
    if(Layout::primaryNode()) res.save();
    Tparts.stop();
    QDPIO::cout << InlineHadSpec2Env::name << ": iog_save time = "
		<< Tparts.getTimeInSeconds() 
		<< " secs" << std::endl;
    Tparts.reset();
    Tparts.start();
    
    if(Layout::primaryNode()&&params.param.save_vec==true)
    {
       std::string str=(params.output+".vec");
       general_data_base res(str.c_str());
       res.load();
       res.sum(res,dim_x);
       res.sum(res,dim_y);
       res.sum(res,dim_z);
       res.print_all();
    }
    Tparts.stop();
    QDPIO::cout << InlineHadSpec2Env::name << ": clear time = "
		<< Tparts.getTimeInSeconds() 
		<< " secs" << std::endl;
    Tparts.reset();
    Tparts.start();
    
    pop(xml_out);  // Wilson_spectroscopy
    pop(xml_out);  // hadspec

    snoop.stop();
    QDPIO::cout << InlineHadSpec2Env::name << ": total time = "
		<< snoop.getTimeInSeconds() 
		<< " secs" << std::endl;

    QDPIO::cout << InlineHadSpec2Env::name << ": ran successfully" << std::endl;

    END_CODE();
  } 

};

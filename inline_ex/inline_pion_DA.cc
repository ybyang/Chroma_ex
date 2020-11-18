// $Id: inline_pion_decay.cc, v1.2 2015-07-05 22:29:00 sdcohen Exp $
/*! \file
 * \brief Inline construction of pion-decay two-point correlators
 *
 * Calculates the two-point correlator of the pseudoscalar meson corresponding
 * to leptonic decay by temporal-axial current.
 *
 * v1.2 adds input of momentum using a list
 * v1.1 adds multiple-source location trickery
 */
 
#include "inline_ex/inline_pion_DA.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "io/param_io.h"
#include "io/qprop_io.h"
#include "meas/inline/make_xml_file.h"
#include "meas/inline/io/named_objmap.h"
#include "meas/smear/no_quark_displacement.h"
#include "io_ex/io_general_class.h"

namespace Chroma
{

//------------------------------------------------------------------------------
// Factory Registration

  // Environment in which the measurement lives (holds params and such)
  namespace InlinePionDAEnv
  {
    namespace
    {
      // Function to register with a factory
      AbsInlineMeasurement* createMeasurement( XMLReader& xml_in,
                                       const std::string& path )
      {
        return new InlinePionDA( InlinePionDAParams(xml_in, path) );
      }

      //! Local registration flag
      bool registered = false;
    }
    
    // The name of the measurement for the XML file
    const std::string name = "PION_DA";

    //! Register all the factories
    bool registerAll()
    {
      bool success = true;
      if (!registered)
      {
        success &= TheInlineMeasurementFactory::Instance()
            .registerObject(name, createMeasurement);
        registered = true;
      }
      return success;
    }
  } // end namespace InlinePionDAEnv

//------------------------------------------------------------------------------
// Parameter reading, writing, handling

  //! Reader for parameters
  void read(             XMLReader& xml,
                 const std::string& path,
    InlinePionDAParams::Param_t& param )
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "mom", param.mom);
    read(paramtop, "links_dir", param.links_dir);
    read(paramtop, "links_max", param.links_max);
    read(paramtop, "file_name", param.file_name);
    read(paramtop, "cfg_serial", param.cfg_serial);
    
      read(paramtop, "Src", param.src);

  }

  //! Writer for parameters
  void write(                  XMLWriter& xml,
                       const std::string& path,
    const InlinePionDAParams::Param_t& param )
  {
    push(xml, path);

    write(xml, "mom", param.mom);
    write(xml, "links_dir", param.links_dir);
    write(xml, "links_max", param.links_max);
    write(xml, "file_name", param.file_name);
    write(xml, "cfg_serial", param.cfg_serial);

      write(xml, "Src", param.src);

    pop(xml);
  }

  //! Named objects (gauge config, prop pairs) input
  void read(                   XMLReader& xml,
                       const std::string& path,
    InlinePionDAParams::NamedObject_t& input )
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "gauge_id", input.gauge_id);
    read(inputtop, "quark_props", input.quark_props);
    read(inputtop, "antiquark_props", input.antiquark_props);
  }

  //! Named objects (gauge config, prop pairs) output
  void write(                        XMLWriter& xml,
                             const std::string& path,
    const InlinePionDAParams::NamedObject_t& input )
  {
    push(xml, path);

    write(xml, "gauge_id", input.gauge_id);
    write(xml, "quark_props", input.quark_props);
    write(xml, "antiquark_props", input.antiquark_props);

    pop(xml);
  }

  // Params default constructor
  InlinePionDAParams::InlinePionDAParams(void)
  {
    frequency = 0;
  }

  // Construct params from XML
  InlinePionDAParams::InlinePionDAParams( XMLReader& xml_in,
                                        const std::string& path )
  {
    try
    {
      XMLReader paramtop(xml_in, path);

      if (paramtop.count("Frequency") == 1)
        read(paramtop, "Frequency", frequency);
      else
        frequency = 1;

      // Read in the parameters
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
      QDPIO::cerr << __func__ << ": Caught Exception reading XML: "
                  << e << std::endl;
      QDP_abort(1);
    }
  }

  // Write out the parameters we constructed
  void InlinePionDAParams::write( XMLWriter& xml_out,
                             const std::string& path )
  {
    push(xml_out, path);

    Chroma::write(xml_out, "Param", param);
    Chroma::write(xml_out, "NamedObject", named_obj);
    QDP::write(xml_out, "xml_file", xml_file);

    pop(xml_out);
  }

//------------------------------------------------------------------------------
// Helper functions to read in pairs of propagators

  // Anonymous namespace for propagator-handling utilities
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


    //! Read a sink prop
    void readSinkProp( SinkPropContainer_t& s,
                         const std::string& id )
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
        std::string xpath;
        read(prop_record_xml, "/SinkSmear", s.prop_header);

        xpath = "/SinkSmear/PropSource/Source/SourceType";
        read(prop_record_xml, xpath, s.source_type);
        xpath = "/SinkSmear/PropSource/Source/Displacement/DisplacementType";
        if (prop_record_xml.count(xpath) != 0)
          read(prop_record_xml, xpath, s.source_disp_type);
        else
          s.source_disp_type = NoQuarkDisplacementEnv::getName();
        
        xpath = "/SinkSmear/PropSink/Sink/SinkType";
        read(prop_record_xml, xpath, s.sink_type);
        xpath = "/SinkSmear/PropSink/Sink/Displacement/DisplacementType";
        if (prop_record_xml.count(xpath) != 0)
          read(prop_record_xml, xpath, s.sink_disp_type);
        else
          s.sink_disp_type = NoQuarkDisplacementEnv::getName();
      }
      catch( std::bad_cast )
      {
        QDPIO::cerr << InlinePionDAEnv::name << ": caught dynamic cast error"
                    << std::endl;
        QDP_abort(1);
      }
      catch (const std::string& e)
      {
        QDPIO::cerr << InlinePionDAEnv::name << ": error message: " << e
                    << std::endl;
        QDP_abort(1);
      }

      // Derived from input prop
      // Hunt around to find the mass
      // NOTE: This may be problematic in the future if actions are used with no
      // clear definition of a mass
      QDPIO::cout << "Try action and mass" << std::endl;
      s.Mass = getMass(s.prop_header.prop_header.fermact);

      // Only baryons care about boundaries
      // Try to find them. If not present, assume Dirichlet.
      // This turns off any attempt to time reverse which is the
      // only thing that the boundary conditions are affecting.
      s.bc.resize(Nd);
      s.bc = 0;

      try
      {
        s.bc = getFermActBoundary(s.prop_header.prop_header.fermact);
      }
      catch (const std::string& e)
      {
        QDPIO::cerr << InlinePionDAEnv::name
                    << ": caught exception. No BC found in these headers. "
                    << "Will assume dirichlet: " << e << std::endl;
      }

      QDPIO::cout << "FermAct = " << s.prop_header.prop_header.fermact.id << std::endl;
      QDPIO::cout << "Mass = " << s.Mass << std::endl;
    }

  } // end namespace

//------------------------------------------------------------------------------
// The real work is done here

  // Set up the XML and invoke func, which does the actual work
  void InlinePionDA::operator()( unsigned long update_no,
                                       XMLWriter& xml_out )
  {
    // If xml file not empty, then use alternate
    if (params.xml_file != "")
    {
      std::string xml_file = makeXMLFileName(params.xml_file, update_no);

      push(xml_out, "PionDA");
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

  // Real work done here
  void InlinePionDA::func( unsigned long update_no,
                                 XMLWriter& xml_out)
  {
    START_CODE();

    StopWatch snoop;
    snoop.reset();
    snoop.start();

    // Grab gauge configuration
    multi1d<LatticeColorMatrix> U;
    XMLBufferWriter gauge_xml;
    try
    {
      // Try to get the gauge field.
      // If it doesn't exist, an exception will be thrown.
      U=TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);
      TheNamedObjMap::Instance().get(params.named_obj.gauge_id).getRecordXML(gauge_xml);
    }
    catch(std::bad_cast)
    {
      // If the object exists and is of the wrong type we will land here
      QDPIO::cerr << InlinePionDAEnv::name << ": caught dynamic cast error"
                  << std::endl;
      QDP_abort(1);
    }
    catch (const std::string& e)
    {
      // If the object is not in the map we will land here
      QDPIO::cerr << InlinePionDAEnv::name << ": map call failed: " << e
                  << std::endl;
      QDP_abort(1);
    }
    // Bind gauge field reference to a local name
    const multi1d<LatticeColorMatrix>& u =
      TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);

    // Boilerplate stuff to the output XML
    push(xml_out, "PionDA");
    write(xml_out, "update_no", update_no);

    // Write some diagnostic junk about the lattice
    QDPIO::cout << " PionDA: Pion decay correlators for Wilson" << std::endl;
    QDPIO::cout << std::endl << "     Gauge group: SU(" << Nc << ")" << std::endl;
    QDPIO::cout << "     volume: " << Layout::lattSize()[0];
    for (int i=1; i<Nd; ++i) {
      QDPIO::cout << " x " << Layout::lattSize()[i];
    }
    QDPIO::cout << std::endl;

    // Write info about the program
    proginfo(xml_out);

    // Write out the input
    params.write(xml_out, "Input");

    // Write out the config info
    write(xml_out, "Config_info", gauge_xml);

    // Calculate the plaquette. Why not? It's always fun!
    // MesPlq(xml_out, "Observables", u);

    // Keep an array of all the xml output buffers
    push(xml_out, "pion_decay_measurements");


    multi1d<SinkPropContainer_t> quark_props(params.named_obj.quark_props.size());
    multi1d<SinkPropContainer_t> antiquark_props(params.named_obj.antiquark_props.size());
    for(int lnum=0; lnum < params.named_obj.quark_props.size(); ++lnum) readSinkProp(quark_props[lnum], params.named_obj.quark_props[lnum]);
    for(int rnum=0; rnum < params.named_obj.antiquark_props.size(); ++rnum) readSinkProp(antiquark_props[rnum], params.named_obj.antiquark_props[rnum]);

    multi1d<int> xsrc
                = quark_props[0].prop_header.source_header.getTSrce();
    int j_decay = quark_props[0].prop_header.source_header.j_decay;
    int t0      = quark_props[0].prop_header.source_header.t_source;

    SftMom *p_phases;  // Pointer for phases

/// set the phase.
    // Keep a copy of the phases with NO momenta
    QDPIO::cout << "J_decay=" << j_decay <<  std::endl;
    SftMom phases_nomom(0, true, j_decay);
    Set timeslice= phases_nomom.getSet();
    int tlen = Layout::lattSize()[j_decay];
    
    LatticeReal p_dot_x;
    p_dot_x = 0.;
    for (int mu = 0, j=0; mu < Nd; ++mu) {
       const Real twopi = 6.283185307179586476925286;
       if (mu == j_decay) continue;
       p_dot_x += LatticeReal(Layout::latticeCoordinate(mu) - params.param.src[mu]) * twopi *
            Real(params.param.mom[j]) / Layout::lattSize()[mu];
       ++j;
    }
    LatticeComplex  phase=cmplx(cos(p_dot_x),sin(p_dot_x));
/// set the phase..

   general_data_base res(params.param.file_name.c_str());
   res.add_dimension(dim_conf, 1, &params.param.cfg_serial);
   res.add_dimension(dim_mass,params.named_obj.quark_props.size());
   res.add_dimension(dim_mass2,params.named_obj.antiquark_props.size());
   res.add_dimension(dim_displacement,params.param.links_max+1);
   res.add_dimension(dim_operator,2);
   res.add_dimension(dim_t,tlen);
   res.add_dimension(dim_complex, 2);
   if(Layout::primaryNode()) res.initialize();
   
    // Now loop over the various fermion pairs
    for(int lnum=0; lnum < params.named_obj.quark_props.size(); ++lnum)
    {
      // Give propagators names for use later
      const LatticePropagator& sink_prop_1 =TheNamedObjMap::Instance().getData<LatticePropagator>(quark_props[lnum].quark_propagator_id);

	// Setting the staring position of data writing	
	int offset=0;
	long data_separation=(params.param.links_max+1)*2*tlen*2;
	long data_start_position=lnum*params.named_obj.antiquark_props.size()*data_separation;
      // Do the contractions
      int LinksDirection=params.param.links_dir;
      if(params.param.links_dir<0 || params.param.links_dir>2 ) LinksDirection=2;
      int gamma_value0=11;
      if(LinksDirection==0) gamma_value0=14;
      if(LinksDirection==1) gamma_value0=13;

      //PS to PS;
      LatticePropagator sink_prop_f=sink_prop_1;
      LatticePropagator sink_prop_b=sink_prop_1;
      LatticePropagator F_tmp;
      for(int is=0;is<params.param.links_max+1;is++)
      for(int icase=0;icase<2;icase++)
      {
        if(is>0)
        {
          if(icase==0)
          {
            F_tmp=U[ LinksDirection ] * shift( sink_prop_f, FORWARD, LinksDirection );
            sink_prop_f=F_tmp;
          }
          if(icase==1)
          {
            F_tmp=shift( adj( U[ LinksDirection ] ) * sink_prop_b, BACKWARD, LinksDirection );
            sink_prop_b=F_tmp;
          }
        }
        int gamma_value=(is==0&&icase==0)?15:gamma_value0;
        for(int rnum=0; rnum < params.named_obj.antiquark_props.size(); ++rnum)
        {
          const LatticePropagator& sink_prop_2 =TheNamedObjMap::Instance().getData<LatticePropagator>(antiquark_props[rnum].quark_propagator_id);
          // Construct the anti-quark propagator from sink_prop_2
          LatticePropagator anti_sink_prop =  Gamma(15) * sink_prop_2 * Gamma(15);          
          if(icase==0&&is==0)
          {
             // Masses
             push(xml_out, "elem");
             write(xml_out, "Mass_1", quark_props[lnum].Mass);
             write(xml_out, "Mass_2", antiquark_props[rnum].Mass);
             write(xml_out, "t0", t0);
       
             // Save prop input
             push(xml_out, "Forward_prop_headers");
             write(xml_out, "First_forward_prop", quark_props[lnum].prop_header);
             write(xml_out, "Second_forward_prop", antiquark_props[rnum].prop_header);
             pop(xml_out);
       
             // Sanity check: write out the norm2 of the forward prop
             push(xml_out, "Forward_prop_correlator");
             {
               const LatticePropagator& sink_prop_1 =
                 TheNamedObjMap::Instance()
                     .getData<LatticePropagator>(quark_props[lnum]
                                                     .quark_propagator_id);
               const LatticePropagator& sink_prop_2 =
                 TheNamedObjMap::Instance()
                     .getData<LatticePropagator>(antiquark_props[rnum]
                                                     .quark_propagator_id);
       
               write(xml_out, "forward_prop_corr_1",
                     sumMulti(localNorm2(sink_prop_1), phases_nomom.getSet()));
               write(xml_out, "forward_prop_corr_2",
                     sumMulti(localNorm2(sink_prop_2), phases_nomom.getSet()));
             }
             pop(xml_out);
       
             // Dump out a bunch of information about the propagators
             push(xml_out, "SourceSinkType");
             {
               QDPIO::cout << "Source_type_1 = "
                           << quark_props[lnum].source_type << std::endl;
               QDPIO::cout << "Sink_type_1 = "
                           << quark_props[lnum].sink_type << std::endl;
               QDPIO::cout << "Source_type_2 = "
                           << antiquark_props[rnum].source_type << std::endl;
               QDPIO::cout << "Sink_type_2 = "
                           << antiquark_props[rnum].sink_type << std::endl;
       
               write(xml_out, "source_type_1", quark_props[lnum].source_type);
               write(xml_out, "source_disp_type_1",
                         quark_props[lnum].source_disp_type);
               write(xml_out, "sink_type_1", quark_props[lnum].sink_type);
               write(xml_out, "sink_disp_type_1",
                         quark_props[lnum].sink_disp_type);
       
               write(xml_out, "source_type_2", antiquark_props[rnum].source_type);
               write(xml_out, "source_disp_type_2",
                         antiquark_props[rnum].source_disp_type);
               write(xml_out, "sink_type_2", antiquark_props[rnum].sink_type);
               write(xml_out, "sink_disp_type_2",
                         antiquark_props[rnum].sink_disp_type);
             }
             pop(xml_out);
       

       
             // Construct group name for output
             std::string src_type;
             if (quark_props[lnum].source_type == "POINT_SOURCE")
               src_type = "Point";
             else if (quark_props[lnum].source_type == "SF_POINT_SOURCE")
               src_type = "Point";
             else if (quark_props[lnum].source_type == "SHELL_SOURCE")
               src_type = "Shell";
             else if (quark_props[lnum].source_type == "NORM_SHELL_SOURCE")
               src_type = "Shell";
             else if (quark_props[lnum].source_type == "SF_SHELL_SOURCE")
               src_type = "Shell";
             else if (quark_props[lnum].source_type == "WALL_SOURCE")
               src_type = "Wall";
             else if (quark_props[lnum].source_type == "SF_WALL_SOURCE")
               src_type = "Wall";
             else if (quark_props[lnum].source_type == "RAND_ZN_WALL_SOURCE")
               src_type = "Wall";
             else
             {
               QDPIO::cerr << "Unsupported source type = "
                           << quark_props[lnum].source_type << std::endl;
               QDP_abort(1);
             }
       
             std::string snk_type;
             if (quark_props[lnum].sink_type == "POINT_SINK")
               snk_type = "Point";
             else if (quark_props[lnum].sink_type == "SHELL_SINK")
               snk_type = "Shell";
             else if (quark_props[lnum].sink_type == "NORM_SHELL_SINK")
               snk_type = "Shell";
             else if (quark_props[lnum].sink_type == "WALL_SINK")
               snk_type = "Wall";
             else
             {
               QDPIO::cerr << "Unsupported sink type = "
                           << quark_props[lnum].sink_type << std::endl;
               QDP_abort(1);
             }
       
             std::string source_sink_type = src_type + "_" + snk_type;
             QDPIO::cout << "Source type = " << src_type << std::endl;
             QDPIO::cout << "Sink type = "   << snk_type << std::endl;
             pop(xml_out); // elem
          }
       
          // Construct the meson correlation function
          LatticeComplex corr;
          if(icase==0)corr = trace(adj(anti_sink_prop) * Gamma(gamma_value) * sink_prop_f * Gamma(15));
          if(icase==1)corr = trace(adj(anti_sink_prop) * Gamma(gamma_value) * sink_prop_b * Gamma(15));
          multi1d<DComplex> hsum = sumMulti(phase*corr,timeslice);
          if(Layout::primaryNode())
          for (int t=0; t < tlen; ++t) 
          {
            int t_eff = (t - t0 + tlen) % tlen;
            res.data[rnum*data_separation+data_start_position+offset*tlen*2+2*t_eff] = hsum[t].elem().elem().elem().real();
            res.data[rnum*data_separation+data_start_position+offset*tlen*2+2*t_eff+1] = hsum[t].elem().elem().elem().imag();
          }
          
        }// end of loop over antiquarks
        offset++;
      }//end of loop over directions   

    } // end of loop over quarks
    if(Layout::primaryNode()) res.save();
    pop(xml_out);  // Wilson_spectroscopy
    pop(xml_out);  // PionDA

    snoop.stop();
    QDPIO::cout << InlinePionDAEnv::name << ": total time = "
                << snoop.getTimeInSeconds()
                << " secs" << std::endl;

    QDPIO::cout << InlinePionDAEnv::name << ": ran successfully" << std::endl;

    END_CODE();
  } // end of InlinePionDA::func

}; // end of namespace Chroma


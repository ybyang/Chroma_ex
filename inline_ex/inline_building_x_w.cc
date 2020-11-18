/*! \file
 * \brief Inline construction of BuildingX
 *
 * Building Blocks on forward and sequential props
 */
#include "inline_ex/inline_building_x_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "meas_ex/sftmom.h"
#include "util/info/proginfo.h"
#include "meas/inline/make_xml_file.h"

#include "meas/inline/io/named_objmap.h"

#include "actions/ferm/fermstates/ferm_createstate_factory_w.h"
#include "actions/ferm/fermstates/ferm_createstate_aggregate_w.h"
#include "io_ex/io_general_class.h"


namespace Chroma
{

//------------------------------------------------------------------------------
// Factory Registration


  // Environment in which the measurement lives (holds params and such)
  namespace InlineBuildingXIOGEnv
  {
    namespace
    {
      // Function to register with a factory
      AbsInlineMeasurement* createMeasurement( XMLReader& xml_in,
                                       const std::string& path )
      {
        return new InlineBuildingXIOG(
                        InlineBuildingXIOGParams(xml_in, path) );
      }

      //! Local registration flag
      bool registered = false;
    }

    // The name of the measurement for the XML file
    const std::string name = "BUILDING_X_IOG";

    //! Register all the factories
    bool registerAll()
    {
      bool success = true;
      if (!registered)
      {
        success &= CreateFermStateEnv::registerAll(); // ??? Why?
        success &= TheInlineMeasurementFactory::Instance().
          registerObject(name, createMeasurement);
        registered = true;
      }
      return success;
    }
  } // end namespace InlineBuildingXIOGEnv

//------------------------------------------------------------------------------
// Parameter reading, writing, handling

  //! Reader for parameters
  void read(                         XMLReader& xml,
                             const std::string& path,
           InlineBuildingXIOGParams::Param_t& param )
  {
    XMLReader paramtop(xml, path);

    int version;
    read(paramtop, "version", version);
    param.use_sink_offset = false;
    param.canonical = false;
    param.time_reverse = false;
    param.translate = false;

    param.cfs = CreateFermStateEnv::nullXMLGroup();

    // Each version adds new parameters and is a superset of all previous versions
    switch (version)
    {
    case 6:
      // Read Multi_Src info
    case 5:
      read(paramtop, "time_reverse", param.time_reverse);
      read(paramtop, "translate", param.translate);

    case 4:
      if (paramtop.count("FermState") > 0)
        param.cfs = readXMLGroup(paramtop, "FermState", "Name");

    case 3:
      read(paramtop, "canonical", param.canonical);

    case 2:
      read(paramtop, "use_sink_offset", param.use_sink_offset);

    case 1:
      break;

    default :
      QDPIO::cerr << InlineBuildingXIOGEnv::name
                  << ": input parameter version " << version << " unsupported."
                  << std::endl;
      QDP_abort(1);
    }

    read(paramtop, "mom_dir", param.mom_dir);

    read(paramtop, "cfg_serial", param.icfg);
    if (paramtop.count("mom2_max") > 0) {
      read(paramtop, "mom2_max", param.mom2_max);
    } else {
      param.mom2_max = -1;
      read(paramtop, "mom2_list", param.mom2_list);
    }

  }

  //! Writer for parameters
  void write(                              XMLWriter& xml,
                                   const std::string& path,
           const InlineBuildingXIOGParams::Param_t& param )
  {
    push(xml, path);

    int version = 6;
    write(xml, "version", version);
    if (param.mom2_max >= 0) {
      write(xml, "mom2_max", param.mom2_max);
    } else {
      write(xml, "mom2_list", param.mom2_list);
    }
    write(xml, "canonical", param.canonical);
    write(xml, "time_reverse", param.time_reverse);
    write(xml, "translate", param.translate);
    xml << param.cfs.xml;

    pop(xml);
  }

  //! BB parameters (gauge id, props) input
  void read(                      XMLReader& xml,
                          const std::string& path,
           InlineBuildingXIOGParams::BB_t& bb )
  {
    XMLReader bbtop(xml, path);

    read(bbtop, "GaugeId", bb.GaugeId);
    read(bbtop, "FrwdPropId", bb.FrwdPropId);
    read(bbtop, "BkwdProps", bb.BkwdProps);
  }

  //! BB parameters (gauge id, props) output
  void write(                           XMLWriter& xml,
                                const std::string& path,
           const InlineBuildingXIOGParams::BB_t& bb )
  {
    push(xml, path);

    write(xml, "GaugeId", bb.GaugeId);
    write(xml, "FrwdPropId", bb.FrwdPropId);
    write(xml, "BkwdProps", bb.BkwdProps);

    pop(xml);
  }

  //! Correlator description "backward prop" (flavor, gamma, props, etc.) input
  void read(                            XMLReader& xml,
                                const std::string& path,
           InlineBuildingXIOGParams::BkwdProp_t& bkwdprop )
  {
    XMLReader bkwdproptop(xml, path);

    read(bkwdproptop, "BkwdPropId", bkwdprop.BkwdPropId);
    read(bkwdproptop, "BkwdPropG5Format", bkwdprop.BkwdPropG5Format);
    read(bkwdproptop, "GammaInsertion", bkwdprop.GammaInsertion);
    read(bkwdproptop, "Flavor", bkwdprop.Flavor);
    read(bkwdproptop, "BBFileNamePattern", bkwdprop.BBFileNamePattern);
  }

  //! Correlator description "backward prop" (flavor, gamma, props, etc.) output
  void write(                                 XMLWriter& xml,
                                      const std::string& path,
           const InlineBuildingXIOGParams::BkwdProp_t& bkwdprop )
  {
    push(xml, path);

    write(xml, "BkwdPropId", bkwdprop.BkwdPropId);
    write(xml, "BkwdPropG5Format", bkwdprop.BkwdPropG5Format);
    write(xml, "GammaInsertion", bkwdprop.GammaInsertion);
    write(xml, "Flavor", bkwdprop.Flavor);
    write(xml, "BBFileNamePattern", bkwdprop.BBFileNamePattern);

    pop(xml);
  }

  // Params default constructor
  InlineBuildingXIOGParams::InlineBuildingXIOGParams()
  {
    frequency = 0;
  }

  // Construct params from XML
  InlineBuildingXIOGParams::InlineBuildingXIOGParams(
            XMLReader& xml_in,
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
      read(paramtop, "BuildingX", bb);

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


  // Write out the parameters we constructed
  void InlineBuildingXIOGParams::write( XMLWriter& xml_out,
                                  const std::string& path )
  {
    push(xml_out, path);
    
    Chroma::write(xml_out, "Param", param);
    Chroma::write(xml_out, "BuildingX", bb);
    QDP::write(xml_out, "xml_file", xml_file);

    pop(xml_out);
  }


//------------------------------------------------------------------------------
// The real work is done here

  void AllLinkPatterns(      bool & DoThisPattern,
                             bool & DoFurtherPatterns,
      multi1d<unsigned short int> & LinkPattern );


  // Set up the XML and invoke func, which does the acual work
  void InlineBuildingXIOG::operator()( unsigned long update_no,
                                            XMLWriter& xml_out )
  {
    // If xml file not empty, then use alternate
    if (params.xml_file != "")
    {
      std::string xml_file = makeXMLFileName(params.xml_file, update_no);

      push(xml_out, "BuildingX");
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


void BuildingXIOG( const multi1d< LatticePropagator > &  B,
                     const LatticePropagator &             F,
                     const multi1d< LatticeColorMatrix > & U,
                     const multi1d< int > &                GammaInsertions,
                     const multi1d< int > &                Flavors,
                     const SftMom_cohen &                        Phases,
                     const SftMom_cohen &                        PhasesCanonical,
                     const multi2d< std::string > &        BinaryDataFileNames,
                     const signed short int                T1,
                     const signed short int                T2,
                     const signed short int                Tsrc,
                     const signed short int                Tsnk,
                     const std::string&                    SeqSourceType, 
                     const multi1d< int >&                 SnkMom, 
                     const signed short int                DecayDir,
                     const bool                            TimeReverse,
                     const bool                            ShiftFlag,
                     int mom_dir,int icfg
);

  // Real work done here
  void InlineBuildingXIOG::func( unsigned long update_no,
                                             XMLWriter& xml_out)
  {
    START_CODE();

    StopWatch snoop;
    snoop.reset();
    snoop.start();

    //--------------------------------------------------------------------------
    // Start building the output XML
    push(xml_out, "BuildingX");
    write(xml_out, "update_no", update_no);

    QDPIO::cout << " BuildingX" << std::endl;
    QDPIO::cout << "     volume: " << QDP::Layout::lattSize()[0];
    for (int i=1; i<Nd; ++i) {
      QDPIO::cout << " x " << QDP::Layout::lattSize()[i];
    }
    QDPIO::cout << std::endl;

    proginfo(xml_out); // Print out basic program info

    push(xml_out, "Output_version");
    write(xml_out, "out_version", 2);
    pop(xml_out);

    //--------------------------------------------------------------------------
    // Grab gauge configuration
    multi1d<LatticeColorMatrix> U;
    XMLBufferWriter gauge_xml;

    try
    {
      // Try to get the gauge field.
      // If it doesn't exist, an exception will be thrown.
      U = TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >
            (params.bb.GaugeId);
      TheNamedObjMap::Instance().get(params.bb.GaugeId).getRecordXML(gauge_xml);

      // Set the construct state and modify the fields
      QDPIO::cout << "cfs=XX" << params.param.cfs.xml << "XX" << std::endl;
      std::istringstream xml_s(params.param.cfs.xml);
      XMLReader fermtop(xml_s);

      Handle< CreateFermState<
        LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > >
          cfs( TheCreateFermStateFactory::Instance().createObject(
                 params.param.cfs.id, fermtop, params.param.cfs.path ) );

      Handle< FermState<
        LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > >
          state( (*cfs)(U) );

      // Pull the gauge fields out from state which may have fermBCs
      U = state->getLinks();
    }
    catch (std::bad_cast)
    {
      // If the object exists and is of the wrong type we will land here
      QDPIO::cerr << InlineBuildingXIOGEnv::name
                  << ": caught dynamic cast error" << std::endl;
      QDP_abort(1);
    }
    catch (const std::string& e)
    {
      // If the object is not in the map we will land here
      QDPIO::cerr << InlineBuildingXIOGEnv::name
                  << ": map call failed: " << e << std::endl;
      QDP_abort(1);
    }

    // Write out the input
    params.write(xml_out, "Input");

    // Write out the config info
    write(xml_out, "Config_info", gauge_xml);

    // Calculate a bunch of plaquettes. Why not? It's always fun!
    Double ave_plaq, ave_spacelike_plaq, ave_timelike_plaq, ave_link_trace;
    MesPlq( U, ave_plaq, ave_spacelike_plaq, ave_timelike_plaq, ave_link_trace );
    push(xml_out, "Observables");
    write(xml_out, "ave_plaq", ave_plaq);
    write(xml_out, "ave_spacelike_plaq", ave_spacelike_plaq);
    write(xml_out, "ave_timelike_plaq", ave_timelike_plaq);
    write(xml_out, "ave_link_trace", ave_link_trace);
    pop(xml_out);

    //--------------------------------------------------------------------------
    // Grab forward propagator

    // Keep a copy of the phases with no momenta
    SftMom phases_nomom(0, true, Nd-1);

    LatticePropagator F;
    ChromaProp_t prop_header;
    PropSourceConst_t source_header;
    QDPIO::cout << "Attempt to parse forward propagator" << std::endl;

    try
    {
      // Grab the forward propagator
      F = TheNamedObjMap::Instance().getData<LatticePropagator>(params.bb.FrwdPropId);
      // Grab its associated XML structures
      XMLReader FrwdPropXML, FrwdPropRecordXML;
      TheNamedObjMap::Instance().get(params.bb.FrwdPropId).getFileXML(FrwdPropXML);
      TheNamedObjMap::Instance().get(params.bb.FrwdPropId).getRecordXML(FrwdPropRecordXML);

      // Read the source and prop headers from the XML
      read(FrwdPropRecordXML, "/Propagator/ForwardProp", prop_header);
      read(FrwdPropRecordXML, "/Propagator/PropSource", source_header);

      // Diagnostic! Write out the zero-mom pion for this prop
      multi1d<Double> FrwdPropCheck =
        sumMulti( localNorm2(F), phases_nomom.getSet() );

      // Write out the forward propagator header
      push(xml_out, "ForwardProp");
      write(xml_out, "FrwdPropId", params.bb.FrwdPropId);
      write(xml_out, "FrwdPropXML", FrwdPropXML);
      write(xml_out, "FrwdPropRecordXML", FrwdPropRecordXML);
      write(xml_out, "FrwdPropCheck", FrwdPropCheck);
      pop(xml_out);
    }
    catch (std::bad_cast)
    {
      QDPIO::cerr << InlineBuildingXIOGEnv::name
                  << ": caught dynamic cast error" << std::endl;
      QDP_abort(1);
    }
    catch (const std::string& e)
    {
      QDPIO::cerr << InlineBuildingXIOGEnv::name
                  << ": forward prop: error message: " << e << std::endl;
      QDP_abort(1);
    }

    // If we made it through the try block, the parsing was successful
    QDPIO::cout << "Forward propagator successfully parsed" << std::endl;

    // Look at the forward prop to figure out which way is time
    int j_decay = source_header.j_decay;
    multi1d<int> tsrc = source_header.getTSrce() ;

    //--------------------------------------------------------------------------
    // Grab backward propagators (that is, sequential propagators)

    StopWatch swatch;
    push(xml_out, "SequentialSource");

    for (int loop = 0; loop < params.bb.BkwdProps.size(); ++loop)
    {
      push(xml_out, "elem");
      write(xml_out, "loop_ctr", loop);

      QDPIO::cout << "Backward propagator #" << loop << std::endl;

      multi1d<LatticePropagator> B(1);
      SeqSource_t seqsource_header;
      QDPIO::cout << "Attempt to parse backward propagator" << std::endl;
      try
      {
        // Grab the backward prop
        B[0] = TheNamedObjMap::Instance().getData<LatticePropagator>(
                  params.bb.BkwdProps[loop].BkwdPropId );

        // Grab its associated XML structures
        XMLReader BkwdPropXML, BkwdPropRecordXML;
        TheNamedObjMap::Instance().get(params.bb.BkwdProps[loop].BkwdPropId).
          getFileXML(BkwdPropXML);
        TheNamedObjMap::Instance().get(params.bb.BkwdProps[loop].BkwdPropId).
          getRecordXML(BkwdPropRecordXML);

        // Read the source and prop headers from the XML
        read(BkwdPropRecordXML, "/SequentialProp/SeqSource", seqsource_header);

        // Diagnostic! Write out the zero-mom pion for this prop
        // ??? Uh, why? The backward prop does not have a sensible time behavior
        multi1d<Double> BkwdPropCheck =
          sumMulti( localNorm2( B[0] ), phases_nomom.getSet() );

        // Write out the forward propagator header
        push(xml_out, "BackwardProp");
        write(xml_out, "BkwdPropId", params.bb.BkwdProps[loop].BkwdPropId);
        write(xml_out, "BkwdPropG5Format", params.bb.BkwdProps[loop].BkwdPropG5Format);
        write(xml_out, "SequentialSourceType", seqsource_header.seqsrc.id);
        write(xml_out, "BkwdPropXML", BkwdPropXML);
        write(xml_out, "BkwdPropRecordXML", BkwdPropRecordXML);
        write(xml_out, "BkwdPropCheck", BkwdPropCheck);
        pop(xml_out);
      }
      catch (std::bad_cast)
      {
        QDPIO::cerr << InlineBuildingXIOGEnv::name
                    << ": forward prop: caught dynamic cast error" << std::endl;
        QDP_abort(1);
      }
      catch (const std::string& e)
      {
        QDPIO::cerr << InlineBuildingXIOGEnv::name
                    << ": forward prop: error message: " << e << std::endl;
        QDP_abort(1);
      }

      QDPIO::cout << "Backward propagator #" << loop << " successfully parsed"
                  << std::endl;

      // Multiply in the specified gamma matrix
      multi1d<int> GammaInsertions(1);
      GammaInsertions[0] = params.bb.BkwdProps[loop].GammaInsertion;

      if (GammaInsertions[0] < 0 || GammaInsertions[0] >= Ns*Ns)
      {
        QDPIO::cerr << "InlineBuildingXIOG: Gamma insertion out of bounds: "
                    << GammaInsertions[0] << std::endl;
        QDP_abort(1);
      }

      // Dump some diagnostics
      QDPIO::cout << "Seqsource name  = "
                  << seqsource_header.seqsrc.id << std::endl;
      QDPIO::cout << "Gamma insertion = "
                  << params.bb.BkwdProps[loop].GammaInsertion << std::endl;
      QDPIO::cout << "Flavor          = "
                  << params.bb.BkwdProps[loop].Flavor << std::endl;

      write(xml_out, "seq_src", seqsource_header.seqsrc.id);
      write(xml_out, "gamma_insertion", GammaInsertions[0]);
      write(xml_out, "flavor", params.bb.BkwdProps[loop].Flavor);

      write(xml_out, "t_source", source_header.t_source);
      write(xml_out, "t_sink", seqsource_header.t_sink);
      write(xml_out, "sink_mom", seqsource_header.sink_mom);

      // Resolve the gamma-5's that may need to be at either end of the
      // backward propagator
      if (params.bb.BkwdProps[loop].BkwdPropG5Format == "G5_B")
      {
        LatticePropagator Bu = B[0];
        B[0] = Gamma(Ns*Ns-1)*Bu;
      }
      else if (params.bb.BkwdProps[loop].BkwdPropG5Format == "B_G5")
      {
        LatticePropagator Bu = B[0];
        B[0] = Bu*Gamma(Ns*Ns-1);
      }
      else if (params.bb.BkwdProps[loop].BkwdPropG5Format == "G5_B_G5")
      {
        LatticePropagator Bu = B[0];
        B[0] = Gamma(Ns*Ns-1)*Bu*Gamma(Ns*Ns-1);
      }

      // Figure out the momentum at the sink
      multi1d<int> SnkMom(Nd-1);
      if (params.param.use_sink_offset)
        SnkMom = seqsource_header.sink_mom;
      else
        SnkMom = 0;


      SftMom_cohen *Phases, *PhasesCanonical;
      {
        if (params.param.mom2_max >= 0) {
          Phases = new SftMom_cohen(params.param.mom2_max, tsrc, SnkMom, false,
                              j_decay);
          PhasesCanonical = new SftMom_cohen( params.param.mom2_max, tsrc, SnkMom,
                                    params.param.canonical, j_decay );
        } else {
          Phases = new SftMom_cohen(params.param.mom2_list, tsrc, SnkMom, false,
                              j_decay);
          PhasesCanonical = new SftMom_cohen( params.param.mom2_list, tsrc, SnkMom,
                                    params.param.canonical, j_decay );
        }
      }

      // Construct file names
      int NumO = PhasesCanonical->numMom();

      multi2d<std::string> Files(1, NumO);

      Files.resize(1, 1);
      Files(0, 0) = params.bb.BkwdProps[loop].BBFileNamePattern;

      // Set flavors
      multi1d<int> Flavors(1);
      switch ((int)(params.bb.BkwdProps[loop].Flavor[0]))
      {
      case (int)'U': Flavors[0] = 0; break;
      case (int)'D': Flavors[0] = 1; break;
      case (int)'S': Flavors[0] = 2; break;
      case (int)'C': Flavors[0] = 3; break;
      case (int)'B': Flavors[0] = 4; break;
      case (int)'T': Flavors[0] = 5; break;
      default:
        QDPIO::cerr << InlineBuildingXIOGEnv::name
                    << ": invalid flavor tag = "
                    << params.bb.BkwdProps[loop].Flavor
                    << ", should be one of U,D,S,C,T,B" << std::endl;
        QDP_abort(1);
      }


      // Diagnostic!
      //! Test a meson sequential source.
      /*!
       *  For the case of a meson, we have evaluated as the sequential source
       *
       *  H(y, 0; tx, p) = \sum exp{ip.x} U(y,x) \gamma_5\Gamma_f^\dag\gamma_5 D(x,0)
       *
       *  H^\dag(y, 0; tx, p) =
       *    \sum_x exp{-ip.x} \gamma_5 D(0,x) \Gamma_f U(x,y) \gamma_5
       *
       *  Thus we can see that
       *
       *  Tr[ \gamma_5 H^\dag(0,0; tx, p)\gamma_5 \Gamma_i] =
       *                         \sum_x exp{-ip.x} Tr[ D(0,x)\Gamma_f U(x,0) \Gamma_i ]
       *
       *  which is the desired meson correlator at momentum p and timslice tx
       */
      {
        push(xml_out, "LoopBackTest");
        write(xml_out, "seq_src", seqsource_header.seqsrc.id);
        write(xml_out, "gamma_insertion", GammaInsertions[0]);
        write(xml_out, "flavor", params.bb.BkwdProps[loop].Flavor);
        write(xml_out, "tsrc", tsrc);
        write(xml_out, "t_source", source_header.t_source);
        write(xml_out, "t_sink", seqsource_header.t_sink);
        write(xml_out, "sink_mom", seqsource_header.sink_mom);
        LatticeComplex tr = trace(adj(B[0]) * Gamma(GammaInsertions[0]));
        Complex seq_src_value = peekSite(tr, tsrc);
        write(xml_out, "source_value", seq_src_value);
        pop(xml_out);
      }

      // Diagnostic! Further diagnostics...
      {
        // assumes any Gamma5 matrices have already been absorbed
        int GammaInsertion = GammaInsertions[0];
        LatticePropagator GFG;
        LatticeComplex tr;
        multi1d<DComplex> pr;

        push(xml_out, "DiagnosticTest");
        write(xml_out, "GammaInsertion", GammaInsertion);

        GFG = Gamma(0) * F * Gamma(GammaInsertion);
        tr = localInnerProduct( B[0], GFG );
        pr = sumMulti(tr, Phases->getSet());
        write(xml_out, "formFactor_G0", pr);

        GFG = Gamma(1<<(Nd-1)) * F * Gamma(GammaInsertion);
        tr = localInnerProduct( B[0], GFG );
        pr = sumMulti(tr, Phases->getSet());
        write(xml_out, "formFactor_G8", pr);

        GFG = Gamma(Ns*Ns-1) * F * Gamma(GammaInsertion);
        tr = localInnerProduct( B[0], GFG );
        pr = sumMulti(tr, Phases->getSet());
        write(xml_out, "formFactor_G15", pr);

        pop(xml_out);
      }

      // Invoke the master building blocks function!
      swatch.reset();
      QDPIO::cout << "calculating building x" << std::endl;

      const signed short int T1 = 0;
      const signed short int T2 = QDP::Layout::lattSize()[j_decay] - 1;
      const signed short int DecayDir = j_decay;
      const signed short int Tsrc = source_header.t_source;
      const signed short int Tsnk = seqsource_header.t_sink;

      swatch.start();
      BuildingXIOG(B, F, U,
		     GammaInsertions, Flavors,
		     *Phases, *PhasesCanonical,
		     Files, T1, T2,
		     Tsrc, Tsnk,
		     seqsource_header.seqsrc.id, seqsource_header.sink_mom, DecayDir,
		     params.param.time_reverse,
		     params.param.translate,
		     params.param.mom_dir,
         params.param.icfg);
      swatch.stop();
     
      delete Phases;
      delete PhasesCanonical;

      QDPIO::cout << "finished building blocks for backward prop #" << loop
                  << "  time= " << swatch.getTimeInSeconds() << " secs"
                  << std::endl;

      pop(xml_out);   // elem
    } // end loop over sequential sources

    pop(xml_out);  // SequentialSource

    pop(xml_out);  // BuildingX

    snoop.stop();
    QDPIO::cout << InlineBuildingXIOGEnv::name << ": total time = "
                << snoop.getTimeInSeconds()
                << " secs" << std::endl;

    QDPIO::cout << InlineBuildingXIOGEnv::name << ": ran successfully"
                << std::endl;

    END_CODE();
  } // end of InlineBuildingXIOG::func

LatticePropagator derivative(const LatticePropagator &F,
                const multi1d< LatticeColorMatrix > & U,int mu)
{
      return 0.5*(U[mu] *shift( F, FORWARD, mu)-shift( adj( U[mu] ) * F, BACKWARD, mu ));
}

LatticeComplex T_munurho(const LatticePropagator &B,
                const LatticePropagator &F,
                const multi1d< LatticeColorMatrix > & U,
                int rho,int mu,int nu)
{
// not the real general case. 
     int igrho=1;for(int ir=0;ir<rho;ir++)igrho*=2;
     LatticePropagator tmp,DDF,DDB;
     DDF=derivative(derivative(F,U,mu),U,nu);
     DDB=derivative(derivative(B,U,nu),U,mu);
     return 0.5*localInnerProduct( B, Gamma(igrho)*DDF)
             +0.5*localInnerProduct( DDB, Gamma(igrho)*F);
}

void BuildingXIOG( const multi1d< LatticePropagator > &  B,
                     const LatticePropagator &             F,
                     const multi1d< LatticeColorMatrix > & U,
                     const multi1d< int > &                GammaInsertions,
                     const multi1d< int > &                Flavors,
                     const SftMom_cohen &                        Phases,
                     const SftMom_cohen &                        PhasesCanonical,
                     const multi2d< std::string > &        BinaryDataFileNames,
                     const signed short int                T1,
                     const signed short int                T2,
                     const signed short int                Tsrc,
                     const signed short int                Tsnk,
                     const std::string&                    SeqSourceType, 
                     const multi1d< int >&                 SnkMom, 
                     const signed short int                DecayDir,
                     const bool                            TimeReverse,
                     const bool                            ShiftFlag,
                     int mom_dir, int icfg)
{
  StopWatch TotalTime;
  TotalTime.reset();
  TotalTime.start();

  StopWatch Timer;

  //##########################################################################//
  // open building blocks data files                                          //
  //##########################################################################//

  Timer.reset();
  Timer.start();

  const int NumF = B.size();
  const int NumO = PhasesCanonical.numMom();
  multi1d< general_data_base > iog_data;

  iog_data.resize( NumF);

  multi1d< int > GBB_NLinkPatterns( NumF );

    int nop=2;
    const int NumQ = Phases.numMom();
    for ( int f = 0; f < NumF; f++ )
    {
      sprintf(iog_data(f).name,"%s",BinaryDataFileNames(f, 0).c_str());
      iog_data(f).add_dimension(dim_conf,1,&icfg);
      iog_data(f).add_dimension(dim_channel,1,&f);
      iog_data(f).add_dimension(dim_operator,nop);
      int mom_list[NumQ];
      int icount=0;
      for ( int idx_o=0; idx_o < NumO; ++idx_o)
      for ( int q = 0; q < NumQ; q++ )
      {
          multi1d< int > Q = Phases.numToMom( q );
          int o = PhasesCanonical.momToNum( Q );
          if (o != idx_o) continue;
          if (o == -1 || DecayDir!=3) {mom_list[icount++]=-1; continue;}
          mom_list[icount++]=(50+Q[0])*10000+(50+Q[1])*100+50+Q[2];
      }
      iog_data(f).add_dimension(dim_momentum,NumQ,mom_list);
      iog_data(f).add_dimension(dim_t,QDP::Layout::lattSize()[DecayDir]);
      iog_data(f).add_dimension(dim_complex,2);
      if(Layout::primaryNode())iog_data(f).initialize();

      GBB_NLinkPatterns[f] = 0;
    }

  Timer.stop();
  QDPIO::cout << __func__ << ": time to open files = "
              << Timer.getTimeInSeconds() 
              << " seconds" << std::endl;

  //##########################################################################//
  // calculate building blocks                                                //
  //##########################################################################//

  Timer.reset();
  Timer.start();

  QDPIO::cout << __func__ << ": start x and x^2" << std::endl;

  const unsigned short int NLinks = 0;
  multi1d< unsigned short int > LinkDirs( 0 );
  const int NT   = Phases.numSubsets();
  

  for ( int f = 0; f < NumF; f++ )
  for ( int iop=0; iop<nop; iop++)
  {
    //build data here
    //operator: X and X2
    //momentum: Q
    //nt
    //complex
/*
     T4i=( conj(D_4 B) g_3 F + conj(D_3 B) g_4 F 
                 + conj(B) g_3 D_4 F + conj(B) g_4 D_3 F )/4
         =( conj(D_4 B) g_3 F 
                 + conj(B) g_3 D_4 F + 2* conj(B) g_4 D_3 F )/4 
     Note that conj(D_3 B) g_4 F= conj(B) g_4 D_3 F when we sum over \vec{x} without momentum
*/   
    QDPIO::cout << "before oper" << std::endl;
    LatticeComplex Trace=zero;
    if(iop==0) 
    {
        int iz=mom_dir;
        int igz=1;for(int ir=0;ir<mom_dir;ir++)igz*=2;
        LatticePropagator DF_T,DF_Z,DB_T;
        DF_T=derivative(F,U,3);
        DF_Z=derivative(F,U,iz);
        DB_T=derivative(B[f],U,3);
        Trace = -0.25*localInnerProduct( DB_T, Gamma(igz)*F)
                    +0.25*localInnerProduct( B[f], Gamma(igz)*DF_T)
                    +0.5*localInnerProduct( B[f], Gamma(8)*DF_Z);
    }
/*
    T4ii_sub=T433-(T411+T422)/2;

*/    
    if(iop==1) 
    {
        for(int idr=0;idr<3;idr++)
        {
            double fac=(idr==mom_dir)?1.0:-0.5;
            Trace=Trace
                  +fac*(T_munurho(B[f],F,U,3,idr,idr)
                       +T_munurho(B[f],F,U,idr,3,idr)
                       +T_munurho(B[f],F,U,idr,idr,3));
        }
        Trace=Trace/3;
    }
    multi2d< DComplex > Projections = Phases.sft( Trace );
    QDPIO::cout << "after oper" << std::endl;

    for ( int idx_o=0; idx_o < NumO; ++idx_o)
    for ( int q = 0; q < NumQ; q++ )
    {
      multi1d< DComplex > Projection = Projections[q];
      multi1d< int > Q = Phases.numToMom( q );

      int o = PhasesCanonical.momToNum( Q );

      // Write outputs in momentum order
      if (o != idx_o) continue;

      if (o == -1)
      {
        QDPIO::cerr << __func__ << ": internal error: failed to find index of ordered momentum" << std::endl;
        QDP_abort(1);
      }
      const signed short int QX = Q[0];
      const signed short int QY = Q[1];
      const signed short int QZ = Q[2];

      float real_part[T2 - T1 + 1];
      float imag_part[T2 - T1 + 1];

      // Fill correlator
      for ( int t = T1; t <= T2; t++ )
      {
        float r = toFloat( real( Projection[t] ) );
        float i = toFloat( imag( Projection[t] ) );

        int t_prime = t;

        if ( TimeReverse == true )
        {
          // QDPIO::cout<<"TimeReversing: " ;
          // shift the time origin to the source
          int t_shifted = ( t - Tsrc + NT ) % NT;
          //time reverse around the source
          int t_reversed = ( NT - t_shifted ) % NT;
          // Undo the shift to put time back where it was.
          // We may not want to do this. It may be better to just shift
          // the time origin to Tsrc as we do in the spectrum
          if ( ShiftFlag == false )
          {
            t_prime = ( t_reversed + Tsrc ) % NT;
          }
          else
          {
            t_prime = t_reversed;
          }
          //QDPIO::cout<<t<<" "<<t_prime<<" [Tsrc="<<Tsrc<<",NT="<<NT<<",tsh="<<t_shifted<<"]"<<endl;
        }
        // When TimeReverse is on, shifting is done differently.
        if ( ( ShiftFlag == true ) && ( TimeReverse == false ) )
          t_prime = ( t - Tsrc + NT ) % NT;

        //if (TimeReverse==false)
        // QDPIO::cout<<t<<" "<<t_prime<<" [Tsrc="<<Tsrc<<",NT="<<NT<<"]"<<endl;

        real_part[t_prime] = r;
        imag_part[t_prime] = i;

        #if _DEBUG_BB_C_ == 1
        {
          QDPIO::cout << "DEBUG: " << __FILE__ << " " << __LINE__ << "\n";
          QDPIO::cout << "q = " << q << "\n";
          QDPIO::cout << "o = " << o << "\n";
          QDPIO::cout << "f = " << f << "\n";
          QDPIO::cout << "t = " << t << "\n";
          QDPIO::cout << "r = " << r << "\n";
          QDPIO::cout << "i = " << i << "\n";
        }
        #endif
      }

      // Write correlator
      if(Layout::primaryNode())
      for ( int t = 0; t < (T2-T1+1); t++ )
      {
        int t0=t;
//        if(ShiftFlag==true)t0=(t+Tsrc)%(T2-T1+1);
        iog_data(f).data[0+2*(t+GBB_NLinkPatterns[f]*(T2-T1+1))]=real_part[t0];
        iog_data(f).data[1+2*(t+GBB_NLinkPatterns[f]*(T2-T1+1))]=imag_part[t0];
      }
      GBB_NLinkPatterns[f]++;
    } // for(idx_o, q)
  }

  Timer.stop();
  QDPIO::cout << __func__ << ": total time for 0 links (single BkwdFrwdTr call) = "
              << Timer.getTimeInSeconds() 
              << " seconds" << std::endl;

  Timer.reset();
  Timer.start();

  QDPIO::cout << __func__ << ": start AddLinks" << std::endl;

            
  Timer.stop();
  QDPIO::cout << __func__ << ": total time for remaining links (outermost AddLinks call) = "
              << Timer.getTimeInSeconds() 
              << " seconds" << std::endl;

  //##########################################################################//
  // add footer and close files                                               //
  //##########################################################################//

  Timer.reset();
  Timer.start();

  const unsigned short int Id = 0;  // indicates building blocks
  unsigned short int Version;  // building blocks version
    Version = 4;

  const unsigned short int Contraction = 0;  // 0 indicates connected diagram
  const signed short int   SeqSourceLen = 64;
  std::string SeqSource = SeqSourceType;
  SeqSource.resize(SeqSourceLen, 0);

  for ( int f = 0; f < NumF; f++ )
  {
      if(Layout::primaryNode())
      {
         iog_data(f).print();
         iog_data(f).save();
      fflush(stdout);
      }
  } // for(f)

  Timer.stop();
  QDPIO::cout << __func__ << ": time to write footer = "
              << Timer.getTimeInSeconds() 
              << " seconds" << std::endl;

  TotalTime.stop();
  QDPIO::cout << __func__ << ": total time = "
              << TotalTime.getTimeInSeconds() 
              << " seconds" << std::endl;

  return;
}

//############################################################################//

#undef _DEBUG_BB_C_

//############################################################################//

}  // end namespace Chroma


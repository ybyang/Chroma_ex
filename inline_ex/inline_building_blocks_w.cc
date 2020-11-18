/*! \file
 * \brief Inline construction of BuildingBlocks
 *
 * Building Blocks on forward and sequential props
 */
#include "inline_ex/inline_building_blocks_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "meas_ex/sftmom.h"
#include "meas/hadron/BuildingBlocks_w.h"
#include "util/info/proginfo.h"
#include "meas/inline/make_xml_file.h"

#include "meas/inline/io/named_objmap.h"

#include "actions/ferm/fermstates/ferm_createstate_factory_w.h"
#include "actions/ferm/fermstates/ferm_createstate_aggregate_w.h"
#include "io_ex/io_general_class.h"

#ifdef BUILD_QUDA
#include <quda.h>
#include <comm_quda.h>
#endif
#include "util/ferm/transf.h"


namespace Chroma
{

//------------------------------------------------------------------------------
// Factory Registration


  // Environment in which the measurement lives (holds params and such)
  namespace InlineBuildingBlocksIOGEnv
  {
    bool conserved =false;
    namespace
    {
      // Function to register with a factory
      AbsInlineMeasurement* createMeasurement( XMLReader& xml_in,
                                       const std::string& path )
      {
        return new InlineBuildingBlocksIOG(
                        InlineBuildingBlocksIOGParams(xml_in, path) );
      }

      //! Local registration flag
      bool registered = false;
    }

    // The name of the measurement for the XML file
    const std::string name = "BUILDING_BLOCKS_IOG";

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
  } // end namespace InlineBuildingBlocksIOGEnv

//------------------------------------------------------------------------------
// Parameter reading, writing, handling

  //! Reader for parameters
  void read(                         XMLReader& xml,
                             const std::string& path,
           InlineBuildingBlocksIOGParams::Param_t& param )
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
      QDPIO::cerr << InlineBuildingBlocksIOGEnv::name
                  << ": input parameter version " << version << " unsupported."
                  << std::endl;
      QDP_abort(1);
    }

    if (paramtop.count("no_gaugelink") > 0) {
      read(paramtop, "no_gaugelink", param.no_gaugelink);
    } else {
      param.no_gaugelink=false;
    }

    if (paramtop.count("use_cpu") > 0) {
      read(paramtop, "use_cpu", param.use_cpu);
    } else {
      param.use_cpu=true;
    }
#ifdef BUILD_QUDA
    if (paramtop.count("use_gpu") > 0) {
      read(paramtop, "use_gpu", param.use_gpu);
    } else {
      param.use_gpu=true;
    }
#endif
    if (paramtop.count("gamma_id") > 0) {
      read(paramtop, "gamma_id", param.gamma_id);
    } else {
      param.gamma_id=-1;
    }


    if (paramtop.count("conserved") > 0) {
      read(paramtop, "conserved", param.conserved);
    } else {
      param.conserved=false;
    }
    InlineBuildingBlocksIOGEnv::conserved=param.conserved;

    read(paramtop, "cfg_serial", param.icfg);
    read(paramtop, "links_max", param.links_max);
    if (paramtop.count("links_dir") > 0) {
      read(paramtop, "links_dir", param.links_dir);
    } else {
      param.links_dir = -1;
    }
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
           const InlineBuildingBlocksIOGParams::Param_t& param )
  {
    push(xml, path);

    int version = 6;
    write(xml, "version", version);
    write(xml, "links_max", param.links_max);
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
           InlineBuildingBlocksIOGParams::BB_t& bb )
  {
    XMLReader bbtop(xml, path);

    read(bbtop, "GaugeId", bb.GaugeId);
    read(bbtop, "FrwdPropId", bb.FrwdPropId);
    read(bbtop, "BkwdProps", bb.BkwdProps);
  }

  //! BB parameters (gauge id, props) output
  void write(                           XMLWriter& xml,
                                const std::string& path,
           const InlineBuildingBlocksIOGParams::BB_t& bb )
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
           InlineBuildingBlocksIOGParams::BkwdProp_t& bkwdprop )
  {
    XMLReader bkwdproptop(xml, path);

    read(bkwdproptop, "BkwdPropId", bkwdprop.BkwdPropId);
    read(bkwdproptop, "BkwdPropG5Format", bkwdprop.BkwdPropG5Format);
    read(bkwdproptop, "GammaInsertion", bkwdprop.GammaInsertion);
    read(bkwdproptop, "Flavor", bkwdprop.Flavor);
    read(bkwdproptop, "BBFileNamePattern", bkwdprop.BBFileNamePattern);
#ifdef BUILD_QUDA
    read(bkwdproptop, "BBFileNamePattern_gpu", bkwdprop.BBFileNamePattern_gpu);
#endif  
  }


  //! Correlator description "backward prop" (flavor, gamma, props, etc.) output
  void write(                                 XMLWriter& xml,
                                      const std::string& path,
           const InlineBuildingBlocksIOGParams::BkwdProp_t& bkwdprop )
  {
    push(xml, path);

    write(xml, "BkwdPropId", bkwdprop.BkwdPropId);
    write(xml, "BkwdPropG5Format", bkwdprop.BkwdPropG5Format);
    write(xml, "GammaInsertion", bkwdprop.GammaInsertion);
    write(xml, "Flavor", bkwdprop.Flavor);   
    write(xml, "BBFileNamePattern", bkwdprop.BBFileNamePattern);
#ifdef BUILD_QUDA
    write(xml, "BBFileNamePattern_gpu", bkwdprop.BBFileNamePattern_gpu);
#endif

    pop(xml);
  }

  // Params default constructor
  InlineBuildingBlocksIOGParams::InlineBuildingBlocksIOGParams()
  {
    frequency = 0;
  }

  // Construct params from XML
  InlineBuildingBlocksIOGParams::InlineBuildingBlocksIOGParams(
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
      read(paramtop, "BuildingBlocks", bb);

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
  void InlineBuildingBlocksIOGParams::write( XMLWriter& xml_out,
                                  const std::string& path )
  {
    push(xml_out, path);
    
    Chroma::write(xml_out, "Param", param);
    Chroma::write(xml_out, "BuildingBlocks", bb);
    QDP::write(xml_out, "xml_file", xml_file);

    pop(xml_out);
  }


//------------------------------------------------------------------------------
// The real work is done here

  void AllLinkPatterns(      bool & DoThisPattern,
                             bool & DoFurtherPatterns,
      multi1d<unsigned short int> & LinkPattern );


  // Set up the XML and invoke func, which does the acual work
  void InlineBuildingBlocksIOG::operator()( unsigned long update_no,
                                            XMLWriter& xml_out )
  {
    // If xml file not empty, then use alternate
    if (params.xml_file != "")
    {
      std::string xml_file = makeXMLFileName(params.xml_file, update_no);

      push(xml_out, "BuildingBlocks");
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


void BuildingBlocksIOG( const multi1d< LatticePropagator > &  B,
                     const LatticePropagator &             F,
                     const multi1d< LatticeColorMatrix > & U,
                     const multi1d< int > &                GammaInsertions,
                     const multi1d< int > &                Flavors,
                     const unsigned short int              MaxNLinks,
                     const BBLinkPattern                   LinkPattern,
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
                     const int	links_idr,int icfg,            
                     const int  gamma_id
);

#ifdef BUILD_QUDA
void BuildingBlocksIOG_GPU( const multi1d< LatticePropagator > &  B,
                     const LatticePropagator &             F,
                     const multi1d< LatticeColorMatrix > & U,
                     const multi1d< int > &                GammaInsertions,
                     const multi1d< int > &                Flavors,
                     const unsigned short int              MaxNLinks,
                     const BBLinkPattern                   LinkPattern,
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
                     const int  links_idr,int icfg,
                     const int gamma_id
);
#endif




  // Real work done here
  void InlineBuildingBlocksIOG::func( unsigned long update_no,
                                             XMLWriter& xml_out)
  {
    START_CODE();

    StopWatch snoop;
    snoop.reset();
    snoop.start();

    //--------------------------------------------------------------------------
    // Start building the output XML
    push(xml_out, "BuildingBlocks");
    write(xml_out, "update_no", update_no);

    QDPIO::cout << " BuildingBlocks" << std::endl;
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
      QDPIO::cerr << InlineBuildingBlocksIOGEnv::name
                  << ": caught dynamic cast error" << std::endl;
      QDP_abort(1);
    }
    catch (const std::string& e)
    {
      // If the object is not in the map we will land here
      QDPIO::cerr << InlineBuildingBlocksIOGEnv::name
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
      QDPIO::cerr << InlineBuildingBlocksIOGEnv::name
                  << ": caught dynamic cast error" << std::endl;
      QDP_abort(1);
    }
    catch (const std::string& e)
    {
      QDPIO::cerr << InlineBuildingBlocksIOGEnv::name
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
        QDPIO::cerr << InlineBuildingBlocksIOGEnv::name
                    << ": forward prop: caught dynamic cast error" << std::endl;
        QDP_abort(1);
      }
      catch (const std::string& e)
      {
        QDPIO::cerr << InlineBuildingBlocksIOGEnv::name
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
        QDPIO::cerr << "InlineBuildingBlocksIOG: Gamma insertion out of bounds: "
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
#ifdef BUILD_QUDA
      multi2d<std::string> Files_gpu(1, NumO);

      Files_gpu.resize(1, 1);
      Files_gpu(0, 0) = params.bb.BkwdProps[loop].BBFileNamePattern_gpu;
 #endif




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
        QDPIO::cerr << InlineBuildingBlocksIOGEnv::name
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
      QDPIO::cout << "calculating building blocks" << std::endl;

      const signed short int T1 = 0;
      const signed short int T2 = QDP::Layout::lattSize()[j_decay] - 1;
      const signed short int DecayDir = j_decay;
      const signed short int Tsrc = source_header.t_source;
      const signed short int Tsnk = seqsource_header.t_sink;

      QDPIO::cout << "params.param.no_gaugelink = "<< params.param.no_gaugelink  << std::endl;
   
      if (  params.param.no_gaugelink    )
      {

    for(int mu=0; mu < Nd; ++mu)
    for(int col=0; col < 3; ++col)
    for(int row=0; row < 3; ++row)
      {
       if(col==row)
        pokeColor(U[mu], 
                  cmplx(LatticeReal(1.0),LatticeReal(0.0)), 
                  row, col);
       else
        pokeColor(U[mu], 
                  cmplx(LatticeReal(0.0),LatticeReal(0.0)), 
                  row, col);
      }


     }
    

      swatch.start();

      if (  params.param.use_cpu   )
      {
      BuildingBlocksIOG(B, F, U,
		     GammaInsertions, Flavors,
		     params.param.links_max, AllLinkPatterns,
		     *Phases, *PhasesCanonical,
		     Files, T1, T2,
		     Tsrc, Tsnk,
		     seqsource_header.seqsrc.id, seqsource_header.sink_mom, DecayDir,
		     params.param.time_reverse,
		     params.param.translate,
		     params.param.links_dir,
                     params.param.icfg,
                     params.param.gamma_id);
       }

#ifdef BUILD_QUDA
       if (   params.param.use_gpu   )
       {
       BuildingBlocksIOG_GPU(B, F, U,
                     GammaInsertions, Flavors,
                     params.param.links_max, AllLinkPatterns,
                     *Phases, *PhasesCanonical,
                     Files_gpu, T1, T2,
                     Tsrc, Tsnk,
                     seqsource_header.seqsrc.id, seqsource_header.sink_mom, DecayDir,
                     params.param.time_reverse,
                     params.param.translate,
                     params.param.links_dir,
                     params.param.icfg,
                     params.param.gamma_id);
        }
#endif
      swatch.stop();
     
      delete Phases;
      delete PhasesCanonical;


      QDPIO::cout << "finished building blocks for backward prop #" << loop
                  << "  time= " << swatch.getTimeInSeconds() << " secs"
                  << std::endl;

      pop(xml_out);   // elem
    } // end loop over sequential sources

    pop(xml_out);  // SequentialSource

    pop(xml_out);  // BuildingBlocks

    snoop.stop();
    QDPIO::cout << InlineBuildingBlocksIOGEnv::name << ": total time = "
                << snoop.getTimeInSeconds()
                << " secs" << std::endl;

    QDPIO::cout << InlineBuildingBlocksIOGEnv::name << ": ran successfully"
                << std::endl;

    END_CODE();
  } // end of InlineBuildingBlocksIOG::func

#ifdef BUILD_QUDA
void GPU_kernel( const LatticePropagator &             B,
                 const LatticePropagator &             F,
                 const multi1d< LatticeColorMatrix > & U,
                 std::vector< LatticeSpinMatrix  > &  out,                
                 const int                             lnk_dir  ,
                 const unsigned short int              lnk_max  
               )
{

           StopWatch TotalTime;
           TotalTime.reset();
           TotalTime.start();

           StopWatch Timer;

           int LGCalls=0;
           int CTCalls=0;
 
           double LGTime=0.0;
           double CTTime=0.0;


           QudaPrecision_s cpu_prec;
           QudaPrecision_s gpu_prec;

           int s = sizeof( WordType<Real>::Type_t );
           if (s == 4) {
             cpu_prec = QUDA_SINGLE_PRECISION;
           }
           else {
             cpu_prec = QUDA_DOUBLE_PRECISION;
           }
/*
           switch( params.cudaPrecision ) {
           case HALF:
             gpu_prec = QUDA_HALF_PRECISION;
             break;
           case SINGLE:
             gpu_prec = QUDA_SINGLE_PRECISION;
             break;
           case DOUBLE:
             gpu_prec = QUDA_DOUBLE_PRECISION;
             break;
           default:
             gpu_prec = cpu_prec;
             break;
           }
*/

           gpu_prec = cpu_prec; 
         
           QudaGaugeParam q_gauge_param = newQudaGaugeParam();

           const multi1d<int>& latdims = Layout::subgridLattSize();
           q_gauge_param.X[0] = latdims[0];
           q_gauge_param.X[1] = latdims[1];
           q_gauge_param.X[2] = latdims[2];
           q_gauge_param.X[3] = latdims[3];
           q_gauge_param.type = QUDA_SU3_LINKS;
           q_gauge_param.gauge_order = QUDA_QDP_GAUGE_ORDER;
           q_gauge_param.t_boundary = QUDA_PERIODIC_T;
           q_gauge_param.cpu_prec = cpu_prec;
           q_gauge_param.cuda_prec = gpu_prec;


           q_gauge_param.reconstruct = QUDA_RECONSTRUCT_12;

           q_gauge_param.cuda_prec_sloppy = gpu_prec;
           q_gauge_param.reconstruct_sloppy = QUDA_RECONSTRUCT_12;

           q_gauge_param.gauge_fix = QUDA_GAUGE_FIXED_NO;
           q_gauge_param.anisotropy = 1.0;

           multi1d<int> face_size(4);
           face_size[0] = latdims[1]*latdims[2]*latdims[3]/2;
           face_size[1] = latdims[0]*latdims[2]*latdims[3]/2;
           face_size[2] = latdims[0]*latdims[1]*latdims[3]/2;
           face_size[3] = latdims[0]*latdims[1]*latdims[2]/2;

           int max_face = face_size[0];
           for(int i=1; i <=3; i++) {
             if ( face_size[i] > max_face ) {
               max_face = face_size[i];
             }
           }

           q_gauge_param.ga_pad = max_face;
           q_gauge_param.cuda_prec_precondition = gpu_prec;
           q_gauge_param.reconstruct_precondition = QUDA_RECONSTRUCT_12;

/*
           multi1d<LatticeColorMatrix> links_single(Nd);
           const Real twopi = 6.283185307179586476925286;
           
           for(int mu=0; mu < Nd; mu++)
           {
           links_single[mu] = (state->getLinks())[mu];    
           }
*/
           void* gauge[4];

           QDPIO::cout << __func__ << " loadGaugeQuda begin " << std::endl;

#ifndef BUILD_QUDA_DEVIFACE_GAUGE
           for(int mu=0; mu < Nd; mu++) {
           gauge[mu] = (void *)&(U[mu].elem(all.start()).elem().elem(0,0).real());
           }
#else
           GetMemoryPtrGauge(gauge,U);
#endif


           Timer.reset();
           Timer.start();

           loadGaugeQuda((void *)gauge, &q_gauge_param);
              
           Timer.stop();
           LGCalls +=1;
           LGTime += Timer.getTimeInSeconds();

           QDPIO::cout << __func__ << " loadGaugeQuda done " << std::endl;


           
           QudaInvertParam  quda_inv_param = newQudaInvertParam();
           quda_inv_param.cpu_prec = cpu_prec;
           quda_inv_param.cuda_prec = gpu_prec;
           quda_inv_param.cuda_prec_sloppy = gpu_prec;
           quda_inv_param.preserve_source = QUDA_PRESERVE_SOURCE_NO;
           quda_inv_param.dirac_order = QUDA_DIRAC_ORDER;
           quda_inv_param.input_location = QUDA_CPU_FIELD_LOCATION;
           quda_inv_param.output_location = QUDA_CPU_FIELD_LOCATION;

           quda_inv_param.dslash_type = QUDA_WILSON_DSLASH;
           quda_inv_param.gamma_basis = QUDA_DEGRAND_ROSSI_GAMMA_BASIS;
           quda_inv_param.Ls = 1;
           quda_inv_param.solution_type = QUDA_MAT_SOLUTION ;

           quda_inv_param.tune = QUDA_TUNE_NO;
           quda_inv_param.sp_pad = 0;
           quda_inv_param.cl_pad = 0;
           quda_inv_param.verbosity = QUDA_DEBUG_VERBOSE;


           std::vector< LatticeSpinMatrix > psi_s(lnk_max+1);
           
           LatticeFermion tmp_F;
           LatticeFermion tmp_B;

          
             shiftPatternQuda_chroma  p_s_q;
             shiftPatternQuda_chroma  p_s_qbar; 

             p_s_q.post_shift.push_back(  lnk_max  );
             p_s_q.post_shift.push_back(  lnk_dir  );
           
           for ( int i=0; i<lnk_max+1  ; i++ )out[i]=0;
           
           for(int is=0;is<4;is++)
           for(int ic=0;ic<3;ic++)
           {
           PropToFerm(F,tmp_F,ic,is);
           PropToFerm(B,tmp_B,ic,is);
          
           std::vector< void *  > spinorOut( lnk_max+1 );


#ifndef BUILD_QUDA_DEVIFACE_SPINOR

           void* spinorIn_q    =(void *)&(tmp_F.elem(0).elem(0).elem(0).real());
           void* spinorIn_qbar =(void *)&(tmp_B.elem(0).elem(0).elem(0).real());
         
           for (int i=0; i<lnk_max+1 ;i++)
           {
            spinorOut[i] =(void *)&(psi_s[i].elem(0).elem(0,0).elem().real());   ;
           }
#else

           void* spinorIn_q    = GetMemoryPtr( tmp_F.getId() );
           void* spinorIn_qbar = GetMemoryPtr( tmp_B.getId() );
           for (int i=0; i< lnk_max+1 ; i++  )
           {
           spinorOut[i] = GetMemoryPtr( psi_s[i].getId() );
           }
#endif

           Timer.reset();
           Timer.start();

           LaMETQuda(   spinorOut  ,  spinorIn_q,  spinorIn_qbar, &quda_inv_param, &p_s_q,  &p_s_qbar);
           
           Timer.stop(); 
           CTCalls +=1;
           CTTime += Timer.getTimeInSeconds();



           for ( int i=0; i<lnk_max+1 ; i++)
           {
            out[i]=out[i]+psi_s[i];          
           }
          }

           QDPIO::cout << __func__ << ":  LG time = " << LGTime  << " seconds LGCalls = " << LGCalls  << std::endl;
           QDPIO::cout << __func__ << ":  CT time = " << CTTime  << " seconds CTCalls = " << CTCalls  << std::endl;
           

           TotalTime.stop();
           QDPIO::cout << __func__ << ": total time = " << TotalTime.getTimeInSeconds() << " seconds" << std::endl;
        }


void BkwdFrwdTrIOG_GPU( const LatticePropagator &             B,
                 const LatticePropagator &             F,
                 const multi1d< LatticeColorMatrix > & U,
                 int                                   GammaInsertion,
                 const SftMom_cohen &                        Phases,
                 const SftMom_cohen &                        PhasesCanonical,
                 multi1d< general_data_base > &        iog_data,
                 multi1d< int > &                      GBB_NLinkPatterns,
                 const int                             f,
                 const multi1d< unsigned short int > & LinkDirs,
                 const signed short int                T1,
                 const signed short int                T2,
                 const signed short int                Tsrc,
                 const signed short int                Tsnk,
                 const bool                            TimeReverse,
                 const bool                            ShiftFlag,
                 const int                             lnk_dir  ,
                 const unsigned short int              lnk_max ,
                 const short int                       gamma_id )
{


  StopWatch TotalTime;
  TotalTime.reset();
  TotalTime.start();

  StopWatch Timer;

  int TRCalls = 0;
  int FTCalls = 0;
  int GFGCalls = 0;
  int IPCalls = 0;
  double TRTime = 0.0;
  double FTTime = 0.0;
  double GFGTime = 0.0;
  double IPTime = 0.0;
  double IOTime = 0.0;

  const unsigned short int NLinks = LinkDirs.size();
  unsigned short int Link;
  const int NumQ = Phases.numMom();
  const int NumO = PhasesCanonical.numMom();
  const int NT   = Phases.numSubsets();  // Length of lattice in decay direction

  //##########################################################################//
  // add a tag to identify the link pattern                                   //
  //##########################################################################//

  Timer.reset();
  Timer.start();

  for ( int o = 0; o < NumO; o++ )
  {
    #if _DEBUG_BB_C_ == 1
    {
      QDPIO::cout << "DEBUG: " << __FILE__ << " " << __LINE__ << "\n";
      QDPIO::cout << "q = " << o << "\n";
      QDPIO::cout << "f = " << f << "\n";
      QDPIO::cout << "NLinks = " << NLinks << "\n";
      for ( Link = 0; Link < NLinks; Link++ )
      {
        QDPIO::cout << "Link = " << Link << "\n";
        QDPIO::cout << "LinkDirs[Link] = " << LinkDirs[Link] << "\n";
      }
    }
    #endif


    // counts number of link patterns per flavor
  }

  Timer.stop();
  IOTime += Timer.getTimeInSeconds();

  for( int dir=0 ; dir <2 ; dir++    )
  {
      std::vector< LatticeSpinMatrix >  out( lnk_max+1  );
     
      QDPIO::cout << __func__ << "  GPU_kernel setup  " << std::endl;
   
      GPU_kernel( B , F ,U ,out , lnk_dir+4*dir  , lnk_max  );
  
      QDPIO::cout << __func__ << " GPU_kernel done  " << std::endl;
  
      for ( int L=0+dir; L< lnk_max+1 ;L++     )
      {
        
        int gamma_ini, gamma_fin;
        
      
        if (gamma_id < 0 || gamma_id >15 ) 
        {        
            gamma_ini=0;
            gamma_fin=16;
        } 
        else
        {
            gamma_ini=gamma_id;
            gamma_fin=gamma_id+1;
        }   
     
    
        for ( int i = gamma_ini; i < gamma_fin  ; i++ )
        {
          Timer.reset();
          Timer.start();
  
          // assumes any Gamma5 matrices have already been absorbed
  
  
          // There is an overall minus sign from interchanging the initial and final
          // states for baryons.  This might not be present for mesons, so we should
          // think about this carefully.  It seems there should be another sign for
          // conjugating the operator, but it appears to be absent.  There is a minus
          // sign for all Dirac structures with a gamma_t.  In the current scheme this
          // is all gamma_i with i = 8, ..., 15.  If the gamma basis changes, then
          // this must change.
  
  
          Timer.stop();
          GFGCalls += 1;
          GFGTime += Timer.getTimeInSeconds();
          Timer.reset();
          Timer.start();
  
          Timer.stop();
          IPCalls += 1;
          IPTime += Timer.getTimeInSeconds();
          Timer.reset();
          Timer.start();
  
  //        LatticeComplex Trace=0;   
  //        for (    int site=all.start(); site <= all.end(); site++      )
  //    {Trace.elem(site).elem().elem()=out[L].elem(site).elem(i/4,i%4).elem();}
          LatticeComplex Trace=peekSpin( out[L] , i/4  ,  i%4  );
          Trace=adj(Trace);
          multi2d< DComplex > Projections = Phases.sft( Trace );
  
          Timer.stop();
          FTTime += Timer.getTimeInSeconds();
          FTCalls += 1;
  
          Timer.reset();
          Timer.start();
  
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
  //               if(ShiftFlag==true)t0=(t+Tsrc)%(T2-T1+1);
                 iog_data(f).data[0+2*(t+GBB_NLinkPatterns[f]*(T2-T1+1))]=real_part[t0];
                 iog_data(f).data[1+2*(t+GBB_NLinkPatterns[f]*(T2-T1+1))]=imag_part[t0];
               }
               GBB_NLinkPatterns[f]++;
          }    // for(idx_o, q)
  
          Timer.stop();
          IOTime += Timer.getTimeInSeconds();
        }//for gamma
  
    }//for length of z

  }//for plus/minus z
  
  QDPIO::cout << __func__ << ":  io time = " << IOTime << " seconds" << std::endl;
  QDPIO::cout << __func__ << ": gfg time = " << GFGTime / (double) GFGCalls << " seconds" << std::endl;
  QDPIO::cout << __func__ << ":  ip time = " << IPTime / (double) IPCalls << " seconds" << std::endl;
  QDPIO::cout << __func__ << ":  ft time = " << FTTime / (double) FTCalls << " seconds" << std::endl;
  TotalTime.stop();
  QDPIO::cout << __func__ << ": total time = " << TotalTime.getTimeInSeconds() << " seconds" << std::endl;

  return;

}
#endif

void BkwdFrwdTrIOG( const LatticePropagator &             B,
                 const LatticePropagator &             F,
                 const multi1d< LatticeColorMatrix > & U,
                 int                                   GammaInsertion,
                 const SftMom_cohen &                        Phases,
                 const SftMom_cohen &                        PhasesCanonical,
                 multi1d< general_data_base > &        iog_data,
                 multi1d< int > &                      GBB_NLinkPatterns,
                 const int                             f,
                 const multi1d< unsigned short int > & LinkDirs,
                 const signed short int                T1,
                 const signed short int                T2,
                 const signed short int                Tsrc,
                 const signed short int                Tsnk,
                 const bool                            TimeReverse,
                 const bool                            ShiftFlag,
                 const short int                       gamma_id    )
{
  StopWatch TotalTime;
  TotalTime.reset();
  TotalTime.start();

  StopWatch Timer;

  int TRCalls = 0;
  int FTCalls = 0;
  int GFGCalls = 0;
  int IPCalls = 0;
  double TRTime = 0.0;
  double FTTime = 0.0;
  double GFGTime = 0.0;
  double IPTime = 0.0;
  double IOTime = 0.0;

  const unsigned short int NLinks = LinkDirs.size();
  unsigned short int Link;
  const int NumQ = Phases.numMom();
  const int NumO = PhasesCanonical.numMom();
  const int NT   = Phases.numSubsets();  // Length of lattice in decay direction

  //##########################################################################//
  // add a tag to identify the link pattern                                   //
  //##########################################################################//

  Timer.reset();
  Timer.start();

  for ( int o = 0; o < NumO; o++ )
  {
    #if _DEBUG_BB_C_ == 1
    {
      QDPIO::cout << "DEBUG: " << __FILE__ << " " << __LINE__ << "\n";
      QDPIO::cout << "q = " << o << "\n";
      QDPIO::cout << "f = " << f << "\n";
      QDPIO::cout << "NLinks = " << NLinks << "\n";
      for ( Link = 0; Link < NLinks; Link++ )
      {
        QDPIO::cout << "Link = " << Link << "\n";
        QDPIO::cout << "LinkDirs[Link] = " << LinkDirs[Link] << "\n";
      }
    }
    #endif


    // counts number of link patterns per flavor
  }

  Timer.stop();
  IOTime += Timer.getTimeInSeconds();

   int gamma_ini, gamma_fin;


  if (gamma_id < 0 || gamma_id >15 )
  {
    gamma_ini=0;
    gamma_fin=16;
       } else
{
    gamma_ini=gamma_id;
    gamma_fin=gamma_id+1;
}


  for ( int i = gamma_ini; i < gamma_fin  ; i++ )
  {
    Timer.reset();
    Timer.start();

    // assumes any Gamma5 matrices have already been absorbed
    LatticePropagator GFG = Gamma(i) * F * Gamma( GammaInsertion );

    // There is an overall minus sign from interchanging the initial and final
    // states for baryons.  This might not be present for mesons, so we should
    // think about this carefully.  It seems there should be another sign for
    // conjugating the operator, but it appears to be absent.  There is a minus
    // sign for all Dirac structures with a gamma_t.  In the current scheme this
    // is all gamma_i with i = 8, ..., 15.  If the gamma basis changes, then
    // this must change.
    if ( ( TimeReverse == true ) & ( i < 8 ) ) GFG *= -1;

    Timer.stop();
    GFGCalls += 1;
    GFGTime += Timer.getTimeInSeconds();
    Timer.reset();
    Timer.start();

    LatticeComplex Trace = localInnerProduct( B, GFG );

    if(i==8 && InlineBuildingBlocksIOGEnv::conserved==true)
    {
        LatticePropagator Bp=shift( adj( U[3] ) * F, BACKWARD, 3 );
        LatticePropagator tmp_p=0.5*(Bp+Gamma(i)*Bp)*Gamma( GammaInsertion );
        Trace = localInnerProduct( B, tmp_p);
        Bp=shift( adj( U[3] ) * B, BACKWARD, 3 ); 
        tmp_p=0.5*(-Bp+Gamma(i)*Bp);
        Bp=F * Gamma( GammaInsertion );
        Trace+= localInnerProduct( tmp_p, Bp);
    }

    Timer.stop();
    IPCalls += 1;
    IPTime += Timer.getTimeInSeconds();
    Timer.reset();
    Timer.start();

    multi2d< DComplex > Projections = Phases.sft( Trace );

    Timer.stop();
    FTTime += Timer.getTimeInSeconds();
    FTCalls += 1;

    Timer.reset();
    Timer.start();

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

    Timer.stop();
    IOTime += Timer.getTimeInSeconds();
  }

  QDPIO::cout << __func__ << ":  io time = " << IOTime << " seconds" << std::endl;
  QDPIO::cout << __func__ << ": gfg time = " << GFGTime / (double) GFGCalls << " seconds" << std::endl;
  QDPIO::cout << __func__ << ":  ip time = " << IPTime / (double) IPCalls << " seconds" << std::endl;
  QDPIO::cout << __func__ << ":  ft time = " << FTTime / (double) FTCalls << " seconds" << std::endl;
  TotalTime.stop();
  QDPIO::cout << __func__ << ": total time = " << TotalTime.getTimeInSeconds() << " seconds" << std::endl;

  return;
}


void AddLinksIOG( const multi1d< LatticePropagator > &  B,
               const LatticePropagator &             F,
               const multi1d< LatticeColorMatrix > & U,
               const multi1d< int > &                GammaInsertions,
               const SftMom_cohen &                        Phases,
               const SftMom_cohen &                        PhasesCanonical,
               multi1d< unsigned short int > &       LinkDirs,
               const unsigned short int              MaxNLinks,
               BBLinkPattern                         LinkPattern,
               const short int                       PreviousDir,
               const short int                       PreviousMu,
               multi1d< general_data_base > &        iog_data,
               multi1d< int > &                      GBB_NLinkPatterns,
               const signed short int                T1,
               const signed short int                T2,
               const signed short int                Tsrc,
               const signed short int                Tsnk,
               const bool                            TimeReverse,
               const bool                            ShiftFlag,
               const int                             DecayDir,
               const int 			     lnk_dir,
               const int                             gamma_id)
{
  StopWatch Timer;
  int ShiftCalls = 0;
  double ShiftTime = 0.0;

  const unsigned short int NLinks = LinkDirs.size();

  if ( NLinks == MaxNLinks )
  {
    return;
  }

  LatticePropagator F_mu;
  const int NF = B.size();
  multi1d< unsigned short int > NextLinkDirs( NLinks + 1 );
  int Link;

  for ( Link = 0; Link < NLinks; Link++ )
  {
    NextLinkDirs[Link] = LinkDirs[Link];
  }

  // add link in forward mu direction
  for ( int mu = 0; mu < Nd; mu++ )
  if((lnk_dir<0&mu!=DecayDir)||(mu==lnk_dir))
  {
    // skip the double back
    if( ( PreviousDir == 0 ) || (( PreviousDir == 1 ) && ( PreviousMu == mu )) )
    {
      bool DoThisPattern = true;
      bool DoFurtherPatterns = true;

      NextLinkDirs[NLinks] = mu;

      LinkPattern( DoThisPattern, DoFurtherPatterns, NextLinkDirs );

      if ( DoFurtherPatterns == true )
      {
        Timer.reset();
        Timer.start();

        // accumulate product of link fields
        F_mu = shift( adj( U[mu] ) * F, BACKWARD, mu );

        Timer.stop();
        ShiftTime += Timer.getTimeInSeconds();
        ShiftCalls += 1;
         QDPIO::cout << "do the direction " << mu << "..." << std::endl;
      }

      if ( DoThisPattern == true )
      {
        // form correlation functions
        for ( int f = 0; f < NF; f++ )
        {
          BkwdFrwdTrIOG( B[f], F_mu, U, GammaInsertions[f], Phases, PhasesCanonical,
                      iog_data, GBB_NLinkPatterns,
                      f, NextLinkDirs, T1, T2, Tsrc, Tsnk,
                      TimeReverse, ShiftFlag, gamma_id);
        }
      }

      if ( DoFurtherPatterns == true )
      {
        // add another link
        AddLinksIOG( B, F_mu, U, GammaInsertions,
                  Phases, PhasesCanonical,
                  NextLinkDirs, MaxNLinks, LinkPattern, 1, mu,
                  iog_data, GBB_NLinkPatterns,
                  T1, T2, Tsrc, Tsnk, TimeReverse, ShiftFlag, DecayDir ,lnk_dir, gamma_id );
      }
    }
  }

  // add link in backward mu direction
  for ( int mu = 0; mu < Nd; mu++ )
  if((lnk_dir<0&mu!=DecayDir)||(mu==lnk_dir))
  {
    // skip the double back
    if( ( PreviousDir == 0 ) || (( PreviousDir == -1 ) && ( PreviousMu == mu )) )
    {
      bool DoThisPattern = true;
      bool DoFurtherPatterns = true;

      NextLinkDirs[NLinks] = mu + Nd;

      LinkPattern( DoThisPattern, DoFurtherPatterns, NextLinkDirs );

      if ( DoFurtherPatterns == true )
      {
        Timer.reset();
        Timer.start();

        // accumulate product of link fields
        F_mu = U[mu] * shift( F, FORWARD, mu );

        Timer.stop();
        ShiftTime += Timer.getTimeInSeconds();
        ShiftCalls += 1;
         QDPIO::cout << "do the direction " << mu+4 << "..." << std::endl;
      }

      if ( DoThisPattern == true )
      {
        // form correlation functions
        for ( int f = 0; f < NF; f++ )
        {
          BkwdFrwdTrIOG( B[f], F_mu, U, GammaInsertions[f], Phases, PhasesCanonical,
                      iog_data, GBB_NLinkPatterns,
                      f, NextLinkDirs, T1, T2, Tsrc, Tsnk,
                      TimeReverse, ShiftFlag, gamma_id );
        }
      }

      if ( DoFurtherPatterns == true )
      {
        // add another link
        AddLinksIOG( B, F_mu, U, GammaInsertions, Phases, PhasesCanonical,
                  NextLinkDirs, MaxNLinks, LinkPattern, -1, mu, iog_data,
                  GBB_NLinkPatterns,
                  T1, T2, Tsrc, Tsnk, TimeReverse, ShiftFlag, DecayDir ,lnk_dir, gamma_id );
      }
    }
  }

  QDPIO::cout << __func__ << ": shift time = "
        << ShiftTime
        << " seconds with shift calls = "
              << ShiftCalls
              << std::endl;

  return;
}

#ifdef BUILD_QUDA
void BuildingBlocksIOG_GPU( const multi1d< LatticePropagator > &  B,
                     const LatticePropagator &             F,
                     const multi1d< LatticeColorMatrix > & U,
                     const multi1d< int > &                GammaInsertions,
                     const multi1d< int > &                Flavors,
                     const unsigned short int              MaxNLinks,
                     const BBLinkPattern                   LinkPattern,
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
                     const int lnk_dir, int icfg,
                     const int                             gamma_id)
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

    int n_dir=(lnk_dir<0)?6:2;
    for ( int f = 0; f < NumF; f++ )
    {
      sprintf(iog_data(f).name,"%s",BinaryDataFileNames(f, 0).c_str());
      iog_data(f).add_dimension(dim_conf,1,&icfg);
      iog_data(f).add_dimension(dim_channel,1,&f);
      int links_id[1+6*MaxNLinks];
      links_id[0]=0;
      int icount=1;
      for(int idr=0;idr<Nd;idr++)
      if((lnk_dir<0&&idr!=DecayDir)||(idr==lnk_dir))
      for(int is=0;is<MaxNLinks;is++)
         links_id[icount++]=idr*100+is+1;
      for(int idr=0;idr<Nd;idr++)
      if((lnk_dir<0&&idr!=DecayDir)||(idr==lnk_dir))
      for(int is=0;is<MaxNLinks;is++)
         links_id[icount++]=(idr+4)*100+is+1;
      iog_data(f).add_dimension(dim_displacement,1+n_dir*MaxNLinks,links_id);
       if (gamma_id < 0 || gamma_id >15 ){iog_data(f).add_dimension(dim_operator,Ns*Ns);}
      else{iog_data(f).add_dimension(dim_operator,1,&gamma_id);}          
      const int NumQ = Phases.numMom();
      int mom_list[NumQ];
      icount=0;
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
      if(Layout::primaryNode()) iog_data(f).initialize();

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

  QDPIO::cout << __func__ << ": start BkwdFrwdTr" << std::endl;

  const unsigned short int NLinks = 0;
  multi1d< unsigned short int > LinkDirs( 0 );

  for ( int f = 0; f < NumF; f++ )
  {

    QDPIO::cout << __func__ << ": start BkwdFrwdTrIOG_GPU " << std::endl;

      BkwdFrwdTrIOG_GPU( B[f], F, U, GammaInsertions[f], Phases, PhasesCanonical,
                   iog_data, GBB_NLinkPatterns, f, LinkDirs,
                   T1, T2, Tsrc, Tsnk,
                   TimeReverse, ShiftFlag, lnk_dir, MaxNLinks, gamma_id);
  }

  Timer.stop();
  QDPIO::cout << __func__ << ": total time for all links (single BkwdFrwdTr call) = "
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
  const unsigned short int NX = Layout::lattSize()[0];
  const unsigned short int NY = Layout::lattSize()[1];
  const unsigned short int NZ = Layout::lattSize()[2];
  const unsigned short int NT = Layout::lattSize()[3];
  const signed short int   PX = SnkMom[0];
  const signed short int   PY = SnkMom[1];
  const signed short int   PZ = SnkMom[2];
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
#endif

void BuildingBlocksIOG( const multi1d< LatticePropagator > &  B,
                     const LatticePropagator &             F,
                     const multi1d< LatticeColorMatrix > & U,
                     const multi1d< int > &                GammaInsertions,
                     const multi1d< int > &                Flavors,
                     const unsigned short int              MaxNLinks,
                     const BBLinkPattern                   LinkPattern,
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
                     const int lnk_dir, int icfg,
                     const int                             gamma_id)
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

    int n_dir=(lnk_dir<0)?6:2;
    for ( int f = 0; f < NumF; f++ )
    {
      sprintf(iog_data(f).name,"%s",BinaryDataFileNames(f, 0).c_str());
      iog_data(f).add_dimension(dim_conf,1,&icfg);
      iog_data(f).add_dimension(dim_channel,1,&f);
      int links_id[1+6*MaxNLinks];
      links_id[0]=0;
      int icount=1;
      for(int idr=0;idr<Nd;idr++)
      if((lnk_dir<0&&idr!=DecayDir)||(idr==lnk_dir))
      for(int is=0;is<MaxNLinks;is++)
         links_id[icount++]=idr*100+is+1;
      for(int idr=0;idr<Nd;idr++)
      if((lnk_dir<0&&idr!=DecayDir)||(idr==lnk_dir))
      for(int is=0;is<MaxNLinks;is++)
         links_id[icount++]=(idr+4)*100+is+1;
      iog_data(f).add_dimension(dim_displacement,1+n_dir*MaxNLinks,links_id);
      if (gamma_id < 0 || gamma_id >15 ){iog_data(f).add_dimension(dim_operator,Ns*Ns);}
      else{iog_data(f).add_dimension(dim_operator,1,&gamma_id);}
      const int NumQ = Phases.numMom();
      int mom_list[NumQ];
      icount=0;
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
      if(Layout::primaryNode()) iog_data(f).initialize();

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

  QDPIO::cout << __func__ << ": start BkwdFrwdTr" << std::endl;

  const unsigned short int NLinks = 0;
  multi1d< unsigned short int > LinkDirs( 0 );

  for ( int f = 0; f < NumF; f++ )
  {
    BkwdFrwdTrIOG( B[f], F, U, GammaInsertions[f], Phases, PhasesCanonical,
                iog_data, GBB_NLinkPatterns, f, LinkDirs, 
                T1, T2, Tsrc, Tsnk,
                TimeReverse, ShiftFlag,gamma_id);
  }

  Timer.stop();
  QDPIO::cout << __func__ << ": total time for 0 links (single BkwdFrwdTr call) = "
              << Timer.getTimeInSeconds() 
              << " seconds" << std::endl;

  Timer.reset();
  Timer.start();

  QDPIO::cout << __func__ << ": start AddLinks" << std::endl;

  AddLinksIOG( B, F, U, GammaInsertions, 
            Phases, PhasesCanonical,
            LinkDirs, MaxNLinks, LinkPattern, 0, -1, 
            iog_data, GBB_NLinkPatterns, 
            T1, T2, Tsrc, Tsnk, TimeReverse, ShiftFlag, DecayDir, lnk_dir, gamma_id);
            
  QDPIO::cout << "total data size is (" << GBB_NLinkPatterns[0] <<","
              << (1+n_dir*MaxNLinks)*Ns*Ns*Phases.numMom() <<std::endl;

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
  const unsigned short int NX = Layout::lattSize()[0];
  const unsigned short int NY = Layout::lattSize()[1];
  const unsigned short int NZ = Layout::lattSize()[2];
  const unsigned short int NT = Layout::lattSize()[3];
  const signed short int   PX = SnkMom[0];
  const signed short int   PY = SnkMom[1];
  const signed short int   PZ = SnkMom[2];
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


/*! \file
 * \brief Inline construction of BuildingBlocks
 *
 * Building Blocks on forward and sequential props
 */
#include "inline_ex/inline_LaMET.h"
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
#endif


namespace Chroma
{

//------------------------------------------------------------------------------
// Factory Registration


  // Environment in which the measurement lives (holds params and such)
  namespace InlineLaMETEnv
  {
    bool conserved =false;
    namespace
    {
      // Function to register with a factory
      AbsInlineMeasurement* createMeasurement( XMLReader& xml_in,
                                       const std::string& path )
      {
        return new InlineLaMET(
                        InlineLaMETParams(xml_in, path) );
      }

      //! Local registration flag
      bool registered = false;
    }

    // The name of the measurement for the XML file
    const std::string name = "LaMET";

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
  } // end namespace InlineLaMETEnv

//------------------------------------------------------------------------------
// Parameter reading, writing, handling

  //! Reader for parameters
  void read(                         XMLReader& xml,
                             const std::string& path,
           InlineLaMETParams::Param_t& param )
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
    case 1:
      read(paramtop, "translate", param.translate);
      break;

    default :
      QDPIO::cerr << InlineLaMETEnv::name
                  << ": input parameter version " << version << " unsupported."
                  << std::endl;
      QDP_abort(1);
    }

    read(paramtop, "cfg_serial", param.icfg);
    multi1d<int> tmp;
    read(paramtop, "pre_shift", tmp);
    param.pre_shift.load(tmp);
    if(param.pre_shift.check("pre_shift")!=0)
    {
      QDPIO::cerr << InlineLaMETEnv::name
                  << ": pre_shift pattern is incorrect."
                  << std::endl;
      QDP_abort(1);
    }

    multi1d< multi1d<int> > tmp2;    
    read(paramtop, "post_shift", tmp2);
    param.post_shift.resize(tmp2.size());
    for(int i=0;i<tmp2.size();i++)
    {
        param.post_shift[i].load(tmp2(i));
        if(param.post_shift[i].check())
        {
           QDPIO::cerr << InlineLaMETEnv::name
                  << ": post_shift pattern (" << i << ") is incorrect."
                  << std::endl;
           QDP_abort(1);
        }
    }

    read(paramtop, "mom_list", param.mom_list);
    read(paramtop, "npr", param.npr);
    if(param.npr==true && mom_list.size()!=1) 
    {
      QDPIO::cerr << InlineLaMETEnv::name
                  << ": momentum of npr should be fixed to given value."
                  << std::endl;
      QDP_abort(1);
    }

    read(paramtop, "corr_t", param.corr_t);
#ifdef BUILD_QUDA    
    read(paramtop, "qudaLaMET", param.qudaLaMET);
#endif
    param.qudaLaMET=false;
#endif        

  }

  //! Writer for parameters
  void write(                              XMLWriter& xml,
                                   const std::string& path,
           const InlineLaMETParams::Param_t& param )
  {
    push(xml, path);

    int version = 1;
    write(xml, "version", version);
    write(xml, "pre_shift", param.pre_shift);
    write(xml, "post_shift", param.post_shift);
    write(xml, "mom_list", param.mom_list);
    write(xml, "translate", param.translate);
    xml << param.cfs.xml;

    pop(xml);
  }

  //! BB parameters (gauge id, props) input
  void read(                      XMLReader& xml,
                          const std::string& path,
           InlineLaMETParams::BB_t& bb )
  {
    XMLReader bbtop(xml, path);

    read(bbtop, "GaugeId", bb.GaugeId);
    read(bbtop, "FrwdPropId", bb.FrwdPropId);
    read(bbtop, "BkwdProps", bb.BkwdProps);
  }

  //! BB parameters (gauge id, props) output
  void write(                           XMLWriter& xml,
                                const std::string& path,
           const InlineLaMETParams::BB_t& bb )
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
           InlineLaMETParams::BkwdProp_t& bkwdprop )
  {
    XMLReader bkwdproptop(xml, path);

    read(bkwdproptop, "BkwdPropId", bkwdprop.BkwdPropId);
    read(bkwdproptop, "BkwdPropG5Format", bkwdprop.BkwdPropG5Format);
    read(bkwdproptop, "BBFileNamePattern", bkwdprop.BBFileNamePattern);
    read(paramtop, "GammaInsertion", bkwdprop.GammaInsertion);
  }

  //! Correlator description "backward prop" (flavor, gamma, props, etc.) output
  void write(                                 XMLWriter& xml,
                                      const std::string& path,
           const InlineLaMETParams::BkwdProp_t& bkwdprop )
  {
    push(xml, path);

    write(xml, "BkwdPropId", bkwdprop.BkwdPropId);
    write(xml, "BkwdPropG5Format", bkwdprop.BkwdPropG5Format);
    write(xml, "BBFileNamePattern", bkwdprop.BBFileNamePattern);

    pop(xml);
  }

  // Params default constructor
  InlineLaMETParams::InlineLaMETParams()
  {
    frequency = 0;
  }

  // Construct params from XML
  InlineLaMETParams::InlineLaMETParams(
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
  void InlineLaMETParams::write( XMLWriter& xml_out,
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


  // Set up the XML and invoke func, which does the acual work
  void InlineLaMET::operator()( unsigned long update_no,
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


  // Real work done here
  void InlineLaMET::func( unsigned long update_no,
                                             XMLWriter& xml_out)
  {
    START_CODE();

    StopWatch snoop;
    snoop.reset();
    snoop.start();

    //--------------------------------------------------------------------------
    // Start building the output XML
    push(xml_out, "LaMET");
    write(xml_out, "update_no", update_no);

    QDPIO::cout << " LaMET" << std::endl;
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
      QDPIO::cerr << InlineLaMETEnv::name
                  << ": caught dynamic cast error" << std::endl;
      QDP_abort(1);
    }
    catch (const std::string& e)
    {
      // If the object is not in the map we will land here
      QDPIO::cerr << InlineLaMETEnv::name
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
      QDPIO::cerr << InlineLaMETEnv::name
                  << ": caught dynamic cast error" << std::endl;
      QDP_abort(1);
    }
    catch (const std::string& e)
    {
      QDPIO::cerr << InlineLaMETEnv::name
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

    int Bsize=params.bb.BkwdProps.size();
    multi1d<LatticePropagator> B(Bsize);
    for (int loop = 0; loop < params.bb.BkwdProps.size(); ++loop)
    {
      push(xml_out, "elem");
      write(xml_out, "loop_ctr", loop);

      QDPIO::cout << "Backward propagator #" << loop << std::endl;

      SeqSource_t seqsource_header;
      QDPIO::cout << "Attempt to parse backward propagator" << std::endl;
      try
      {
        // Grab the backward prop
        B[loop] = TheNamedObjMap::Instance().getData<LatticePropagator>(
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
        multi1d<Double> BkwdPropCheck =
          sumMulti( localNorm2( B[loop] ), phases_nomom.getSet() );

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
        QDPIO::cerr << InlineLaMETEnv::name
                    << ": forward prop: caught dynamic cast error" << std::endl;
        QDP_abort(1);
      }
      catch (const std::string& e)
      {
        QDPIO::cerr << InlineLaMETEnv::name
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
        QDPIO::cerr << "InlineLaMET: Gamma insertion out of bounds: "
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
        LatticePropagator Bu = B[loop];
        B[loop] = Gamma(Ns*Ns-1)*Bu;
      }
      else if (params.bb.BkwdProps[loop].BkwdPropG5Format == "B_G5")
      {
        LatticePropagator Bu = B[loop];
        B[loop] = Bu*Gamma(Ns*Ns-1);
      }
      else if (params.bb.BkwdProps[loop].BkwdPropG5Format == "G5_B_G5")
      {
        LatticePropagator Bu = B[loop];
        B[loop] = Gamma(Ns*Ns-1)*Bu*Gamma(Ns*Ns-1);
      }
      
      LatticePropagator Bu = B[loop];
      B[loop] = Bu*adj(Gamma(params.bb.BkwdProps[loop].GammaInsertion));
   }//end loop of B

      // Invoke the master building blocks function!
      swatch.reset();
      QDPIO::cout << "calculating building blocks" << std::endl;

      const signed short int DecayDir = j_decay;
      const signed short int Tsrc = source_header.t_source;
      const signed short int Tsnk = seqsource_header.t_sink;

      // generate the momentum list;
      
      multi2d<std::string> Files(Bsize);
      for(int loop=0;loop<Bsize;loop++)
          Files(loop) = params.bb.BkwdProps[loop].BBFileNamePattern;
          
      multi1d< general_data_base > iog(Bsize);
      
      swatch.start();
      if(params.param.npr==true)   
      {
          //save forward prop;, apply the momentum;
          
          //determine s_q and  s_qbar
      
          LaMET_kernel(B, F, U, iog, 
                   s_q, s_qbar, 
                   params.param.mom_list, 
                   params.param.corr_t,
                   params.param.qudaLaMET,
                   params.param.translate);
                   params.param.icfg);
                   
          // save the data;
      }
      else
      {
          LaMET_kernel(B, F, U, iog, 
                   s_q, s_qbar, 
                   params.param.mom_list, 
                   params.param.corr_t,
                   params.param.qudaLaMET,
                   params.param.translate);
                   params.param.icfg);
      
      
      }
      swatch.stop();
     
    pop(xml_out);  // SequentialSource

    pop(xml_out);  // BuildingBlocks

    snoop.stop();
    QDPIO::cout << InlineLaMETEnv::name << ": total time = "
                << snoop.getTimeInSeconds()
                << " secs" << std::endl;

    QDPIO::cout << InlineLaMETEnv::name << ": ran successfully"
                << std::endl;

    END_CODE();
  } // end of InlineLaMET::func

  void LaMETChroma(
                   multi1d< LatticeSpinMatrix > &psi_s,
                   multi1d< Real > corr,
                   const LatticePropagator&  B,
                   const LatticePropagator&  F,
                   multi1d< LatticeColorMatrix > & U,
                   shiftPatternChroma &s_q,
                   shiftPatternChroma &s_qbar)
  {
        
     if(s_qbar.post_shift.size()>0) if(Layout::primaryNode()) printf("Only the post_shift of quark is supported at present\n");
     s_q.check();
     s_qbar.check();
     
     LatticePropagator qbar=B;
     LatticePropagator q=F;
     LatticePropagator tmp;
     
     for(int i=0;i<s_q.pre_shift.size();i+=2)
     for(int j=0;j<s_q.pre_shift[i];j++)
     {
         int LinksDirection=s_q.pre_shift[i+1];
         if(LinksDirection<Nd)
            tmp=shift( adj( U[ LinksDirection ] ) * q, BACKWARD, LinksDirection );
         else
            tmp=U[ LinksDirection-Nd ] * shift( q, FORWARD, LinksDirection-Nd );
         q = tmp;
     }

     for(int i=0;i<s_qbar.pre_shift.size();i+=2)
     for(int j=0;j<s_qbar.pre_shift[i];j++)
     {
         int LinksDirection=s_qbar.pre_shift[i+1];
         if(LinksDirection<Nd)
            tmp=shift( adj( U[ LinksDirection ] ) * qbar, BACKWARD, LinksDirection );
         else
            tmp=U[ LinksDirection-Nd ] * shift( qbar, FORWARD, LinksDirection-Nd );
         qbar = tmp;
     }
     
     for(int i=0;i<s_q.post_shift.size();i+=2)
     for(int j=0;j<s_q.post_shift[i];j++)
     {
         int LinksDirection=s_q.post_shift[i+1];
         if(LinksDirection<Nd)
            tmp=shift( adj( U[ LinksDirection ] ) * q, BACKWARD, LinksDirection );
         else
            tmp=U[ LinksDirection-Nd ] * shift( q, FORWARD, LinksDirection-Nd );
         q = tmp;
     }

  }

void LaMET_kernel( const LatticePropagator &             F,
                 const multi1d< LatticePropagator > &    B,
                 const multi1d< LatticeColorMatrix > &   U,
                 const multi1d< general_data_base > &   iog,
                 shiftPatternChroma &s_q,
                 shiftPatternChroma &s_qbar,
                 multi1d<int> &mom_list,
                 bool npr, bool corr_t, bool qudaLaMET,
                 bool translate, int icfg)
{

   //check input; 
   StopWatch TotalTime;
   TotalTime.reset();
   TotalTime.start();

   StopWatch Timer;

   double setupTime=0.0;
   double QudaTime=0.0;

   int link_max=1;
   for(int i=0;i<s_q.post_shift.size();i+=2)
         link_max+=s_q.post_shift[i];

   int Bsize=B.size();
   multi1d<LatticeFermion> tmp_F(12);
   multi1d<LatticeFermion> tmp_B(12*Bsize);
   std::vector< LatticeSpinMatrix > psi_s(lnk_max+1); //for corr
   std::vector< Real > corr((link_max+1)*16) // for 
   std::vector< void *  > spinorin_q(12);
   std::vector< void *  > spinorin_qbar(12*Bsize);
   
   for(int is=0;id<4;id++)
   for(int ic=0;ic<3;ic++)
   {
       PropToFerm(F,tmp_F[is*3+ic],ic,is);
       
#ifndef BUILD_QUDA_DEVIFACE_SPINOR
       spinorIn_q[is*3+ic]    =(void *)&(tmp_F.elem(0).elem(0).elem(0).real());
#else
       spinorIn_q[is*3+ic]    = GetMemoryPtr( tmp_F.getId() );
#endif       
       for(int ib=0;ib<Bsize;ib++)
       {
          PropToFerm(B,tmp_B[is*3+ic],ic,is);
#ifndef BUILD_QUDA_DEVIFACE_SPINOR
          spinorIn_qbar[ib*12+is*3+ic] =(void *)&(tmp_B.elem(0).elem(0).elem(0).real());
#else
          spinorIn_qbar[ib*12+is*3+ic] = GetMemoryPtr( tmp_B.getId() );
#endif      
       }
   }

   std::vector< void *  > spinorOut( lnk_max+1 );
   if(npr==false)
   for (int i=0; i<lnk_max+1 ;i++)
   {
#ifndef BUILD_QUDA_DEVIFACE_SPINOR
       spinorOut[i] =(void *)&(psi_s[i].elem(0).elem(0,0).elem().real());   ;
#else
       spinorOut[i] = GetMemoryPtr( psi_s[i].getId() );
#endif
   }
   else
       spinorOut[i] = corr.data()+i*16;

   Timer.reset();
   Timer.start();
   
   if(qudaLaMET==false)
       LaMETChroma(psi_s, corr ,F  ,B, U, &s_q,  &s_qbar);   
   else
   {
#ifdef BUILD_QUDA
       QudaPrecision_s cpu_prec;
       QudaPrecision_s gpu_prec;

       int s = sizeof( WordType<Real>::Type_t );
       if (s == 4) {
         cpu_prec = QUDA_SINGLE_PRECISION;
       }
       else {
         cpu_prec = QUDA_DOUBLE_PRECISION;
       }

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
       setupTime += Timer.getTimeInSeconds();

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
           
       LaMETQuda(   spinorOut  ,  spinorIn_q,  spinorIn_qbar, &quda_inv_param, &s_q,  &s_qbar);
#endif // endif BUILD_QUDA
       
   }
       Timer.stop(); 
       CTTime = Timer.getTimeInSeconds();
      
   if(npr==false)
   {
      // get 2pt
      
      
      
   }

   QDPIO::cout << __func__ << ": Calculation time = " << CTTime  << " seconds "  << std::endl;
           
   TotalTime.stop();
   QDPIO::cout << __func__ << ": total time = " << TotalTime.getTimeInSeconds() << " seconds" << std::endl;
}


//############################################################################//

}  // end namespace Chroma

    

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
 
#include "inline_ex/inline_vev.h"
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
  namespace InlineVEVEnv
  {
    namespace
    {
      // Function to register with a factory
      AbsInlineMeasurement* createMeasurement( XMLReader& xml_in,
                                       const std::string& path )
      {
        return new InlineVEV( InlineVEVParams(xml_in, path) );
      }

      //! Local registration flag
      bool registered = false;
    }
    
    // The name of the measurement for the XML file
    const std::string name = "VEV";

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
  } // end namespace InlineVEVEnv

//------------------------------------------------------------------------------
// Parameter reading, writing, handling

  //! Reader for parameters
  void read(             XMLReader& xml,
                 const std::string& path,
    InlineVEVParams::Param_t& param )
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "links_dir", param.links_dir);
    read(paramtop, "links_max", param.links_max);
    read(paramtop, "t_src", param.t_src);
    read(paramtop, "t_max", param.t_max);
    read(paramtop, "file_name", param.file_name);
    read(paramtop, "cfg_serial", param.cfg_serial);
  }

  //! Writer for parameters
  void write(                  XMLWriter& xml,
                       const std::string& path,
    const InlineVEVParams::Param_t& param )
  {
    push(xml, path);

    write(xml, "links_dir", param.links_dir);
    write(xml, "links_max", param.links_max);
    write(xml, "t_src", param.t_src);
    write(xml, "t_max", param.t_max);
    write(xml, "file_name", param.file_name);
    write(xml, "cfg_serial", param.cfg_serial);
    pop(xml);
  }

  //! Named objects (gauge config, prop pairs) input
  void read(                   XMLReader& xml,
                       const std::string& path,
    InlineVEVParams::NamedObject_t& input )
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "gauge_id", input.gauge_id);
    read(inputtop, "quark_prop", input.quark_prop);
  }

  //! Named objects (gauge config, prop pairs) output
  void write(                        XMLWriter& xml,
                             const std::string& path,
    const InlineVEVParams::NamedObject_t& input )
  {
    push(xml, path);

    write(xml, "gauge_id", input.gauge_id);
    write(xml, "quark_prop", input.quark_prop);

    pop(xml);
  }

  // Params default constructor
  InlineVEVParams::InlineVEVParams(void)
  {
    frequency = 0;
  }

  // Construct params from XML
  InlineVEVParams::InlineVEVParams( XMLReader& xml_in,
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
  void InlineVEVParams::write( XMLWriter& xml_out,
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
//------------------------------------------------------------------------------
// The real work is done here

  // Set up the XML and invoke func, which does the actual work
  void InlineVEV::operator()( unsigned long update_no,
                                       XMLWriter& xml_out )
  {
    // If xml file not empty, then use alternate
    if (params.xml_file != "")
    {
      std::string xml_file = makeXMLFileName(params.xml_file, update_no);

      push(xml_out, "VEV");
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
  void InlineVEV::func( unsigned long update_no,
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
      QDPIO::cerr << InlineVEVEnv::name << ": caught dynamic cast error"
                  << std::endl;
      QDP_abort(1);
    }
    catch (const std::string& e)
    {
      // If the object is not in the map we will land here
      QDPIO::cerr << InlineVEVEnv::name << ": map call failed: " << e
                  << std::endl;
      QDP_abort(1);
    }
    // Bind gauge field reference to a local name
    const multi1d<LatticeColorMatrix>& u =
      TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);

    // Boilerplate stuff to the output XML
    push(xml_out, "VEV");
    write(xml_out, "update_no", update_no);

    // Write some diagnostic junk about the lattice
    QDPIO::cout << " VEV: Vaccum expectation values of quasi-PDF operator" << std::endl;
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
    push(xml_out, "vev_measurements");

    LatticePropagator quark_prop;
    quark_prop= TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.quark_prop);
   //Create phase;
   LatticeReal phase=where(Layout::latticeCoordinate(3)==params.param.t_src,LatticeReal(1.0),LatticeReal(zero)); 
   
    int links_max=params.param.links_max;
   std::string name1=params.param.file_name+".vev1";
   general_data_base res(name1.c_str());
   res.add_dimension(dim_conf, 1, &params.param.cfg_serial);
   res.add_dimension(dim_displacement,2*links_max+1);
   for(int i=0;i<2*links_max+1;i++)
       res.dim[1].indices[i]=-links_max+i;
   res.add_dimension(dim_operator,2);
   res.dim[2].indices[0]=0;res.dim[2].indices[1]=8;
   res.add_dimension(dim_complex, 2);
   if(Layout::primaryNode()) res.initialize();
   
    // Now loop over the various fermion pairs
    for(int idr=0;idr<2;idr++)
    {
        LatticePropagator F=quark_prop,F_tmp;
        for(int ir=0;ir<links_max+1;ir++)
        if(ir!=0||idr==0)
        {
           int off=(idr==0)?(links_max+ir):(links_max-ir);
           for(int ig=0;ig<2;ig++)
           {
                DComplex value=sum(trace(F*Gamma(0+8*ig))*phase);
                if(Layout::primaryNode())
                {
                   res.data[4*off+0+2*ig]=value.elem().elem().elem().real();
                   res.data[4*off+1+2*ig]=value.elem().elem().elem().imag();
                }
           }
           if(idr==0&&ir!=links_max)
           {
                F_tmp=U[3] * shift( F, FORWARD, 3);
                F=F_tmp;
           }
           if(idr==1&&ir!=links_max)
           {
                F_tmp=shift( adj(U[3]) * F, BACKWARD, 3);
                F=F_tmp;
           }
        }
    }
    if(Layout::primaryNode()) res.save();

   std::string name2=params.param.file_name+".vev2";
   int nlnk=(params.param.links_dir>=0&&params.param.links_dir<3)?1:3;
   general_data_base res2(name2.c_str());
   res2.add_dimension(dim_conf, 1, &params.param.cfg_serial);
   int nt = Layout::lattSize()[3];
   res2.add_dimension(dim_channel,nlnk*2);
   res2.add_dimension(dim_displacement,links_max);
   res2.add_dimension(dim_operator,16);
   res2.add_dimension(dim_t,nt);
   res2.add_dimension(dim_complex, 2);
   if(Layout::primaryNode()) res2.initialize();

   int dr_st=(nlnk==1)?params.param.links_dir:0;
   int dr_ed=(nlnk==1)?(params.param.links_dir+1):3;
   
      SftMom phases_nomom(0, true, 3);
      Set timeslice= phases_nomom.getSet();
   
   for(int dir=dr_st;dir<dr_ed;dir++)
   {
       //create the tempoary link;
       LatticeColorMatrix lnk=U[dir]*phase;
       for(int idr=0;idr<2;idr++)
       {
            LatticeColorMatrix lnk_tmp=U[dir]*phase,lnk_tmp2;
            for(int i=1;i<=params.param.t_max;i++)
            {
                   lnk_tmp2=shift( lnk_tmp, (idr==0)?FORWARD:BACKWARD, 3);
                   lnk_tmp = lnk_tmp2;
                   lnk_tmp2 = lnk_tmp + lnk;
                   lnk = lnk_tmp2;
            }
       }
       for(int idr=0;idr<2;idr++)
       {
          LatticePropagator F=Gamma(15)*quark_prop*Gamma(15),F_tmp;
          for(int iz=0;iz<links_max;iz++)
          {
            int off=((dir-dr_st)*2+idr)*links_max+iz;
            if(idr==0)
            {
                  F_tmp=U[dir] * shift( F, FORWARD, dir) * adj(lnk);
                  F=F_tmp;
            }
            if(idr==1)
            {
                  F_tmp=shift( adj(U[dir])* F * lnk, BACKWARD, dir);
                  F=F_tmp;
            }
            for(int ig=0;ig<16;ig++)
            {
                 LatticeComplex corr=trace(adj(F)*Gamma(ig)*quark_prop*Gamma(ig));
                 multi1d<DComplex> hsum = sumMulti(corr,timeslice);
                 if(Layout::primaryNode())
                 for(int it=0;it<nt;it++)
                 {
                      res2.data[(off*16+ig)*nt*2+it*2+0]=hsum[it].elem().elem().elem().real();
                      res2.data[(off*16+ig)*nt*2+it*2+1]=hsum[it].elem().elem().elem().imag();
                 }
            }
          }
        }
    }
       

    if(Layout::primaryNode()) res2.save();
    pop(xml_out);  // VEV

    snoop.stop();
    QDPIO::cout << InlineVEVEnv::name << ": total time = "
                << snoop.getTimeInSeconds()
                << " secs" << std::endl;

    QDPIO::cout << InlineVEVEnv::name << ": ran successfully" << std::endl;

    END_CODE();
  } // end of InlineVEV::func

}; // end of namespace Chroma


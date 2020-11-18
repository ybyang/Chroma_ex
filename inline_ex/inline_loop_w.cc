// $Id: inline_propagator_w.cc,v 3.11 2007-10-13 20:46:29 edwards Exp $
/*! \file
 * \brief Inline construction of propagator
 *
 * Propagator calculations
 */

#include "fermact.h"
#include "inline_ex/inline_loop_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "util/info/unique_id.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"
#include "meas/inline/make_xml_file.h"

#include "meas/inline/io/named_objmap.h"

#include "util/ferm/transf.h"
#include<time.h>
#include "io_ex/io_general_class.h"

namespace Chroma 
{ 
  namespace InlineLoopEnv 
  { 
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineLoop(InlineLoopParams(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "LOOP";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= WilsonTypeFermActsEnv::registerAll();
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	registered = true;
      }
      return success;
    }
  } // end namespace


  //! Propagator input
  void read(XMLReader& xml, const std::string& path, InlineLoopParams::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "gauge_id", input.gauge_id);
  }

  //! Propagator output
  void write(XMLWriter& xml, const std::string& path, const InlineLoopParams::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "gauge_id", input.gauge_id);

    pop(xml);
  }


  // Param stuff
  InlineLoopParams::InlineLoopParams() { frequency = 0; }

  InlineLoopParams::InlineLoopParams(XMLReader& xml_in, const std::string& path) 
  {
    try 
    {
      XMLReader paramtop(xml_in, path);

      if (paramtop.count("Frequency") == 1)
	read(paramtop, "Frequency", frequency);
      else
	frequency = 1;

      // Parameters for source construction
      read(paramtop, "Grid", grid);
      if (grid.size() != 3) {
        QDPIO::cerr << InlineLoopEnv::name << ": wrong size of grid array: expected length=" << 3
        << std::endl;
        QDP_abort(1);
      }

      read(paramtop, "link_lengths", link_lengths);
      if (link_lengths.size() != 3) {
        QDPIO::cerr << InlineLoopEnv::name << ": wrong size of link_length array: expected length=" << 3
        << std::endl;
        QDP_abort(1);
      }

      if (paramtop.count("use_ckpoint") != 0)
      {
         read(paramtop, "use_ckpoint", use_ckpoint);
      }
      else
      {
          use_ckpoint=true;
      }
      read(paramtop, "nsets", nsets);
      read(paramtop, "HL_ratio", HL_ratio);
      if(HL_ratio<=0)
      {
         QDPIO::cerr << InlineLoopEnv::name << ": set HL_ratio<=0 is pointless, should be larger than 0."
         << std::endl;
         QDP_abort(1);
      }
      if(nsets%HL_ratio!=0)
      {
         QDPIO::cerr << InlineLoopEnv::name << ": HL_ratio (" <<HL_ratio<<") should be a factor of nsets ("<<nsets<<")."
         << std::endl;
         QDP_abort(1);
      }
      read(paramtop, "Param", param);
      inv_param_h=readXMLGroup(paramtop, "InvParamH", "invType");
      inv_param_l=readXMLGroup(paramtop, "InvParamL", "invType");
      read(paramtop, "cfg_serial", cfg_serial);

      // Read in the output propagator/source configuration info
      read(paramtop, "NamedObject", named_obj);

      read(paramtop, "iog_file", iog_file);
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
  InlineLoopParams::writeXML(XMLWriter& xml_out, const std::string& path) 
  {
    push(xml_out, path);
    
    write(xml_out, "Param", param);
    write(xml_out, "NamedObject", named_obj);

    pop(xml_out);
  }

  int get_ckpoint_info(char *name)
  {
      int info=-1;
      if (Layout::primaryNode())
      {
         FILE *pfile=fopen(name,"r");
         if(pfile!=NULL)
         {
             fclose(pfile);
             general_data_base tmp(name);
             tmp.load_type();
             info=tmp.dim[1].indices[0];
         }
      }
      QDPInternal::broadcast(info);
      return info;
  }


  // Function call
  void 
  InlineLoop::operator()(unsigned long update_no,
			       XMLWriter& xml_out) 
  {
    // If xml file not empty, then use alternate
    if (params.xml_file != "")
    {
      std::string xml_file = makeXMLFileName(params.xml_file, update_no);

      push(xml_out, "propagator");
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
  void 
  InlineLoop::func(unsigned long update_no,
			 XMLWriter& xml_out) 
  {
    START_CODE();

    QDPIO::cout << InlineLoopEnv::name << ": loop calculation" << std::endl;

    StopWatch snoop;
    snoop.reset();
    snoop.start();

    // Test and grab a reference to the gauge field
    XMLBufferWriter gauge_xml;
    try
    {
      TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);
      TheNamedObjMap::Instance().get(params.named_obj.gauge_id).getRecordXML(gauge_xml);
    }
    catch( std::bad_cast ) 
    {
      QDPIO::cerr << InlineLoopEnv::name << ": caught dynamic cast error" 
		  << std::endl;
      QDP_abort(1);
    }
    catch (const std::string& e) 
    {
      QDPIO::cerr << InlineLoopEnv::name << ": std::map call failed: " << e 
		  << std::endl;
      QDP_abort(1);
    }
    const multi1d<LatticeColorMatrix>& u = 
      TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);

    push(xml_out, "propagator");
    write(xml_out, "update_no", update_no);

    proginfo(xml_out);    // Print out basic program info

    // Write out the input
    params.writeXML(xml_out, "Input");

    // Write out the config header
    write(xml_out, "Config_info", gauge_xml);

    push(xml_out, "Output_version");
    write(xml_out, "out_version", 1);
    pop(xml_out);

    // Calculate some gauge invariant observables just for info.
    MesPlq(xml_out, "Observables", u);

    // These pesky variables are needed in the quarkprop call - only chiral dudes
    // need this stuff, but it must be there for the bleeping virtual function
    // to live at the base class
    int t0;
    int j_decay;
    SftMom phases(0, true, Nd-1);

    int op_size=1;
    for(int idr=0;idr<3;idr++)
       op_size+=2*params.link_lengths[idr];
    multi1d<LatticeComplex> quark_loop(16*op_size);

        LatticeComplex c;
        LatticeReal rnd1, theta;
        random(rnd1);
        Real twopiN = Chroma::twopi / 4; //Z_4 grid for loop;
        theta = twopiN * floor(4*rnd1);
        c = cmplx(cos(theta),sin(theta));

     srand(time(0));
     multi1d<int> offset(4);
     for(int idr=0;idr<3;idr++) offset[idr]=(params.grid[idr]>0)?params.grid[idr]:(Layout::lattSize()[idr]);
     offset[3]=2;
     
     int nt=Layout::lattSize()[3];
     
     general_data_base result;
     result.add_dimension(dim_conf,1,&params.cfg_serial);
     result.add_dimension(dim_channel,1);
     result.add_dimension(dim_displacement,op_size);
     result.dim[2].indices[0]=0;
     {
        int count=1;
        for(int idr=0;idr<3;idr++)
        {
           for(int is=0;is<params.link_lengths[idr];is++)
               result.dim[2].indices[count++]=idr*100+is+1;
           for(int is=0;is<params.link_lengths[idr];is++)
               result.dim[2].indices[count++]=(idr+4)*100+is+1;
         }
     }
     result.add_dimension(dim_operator,16);
     result.add_dimension(dim_t,nt);
     result.add_dimension(dim_complex,2);
     if(Layout::primaryNode())
        result.initialize();
        
      QDPIO::cout << InlineLoopEnv::name << ": initialized " << std::endl;

     bool mg_setup_flag=false;
     
     std::vector<int> pos_lst;
     for(int iset=0;iset<params.nsets;iset++)
     {
        sprintf(result.name,"%s.%04d",params.iog_file.c_str(),iset);
        int pos_l=get_ckpoint_info(result.name);
        if(pos_l>=0 && params.use_ckpoint==true)
        {
             pos_lst.push_back(pos_l);
             QDPIO::cout << InlineLoopEnv::name << ": the" << iset <<"-th set with pos_l="
                << pos_l <<", skipped" << std::endl;
             continue;
        }
     
     // create the sources; 
        multi1d<int> pos(3);int flag=0;
         while(flag==0)
         {
           pos_l=0;
           for(int i=0,order=1;i<3;i++,order*=1000)
           {
                pos[i]=rand()%params.grid[i];
                pos_l+=order*pos[i];
           }
           flag=1;
           for(int i=0;i<iset;i++)
             if(pos_lst[i]==pos_l){flag=0;break;}
         }
         pos_lst.push_back(pos_l);
         result.clear();
         result.dim[1].indices[0]=pos_l;

          QDPIO::cout << InlineLoopEnv::name << ": the " << iset <<"-th set with pos_l="
                << pos_l <<", start the calculation" << std::endl;

         for(int is=0;is<16*op_size;is++) quark_loop[is]=zero;

         for(int it0=0;it0<2;it0++)
         for(int ieo=0;ieo<2;ieo++)
         {
             LatticeComplex c1=where(
                ((((Layout::latticeCoordinate(0))-pos[0])%offset[0]==0)&&
                 (((Layout::latticeCoordinate(1))-pos[1])%offset[1]==0)&&
                 (((Layout::latticeCoordinate(2))-pos[2])%offset[2]==0)&&
                 (((Layout::latticeCoordinate(3))-it0)%offset[3]==0)),
                c, LatticeComplex(zero));
              LatticeComplex c2=where(
                 (((Layout::latticeCoordinate(0))-pos[0])/offset[0]
                 +((Layout::latticeCoordinate(1))-pos[1])/offset[1]
                 +((Layout::latticeCoordinate(2))-pos[2])/offset[2]
                 +((Layout::latticeCoordinate(3))-it0)/offset[3])%2==ieo,
                 c1,LatticeComplex(zero));
            
                LatticePropagator quark_prop_source;
     
             for(int color_source = 0; color_source < Nc; ++color_source)
             for(int spin_source = 0; spin_source < Ns; ++spin_source)
             {
                 LatticeColorVector colorvec = zero;
                 LatticeFermion chi = zero;
                 pokeSpin(chi,pokeColor(colorvec,c2,color_source),spin_source);
                 FermToProp(chi, quark_prop_source, color_source, spin_source);
             } 

           // original propagator generation starts here
           
               // Sanity check - write out the norm2 of the source in the Nd-1 direction
               // Use this for any possible verification
               {
                 // Initialize the slow Fourier transform phases
           
                 multi1d<DComplex> source_corr = sumMulti(trace(quark_prop_source*adj(quark_prop_source)),phases.getSet());
                 if(Layout::primaryNode())
                 {
                     for(int it=0;it<Layout::lattSize()[3];it++)
                        printf("srcNROM:%4d%13.5f\n",it,source_corr[it].elem().elem().elem().real());
                     fflush(stdout);
                 }
               }
           
               LatticePropagator quark_propagator,quark_propagator_H;
               int ncg_had = 0;
           
               //
               // Initialize fermion action
               //
               std::istringstream  xml_s(params.param.fermact.xml);
               XMLReader  fermacttop(xml_s);
               QDPIO::cout << "FermAct = " << params.param.fermact.id << std::endl;
           
           
               //
               // Try the factories
               //
               bool success = false;
           
               if (! success)
               {
                 try
                 {
           	StopWatch swatch;
           	swatch.reset();
           	QDPIO::cout << "Try the various factories" << std::endl;
           
           	// Typedefs to save typing
           	typedef LatticeFermion               T;
           	typedef multi1d<LatticeColorMatrix>  P;
           	typedef multi1d<LatticeColorMatrix>  Q;
           
           	// Generic Wilson-Type stuff
           	Handle< FermionAction<T,P,Q> >
           	  S_f(TheFermionActionFactory::Instance().createObject(params.param.fermact.id,
           							       fermacttop,
           							       params.param.fermact.path));
           
           	Handle< FermState<T,P,Q> > state(S_f->createState(u));
           
           	QDPIO::cout << "Suitable factory found: compute the quark prop" << std::endl;
           	swatch.start();
           	
           	QDPIO::cout << "iters: " << iset << "(" <<  (iset+1)%params.HL_ratio <<") " << it0*2+ieo << std::endl; 
           	
           	GroupXML_t invParamR=(mg_setup_flag==false)?params.param.invParam:params.inv_param_l;
//           	if(params.HL_ratio==1)invParamR=params.inv_param_h;

           	QDPIO::cout << "Calling quarkProp" << std::endl;
           	S_f->quarkProp(quark_propagator, 
           		       xml_out, 
           		       quark_prop_source,
           		       t0, j_decay,
           		       state, 
           		       invParamR, 
           		       params.param.quarkSpinType,
           		       params.param.obsvP,
           		       ncg_had);
           	swatch.stop();
           	QDPIO::cout << "Propagator computed: time= " 
           		    << swatch.getTimeInSeconds() 
           		    << " secs" << std::endl;
           		    
                 if(mg_setup_flag==false) mg_setup_flag=true;
                 
                 if(params.HL_ratio>1&&(iset+1)%params.HL_ratio==0)//make the propagator needed by the AMA correction
                 {
                      GroupXML_t invParamR2=params.inv_param_h;
                      
               	      QDPIO::cout << "Calling quarkProp" << std::endl;
           	      S_f->quarkProp(quark_propagator_H, 
           	      	       xml_out, 
           	      	       quark_prop_source,
           	      	       t0, j_decay,
           	      	       state, 
           	      	       invParamR2, 
           	      	       params.param.quarkSpinType,
           	      	       params.param.obsvP,
           	      	       ncg_had);
           	      swatch.stop();
           	      QDPIO::cout << "Propagator computed: time= " 
           	      	    << swatch.getTimeInSeconds() 
           	      	    << " secs" << std::endl;                       
                 }
                 
           	success = true;
                 }
                 catch (const std::string& e) 
                 {
           	QDPIO::cout << InlineLoopEnv::name << ": caught exception around quarkprop: " << e << std::endl;
                 }
               }
             
               if (! success)
               {
                 QDPIO::cerr << "Error: no fermact found" << std::endl;
                 QDP_abort(1);
               }
         
               // handle the loop contraction.   
               LatticeComplex c2_dagger=adj(c2);
               LatticePropagator F,F_mu;
               LatticeComplex tmp,tmp2;
               int icount=1;
               for(int mu=0;mu<3;mu++)
               {
                   F=quark_propagator;
                   if(mu==0)
                   for(int ig=0;ig<16;ig++)
                      quark_loop[ig]+=LatticeComplex(trace(F*c2_dagger*Gamma(ig)));
                   for(int is=0;is<params.link_lengths[mu];is++)
                   {
                       F_mu = shift( adj( u[mu] ) * F, BACKWARD, mu );
                       for(int ig=0;ig<16;ig++)
                           quark_loop[ig+16*icount]+=LatticeComplex(trace(F_mu*c2_dagger*Gamma(ig)));  
                       icount++;
                       F=F_mu;
                   }
                   F=quark_propagator;
                   for(int is=0;is<params.link_lengths[mu];is++)
                   {
                       F_mu = u[mu] * shift( F, FORWARD, mu );
                       for(int ig=0;ig<16;ig++)
                           quark_loop[ig+16*icount]+=LatticeComplex(trace(F_mu*c2_dagger*Gamma(ig)));  
                       icount++;
                       F=F_mu;
                   }
               }
               if(params.HL_ratio>1&&(iset+1)%params.HL_ratio==0)//make the AMA correction
               {
                    icount=1;
                    quark_propagator_H-=quark_propagator;
                    for(int mu=0;mu<3;mu++)
                    {
                        F=quark_propagator_H;
                        if(mu==0)
                        for(int ig=0;ig<16;ig++)
                           quark_loop[ig]+=LatticeComplex(trace(F*c2_dagger*Gamma(ig)));
                        for(int is=0;is<params.link_lengths[mu];is++)
                        {
                            F_mu = shift( adj( u[mu] ) * F, BACKWARD, mu );
                            for(int ig=0;ig<16;ig++)
                                quark_loop[ig+16*icount]+=LatticeComplex(trace(F_mu*c2_dagger*Gamma(ig)))*params.HL_ratio;  
                            icount++;
                            F=F_mu;
                        }
                        F=quark_propagator_H;
                        for(int is=0;is<params.link_lengths[mu];is++)
                        {
                            F_mu = u[mu] * shift( F, FORWARD, mu );
                            for(int ig=0;ig<16;ig++)
                                quark_loop[ig+16*icount]+=LatticeComplex(trace(F_mu*c2_dagger*Gamma(ig)))*params.HL_ratio;  
                            icount++;
                            F=F_mu;
                        }
                    }
               }
               
               QDPIO::cout << InlineLoopEnv::name << "quark loop done" << std::endl;
              
          }
          
            for(int is=0;is<op_size;is++) 
            {
               for(int ig=0;ig<16;ig++)
               {
                  multi1d<DComplex> hsum =sumMulti(quark_loop[is*16+ig],phases.getSet());
                  if(Layout::primaryNode())
                  for(int t=0;t<nt;t++)//addtional minus sign for the disconnect quark loop
                  {
                     result.data[(is*16+ig)*2*nt+2*t]-=hsum[t].elem().elem().elem().real();
                     result.data[(is*16+ig)*2*nt+2*t+1]-=hsum[t].elem().elem().elem().imag();
                  }
               }
            }
               for(int is=0;is<op_size;is++)
               if(Layout::primaryNode())
               {
                  for(int t=0;t<nt;t++)
                  printf("scaLOOP:%4d%4d%15.5f%13.5f%15.5f%13.5f%15.5f%13.5f%15.5f%13.5f\n",result.dim[2].indices[is],t,
                        result.data[(is*16+0)*2*nt+2*t],result.data[(is*16+0)*2*nt+2*t+1],
                        result.data[(is*16+1)*2*nt+2*t],result.data[(is*16+1)*2*nt+2*t+1],
                        result.data[(is*16+2)*2*nt+2*t],result.data[(is*16+2)*2*nt+2*t+1],
                        result.data[(is*16+4)*2*nt+2*t],result.data[(is*16+4)*2*nt+2*t+1]);
               fflush(stdout);
               }
           if(Layout::primaryNode()) result.save();
      }//end of the loops
      
      //merge the data
      if(Layout::primaryNode())
      {
          for(int i=0;i<params.nsets;i++) pos_lst[i]=i;
          std::string name=params.iog_file+".%04d";
          result.gather_file(name.c_str(),pos_lst);
          result.aver(result,dim_channel,0);
          result.dim[1].indices[0]=1000*(1000*params.grid[0]+params.grid[1])+params.grid[2];
          int fac=params.grid[0]*params.grid[1]*params.grid[2];
          for(int i=0;i<result.size;i++)
           result.data[i]*=fac;
          sprintf(result.name,params.iog_file.c_str());
          result.save();
      }
     

    snoop.stop();
    QDPIO::cout << InlineLoopEnv::name << ": total time = "
		<< snoop.getTimeInSeconds() 
		<< " secs" << std::endl;

    QDPIO::cout << InlineLoopEnv::name << ": ran successfully" << std::endl;

    END_CODE();
  } 

}

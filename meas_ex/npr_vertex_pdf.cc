//  $Id: npr_vertex_pdf.cc,v 1.0 2017-8-14 R. Zhang & Y. Yang, Exp $
/*! \file
 *  \brief NPR vertex calculations
 */

#include "util/ft/sftmom.h"
#include "meas_ex/npr_vertex_pdf.h"
#include <string>
using std::string;

namespace Chroma 
{


  QDP::StandardOutputStream& operator<<(QDP::StandardOutputStream& s, const multi1d<int>& d);
/*  {
    if (d.size() > 0)
    {
      s << d[0];
      for(int i=1; i < d.size(); ++i)
	s << " " << d[i];
    }

    return s;
  }*/


  void BkwdFrwd(const LatticePropagator&  B,
		const LatticePropagator&  F,
		general_data_base& io_output,
		int 		LinkNumber,
		int& GBB_NLinkPatterns,
		const multi1d< int > & LinkDirs)
  {
    StopWatch TotalTime;
    TotalTime.reset();
    TotalTime.start();

    for( int i = 0; i < Ns * Ns; i ++ )
    {
      XMLBufferWriter record_xml;
      push(record_xml, "Vertex");

      QDPIO::cout << __func__ << ": LinkDirs = " << LinkDirs 
		  << "  gamma = " << i << std::endl;

      write(record_xml, "linkDirs", LinkDirs);   // link pattern
      write(record_xml, "gamma", i);

      // assumes any Gamma5 matrices have already been absorbed
      int G5 = Ns*Ns-1;
      
      // Compute the single site propagator and write it
      DPropagator prop;
      {
	// assumes any Gamma5 matrices have already been absorbed into B
	LatticePropagator tmp = B * Gamma(i) * F;
	// The site's worth of data of interest
	prop = sum(tmp)/Double(Layout::vol()); // and normalize by the volume
	//QDPIO::cout<<"  The 1/12*trace is: "<<trace(Gamma(i)*prop)/12.0;
	//QDPIO::cout<<std::endl ;
      }
    
    int data_index=GBB_NLinkPatterns*288;
    if(io_output.data!=NULL) if(Layout::primaryNode())
    {
      SpinMatrixD U = KYToDRMat();
      
      prop = adj(U) * prop * U;
      for(int j=0; j!=Ns; j++)
      for(int l=0; l!=Nc; l++)
      for(int y=0; y!=Ns; y++)
      for(int k=0; k!=Nc; k++)
      {
		io_output.data[data_index]=prop.elem().elem(y,j).elem(k,l).real();
		data_index++;
		io_output.data[data_index]=prop.elem().elem(y,j).elem(k,l).imag();
		data_index++;
      }
      pop(record_xml);
    }

      // counts number of link patterns
      GBB_NLinkPatterns++;
    }

    TotalTime.stop();
    QDPIO::cout << __func__ << ": total time = " << TotalTime.getTimeInSeconds() << " seconds" << std::endl;

    return;
  }

//###################################################################################//
// accumulate link operators                                                         //
//###################################################################################//

  void AddLinks(const LatticePropagator&  B,
		const LatticePropagator&  F,
		const multi1d< LatticeColorMatrix > & U,
		multi1d< int >&    LinkDirs,
		const int          MaxNLinks,
		const int	   LinksDirection,
		BBLinkPattern      LinkPattern,
		general_data_base& io_output,
		int&               GBB_NLinkPatterns)
  {
    StopWatch Timer;
    Timer.reset();
    Timer.start();

    const int NLinks = LinkDirs.size();

    if( NLinks == MaxNLinks )
    {
      return;
    }

    LatticePropagator F_mu;
    multi1d< int > NextLinkDirs( NLinks + 1 );

    for(int Link = 0; Link < NLinks; Link ++)
    {
      NextLinkDirs[ Link ] = LinkDirs[ Link ];
    }

		bool DoThisPattern = true;
		bool DoFurtherPatterns = true;
	
		NextLinkDirs[ NLinks ] = LinksDirection;
	
		LinkPattern( DoThisPattern, DoFurtherPatterns, NextLinkDirs );
	
		if( DoFurtherPatterns == true )
		{
		  // accumulate product of link fields
		  if( LinksDirection >= 0 && LinksDirection < Nd )
          		  F_mu = shift( adj( U[ LinksDirection ] ) * F, BACKWARD, LinksDirection );
                  else if( LinksDirection >= Nd && LinksDirection <= Nd*2-1 )
                          F_mu = U[ LinksDirection-Nd ] * shift( F, FORWARD, LinksDirection-Nd );
                  else throw "invalid linkdirection input";
		}
	
		if( DoThisPattern == true )
		{
		  BkwdFrwd(B, F_mu, io_output, NLinks+1, GBB_NLinkPatterns, NextLinkDirs);
		}
	
		if( DoFurtherPatterns == true )
		{
		  // add another link
		  AddLinks(B, F_mu, U,
			   NextLinkDirs, MaxNLinks, LinksDirection, LinkPattern, 
			   io_output, GBB_NLinkPatterns);
		}
	      
    Timer.stop();
    QDPIO::cout << __func__ << ": total time = " << Timer.getTimeInSeconds() << " seconds" << std::endl;

    return;
  }


  //! NPR vertices
  void NprVertexPDF(const LatticePropagator &             F0,
		 const multi1d< LatticeColorMatrix > & U,
		 const unsigned short int              MaxNLinks,
		 const unsigned short int              MaxNLinks_prep,
		 const BBLinkPattern                   LinkPattern,
		 const char* name,
		 int  cfg_serial)
  {
    StopWatch TotalTime;
    TotalTime.reset();
    TotalTime.start();
	
    StopWatch Timer;



    general_data_base output(name);

	output.add_dimension(dim_conf, 1, &cfg_serial);
	if(MaxNLinks_prep>0)
	{
       	  output.add_dimension(dim_channel, 1+MaxNLinks_prep*4);
       	  output.dim[1].indices[0]=0;
       	  for(int mu=0;mu<4;mu++)
       	  for(int is=0;is<MaxNLinks_prep;is++)
       	     output.dim[1].indices[1+is+MaxNLinks_prep*mu]=mu*100+is+1;
        }
	output.add_dimension(dim_direction, 4);
	output.add_dimension(dim_displacement, 1+MaxNLinks);
	output.add_dimension(dim_operator, 16);
	output.add_dimension(dim_temporary, 144);
	output.add_dimension(dim_complex, 2);
	if(Layout::primaryNode()) output.initialize();

    int GBB_NLinkPatterns=0;

    //#################################################################################//
    // open building blocks data files                                                 //
    //#################################################################################//
    LatticePropagator F=F0;

    Timer.reset();
    Timer.start();

    //#################################################################################//
    // calculate building blocks                                                       //
    //#################################################################################//
    for(int ilnk=-1;ilnk<MaxNLinks_prep*4;ilnk++)
    {
          if(ilnk>=0)
          {
              int mu=ilnk/MaxNLinks_prep;
              if(ilnk%MaxNLinks_prep==0) F=F0;
              LatticePropagator F_tmp;
            for(int iter=0;iter<3;iter++)
            {
              F_tmp = shift( adj( U[ mu ] ) * F, BACKWARD, mu );
              F=F_tmp;
            }
          }
          for(int idir=0;idir<4;idir++)
          {
            QDPIO::cout << __func__ << ": start BkwdFrwd" << std::endl;
        
            const int NLinks = 0;
            multi1d< int > LinkDirs( 0 );
        
            LatticePropagator B = Gamma(15)*adj(F)*Gamma(15);
            BkwdFrwd(B, F, output, 0, GBB_NLinkPatterns, LinkDirs);
            
        
            Timer.stop();
            QDPIO::cout << __func__ << ": total time for 0 links (single BkwdFrwdTr call) = "
        		<< Timer.getTimeInSeconds() 
        		<< " seconds" << std::endl;
        
            Timer.reset();
            Timer.start();
        
            QDPIO::cout << __func__ << ": start AddLinks" << std::endl;
        
            AddLinks(B, F, U, 
        	     LinkDirs, MaxNLinks, idir, LinkPattern, 
        	      output, GBB_NLinkPatterns);
        
            Timer.stop();
            QDPIO::cout << __func__ << ": total time for remaining links (outermost AddLinks call) = "
        		<< Timer.getTimeInSeconds() 
        		<< " seconds" << std::endl;
        
            QDPIO::cout << __func__ << "saving io_general output" << std::endl;
          }
    }
    if(output.data!=NULL) if(Layout::primaryNode()) output.save();
    TotalTime.stop();
    QDPIO::cout << __func__ << ": total time = "
		<< TotalTime.getTimeInSeconds() 
		<< " seconds" << std::endl;


    return;
  }
}  // end namespace Chroma

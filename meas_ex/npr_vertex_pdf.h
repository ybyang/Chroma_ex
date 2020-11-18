// -*- C++ -*-
//  $Id: npr_vertex_pdf.h,v 1.0 2017-8-14 R. Zhang & Y. Yang, Exp $
/*! \file
 *  \brief NPR vertex calculations
 */

#ifndef __npr_vertex_pdf_h__
#define __npr_vertex_pdf_h__

#include "chromabase.h"
#include "io_ex/io_general_class.h"

namespace Chroma 
{
  //! Used to Set Requested Link Patterns
  /*! \ingroup hadron */
  typedef void (*BBLinkPattern)(bool &                          DoThisPattern,
				bool &                          DoFurtherPatterns,
				multi1d< int > & LinkPattern);


  //! NPR vertices
  /*! \ingroup hadron */
  void NprVertexPDF(const LatticePropagator &             F,
		 const multi1d< LatticeColorMatrix > & U,
		 const unsigned short int              MaxNLinks,
		 const unsigned short int              MaxNLinks_prep,
		 const BBLinkPattern                   LinkPattern,
		 const char* name,
		 const int   cfg_serial);

SpinMatrixD KYToDRMat();
void writeKYUQprop2(LatticePropagator& q, const std::string& file);
void readKYUQprop2(LatticePropagator& q, const std::string& file);
}  // end namespace Chroma

#endif

//###################################################################################//
//###################################################################################//

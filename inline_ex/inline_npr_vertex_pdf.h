// -*- C++ -*-
/*! \file
 * \brief Inline construction of NPR vertices for quasi PDF
 *
 * NPR vertices on  props
 */

#ifndef __inline_npr_vertex_pdf_h__
#define __inline_npr_vertex_pdf_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineNprVertexPDFEnv 
  {
    extern const std::string name;
    bool registerAll();
  }

  //! Parameter structure
  /*! \ingroup inlinehadron */
  struct InlineNprVertexPDFParams 
  {
    InlineNprVertexPDFParams();
    InlineNprVertexPDFParams(XMLReader& xml_in, const std::string& path);
    void write(XMLWriter& xml_out, const std::string& path);

    unsigned long frequency;

    //! Parameters
    struct Param_t
    {
      int          links_max;	/*!< maximum number of links */
      int          links_prep_max;	/*!< maximum number of prep links */
      bool	   save_prop;
      std::string  file_name;          /*!< bb output file name pattern */
      int          cfg_serial;
      GroupXML_t   cfs;                /*!< Fermion state */
    } param;

    //! Propagators
    struct NamedObject_t
    {
      std::string       gauge_id;        /*!< Input Gauge id */
      std::string       prop_id;         /*!< Input forward prop */
    } named_obj;

    std::string xml_file;  // Alternate XML file pattern
  };


  //! Inline measurement of NPR vertices
  /*! \ingroup inlinehadron */
  class InlineNprVertexPDF : public AbsInlineMeasurement 
  {
  public:
    ~InlineNprVertexPDF() {}
    InlineNprVertexPDF(const InlineNprVertexPDFParams& p) : params(p) {}
    InlineNprVertexPDF(const InlineNprVertexPDF& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const unsigned long update_no,
		    XMLWriter& xml_out); 

  protected:
    //! Do the measurement
    void func(const unsigned long update_no,
	      XMLWriter& xml_out); 

  private:
    InlineNprVertexPDFParams params;
  };

};

#endif

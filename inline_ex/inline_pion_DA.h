// -*- C++ -*-
// $Id: inline_pion_decay.h, v1.1 2015-06-22 19:43:00 sdcohen Exp $
/*! \file
 * \brief Inline construction of pion-decay two-point correlators
 *
 * Calculates the two-point correlator of the pseudoscalar meson corresponding
 * to leptonic decay by temporal-axial current.
 */

#ifndef __inline_PionDA_h__
#define __inline_PionDA_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "util/ft/sftmom.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlinePionDAEnv 
  {
    extern const std::string name;
    bool registerAll();
  }

  //! Parameter structure
  /*! \ingroup inlinehadron */
  struct InlinePionDAParams 
  {
    // Default constructor
    InlinePionDAParams(void);
    // Construct from XML
    InlinePionDAParams(XMLReader& xml_in, const std::string& path);
    // Write out the configuration of the parameters
    void write(XMLWriter& xml_out, const std::string& path);

    // How often should I be called during a gauge evolution?
    unsigned long frequency;

    // Holds the non-lattice parameters
    struct Param_t
    {
      multi1d<int>     mom ;   /*<! momentum at sink*/
      int  links_dir;   /*<! the direction of the wilson link*/
      int  links_max;   /*<! the maxinum length of the wilson link*/
      std::string  file_name;          /*!< bb output file name pattern */
      int  cfg_serial;          /*!< the configuration serial number */
      
      multi1d<int> src;  /*!< source position */
    } param;

    // Holds the names of the lattice objects
    struct NamedObject_t
    {
      std::string  gauge_id;           /*!< Input gauge field */
      multi1d<std::string> quark_props;  /*!< List of propagator pairs to contract */
      multi1d<std::string> antiquark_props;
    } named_obj;

    std::string xml_file;  // Alternate XML file pattern
  }; // end of struct InlinePionDAParams


  //! Inline measurement of hadron spectrum
  /*! \ingroup inlinehadron */
  class InlinePionDA : public AbsInlineMeasurement 
  {
  public:
    // Default destructor
    ~InlinePionDA() {}
    // Constructor from param struct: copies param struct
    InlinePionDA(const InlinePionDAParams& p) : params(p) {}
    // Copy constructor
    InlinePionDA(const InlinePionDA& p) : params(p.params) {}

    // Getter for measurement frequency
    unsigned long getFrequency() const {return params.frequency;}

    //! Sets up the XML and invokes func, which does the acual work
    void operator()( const unsigned long update_no, XMLWriter& xml_out); 

  protected:
    //! Does the actual work
    void func(const unsigned long update_no, XMLWriter& xml_out); 

  private:
    //! The parameter structure; holds names of props, gauge field, XML, etc.
    InlinePionDAParams params;
  }; // end of class InlinePionDA

}; // end of namespace Chroma

#endif


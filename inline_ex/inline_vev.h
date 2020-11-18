// -*- C++ -*-
// $Id: inline_pion_decay.h, v1.1 2015-06-22 19:43:00 sdcohen Exp $
/*! \file
 * \brief Inline construction of pion-decay two-point correlators
 *
 * Calculates the two-point correlator of the pseudoscalar meson corresponding
 * to leptonic decay by temporal-axial current.
 */

#ifndef __inline_vev_h__
#define __inline_vev_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "util/ft/sftmom.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineVEVEnv 
  {
    extern const std::string name;
    bool registerAll();
  }

  //! Parameter structure
  /*! \ingroup inlinehadron */
  struct InlineVEVParams 
  {
    // Default constructor
    InlineVEVParams(void);
    // Construct from XML
    InlineVEVParams(XMLReader& xml_in, const std::string& path);
    // Write out the configuration of the parameters
    void write(XMLWriter& xml_out, const std::string& path);

    // How often should I be called during a gauge evolution?
    unsigned long frequency;

    // Holds the non-lattice parameters
    struct Param_t
    {
      int  links_dir;   /*<! the direction of the wilson link*/
      int  links_max;   /*<! the maxinum length of the wilson link*/
      int  t_src;   /*<! the source t*/
      int  t_max;   /*<! the maxinum length of t*/
      std::string  file_name;          /*!< bb output file name pattern */
      int  cfg_serial;          /*!< the configuration serial number */
      
    } param;

    // Holds the names of the lattice objects
    struct NamedObject_t
    {
      std::string quark_prop;  /*!< List of propagator pairs to contract */
      std::string gauge_id;     /*!< input gauge id */
    } named_obj;

    std::string xml_file;  // Alternate XML file pattern
  }; // end of struct InlineVEVParams


  //! Inline measurement of hadron spectrum
  /*! \ingroup inlinehadron */
  class InlineVEV : public AbsInlineMeasurement 
  {
  public:
    // Default destructor
    ~InlineVEV() {}
    // Constructor from param struct: copies param struct
    InlineVEV(const InlineVEVParams& p) : params(p) {}
    // Copy constructor
    InlineVEV(const InlineVEV& p) : params(p.params) {}

    // Getter for measurement frequency
    unsigned long getFrequency() const {return params.frequency;}

    //! Sets up the XML and invokes func, which does the acual work
    void operator()( const unsigned long update_no, XMLWriter& xml_out); 

  protected:
    //! Does the actual work
    void func(const unsigned long update_no, XMLWriter& xml_out); 

  private:
    //! The parameter structure; holds names of props, gauge field, XML, etc.
    InlineVEVParams params;
  }; // end of class InlineVEV

}; // end of namespace Chroma

#endif


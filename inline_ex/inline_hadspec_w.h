// -*- C++ -*-
/*! \file
 * \brief Inline hadron spectrum calculations
 *
 * Hadron spectrum calculations
 */

#ifndef __inline_hadspec_h_v2__
#define __inline_hadspec_h_v2__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineHadSpec2Env 
  {
    extern const std::string name;
    bool registerAll();
  }

  //! Parameter structure
  /*! \ingroup inlinehadron */
  struct InlineHadSpec2Params 
  {
    InlineHadSpec2Params();
    InlineHadSpec2Params(XMLReader& xml_in, const std::string& path);
    void write(XMLWriter& xml_out, const std::string& path);

    unsigned long frequency;

    struct Param_t
    {
//      bool time_rev;           // Use time reversal in baryon spectroscopy
      bool translate;           // Shift the data to the original point if translate=true
      bool save_vec;
      multi1d<int> src;
      int boundary_fac;
      
      int prj_type; // the type of the projection. =0 by default. non-zero value case will be used for projecting proton in both forward and backward cases.

      int mom2_max;            // (mom)^2 <= mom2_max. mom2_max=7 in szin.
      multi1d<int> mom_list;
      multi1d<int> hadron_list;
      bool avg_equiv_mom;      // average over equivalent momenta
      int cfg_serial;
    } param;

    struct NamedObject_t
    {
      std::string  gauge_id;           /*!< Input gauge field */

      struct Props_t
      {
	std::string  first_id;
	std::string  second_id;
      };

      multi1d<Props_t> sink_pairs;
    } named_obj;

    std::string output;
    std::string xml_file;  // Alternate XML file pattern
  };


  //! Inline measurement of hadron spectrum
  /*! \ingroup inlinehadron */
  class InlineHadSpec2 : public AbsInlineMeasurement 
  {
  public:
    ~InlineHadSpec2() {}
    InlineHadSpec2(const InlineHadSpec2Params& p) : params(p) {}
    InlineHadSpec2(const InlineHadSpec2& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const unsigned long update_no,
		    XMLWriter& xml_out); 

  protected:
    //! Do the measurement
    void func(const unsigned long update_no,
	      XMLWriter& xml_out); 

  private:
    InlineHadSpec2Params params;
  };

};

#endif

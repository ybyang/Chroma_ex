// -*- C++ -*-
// $Id: inline_building_blocks_w.h, v3.6 2014-04-05 17:56 sdcohen Exp $
/*! \file
 * \brief Inline construction of BuildingX
 *
 * Building Blocks on forward and sequential props
 */

#ifndef __inline_BuildingX_iog_h__
#define __inline_BuildingX_iog_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"
#include "util/ft/sftmom.h"

namespace Chroma 
{ 
  // Environment in which the measurement lives (holds params and such)
  namespace InlineBuildingXIOGEnv 
  {
    extern const std::string name;
    bool registerAll();
  }

  //! Parameter structure
  /*! \ingroup inlinehadron */
  struct InlineBuildingXIOGParams 
  {
    // Default constructor
    InlineBuildingXIOGParams();
    // Construct from XML
    InlineBuildingXIOGParams(XMLReader& xml_in, const std::string& path);
    // Write out the configuration of the parameters
    void write(XMLWriter& xml_out, const std::string& path);

    unsigned long frequency;

    // Holds the non-lattice parameters
    struct Param_t
    {
      bool       use_sink_offset; /*!< should insertion origin be sink_mom */
      int        mom2_max;        /*!< (mom)^2 <= mom2_max */
      multi1d<int> mom2_list;     /*!< (mom)^2 in mom2_list */
      bool       canonical;       /*!< true if mom in BB filenames is canonicalized */
      bool       time_reverse;    /*!< time reverse the building blocks */
      bool       translate;       /*!< shift BB to start at tsrc 0 */
      int        icfg;            /*!< The configuration serial number*/
      GroupXML_t cfs;             /*!< fermion state */
      int 	mom_dir;

      multi1d<int> src;  /*!< list of source positions */
      // SftMomSrcPos_t is defined in lib/util/ft/sftmom.h
    } param;

    // Holds the description of a three-point correlator (except forward prop)
    struct BkwdProp_t
    {
      std::string BkwdPropId;        /*!< backward propagator */
      std::string BkwdPropG5Format;  /*!< Where are gamma5's multiplied in? */
      int         GammaInsertion;    /*!< second gamma insertion */
      std::string Flavor;            /*!< flavor in {U,D,S,C,B,T} */
      std::string BBFileNamePattern; /*!< BB output file name pattern */
    };

    //! BB output
    struct BB_t
    {
      std::string         GaugeId;     /*!< input gauge id */
      std::string         FrwdPropId;  /*!< input forward prop */
      multi1d<BkwdProp_t> BkwdProps;   /*!< array of backward prop structs */
    } bb;

    std::string xml_file; /*!< Alternate XML file pattern */
  }; // end of struct InlineBuildingXIOGParams


  //! Inline measurement of 3-pt functions writing building-blocks
  /*! \ingroup inlinehadron */
  class InlineBuildingXIOG : public AbsInlineMeasurement 
  {
  public:
    // Default destructor
    ~InlineBuildingXIOG() {}
    // Constructor from param struct: copies param struct
    InlineBuildingXIOG(const InlineBuildingXIOGParams& p) : params(p) {}
    // Copy constructor
    InlineBuildingXIOG(const InlineBuildingXIOG& p) : params(p.params) {}

    // Getter for measurement frequency
    unsigned long getFrequency(void) const {return params.frequency;}

    //! Sets up the XML and invokes func, which does the acual work
    void operator()( const unsigned long update_no, XMLWriter& xml_out); 

  protected:
    //! Does the actual work
    void func(const unsigned long update_no, XMLWriter& xml_out); 

  private:
    //! The parameter structure; holds names of props, gauge field, XML, etc.
    InlineBuildingXIOGParams params;
  }; // end of class InlineBuildingXIOG

}; // end of namespace Chroma

#endif

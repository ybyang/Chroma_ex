// -*- C++ -*-
// $Id: inline_building_blocks_w.h, v3.6 2014-04-05 17:56 sdcohen Exp $
/*! \file
 * \brief Inline construction of BuildingBlocks
 *
 * Building Blocks on forward and sequential props
 */

#ifndef __inline_BuildingBlocks_iog_h__
#define __inline_BuildingBlocks_iog_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"
#include "util/ft/sftmom.h"

#include "chroma_config.h"

namespace Chroma 
{ 
  // Environment in which the measurement lives (holds params and such)
  namespace InlineBuildingBlocksIOGEnv 
  {
    extern const std::string name;
    bool registerAll();
  }

  //! Parameter structure
  /*! \ingroup inlinehadron */
  struct InlineBuildingBlocksIOGParams 
  {
    // Default constructor
    InlineBuildingBlocksIOGParams();
    // Construct from XML
    InlineBuildingBlocksIOGParams(XMLReader& xml_in, const std::string& path);
    // Write out the configuration of the parameters
    void write(XMLWriter& xml_out, const std::string& path);

    unsigned long frequency;

    // Holds the non-lattice parameters
    struct Param_t
    {
      bool       use_sink_offset; /*!< should insertion origin be sink_mom */
      int        mom2_max;        /*!< (mom)^2 <= mom2_max */
      multi1d<int> mom2_list;     /*!< (mom)^2 in mom2_list */
      int        links_max;       /*!< maximum number of links */
      int        links_dir;       /*!< direction of links */
      bool       canonical;       /*!< true if mom in BB filenames is canonicalized */
      bool       time_reverse;    /*!< time reverse the building blocks */
      bool       translate;       /*!< shift BB to start at tsrc 0 */
      int        icfg;            /*!< The configuration serial number*/
      GroupXML_t cfs;             /*!< fermion state */
      bool       conserved;       /*replace the gamma_t case with conserved one*/

      bool       no_gaugelink;       /*replace gauge link with one*/
      bool       use_gpu;         /* computed by gpu */       
      bool       use_cpu;         /* computed by cpu */ 
      int        gamma_id;         /* gamma matrix inserted */
  

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
      std::string BBFileNamePattern_gpu; /*!< BB output file name pattern for gpu*/
    };

    //! BB output
    struct BB_t
    {
      std::string         GaugeId;     /*!< input gauge id */
      std::string         FrwdPropId;  /*!< input forward prop */
      multi1d<BkwdProp_t> BkwdProps;   /*!< array of backward prop structs */
    } bb;

    std::string xml_file; /*!< Alternate XML file pattern */
  }; // end of struct InlineBuildingBlocksIOGParams


  //! Inline measurement of 3-pt functions writing building-blocks
  /*! \ingroup inlinehadron */
  class InlineBuildingBlocksIOG : public AbsInlineMeasurement 
  {
  public:
    // Default destructor
    ~InlineBuildingBlocksIOG() {}
    // Constructor from param struct: copies param struct
    InlineBuildingBlocksIOG(const InlineBuildingBlocksIOGParams& p) : params(p) {}
    // Copy constructor
    InlineBuildingBlocksIOG(const InlineBuildingBlocksIOG& p) : params(p.params) {}

    // Getter for measurement frequency
    unsigned long getFrequency(void) const {return params.frequency;}

    //! Sets up the XML and invokes func, which does the acual work
    void operator()( const unsigned long update_no, XMLWriter& xml_out); 

  protected:
    //! Does the actual work
    void func(const unsigned long update_no, XMLWriter& xml_out); 

  private:
    //! The parameter structure; holds names of props, gauge field, XML, etc.
    InlineBuildingBlocksIOGParams params;
  }; // end of class InlineBuildingBlocksIOG


    class shiftPatternQuda_chroma
  {
      public:
      std::vector<int> pre_shift;
      std::vector<int> post_shift;
/*
      void check()
      {
           if(pre_shift.size()%2!=0||post_shift.size()%2!=0)
                errorQuda("shift vector length incorrect (%5d, %5d)\n",pre_shift.size(),post_shift.size());
           for(int i=0;i<pre_shift.size();i+=2)
           if(pre_shift[i]<0||pre_shift[i+1]<0||pre_shift[i+1]>7)
                errorQuda("pre_shift patterm %d incorrect, (%6d,%6d)\n",i/2,pre_shift[i],pre_shift[i+1]);
           for(int i=0;i<post_shift.size();i+=2)
           if(post_shift[i]<0||post_shift[i+1]<0||post_shift[i+1]>7)
                errorQuda("post_shift patterm %d incorrect, (%6d,%6d)\n",i/2,post_shift[i],post_shift[i+1]);
      }
*/
  };



}; // end of namespace Chroma

#endif

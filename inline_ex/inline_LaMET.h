// -*- C++ -*-
// $Id: inline_building_blocks_w.h, v3.6 2014-04-05 17:56 sdcohen Exp $
/*! \file
 * \brief Inline construction of BuildingBlocks
 *
 * Building Blocks on forward and sequential props
 */

#ifndef __inline_LaMET_h__
#define __inline_LaMET_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"
#include "util/ft/sftmom.h"

namespace Chroma 
{ 
  // Environment in which the measurement lives (holds params and such)
  namespace InlineLaMETEnv 
  {
    extern const std::string name;
    bool registerAll();
  }

  class shift_pattern:: public std::vector<int>
  {
      void load(multi1d<int> &res)
      {
          resize(res.size());
          for(int i=0;i<size();i++)
              data()[i]=res[i];
      }
  
      int check(std::string str="")
      {
           if(size()%2!=0) 
                if(Layout::primaryNode()) { printf("%s: shift vector length incorrect (%5d)\n",str.c_str(),size()); return 1;}
           for(int i=0;i<pre_shift.size();i+=2)
           if(data()[i]<0||data()[i+1]<0||data()[i+1]>11)
                if(Layout::primaryNode()) { printf("%s: pattern %d incorrect, (%6d,%6d)\n",str.c_str(),i/2,str.c_str()[i],str.c_str()[i+1]); return 2;}
      }      
  };

  class shiftPatternChroma
  {
      public:
      shift_pattern pre_shift;
      shift_pattern post_shift;

      int check()
      {
           return pre_shift.check("pre_shfit")+post_shift.check("post_shift");
      }
  };

  //! Parameter structure
  /*! \ingroup inlinehadron */
  struct InlineLaMETParams 
  {
    // Default constructor
    InlineLaMETParams();
    // Construct from XML
    InlineLaMETParams(XMLReader& xml_in, const std::string& path);
    // Write out the configuration of the parameters
    void write(XMLWriter& xml_out, const std::string& path);

    unsigned long frequency;

    // Holds the non-lattice parameters
    struct Param_t
    {
      int        mom2_max;        /*!< (mom)^2 <= mom2_max */
      multi1d<int> mom_list;     /*!< (mom)^2 in mom2_list */
      shift_pattern pre_shift;
      std::vector< shift_pattern > post_shift;
      bool       translate;       /*!< shift BB to start at tsrc 0 */
      int        icfg;            /*!< The configuration serial number*/
      GroupXML_t cfs;             /*!< fermion state */
      bool	 npr;             /* src contraction type, true=12x12, false=1;*/
      bool	 corr_t;	  /* sink/current contraction type, true=nt, false=1;*/
      bool	 qudaLaMET;	  /* Use quda to do the LaMET contraction */

    } param;

    // Holds the description of a three-point correlator (except forward prop)
    struct BkwdProp_t
    {
      std::string BkwdPropId;        /*!< backward propagator */
      std::string BkwdPropG5Format;  /*!< Where are gamma5's multiplied in? */
      std::string BBFileNamePattern; /*!< BB output file name pattern */
      int        GammaInsertion;     /* gamma insertion of source */
    };

    //! BB output
    struct BB_t
    {
      std::string         GaugeId;     /*!< input gauge id */
      std::string         FrwdPropId;  /*!< input forward prop */
      multi1d<BkwdProp_t> BkwdProps;   /*!< array of backward prop structs */
    } bb;

    std::string xml_file; /*!< Alternate XML file pattern */
  }; // end of struct InlineLaMETParams


  //! Inline measurement of 3-pt functions writing building-blocks
  /*! \ingroup inlinehadron */
  class InlineLaMET : public AbsInlineMeasurement 
  {
  public:
    // Default destructor
    ~InlineLaMET() {}
    // Constructor from param struct: copies param struct
    InlineLaMET(const InlineLaMETParams& p) : params(p) {}
    // Copy constructor
    InlineLaMET(const InlineLaMET& p) : params(p.params) {}

    // Getter for measurement frequency
    unsigned long getFrequency(void) const {return params.frequency;}

    //! Sets up the XML and invokes func, which does the acual work
    void operator()( const unsigned long update_no, XMLWriter& xml_out); 

  protected:
    //! Does the actual work
    void func(const unsigned long update_no, XMLWriter& xml_out); 

  private:
    //! The parameter structure; holds names of props, gauge field, XML, etc.
    InlineLaMETParams params;
  }; // end of class InlineLaMET



}; // end of namespace Chroma

#endif

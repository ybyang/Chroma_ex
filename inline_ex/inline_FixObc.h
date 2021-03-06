// -*- C++ -*-
/*! \file
 *  \brief Inline plaquette
 */

#ifndef __inline_test_cls_h__
#define __inline_test_cls_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/xml_group_reader.h"

namespace Chroma 
{ 
  /*! \ingroup inlineglue */
  namespace InlineOBCEnv 
  {
    extern const std::string name;
    bool registerAll();

    /*! \ingroup inlineglue */
    struct Params 
    {
      Params();
      Params(XMLReader& xml_in, const std::string& path);

      unsigned long frequency;

      struct NamedObject_t
      {
	std::string   gauge_id;
      } named_obj;

      std::string xml_file;
    };


    /*! \ingroup inlineglue */
    class InlineMeas : public AbsInlineMeasurement 
    {
    public:
      ~InlineMeas() {}
      InlineMeas(const Params& p) : params(p) {}
      InlineMeas(const InlineMeas& p) : params(p.params) {}

      unsigned long getFrequency(void) const {return params.frequency;}

      void operator()(unsigned long update_no,
		      XMLWriter& xml_out); 

    protected:
      void func(const unsigned long update_no, 
		XMLWriter& xml_out);

    private:
      Params params;
    };

  }

}

#endif

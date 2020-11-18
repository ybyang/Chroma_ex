// -*- C++ -*-
// $Id: inline_eigen_maker.h, v0.9 2020-07-22$
/*! \file
 * \brief Inline construction of the H-Wilson and Overlap eigensystem, or just load the existed one.
 *
 * Calculates the eigensystem 
 */

#ifndef __inline_EigenMaker_h__
#define __inline_EigenMaker_h__

#include "actions_ex/eigenop.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineEigenMakerEnv 
  {
    extern const std::string name;
    bool registerAll();
  }

  //! Parameter structure
  /*! \ingroup inlinehadron */
  struct InlineEigenMakerParams 
  {
    // Default constructor
    InlineEigenMakerParams(void);
    // Construct from XML
    InlineEigenMakerParams(XMLReader& xml_in, const std::string& path);
    // Write out the configuration of the parameters
    void write(XMLWriter& xml_out, const std::string& path);

    // How often should I be called during a gauge evolution?
    unsigned long frequency;

    // Holds the non-lattice parameters
    struct Param_t
    {
      GroupXML_t      fermact; 
      int  Noeigen;   /*<! the amount of the eigenvectors*/
      int gpu_level; /* 0 for pure CPU, 1 for CPU(vector)+GPU(calc), 2 for pure GPU*/
      bool check_residual; /* flag to check the residual of the data*/
      std::string  filename;          /*!< bb output file name */

      int  extra_space;   /*<! the extra space needed by this maker*/
      double  chebyshev_cut;   /*<! the upper band of the eigenvalues*/
      int  chebyshev_order;   /*<! the maxinum length of t*/
      bool use_ckpoint; /* flag to enable the checkpointing*/
  
      bool file_exist;
    } param;
    

    // Holds the names of the lattice objects
    struct NamedObject_t
    {
      std::string eigen_id;  /*!< output eigen id */
      std::string gauge_id;     /*!< input gauge id */
    } named_obj;

    std::string xml_file;  // Alternate XML file pattern
  }; // end of struct InlineEigenMakerParams


  //! Inline measurement of hadron spectrum
  /*! \ingroup inlinehadron */
  class InlineEigenMaker : public AbsInlineMeasurement 
  {
  public:
    // Default destructor
    ~InlineEigenMaker() {}
    // Constructor from param struct: copies param struct
    InlineEigenMaker(const InlineEigenMakerParams& p) : params(p) {}
    // Copy constructor
    InlineEigenMaker(const InlineEigenMaker& p) : params(p.params) {}

    // Getter for measurement frequency
    unsigned long getFrequency() const {return params.frequency;}

    //! Sets up the XML and invokes func, which does the acual work
    void operator()( const unsigned long update_no, XMLWriter& xml_out); 

  protected:
    //! Does the actual work
    void func(const unsigned long update_no, XMLWriter& xml_out); 

  private:
    //! The parameter structure; holds names of props, gauge field, XML, etc.
    InlineEigenMakerParams params;
  }; // end of class InlineEigenMaker

//-----------------------------------//
//-----------------------------------//

//-----------------------------------//
//-----------------------------------//

// Begin of class InlineEigenEraser

  struct InlineEigenEraserParams 
  {
      InlineEigenEraserParams() {};
      InlineEigenEraserParams(XMLReader& xml_in, const std::string& path)
      {
        XMLReader paramtop(xml_in, path);
        read(paramtop,"frequency",frequency);
        XMLReader input(paramtop, "NamedObject");
        read(input, "eigen_id", named_obj.eigen_id);
        frequency=1;
      }
      void writeXML(XMLWriter& xml_out, const std::string& path)
      {
        push(xml_out, path);
        write(xml_out,"frequency",frequency);
        push(xml_out, "NamedObject");
        write(xml_out, "eigen_id", named_obj.eigen_id);
        pop(xml_out);
        pop(xml_out);
      }       

      struct NamedObject_t
      {
        std::string  eigen_id;
      } named_obj;
      
      unsigned long frequency;
  };
  
  class InlineEigenEraser : public AbsInlineMeasurement 
  {
  public:
      ~InlineEigenEraser() {}
      InlineEigenEraser(const InlineEigenEraserParams& p) : params(p) {}

      unsigned long getFrequency(void) const {return params.frequency;}

      //! Do the writing
      void operator()(const unsigned long update_no,
                      XMLWriter& xml_out)
      {
        START_CODE();

        push(xml_out, "erase_eigen");
        write(xml_out, "update_no", update_no);

        QDPIO::cout << "EIGENERASER: Eigen erase" << std::endl;

        // Erase the object
        QDPIO::cout << "Attempt to erase object name = "
                        << params.named_obj.eigen_id << std::endl;
        write(xml_out, "object_id", params.named_obj.eigen_id);
        if (TheNamedObjMap::Instance().check(params.named_obj.eigen_id)) {
               EigenOperator<LatticeFermion> *eigen=TheNamedObjMap::Instance().getData<EigenOperator<LatticeFermion>* >
                          (params.named_obj.eigen_id);
               delete eigen;
               TheNamedObjMap::Instance().erase(params.named_obj.eigen_id);
        } else {
                QDPIO::cout << "Eigen system: " << params.named_obj.eigen_id
                                << " is not in the map. Cannot delete" << std::endl;
        }
        QDPIO::cout << "EIGENERASER: ran successfully" << std::endl;

        pop(xml_out);  // erase_named_obj

        END_CODE();
      }

    private:
      InlineEigenEraserParams params;
    }; //end of class InlineEigenEraser   
    

}; // end of namespace Chroma

#endif


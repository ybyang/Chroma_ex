#ifndef __eigenop_h__
#define __eigenop_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "linearop.h"
#include "meas/inline/io/named_objmap.h"
#include "io/xml_group_reader.h"
#include "actions/ferm/fermstates/fermstates.h"
#include "actions/ferm/fermbcs/fermbc_factory_w.h"

namespace Chroma 
{ 

    struct InlineEigenMakerParams;

    template <typename T>
    struct EigenPair{
        T vec;
        Complex val;
    };

    template <typename T>
    class EigenOperator: public UnprecLinearOperator<T, 
           multi1d<LatticeColorMatrix>,
           multi1d<LatticeColorMatrix> >, public multi1d<EigenPair<T> >
    {
    public:
        typedef multi1d<LatticeColorMatrix>  P;
        typedef multi1d<LatticeColorMatrix>  Q;

        EigenOperator(GroupXML_t &fermact,multi1d<LatticeColorMatrix> &u,int neig):multi1d<EigenPair<T> >(neig)
        {
             std::istringstream  is(fermact.xml);
             XMLReader  paramtop(is);
             read(paramtop, "FermAct", fermact_id);
             bc=new SimpleFermBC<T,P,Q>(SimpleFermBCParams(paramtop, "FermionBC"));
             fs=new SimpleFermState<T,P,Q>(bc,u);
        }

        const std::string &get_fermact()
        {
               return fermact_id;
        }

        int check_residual(Real residual_goal)
        {
              LatticeFermion tmp;
              EigenOperator<T> &op=*this;
              for(int ind=0;ind<this->size();ind++)
              {
//                   double res=sqrt(norm2(op[ind])).elem().elem().elem();
                   Real res=sqrt(norm2(op[ind].vec));
                   if(Layout::primaryNode())  printf("Norm evec[%4d] bias = %13.5e", ind, 1.0-res.elem().elem().elem().elem());
                   (*this)(tmp,op[ind].vec,PLUS);
                   tmp=tmp-op[ind].val*op[ind].vec;
                   res=sqrt(norm2(tmp));
                   if(Layout::primaryNode()) 
                       printf(" eval[%4d] = %15.12f + I %15.12f, ||lambda vec - mat vec|| = %10.5e\n",
                                ind,op[ind].val.elem().elem().elem().real(),
                                op[ind].val.elem().elem().elem().imag(),
                                res.elem().elem().elem());
              }
        }

        virtual void operator()(T& chi, const T& psi, enum PlusMinus isign) const
        {
             QDPIO::cerr << "The operator in this base class should not be used" <<std::endl;
             QDP_abort(1);
        }

        void save(std::string, bool to_single=true)
        {


        }

        void load(std::string, bool from_single=true)
        {


        }

        const FermBC<T,P,Q>& getFermBC() const {return *bc;}

        const FermState<T,P,Q>& getFermState() const {return *fs;}

        virtual int create_eigen(InlineEigenMakerParams &Params){};

    protected:
        Handle<FermState<T,P,Q> > fs;
        Handle<FermBC<T,P,Q> > bc;
        std::string fermact_id;
    };


}; // end of namespace Chroma

#endif

// -*- C++ -*-
/*! \file
 *  \brief Clover term linear operator
 */

#ifndef __clover_utils_exp_w_h__
#define __clover_utils_exp_w_h__

#include "state.h"
#include "qdp_allocator.h"
#include "actions/ferm/fermacts/clover_fermact_params_w.h"
#include "actions/ferm/linop/clover_term_base_w.h"
#include "actions/ferm/linop/clover_term_qdp_w.h"
#include "meas/glue/mesfield.h"
#include <complex>
namespace Chroma 
{ 

  template<typename R>
  inline void pauli_extract(PrimitiveClovTriang<R> &A, RComplex<R>* res,int ib)
  {
        for(int i=0;i<6;i++)
             res[7*i]=A.diag[ib][i];
        int off=0;
        for(int i=1;i<6;i++)
        for(int j=0;j<i;j++)
        {
             res[i*6+j]=A.offd[ib][off];
             res[j*6+i]=conj(A.offd[ib][off]);
             off++;
        }
  }

  template<typename R>
  inline void pauli_compress(RComplex<R>* res, PrimitiveClovTriang<R> &A,int ib)
  {
        for(int i=0;i<6;i++)
             A.diag[ib][i]=res[7*i].real();
        int off=0;
        for(int i=1;i<6;i++)
        for(int j=0;j<i;j++)
        {
             A.offd[ib][off]=res[i*6+j];
             off++;
        }  
  }

  template<typename R>
  inline void pauli_prod(PrimitiveClovTriang<R> &A,PrimitiveClovTriang<R> &B,PrimitiveClovTriang<R> &C)
  {
       std::vector<RComplex<R>> tmpA(36),tmpB(36),tmpC(36);
       for(int ib=0;ib<2;ib++)
       {
             pauli_extract<R>(A,tmpA.data(),ib);
             pauli_extract<R>(B,tmpB.data(),ib);
             for(int i=0;i<6;i++)
             for(int j=0;j<=i;j++)
             {
                tmpC[i*6+j]=tmpA[i*6]*tmpB[j];
             for(int k=1;k<6;k++)
                tmpC[i*6+j]+=tmpA[i*6+k]*tmpB[k*6+j];
             }
             pauli_compress<R>(tmpC.data(),C,ib);
       }
  }

   template<typename R>
  inline void pauli_prod(R *fac,PrimitiveClovTriang<R> &A,PrimitiveClovTriang<R> &C)
  {
          for(int ib=0;ib<2;ib++)
       {
            for(int ic=0;ic<6;ic++)
                  C.diag[ib][ic]=RScalar<R>(fac[ib])*A.diag[ib][ic];
            for(int ic=0;ic<15;ic++)
                  C.offd[ib][ic]=RScalar<R>(fac[ib])*A.offd[ib][ic];
       }
    
  }

  template<typename R>
  inline void pauli_add_fac(R *fac, PrimitiveClovTriang<R> &A,PrimitiveClovTriang<R> &B,PrimitiveClovTriang<R> &C)
  {
       for(int ib=0;ib<2;ib++)
       {
            for(int ic=0;ic<6;ic++)
                  C.diag[ib][ic]=RScalar<R>(fac[ib])+RScalar<R>(fac[2+ib])*A.diag[ib][ic]+RScalar<R>(fac[4+ib])*B.diag[ib][ic];
            for(int ic=0;ic<15;ic++)
                  C.offd[ib][ic]=RScalar<R>(fac[2+ib])*A.offd[ib][ic]+RScalar<R>(fac[4+ib])*B.offd[ib][ic];
       }
  }   

  template<typename R>
  inline void pauli_add_fac(R *fac, PrimitiveClovTriang<R> &A,PrimitiveClovTriang<R> &C)
  {
       for(int ib=0;ib<2;ib++)
       {
            for(int ic=0;ic<6;ic++)
                  C.diag[ib][ic]=RScalar<R>(fac[ib])+RScalar<R>(fac[2+ib])*A.diag[ib][ic];
            for(int ic=0;ic<15;ic++)
                  C.offd[ib][ic]=RScalar<R>(fac[2+ib])*A.offd[ib][ic];
       }
  }

  

  template<typename R>
  inline void pauli_add(PrimitiveClovTriang<R> &A,PrimitiveClovTriang<R> &B,PrimitiveClovTriang<R> &C)
  {
       R fac[6]={0.0,0.0,1.0,1.0,1.0,1.0};
       pauli_add_fac(fac,A,B,C);
  } 

  template<typename R>
  inline R pauli_tr(PrimitiveClovTriang<R> &A, int ib)
  {
       RScalar<R> tmp=A.diag[ib][0]+A.diag[ib][1]+A.diag[ib][2]+A.diag[ib][3]+A.diag[ib][4]+A.diag[ib][5];
       return tmp.elem();
  }

  template<typename R>
  inline R pauli_tr2(PrimitiveClovTriang<R> &A,PrimitiveClovTriang<R> &B, int ib)
  {
            RScalar<R> res=0.0;
            for(int ic=0;ic<6;ic++)
                res+=A.diag[ib][ic]*B.diag[ib][ic];
            R tmp=0.0;
            for(int ic=0;ic<15;ic++)
                tmp+=A.offd[ib][ic].real()*B.offd[ib][ic].real()+A.offd[ib][ic].imag()*B.offd[ib][ic].imag();
            return res.elem()+2*tmp;
  }


  //! Clover term
  /*!
   * \ingroup linop
   *
   */
  template<typename T, typename U>
  class ExpCloverTermT : public CloverTermBase<T, U>
  {
  public:
    // Typedefs to save typing
    typedef typename WordType<T>::Type_t REALT;

    typedef OLattice< PScalar< PScalar< RScalar< typename WordType<T>::Type_t> > > > LatticeREAL;
    typedef OScalar< PScalar< PScalar< RScalar<REALT> > > > RealT;

    //! Empty constructor. Must use create later
    ExpCloverTermT();

    //! No real need for cleanup here
    ~ExpCloverTermT() {
    	if ( tri != nullptr ) {
    		QDP::Allocator::theQDPAllocator::Instance().free(tri);
    	}
    	if ( tri_ori != nullptr ) {
    		QDP::Allocator::theQDPAllocator::Instance().free(tri_ori);
    	}
 /*      
        while (!ksv.empty()) {
		REALT* p = ksv.back();
		delete p; p = NULL;
		ksv.pop_back();
	}      

        for (auto it = ksv.begin(); it != ksv.end(); it++)
        {
	if (*it != NULL)
	{
		delete *it;
 		*it = NULL;
		}
	} 


       for(int i=0 ; i<ksv.size() ; i++){

        delete ksv[i];
                                        }
        ksv.clear(); 

*/
    }

    //internal functions
    int exp_order();
    void make_exp(int N_, int s, RealT fac, int cb);

    //! Creation routine
    void create(Handle< FermState<T, multi1d<U>, multi1d<U> > > fs,
		const CloverFermActParams& param_,bool is_exp);

    virtual void create(Handle< FermState<T, multi1d<U>, multi1d<U> > > fs,
			const CloverFermActParams& param_,
			const ExpCloverTermT<T,U>& from_,bool is_exp);

    //! Computes the inverse of the term on cb using Cholesky
    /*!
     * \param cb   checkerboard of work (Read)
     */
    void choles(int cb);

    //! Computes the inverse of the term on cb using Cholesky
    /*!
     * \param cb   checkerboard of work (Read)
     * \return logarithm of the determinant  
     */
    Double cholesDet(int cb) const ;

    /**
     * Apply a dslash
     *
     * Performs the operation
     *
     *  chi <-   (L + D + L^dag) . psi
     *
     * where
     *   L       is a lower triangular matrix
     *   D       is the real diagonal. (stored together in type TRIANG)
     *
     * Arguments:
     * \param chi     result                                      (Write)
     * \param psi     source                                      (Read)
     * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
     * \param cb      Checkerboard of OUTPUT std::vector               (Read) 
     */
    void apply (T& chi, const T& psi, enum PlusMinus isign, int cb) const;

 
    void apply2( int N_, int s, RealT fac,  T (&chi)[8], const T& psi1, const T& psi2,  enum PlusMinus isign, int cb) ;


    void applySite(T& chi, const T& psi, enum PlusMinus isign, int site) const;

    const PrimitiveClovTriang<REALT>* getTriBuffer() const {
      return tri;
    }

    //! Calculates Tr_D ( Gamma_mat L )
    void triacntr(U& B, int mat, int cb) const;

    //! Return the fermion BC object for this linear operator
    const FermBC<T, multi1d<U>, multi1d<U> >& getFermBC() const {return *fbc;}

    //! PACK UP the Clover term for QUDA library:
    void packForQUDA(multi1d<QUDAPackedClovSite<REALT> >& quda_pack, int cb) const; 

    void cholesExp(int cb);

    void derivExp(multi1d<U>& ds_u, 
               const T& chi, const T& psi, 
               enum PlusMinus isign, int cb) const;

    void derivExp(multi1d<U>& ds_u,
               const T& chi, const T& psi,
               enum PlusMinus isign) const;


    void derivExpMultipole(multi1d<U>& ds_u,
                        const multi1d<T>& chi, const multi1d<T>& psi,
                        enum PlusMinus isign) const ;

    void derivExpMultipole(multi1d<U>& ds_u,
                        const multi1d<T>& chi, const multi1d<T>& psi,
                        enum PlusMinus isign, int cb) const ;

    void derivExpTrLn(multi1d<U>& ds_u, 
                   enum PlusMinus isign, int cb) const ;

    void derivExp_loops(const int u, const int mu, const int cb,
                     U& ds_u_mu,
                     U& ds_u_nu,
                     const U& Lambda) const {};
                     
    bool Is_Exp() const{return is_exp;};

      
  protected:
    //! Create the clover term on cb
    /*!
     *  \param f         field strength tensor F(mu,nu)        (Read)
     *  \param cb        checkerboard                          (Read)
     */
    void makeClov(const multi1d<U>& f, const RealT& diag_mass);

    //! Invert the clover term on cb
    void chlclovms(LatticeREAL& log_diag, int cb);
    void ldagdlinv(LatticeREAL& tr_log_diag, int cb);

    //! Get the u field
    const multi1d<U>& getU() const {return u;}

    //! Calculates Tr_D ( Gamma_mat L )
    Real getCloverCoeff(int mu, int nu) const;

  private:
    Handle< FermBC<T,multi1d<U>,multi1d<U> > >      fbc;
    multi1d<U>  u;
    CloverFermActParams          param;
    LatticeREAL                  tr_log_diag_; // Fill this out during create
                                                  // but save the global sum until needed.
    multi1d<bool> choles_done;   // Keep note of whether the decomposition has been done
                                 // on a particular checkerboard. 

    PrimitiveClovTriang<REALT>*  tri; // is_exp=true: (4+m)exp(clv/(4+m)); false: 4+m+clv;
    PrimitiveClovTriang<REALT>*  tri_ori; // clv/(4+m)
    std::vector<REALT> csv[2];
    std::vector<REALT> ksv0;
    std::vector< std::vector<REALT> >ksv;
    int N;
    RealT diag_mass;    
    bool is_exp;    
  };


   // Empty constructor. Must use create later
  template<typename T, typename U>
  ExpCloverTermT<T,U>::ExpCloverTermT() {
	  // Always allocate on construction
	  int nodeSites = Layout::sitesOnNode();
	  tri = (PrimitiveClovTriang<REALT>*)QDP::Allocator::theQDPAllocator::Instance().allocate(nodeSites*sizeof(PrimitiveClovTriang<REALT>),
			  	  QDP::Allocator::DEFAULT);

	  tri_ori = (PrimitiveClovTriang<REALT>*)QDP::Allocator::theQDPAllocator::Instance().allocate(nodeSites*sizeof(PrimitiveClovTriang<REALT>),
			  	  QDP::Allocator::DEFAULT);

  }

  // Now copy
  template<typename T, typename U>
  void ExpCloverTermT<T,U>::create(Handle< FermState<T,multi1d<U>,multi1d<U> > > fs,
				   const CloverFermActParams& param_,
				   const ExpCloverTermT<T,U>& from,bool is_exp_)
  {
#ifndef QDP_IS_QDPJIT
    START_CODE();
    for(int i=0; i < rb.numSubsets(); i++) 
    if( from.choles_done[i]==true&&is_exp==true)
    {
        QDPIO::cerr << "ExpCloverTerm: error: ExpClover Should be created from original clover term, not the inversed one." << std::endl;
        QDP_abort(1);
    }
    
    is_exp=is_exp_;
    
    u.resize(Nd);

    u = fs->getLinks();
    fbc = fs->getFermBC();
    param = param_;
    
    // Sanity check
    if (fbc.operator->() == 0) {
      QDPIO::cerr << "ExpCloverTerm: error: fbc is null" << std::endl;
      QDP_abort(1);
    }
    
    {
      RealT ff = where(param.anisoParam.anisoP, Real(1) / param.anisoParam.xi_0, Real(1));
      param.clovCoeffR *= Real(0.5) * ff;
      param.clovCoeffT *= Real(0.5);
    }
    
    //
    // Yuk. Some bits of knowledge of the dslash term are buried in the 
    // effective mass term. They show up here. If I wanted some more 
    // complicated dslash then this will have to be fixed/adjusted.
    //
    {
      RealT ff = where(param.anisoParam.anisoP, param.anisoParam.nu / param.anisoParam.xi_0, Real(1));
      diag_mass = 1 + (Nd-1)*ff + param.Mass;
    }
    RealT inv_diag_mass=1.0/diag_mass; 

    int nodeSites = Layout::sitesOnNode();
    if(&from==this)
    {
       multi1d<LatticeColorMatrix> f;
       mesField(f, u);
       makeClov(f, diag_mass);
    }
    else
    {
       choles_done.resize(rb.numSubsets());
       for(int i=0; i < rb.numSubsets(); i++) {
         choles_done[i] = from.choles_done[i];
       }
    
    // This is for the whole lattice (LatticeReal)
       tr_log_diag_ = from.tr_log_diag_;
   
    // Deep copy.
#pragma omp parallel for
        for(int site=0; site < nodeSites;++site)
    	for(int block=0; block < 2; ++block) {
    		for(int d=0; d < 6; ++d) {
    			tri[site].diag[block][d] = from.tri[site].diag[block][d];
    		}
    		for(int od=0; od < 15; ++od) {
    			tri[site].offd[block][od] = from.tri[site].offd[block][od];
    		}
    	}
    }
    
    if(is_exp==false)
    {
      if(from.Is_Exp()==true&&(&from)!=this) 
      {
#pragma omp parallel for
        for(int site=0; site < nodeSites;++site) 
        for(int block=0; block < 2; ++block) {
           for(int d=0; d < 6; ++d) {
                        tri[site].diag[block][d] = (from.tri_ori[site].diag[block][d]+RScalar<REALT>(1.0))*diag_mass.elem().elem().elem();
           }
           for(int od=0; od < 15; ++od) {
                        tri[site].offd[block][od] = from.tri_ori[site].offd[block][od]*diag_mass.elem().elem().elem();
           }
        }
        // if &from==this, then from.Is_Exp()==false;
        // do nothing for (&from)!=this && from.Is_Exp()=false.
      }
    }
    else
    {
      if(from.Is_Exp()==true&&(&from)!=this)
      {
    // Deep copy.
#pragma omp parallel for
        for(int site=0; site < nodeSites;++site)
    	for(int block=0; block < 2; ++block) {
    		for(int d=0; d < 6; ++d) {
    			tri_ori[site].diag[block][d] = from.tri_ori[site].diag[block][d];
    		}
    		for(int od=0; od < 15; ++od) {
    			tri_ori[site].offd[block][od] = from.tri_ori[site].offd[block][od];
    		}
    	}
      }
      else
      {
        // if &from==this, then from.Is_Exp()==true;
        // make_exp for both for &from==this and ((&from)!=this && from.Is_Exp()=false)
#pragma omp parallel for
    	for(int site=0; site < nodeSites;++site) 
    	for(int block=0; block < 2; ++block) {
           for(int d=0; d < 6; ++d) {
                        tri_ori[site].diag[block][d] = tri[site].diag[block][d]*inv_diag_mass.elem().elem().elem()-RScalar<REALT>(1.0);
                        tri[site].diag[block][d] = diag_mass.elem().elem().elem();
           }
           for(int od=0; od < 15; ++od) {
                        tri_ori[site].offd[block][od] = tri[site].offd[block][od]*inv_diag_mass.elem().elem().elem();
                        tri[site].offd[block][od] = RComplex<REALT>(0.0,0.0);
           }
        }

        N=exp_order();

        END_CODE();

        for(int cb=0; cb < 2; ++cb)
	   make_exp(N,0,diag_mass,cb);
      }
    }

#endif
  }


  //! Creation routine
  template<typename T, typename U>
  void ExpCloverTermT<T,U>::create(Handle< FermState<T,multi1d<U>,multi1d<U> > > fs,
				   const CloverFermActParams& param_,bool is_exp)
  {
    choles_done.resize(rb.numSubsets());
    for(int i=0; i < rb.numSubsets(); i++)
      choles_done[i] = false;
    create(fs,param_,*this,is_exp);
  }

  /* This now just sets up and dispatches... */
  template<typename T, typename U>
  void ExpCloverTermT<T,U>::makeClov(const multi1d<U>& f, const RealT& diag_mass)
  {
    START_CODE();
    
    if ( Nd != 4 ){
      QDPIO::cerr << __func__ << ": expecting Nd==4" << std::endl;
      QDP_abort(1);
    }
    
    if ( Ns != 4 ){
      QDPIO::cerr << __func__ << ": expecting Ns==4" << std::endl;
      QDP_abort(1);
    }
  
    U f0 = f[0] * getCloverCoeff(0,1);
    U f1 = f[1] * getCloverCoeff(0,2);
    U f2 = f[2] * getCloverCoeff(0,3);
    U f3 = f[3] * getCloverCoeff(1,2);
    U f4 = f[4] * getCloverCoeff(1,3);
    U f5 = f[5] * getCloverCoeff(2,3);    

    const int nodeSites = QDP::Layout::sitesOnNode();
    QDPCloverEnv::QDPCloverMakeClovArg<U> arg = {diag_mass, f0,f1,f2,f3,f4,f5,tri };
    dispatch_to_threads(nodeSites, arg, QDPCloverEnv::makeClovSiteLoop<U>);

    END_CODE();
  }
  
  template<typename T, typename U>
  int ExpCloverTermT<T,U>::exp_order()
  {
      int n=0;

      RealT c0=3.0*2*this->param.clovCoeffT/diag_mass; // factor 2 for the defintion used in Chroma;
      REALT c=c0.elem().elem().elem().elem();
      REALT a=c*exp(c);
      REALT b=fuzz.elem().elem().elem().elem()*0.01; // epsl

      for (n=1;n<=100;n++)
      {
         a*=c;
         b*=(double)(n+1);

         if (a<b)
            break;
      }
      if(n==100)
           QDPIO::cout << "ExpCloverTerm: warning: Clover order is higher than 100" << std::endl;
      for(int i=0;i<2;i++) csv[i].resize(n+1);
      csv[0][0]=1.0;csv[1][0]=1.0;
      
      ksv0.resize((n+1)*n/2);
      ksv.resize(n+1);
//      REALT *p=ksv0.data();
//       REALT *p= new REALT;
//       p=ksv0.data();


      for (int k=0;k<=n;k++)
      {
          if(k>0)
          {
              csv[0][k]=csv[0][k-1]/((REALT)(k));
              if (k&0x1)
                  csv[1][k]=-csv[0][k];
              else
                  csv[1][k]=csv[0][k];
          }
          
//          ksv[k]=p;p+=n-k+1;
          ksv[k].resize(n-k+1);          
          ksv[k][0]=csv[0][k]/((REALT)(k+1));
          for (int l=1;l<=(n-k);l++)
               ksv[k][l]=ksv[k][l-1]/((REALT)(k+l+1));
      }
      
      return n;
  }
  
  namespace ExpCloverEnv {   

    template<typename R>
    inline
    void clv_cpoly(PrimitiveClovTriang<R> &A, PrimitiveClovTriang<R> &A2, PrimitiveClovTriang<R> &A3, R *psv)
    {
         pauli_prod<R>(A,A,A2);
         pauli_prod<R>(A,A2,A3); 
         
       for(int ib=0;ib<2;ib++)
       {  
         R tr[5]={pauli_tr<R>(A2,ib),pauli_tr<R>(A3,ib),pauli_tr2<R>(A2,A2,ib),pauli_tr2<R>(A2,A3,ib),pauli_tr2<R>(A3,A3,ib)};
         
         psv[0*2+ib]=(1.0/144.0)*(8.0*tr[1]*tr[1]-24.0*tr[4]+tr[0]*(18.0*tr[2]-3.0*tr[0]*tr[0]));
         psv[1*2+ib]=(1.0/30.0)*(5.0*tr[0]*tr[1]-6.0*tr[3]);
         psv[2*2+ib]=(1.0/8.0)*(tr[0]*tr[0]-2.0*tr[2]);
         psv[3*2+ib]=(-1.0/3.0)*tr[1];
         psv[4*2+ib]=-0.5*tr[0];         
       }
    }
    
    template<typename R>
    inline 
    void clv_cayley(int N,std::vector<R> &c,R* psv, R *q)
    {
        R q5;
        for(int ib=0;ib<2;ib++)
        if(N>5)
        {
            for(int i=0;i<6;i++)
                 q[ib+2*i]=c[N-5+i];
            
            for(int i=N-6;i>=0;i--)
            {
                 q5=q[ib+2*5];
                 q[ib+2*5]=q[ib+2*4];
                 q[ib+2*4]=q[ib+2*3]-q5*psv[ib+2*4];
                 q[ib+2*3]=q[ib+2*2]-q5*psv[ib+2*3];
                 q[ib+2*2]=q[ib+2*1]-q5*psv[ib+2*2];
                 q[ib+2*1]=q[ib+2*0]-q5*psv[ib+2*1];
                 q[ib+2*0]=c[i]-q5*psv[ib+2*0];
            }
        }
        else
        {
            for (int k=0;k<=5;k++)
            {
               if (k<=N)
                  q[ib+2*k]=c[k];
               else
                  q[ib+2*k]=0.0;
            }
        }
    }
  
    template<typename REALT>
    struct ExpCloverMakeExpArg {
      const REALT& fac;
      const int N;
      const int cb;
      PrimitiveClovTriang < REALT >* tri;
      PrimitiveClovTriang < REALT >* tri_ori;
      std::vector<REALT> csv;
    };

    template<typename REALT>
    inline 
    void MakeExp(int lo, int hi, int myId, ExpCloverMakeExpArg<REALT> *a)
    {
#ifndef QDP_IS_QDPJIT    
      const REALT& fac=a->fac;
      int N=a->N;
      int cb=a->cb;
      PrimitiveClovTriang < REALT >* tri=a->tri;
      PrimitiveClovTriang < REALT >* tri_ori=a->tri_ori;
      std::vector<REALT> &csv=a->csv;

      REALT psv[10],q[12];
      PrimitiveClovTriang < REALT > A2,A3,tmp;

      int test;

      for(int ssite=lo; ssite < hi; ++ssite)  {
          int site = rb[cb].siteTable()[ssite];
          test =site; 
          clv_cpoly<REALT>(tri_ori[site],A2,A3,psv);
          clv_cayley<REALT>(N,csv,psv,q);

          for(int k=0;k<12;k++)
                q[k]*=fac;

          pauli_add_fac<REALT>(q+6,tri_ori[site],A2,tmp); // tmp=c3+c4*A+c5*A^2
          pauli_prod<REALT>(A3,tmp,tri[site]); // tri[site]=A^3*(c3+c4*A+c5*A^2)
          pauli_add_fac<REALT>(q,tri_ori[site],A2,tmp); // tmp=c0+c1*A+c2*A^2
          pauli_add<REALT>(tri[site],tmp,tri[site]); //tri[site]=tri[site]+tmp
        

      }/* End Site Loop */
/*
          int site= test;
          QDPIO::cout << "     tri_ori   "  <<   "factor=" << fac  << std::endl;
          int off=0;        
for(int i=0; i<6 ;i++) QDPIO::cout << "     A["<<i<<"]["<<i<<"]=   "    <<  tri_ori[site].diag[0][i]<< std::endl;
for (int i=0; i<6 ;i++) QDPIO::cout << "    B["<<i<<"]["<<i<<"]=   "    <<  tri[site].diag[0][i]<< std::endl;
for(int i=1; i<6 ;i++)
         {for(int j=0; j<i ;j++)
         {
           QDPIO::cout <<"A["<<j<<"]["<<i<<"]=" <<tri_ori[site].offd[0][off].real()<<"+"<<tri_ori[site].offd[0][off].imag()<<"*I"<< std::endl;
           QDPIO::cout <<"B["<<j<<"]["<<i<<"]=" <<tri[site].offd[0][off].real()<<"+"<<tri[site].offd[0][off].imag()<<"*I"<< std::endl;
           off+=1;
         }}
*/

#endif
   }/* End Function */
  }/* End Namespace */
  

    template<typename T, typename U>
    void ExpCloverTermT<T,U>::make_exp(int N_, int s, RealT fac, int cb)
    {
      START_CODE();
      typedef typename WordType<U>::Type_t REALT;
      const int nodeSites = rb[cb].numSiteTable();
      ExpCloverEnv::ExpCloverMakeExpArg<REALT> arg={fac.elem().elem().elem().elem(),N,cb,tri,tri_ori,csv[s]};
      dispatch_to_threads(nodeSites, arg, ExpCloverEnv::MakeExp<REALT>);
      END_CODE();
    }

  template<typename T, typename U> 
  void ExpCloverTermT<T,U>::cholesExp(int cb)  // exp(-clv/(4+m))/(4+m)
  {
     if(cb<0||cb>=rb.numSubsets())
     {
        QDPIO::cerr << "ExpCloverTerm: error: cb = " << cb << "is not allowed." << std::endl;
        QDP_abort(1);
     }
     QDPIO::cout<< "N=" << N << "       begin make_exp(N,1,1.0/diag_mass,cb)  " << std::endl;
     make_exp(N,1,1.0/diag_mass,cb);
     QDPIO::cout << "end make_exp(N,1,1.0/diag_mass,cb)  " << std::endl;
     choles_done[cb] = true;
  }
  
  //! Invert
  /*!
   * Computes the inverse of the term on cb using Cholesky
   */
  template<typename T, typename U>
  void ExpCloverTermT<T,U>::choles(int cb)
  {
    START_CODE();

    // When you are doing the cholesky - also fill out the trace_log_diag piece)
    // chlclovms(tr_log_diag_, cb);
    // Switch to LDL^\dag inversion
    ldagdlinv(tr_log_diag_,cb);

    END_CODE();
  }


  //! Invert
  /*!
   * Computes the inverse of the term on cb using Cholesky
   *
   * \return logarithm of the determinant  
   */
  template<typename T, typename U>
  Double ExpCloverTermT<T,U>::cholesDet(int cb) const
  {
#ifndef QDP_IS_QDPJIT
    START_CODE();

    if( choles_done[cb] == false ) 
    {
      QDPIO::cout << __func__ << ": Error: you have not done the Cholesky.on this operator on this subset" << std::endl;
      QDPIO::cout << "You sure you should not be asking invclov?" << std::endl;
      QDP_abort(1);
    }


    END_CODE();

    // Need to thread generic sums in QDP++?
    // Need to thread generic norm2() in QDP++?


    return sum(tr_log_diag_, rb[cb]);
#else
    assert(!"ni");
    Double ret=0.;
    return ret;
#endif
  }

  /*! An LDL^\dag decomposition and inversion? */
  template<typename T, typename U>
  void ExpCloverTermT<T,U>::ldagdlinv(LatticeREAL& tr_log_diag, int cb)
  {
#ifndef QDP_IS_QDPJIT
    START_CODE();

    if ( 2*Nc < 3 )
    {
      QDPIO::cerr << __func__ << ": Matrix is too small" << std::endl;
      QDP_abort(1);
    }

    // Zero trace log
    tr_log_diag[rb[cb]] = zero;

    QDPCloverEnv::LDagDLInvArgs<U> a = { tr_log_diag, tri, cb };
    int num_site_table = rb[cb].numSiteTable();
    dispatch_to_threads(num_site_table, a, QDPCloverEnv::LDagDLInvSiteLoop<U>);

    
    // This comes from the days when we used to do Cholesky
    choles_done[cb] = true;

    END_CODE();
#endif
  }
 
  template<typename T, typename U>
  void ExpCloverTermT<T,U>::chlclovms(LatticeREAL& tr_log_diag, int cb)
  {
#ifndef QDP_IS_QDPJIT
    START_CODE();

    if ( 2*Nc < 3 )
    {
      QDPIO::cerr << __func__ << ": Matrix is too small" << std::endl;
      QDP_abort(1);
    }
  
    tr_log_diag = zero;
    QDPCloverEnv::LDagDLInvArgs<U> a = { tr_log_diag, tri, cb};
    dispatch_to_threads(rb[cb].numSiteTable(), a, QDPCloverEnv::cholesSiteLoop<U>);
    
    choles_done[cb] = true;
    END_CODE();
#endif
  }


  /**
   * Apply a dslash
   *
   * Performs the operation
   *
   *  chi <-   (L + D + L^dag) . psi
   *
   * where
   *   L       is a lower triangular matrix
   *   D       is the real diagonal. (stored together in type TRIANG)
   *
   * Arguments:
   * \param chi     result                                      (Write)
   * \param psi     source                                      (Read)
   * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
   * \param cb      Checkerboard of OUTPUT std::vector               (Read) 
   */
  template<typename T, typename U>
  void ExpCloverTermT<T,U>::applySite(T& chi, const T& psi, 
			    enum PlusMinus isign, int site) const
  {
#ifndef QDP_IS_QDPJIT
    START_CODE();

    if ( Ns != 4 )
    {
      QDPIO::cerr << __func__ << ": CloverTerm::applySite requires Ns==4" << std::endl;
      QDP_abort(1);
    }

    int n = 2*Nc;

    RComplex<REALT>* cchi = (RComplex<REALT>*)&(chi.elem(site).elem(0).elem(0));
    const RComplex<REALT>* ppsi = (const RComplex<REALT>*)&(psi.elem(site).elem(0).elem(0));


    cchi[ 0] = tri[site].diag[0][ 0]  * ppsi[ 0]
      +   conj(tri[site].offd[0][ 0]) * ppsi[ 1]
      +   conj(tri[site].offd[0][ 1]) * ppsi[ 2]
      +   conj(tri[site].offd[0][ 3]) * ppsi[ 3]
      +   conj(tri[site].offd[0][ 6]) * ppsi[ 4]
      +   conj(tri[site].offd[0][10]) * ppsi[ 5];
    
    cchi[ 1] = tri[site].diag[0][ 1]  * ppsi[ 1]
      +        tri[site].offd[0][ 0]  * ppsi[ 0]
      +   conj(tri[site].offd[0][ 2]) * ppsi[ 2]
      +   conj(tri[site].offd[0][ 4]) * ppsi[ 3]
      +   conj(tri[site].offd[0][ 7]) * ppsi[ 4]
      +   conj(tri[site].offd[0][11]) * ppsi[ 5];
    
    cchi[ 2] = tri[site].diag[0][ 2]  * ppsi[ 2]
      +        tri[site].offd[0][ 1]  * ppsi[ 0]
      +        tri[site].offd[0][ 2]  * ppsi[ 1]
      +   conj(tri[site].offd[0][ 5]) * ppsi[ 3]
      +   conj(tri[site].offd[0][ 8]) * ppsi[ 4]
      +   conj(tri[site].offd[0][12]) * ppsi[ 5];
    
    cchi[ 3] = tri[site].diag[0][ 3]  * ppsi[ 3]
      +        tri[site].offd[0][ 3]  * ppsi[ 0]
      +        tri[site].offd[0][ 4]  * ppsi[ 1]
      +        tri[site].offd[0][ 5]  * ppsi[ 2]
      +   conj(tri[site].offd[0][ 9]) * ppsi[ 4]
      +   conj(tri[site].offd[0][13]) * ppsi[ 5];
    
    cchi[ 4] = tri[site].diag[0][ 4]  * ppsi[ 4]
      +        tri[site].offd[0][ 6]  * ppsi[ 0]
      +        tri[site].offd[0][ 7]  * ppsi[ 1]
      +        tri[site].offd[0][ 8]  * ppsi[ 2]
      +        tri[site].offd[0][ 9]  * ppsi[ 3]
      +   conj(tri[site].offd[0][14]) * ppsi[ 5];
    
    cchi[ 5] = tri[site].diag[0][ 5]  * ppsi[ 5]
      +        tri[site].offd[0][10]  * ppsi[ 0]
      +        tri[site].offd[0][11]  * ppsi[ 1]
      +        tri[site].offd[0][12]  * ppsi[ 2]
      +        tri[site].offd[0][13]  * ppsi[ 3]
      +        tri[site].offd[0][14]  * ppsi[ 4];
    
    cchi[ 6] = tri[site].diag[1][ 0]  * ppsi[ 6]
      +   conj(tri[site].offd[1][ 0]) * ppsi[ 7]
      +   conj(tri[site].offd[1][ 1]) * ppsi[ 8]
      +   conj(tri[site].offd[1][ 3]) * ppsi[ 9]
      +   conj(tri[site].offd[1][ 6]) * ppsi[10]
      +   conj(tri[site].offd[1][10]) * ppsi[11];
    
    cchi[ 7] = tri[site].diag[1][ 1]  * ppsi[ 7]
      +        tri[site].offd[1][ 0]  * ppsi[ 6]
      +   conj(tri[site].offd[1][ 2]) * ppsi[ 8]
      +   conj(tri[site].offd[1][ 4]) * ppsi[ 9]
      +   conj(tri[site].offd[1][ 7]) * ppsi[10]
      +   conj(tri[site].offd[1][11]) * ppsi[11];
    
    cchi[ 8] = tri[site].diag[1][ 2]  * ppsi[ 8]
      +        tri[site].offd[1][ 1]  * ppsi[ 6]
      +        tri[site].offd[1][ 2]  * ppsi[ 7]
      +   conj(tri[site].offd[1][ 5]) * ppsi[ 9]
      +   conj(tri[site].offd[1][ 8]) * ppsi[10]
      +   conj(tri[site].offd[1][12]) * ppsi[11];
    
    cchi[ 9] = tri[site].diag[1][ 3]  * ppsi[ 9]
      +        tri[site].offd[1][ 3]  * ppsi[ 6]
      +        tri[site].offd[1][ 4]  * ppsi[ 7]
      +        tri[site].offd[1][ 5]  * ppsi[ 8]
      +   conj(tri[site].offd[1][ 9]) * ppsi[10]
      +   conj(tri[site].offd[1][13]) * ppsi[11];
    
    cchi[10] = tri[site].diag[1][ 4]  * ppsi[10]
      +        tri[site].offd[1][ 6]  * ppsi[ 6]
      +        tri[site].offd[1][ 7]  * ppsi[ 7]
      +        tri[site].offd[1][ 8]  * ppsi[ 8]
      +        tri[site].offd[1][ 9]  * ppsi[ 9]
      +   conj(tri[site].offd[1][14]) * ppsi[11];
    
    cchi[11] = tri[site].diag[1][ 5]  * ppsi[11]
      +        tri[site].offd[1][10]  * ppsi[ 6]
      +        tri[site].offd[1][11]  * ppsi[ 7]
      +        tri[site].offd[1][12]  * ppsi[ 8]
      +        tri[site].offd[1][13]  * ppsi[ 9]
      +        tri[site].offd[1][14]  * ppsi[10];


    END_CODE();
#endif
  }

  template<typename T, typename U>
  void ExpCloverTermT<T,U>::triacntr(U& B, int mat, int cb) const
  {
#ifndef QDP_IS_QDPJIT
    START_CODE();

    B = zero;

    if ( mat < 0  ||  mat > 15 )
    {
      QDPIO::cerr << __func__ << ": Gamma out of range: mat = " << mat << std::endl;
      QDP_abort(1);
    }

    QDPCloverEnv::TriaCntrArgs<U> a = { B, tri, mat, cb };
    dispatch_to_threads(rb[cb].numSiteTable(), a, 
			QDPCloverEnv::triaCntrSiteLoop<U>);

    END_CODE();
#endif
  }

  //! Returns the appropriate clover coefficient for indices mu and nu
  template<typename T, typename U>
  Real
  ExpCloverTermT<T,U>::getCloverCoeff(int mu, int nu) const 
  { 
    START_CODE();

    if( param.anisoParam.anisoP )  {
      if (mu==param.anisoParam.t_dir || nu == param.anisoParam.t_dir) { 
	return param.clovCoeffT;
      }
      else { 
	// Otherwise return the spatial coeff
	return param.clovCoeffR;
      }
    }
    else { 
      // If there is no anisotropy just return the spatial one, it will
      // be the same as the temporal one
      return param.clovCoeffR; 
    } 
    
    END_CODE();
  }

  /**
   * Apply a dslash
   *
   * Performs the operation
   *
   *  chi <-   (L + D + L^dag) . psi
   *
   * where
   *   L       is a lower triangular matrix
   *   D       is the real diagonal. (stored together in type TRIANG)
   *
   * Arguments:
   * \param chi     result                                      (Write)
   * \param psi     source                                      (Read)
   * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
   * \param cb      Checkerboard of OUTPUT std::vector               (Read) 
   */
  template<typename T, typename U>
  void ExpCloverTermT<T,U>::apply(T& chi, const T& psi, 
			    enum PlusMinus isign, int cb) const
  {
#ifndef QDP_IS_QDPJIT
    START_CODE();
    
    if ( Ns != 4 ) {
      QDPIO::cerr << __func__ << ": CloverTerm::apply requires Ns==4" << std::endl;
      QDP_abort(1);
    }

    QDPCloverEnv::ApplyArgs<T> arg = { chi,psi,tri,cb };
    int num_sites = rb[cb].siteTable().size();

    // The dispatch function is at the end of the file
    // ought to work for non-threaded targets too...
    dispatch_to_threads(num_sites, arg, QDPCloverEnv::applySiteLoop<T>);
    (*this).getFermBC().modifyF(chi, QDP::rb[cb]);

    END_CODE();
#endif
  }

  template<typename T, typename U>
  void ExpCloverTermT<T,U>::packForQUDA(multi1d<QUDAPackedClovSite<typename WordType<T>::Type_t> >& quda_array, int cb) const
    {
      typedef typename WordType<T>::Type_t REALT;
      int num_sites = rb[cb].siteTable().size();

      QDPCloverEnv::QUDAPackArgs<REALT> args = { cb, quda_array,tri };
      dispatch_to_threads(num_sites, args, QDPCloverEnv::qudaPackSiteLoop<REALT>);
      
    }  
///////////////////////////////////////////////////     sunpeng begin /////////////////////////////////////////////////////
 namespace ExpCloverEnv {

   template<typename REALT, typename T>
    struct ExpCloverMakeExpArg2 {
      const REALT& fac;
      const int N;
      const int cb;
      PrimitiveClovTriang < REALT >* tri;
      PrimitiveClovTriang < REALT >* tri_ori;
      std::vector<REALT> csv;
/*     
      T &chi15;
      T &chi14;
      T &chi13;
      T &chi12;
            
      T &chi21;
      T &chi22;
      T &chi23;
      T &chi24;     
*/
      T (&chi)[8];     
 
      const T &psi1;
      const T &psi2;
      
      
    };

   
  template<typename REALT, typename T, typename U>
  void applySite2( REALT *R0, T& chi, const T& psi  , const T& psi2,  PrimitiveClovTriang < REALT >* tri, int site) 
  {
#ifndef QDP_IS_QDPJIT
    START_CODE();

    if ( Ns != 4 )
    {
      QDPIO::cerr << __func__ << ": CloverTerm::applySite requires Ns==4" << std::endl;
      QDP_abort(1);
    }

    int n = 2*Nc;

    RComplex<REALT>* cchi = (RComplex<REALT>*)&(chi.elem(site).elem(0).elem(0));
    const RComplex<REALT>* ppsi = (const RComplex<REALT>*)&(psi.elem(site).elem(0).elem(0));
    const RComplex<REALT>* ppsi2 = (const RComplex<REALT>*)&(psi2.elem(site).elem(0).elem(0));

    cchi[ 0] = tri[site].diag[0][ 0]  * ppsi[ 0]  +  RScalar<REALT>(R0[0]) * ppsi2[ 0]
      +   conj(tri[site].offd[0][ 0]) * ppsi[ 1]
      +   conj(tri[site].offd[0][ 1]) * ppsi[ 2]
      +   conj(tri[site].offd[0][ 3]) * ppsi[ 3]
      +   conj(tri[site].offd[0][ 6]) * ppsi[ 4]
      +   conj(tri[site].offd[0][10]) * ppsi[ 5];
    
    cchi[ 1] = tri[site].diag[0][ 1]  * ppsi[ 1]  +  RScalar<REALT>(R0[0]) * ppsi2[ 1]
      +        tri[site].offd[0][ 0]  * ppsi[ 0]
      +   conj(tri[site].offd[0][ 2]) * ppsi[ 2]
      +   conj(tri[site].offd[0][ 4]) * ppsi[ 3]
      +   conj(tri[site].offd[0][ 7]) * ppsi[ 4]
      +   conj(tri[site].offd[0][11]) * ppsi[ 5];
    
    cchi[ 2] = tri[site].diag[0][ 2]  * ppsi[ 2]  +  RScalar<REALT>(R0[0]) * ppsi2[ 2]
      +        tri[site].offd[0][ 1]  * ppsi[ 0]
      +        tri[site].offd[0][ 2]  * ppsi[ 1]
      +   conj(tri[site].offd[0][ 5]) * ppsi[ 3]
      +   conj(tri[site].offd[0][ 8]) * ppsi[ 4]
      +   conj(tri[site].offd[0][12]) * ppsi[ 5];
    
    cchi[ 3] = tri[site].diag[0][ 3]  * ppsi[ 3]  +  RScalar<REALT>(R0[0]) * ppsi2[ 3]
      +        tri[site].offd[0][ 3]  * ppsi[ 0]
      +        tri[site].offd[0][ 4]  * ppsi[ 1]
      +        tri[site].offd[0][ 5]  * ppsi[ 2]
      +   conj(tri[site].offd[0][ 9]) * ppsi[ 4]
      +   conj(tri[site].offd[0][13]) * ppsi[ 5];
    
    cchi[ 4] = tri[site].diag[0][ 4]  * ppsi[ 4]  +  RScalar<REALT>(R0[0]) * ppsi2[ 4]
      +        tri[site].offd[0][ 6]  * ppsi[ 0]
      +        tri[site].offd[0][ 7]  * ppsi[ 1]
      +        tri[site].offd[0][ 8]  * ppsi[ 2]
      +        tri[site].offd[0][ 9]  * ppsi[ 3]
      +   conj(tri[site].offd[0][14]) * ppsi[ 5];
    
    cchi[ 5] = tri[site].diag[0][ 5]  * ppsi[ 5]  +  RScalar<REALT>(R0[0]) * ppsi2[ 5] 
      +        tri[site].offd[0][10]  * ppsi[ 0]
      +        tri[site].offd[0][11]  * ppsi[ 1]
      +        tri[site].offd[0][12]  * ppsi[ 2]
      +        tri[site].offd[0][13]  * ppsi[ 3]
      +        tri[site].offd[0][14]  * ppsi[ 4];
    
    cchi[ 6] = tri[site].diag[1][ 0]  * ppsi[ 6]  +  RScalar<REALT>(R0[1]) * ppsi2[ 6]
      +   conj(tri[site].offd[1][ 0]) * ppsi[ 7]
      +   conj(tri[site].offd[1][ 1]) * ppsi[ 8]
      +   conj(tri[site].offd[1][ 3]) * ppsi[ 9]
      +   conj(tri[site].offd[1][ 6]) * ppsi[10]
      +   conj(tri[site].offd[1][10]) * ppsi[11];
    
    cchi[ 7] = tri[site].diag[1][ 1]  * ppsi[ 7]  +  RScalar<REALT>(R0[1]) * ppsi2[ 7]
      +        tri[site].offd[1][ 0]  * ppsi[ 6]
      +   conj(tri[site].offd[1][ 2]) * ppsi[ 8]
      +   conj(tri[site].offd[1][ 4]) * ppsi[ 9]
      +   conj(tri[site].offd[1][ 7]) * ppsi[10]
      +   conj(tri[site].offd[1][11]) * ppsi[11];
    
    cchi[ 8] = tri[site].diag[1][ 2]  * ppsi[ 8]  +  RScalar<REALT>(R0[1]) * ppsi2[ 8]
      +        tri[site].offd[1][ 1]  * ppsi[ 6]
      +        tri[site].offd[1][ 2]  * ppsi[ 7]
      +   conj(tri[site].offd[1][ 5]) * ppsi[ 9]
      +   conj(tri[site].offd[1][ 8]) * ppsi[10]
      +   conj(tri[site].offd[1][12]) * ppsi[11];
    
    cchi[ 9] = tri[site].diag[1][ 3]  * ppsi[ 9]  +  RScalar<REALT>(R0[1]) * ppsi2[ 9]
      +        tri[site].offd[1][ 3]  * ppsi[ 6]
      +        tri[site].offd[1][ 4]  * ppsi[ 7]
      +        tri[site].offd[1][ 5]  * ppsi[ 8]
      +   conj(tri[site].offd[1][ 9]) * ppsi[10]
      +   conj(tri[site].offd[1][13]) * ppsi[11];
    
    cchi[10] = tri[site].diag[1][ 4]  * ppsi[10]  +  RScalar<REALT>(R0[1]) * ppsi2[10]
      +        tri[site].offd[1][ 6]  * ppsi[ 6]
      +        tri[site].offd[1][ 7]  * ppsi[ 7]
      +        tri[site].offd[1][ 8]  * ppsi[ 8]
      +        tri[site].offd[1][ 9]  * ppsi[ 9]
      +   conj(tri[site].offd[1][14]) * ppsi[11];
    
    cchi[11] = tri[site].diag[1][ 5]  * ppsi[11]  +  RScalar<REALT>(R0[1]) * ppsi2[11]
      +        tri[site].offd[1][10]  * ppsi[ 6]
      +        tri[site].offd[1][11]  * ppsi[ 7]
      +        tri[site].offd[1][12]  * ppsi[ 8]
      +        tri[site].offd[1][13]  * ppsi[ 9]
      +        tri[site].offd[1][14]  * ppsi[10];


    END_CODE();
#endif
  }


  
    template<typename REALT,typename T, typename U>
    inline 
    void MakeExp2(int lo, int hi, int myId, ExpCloverMakeExpArg2<REALT,T> *a)
    {
#ifndef QDP_IS_QDPJIT    
      const REALT& fac=a->fac;
      int N=a->N;
      int cb=a->cb;
//      PrimitiveClovTriang < REALT >* tri=a->tri;
      PrimitiveClovTriang < REALT >* tri_ori=a->tri_ori;
      PrimitiveClovTriang < REALT >* tri_tmp;             
      std::vector<REALT> &csv=a->csv;
/*
      T &chi12= a->chi12;
      T &chi13= a->chi13;
      T &chi14= a->chi14;
      T &chi15= a->chi15;
       
      T &chi21= a->chi21;
      T &chi22= a->chi22;
      T &chi23= a->chi23;
      T &chi24= a->chi24;
 */     
      T  chi_0 = zero;
      T  chi_tmp; 
      T (&chi)[8]=a->chi;
      
    
      const T &psi1= a->psi1;
      const T &psi2= a->psi2;


      REALT psv[10], q[12];
           
      PrimitiveClovTriang < REALT > A2,A3, tmp ,tmp2 ;
                  
      for(int ssite=lo; ssite < hi; ++ssite)  {
          int site = rb[cb].siteTable()[ssite];

          clv_cpoly<REALT>(tri_ori[site],A2,A3,psv);
          clv_cayley<REALT>(N,csv,psv,q);

          for(int k=0;k<12;k++)
                q[k]*=fac;

        //  1                                        chi   X  psi         c5*A^4 
        //                                                                chi[7]
        // 
        //  c4 + c5*A                                chi   X  psi          A^3
        //   chi[0]                                                       chi[6] 
        //
        //  c3 + c4*A + c5*A^2                       chi   X  psi          A^2  
        //   chi[1]                                                       chi[5]
        //
        //  c2 + c3*A + c4*A^2 + c5*A^3              chi   X  psi          A
        //   chi[2]                                                       chi[4]
        //
        //  c1 + c2*A + c3*A^2 + c4*A^3 +c5*A^4      chi   X  psi          1
        //   chi[3]


         applySite2<REALT, T, U>( q+10  ,chi_tmp, chi_0   , psi1, tri_ori, site);
         applySite2<REALT, T, U>( q+8   ,chi[0] , chi_tmp , psi1, tri_ori, site);
         applySite2<REALT, T, U>( q+6   ,chi[1] , chi[0]  , psi1, tri_ori, site);
         applySite2<REALT, T, U>( q+4   ,chi[2] , chi[1]  , psi1, tri_ori, site);
         applySite2<REALT, T, U>( q+2   ,chi[3] , chi[2]  , psi1, tri_ori, site);

         applySite2<REALT, T, U>( q     ,chi[4]  , psi2     , chi_0      , tri_ori, site); 
         applySite2<REALT, T, U>( q     ,chi[5]  , chi[4]   , chi_0      , tri_ori, site);
         applySite2<REALT, T, U>( q     ,chi[6]  , chi[5]   , chi_0      , tri_ori, site);
         applySite2<REALT, T, U>( q     ,chi_tmp , chi[6]   , chi_0      , tri_ori, site);         
         applySite2<REALT, T, U>( q+10  ,chi[7]  , chi_0    , chi_tmp    , tri_ori, site);
/*

          pauli_add_fac<REALT>(q+8,tri_ori[site],tri_tmp[site]); //    tri=c4+c5*A        
          applySite2<REALT, T, U>( chi[0], psi1, tri_tmp, site);           
           
          pauli_add_fac<REALT>(q+6,tri_ori[site],A2,tri_tmp[site]); // tri=c3+c4*A+c5*A^2
          applySite2<REALT, T, U>( chi[1], psi1, tri_tmp, site);


          pauli_add_fac<REALT>(q+8,tri_ori[site],tmp); // tmp=c4+c5*A
          pauli_prod<REALT>(A2,tmp,tri_tmp[site]);  // tri[site]=tmp*A2=A^2*(c4+c5*A)
          pauli_add_fac<REALT>(q+4,tri_ori[site],tmp2); // tmp2=c2+c3*A
          pauli_add<REALT>(tri_tmp[site],tmp2,tri_tmp[site]); //tri[site]=tri[site]+tmp2=c2+c3*A+c4*A^2+c5*A^3 
          applySite2<REALT, T, U>( chi[2], psi1, tri_tmp, site);


          pauli_add_fac<REALT>(q+2,tri_ori[site],A2,tmp2);  // tmp2=c1+c2*A+c3*A^2
          pauli_prod<REALT>(A3,tmp,tri_tmp[site]); // tri[site]=tmp*A3=A^3*(c4+c5*A)   
          pauli_add<REALT>(tri_tmp[site],tmp2,tri_tmp[site]); // tri[site]=c1+c2*A+c3*A^2+c4*A^3+c5*A^4
          applySite2<REALT, T, U>( chi[3], psi1, tri_tmp, site);



          pauli_prod<REALT>(A2,A2,tmp2);     //tmp2=A^4
          pauli_prod<REALT>(q+10,tmp2,tri_tmp[site]); //  tri[site]=c5*A^2      
          applySite2<REALT, T, U>( chi[7], psi2, tri_tmp, site);
  
          pauli_prod<REALT>(tri_ori[site],A2,tri_tmp[site]);    // tri[site]=A^3      
          applySite2<REALT, T, U>( chi[6], psi2, tri_tmp, site); 
        
          pauli_prod<REALT>(tri_ori[site],tri_ori[site],tri_tmp[site]); //tri[site]=A^2
          applySite2<REALT, T, U>( chi[5], psi2, tri_tmp, site);
           
          applySite2<REALT, T, U>( chi[4], psi2, tri_ori, site);    //tri[site]=A 
*/

#endif
      }/* End Site Loop */
    }/* End Function */
  }/* End Namespace */
 
  
  template<typename T, typename U>
  void ExpCloverTermT<T,U>::apply2( int N_, int s, RealT fac,  T (&chi)[8], const T& psi1, const T& psi2, 
			    enum PlusMinus isign, int cb) 
  {
#ifndef QDP_IS_QDPJIT
    START_CODE();
    
    if ( Ns != 4 ) {
      QDPIO::cerr << __func__ << ": CloverTerm::apply requires Ns==4" << std::endl;
      QDP_abort(1);
    }

    ExpCloverEnv::ExpCloverMakeExpArg2<REALT,T> arg = { fac.elem().elem().elem().elem(),N,cb,tri,tri_ori,csv[s], chi, psi1,psi2 };
    int num_sites = rb[cb].siteTable().size();

    // The dispatch function is at the end of the file
    // ought to work for non-threaded targets too...
    dispatch_to_threads(num_sites, arg, ExpCloverEnv::MakeExp2< REALT, T,  U  >);
    for ( int i ; i< 8 ; i++   )(*this).getFermBC().modifyF(chi[i], QDP::rb[cb]);

    END_CODE();
#endif
  }

  
  template<typename T, typename U>
  void ExpCloverTermT<T,U>::derivExpMultipole(multi1d<U>& ds_u, 
			     const multi1d<T>& chi, const multi1d<T>& psi, 
			     enum PlusMinus isign, int cb) const
  {
    START_CODE();


    // Do I still need to do this?
    if( ds_u.size() != Nd ) { 
      ds_u.resize(Nd);
    }

    ds_u = zero;

    // Get the links
    //const multi1d<U>& u = getU();


    // Now compute the insertions
    for(int mu=0; mu < Nd; mu++) {
      for(int nu = mu+1; nu < Nd; nu++) {
	
	// These will be appropriately overwritten - no need to zero them.
	// Contributions to mu links from mu-nu clover piece
	U ds_tmp_mu; 

	// -ve contribs  to the nu_links from the mu-nu clover piece 
	// -ve because of the exchange of gamma_mu gamma_nu <-> gamma_nu gamma_mu
	U ds_tmp_nu;

	// The weight for the terms
//	Real factor = (Real(-1)/Real(8))*getCloverCoeff(mu,nu);

	// Get gamma_mu gamma_nu psi -- no saving here, from storing shifts because
	// I now only do every mu, nu pair only once.

        if(cb<0||cb>=rb.numSubsets())
     {
        QDPIO::cerr << "ExpCloverTerm: error: cb = " << cb << "is not allowed." << std::endl;
        QDP_abort(1);
     }
 
     U s_xy_dag=zero;
 
     for ( int i=0; i < chi.size(); i++  )
     {
     T chi2[8];

     apply2(N,0,diag_mass, chi2, psi[i], chi[i], isign, cb);


	int mu_nu_index = (1 << mu) + (1 << nu); // 2^{mu} 2^{nui}
       
//	T ferm_tmp = Gamma(mu_nu_index)*psi;
	s_xy_dag += traceSpin( outerProduct(   Gamma(mu_nu_index)*psi[i]    ,chi2[7]));
        s_xy_dag += traceSpin( outerProduct(   Gamma(mu_nu_index)*chi2[0]   ,chi2[6]));
        s_xy_dag += traceSpin( outerProduct(   Gamma(mu_nu_index)*chi2[1]   ,chi2[5]));
        s_xy_dag += traceSpin( outerProduct(   Gamma(mu_nu_index)*chi2[2]   ,chi2[4]));      
        s_xy_dag += traceSpin( outerProduct(   Gamma(mu_nu_index)*chi2[3]   ,chi[i] ));
     } 

//	s_xy_dag *= Real(factor);


        // Compute contributions
	deriv_loops(mu, nu, cb, ds_tmp_mu, ds_tmp_nu, s_xy_dag);

	// Accumulate them
	ds_u[mu] += ds_tmp_mu;
	ds_u[nu] -= ds_tmp_nu;


      }
    }


    // Clear out the deriv on any fixed links
    (*this).getFermBC().zero(ds_u);
    END_CODE();
  }



  template<typename T, typename U>
  void ExpCloverTermT<T,U>::derivExp(multi1d<U>& ds_u, 
			     const T& chi, const T& psi, 
			     enum PlusMinus isign, int cb) const
  {
    START_CODE();

    multi1d<T> chi_v(1);
    multi1d<T> psi_v(1);
    chi_v[0]=chi;
    psi_v[0]=psi;

     derivExpMultipole(  ds_u,  chi_v,  psi_v, isign,  cb   );


/*
    // Do I still need to do this?
    if( ds_u.size() != Nd ) { 
      ds_u.resize(Nd);
    }

    ds_u = zero;

    // Get the links
    //const multi1d<U>& u = getU();


    // Now compute the insertions
    for(int mu=0; mu < Nd; mu++) {
      for(int nu = mu+1; nu < Nd; nu++) {
	
	// These will be appropriately overwritten - no need to zero them.
	// Contributions to mu links from mu-nu clover piece
	U ds_tmp_mu; 

	// -ve contribs  to the nu_links from the mu-nu clover piece 
	// -ve because of the exchange of gamma_mu gamma_nu <-> gamma_nu gamma_mu
	U ds_tmp_nu;

	// The weight for the terms
//	Real factor = (Real(-1)/Real(8))*getCloverCoeff(mu,nu);

	// Get gamma_mu gamma_nu psi -- no saving here, from storing shifts because
	// I now only do every mu, nu pair only once.

        if(cb<0||cb>=rb.numSubsets())
     {
        QDPIO::cerr << "ExpCloverTerm: error: cb = " << cb << "is not allowed." << std::endl;
        QDP_abort(1);
     }
     T chi2[8];

     apply2(N,0,diag_mass, chi2, psi, chi, isign, cb);


	int mu_nu_index = (1 << mu) + (1 << nu); // 2^{mu} 2^{nui}
        U s_xy_dag=zero;
//	T ferm_tmp = Gamma(mu_nu_index)*psi;
	s_xy_dag += traceSpin( outerProduct(   Gamma(mu_nu_index)*psi       ,chi2[7]));
        s_xy_dag += traceSpin( outerProduct(   Gamma(mu_nu_index)*chi2[0]   ,chi2[6]));
        s_xy_dag += traceSpin( outerProduct(   Gamma(mu_nu_index)*chi2[1]   ,chi2[5]));
        s_xy_dag += traceSpin( outerProduct(   Gamma(mu_nu_index)*chi2[2]   ,chi2[4]));      
        s_xy_dag += traceSpin( outerProduct(   Gamma(mu_nu_index)*chi2[3]   ,chi    ));


//	s_xy_dag *= Real(factor);


        // Compute contributions
	derivExp_loops(mu, nu, cb, ds_tmp_mu, ds_tmp_nu, s_xy_dag);

	// Accumulate them
	ds_u[mu] += ds_tmp_mu;
	ds_u[nu] -= ds_tmp_nu;


      }
    }

*/
    // Clear out the deriv on any fixed links
    (*this).getFermBC().zero(ds_u);
    END_CODE();
  }


  template<typename T, typename U>
  void ExpCloverTermT<T,U>::derivExpTrLn(multi1d<U>& ds_u, 
				 enum PlusMinus isign, int cb) const
  {
    START_CODE();
    
    // Do I still need to do this?
    if( ds_u.size() != Nd ) { 
      ds_u.resize(Nd);
    }
    
    ds_u = zero;
/*
    for(int mu=0; mu < Nd; mu++) {
      for(int nu = mu+1; nu < Nd; nu++) { 

	  // Index 
	  int mu_nu_index = (1 << mu) + (1 << nu); // 2^{mu} 2^{nu}

	  // The actual coefficient factor
//	  Real factor = Real(-1)*getCloverCoeff(mu,nu)/Real(8);
	  
	  U sigma_XY_dag=zero;

	  // Get  weight*Tr_spin gamma_mu gamma_nu A^{-1} piece
	  triacntr(sigma_XY_dag, mu_nu_index, cb);
//	  sigma_XY_dag[rb[cb]] *= factor;

	  // These will be overwritten so no need to initialize to zero
	  U ds_tmp_mu;
	  U ds_tmp_nu;

	  // Get contributions from the loops and insersions
	  derivExp_loops(mu, nu, cb, ds_tmp_mu, ds_tmp_nu, sigma_XY_dag);

	  // Accumulate
	  ds_u[mu] += ds_tmp_mu;
	  // -ve weight for nu from gamma_mu gamma_nu -> gamma_nu gamma_mu
	  // commutation.
	  ds_u[nu] -= ds_tmp_nu;

      } // End loop over nu

    } // end of loop over mu
    
*/
    // Not sure this is needed here, but will be sure
    (*this).getFermBC().zero(ds_u); // ???   
    END_CODE();

  }


  template<typename T, typename U>
  void ExpCloverTermT<T,U>::derivExp(multi1d<U>& ds_u, 
			     const T& chi, const T& psi, 
			     enum PlusMinus isign) const
  {
    START_CODE();

    // base deriv resizes.
    // Even even checkerboard
    deriv(ds_u, chi, psi, isign,0);
    
    // Odd Odd checkerboard
    multi1d<U> ds_tmp;
    deriv(ds_tmp, chi, psi, isign,1);
    
    ds_u += ds_tmp;
    
    END_CODE();
  }

  template<typename T, typename U>
  void ExpCloverTermT<T,U>::derivExpMultipole(multi1d<U>& ds_u, 
			     const multi1d<T>& chi, const multi1d<T>& psi, 
			     enum PlusMinus isign) const
  {
    START_CODE();

    // base deriv resizes.
    // Even even checkerboard
    derivMultipole(ds_u, chi, psi, isign,0);
    
    // Odd Odd checkerboard
    multi1d<U> ds_tmp;
    derivMultipole(ds_tmp, chi, psi, isign,1);
    
    ds_u += ds_tmp;
    
    END_CODE();
  }


////////////////////////////////////////////////////////////////////////  sunpeng end  ///////////////////////////////////////////////////////////////////////////


  typedef ExpCloverTermT<LatticeFermion, LatticeColorMatrix> ExpCloverTerm;
  typedef ExpCloverTermT<LatticeFermionF, LatticeColorMatrixF> ExpCloverTermF;
  typedef ExpCloverTermT<LatticeFermionD, LatticeColorMatrixD> ExpCloverTermD;
} // End Namespace Chroma


#endif

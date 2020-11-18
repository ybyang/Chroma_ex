/*! \file
 *  \a test code for R_xi gauge (and Landau) gauge fixing by ChunJiang Shi
 */

#include "chromabase.h"
#include "util/gauge/reunit.h"
#include "util/gauge/su2extract.h"
#include "util/gauge/sunfill.h"

#include "meas_ex/Xigauge.h"

namespace Chroma {
//using namespace QDP;

/*

// Hooks for various things
#if defined(QDP_DEBUG_MEMORY)
#define START_CODE() QDP::Allocator::theQDPAllocator::Instance().pushFunc(__func__, __LINE__)
#define END_CODE()   QDP::Allocator::theQDPAllocator::Instance().popFunc()

#else
#define START_CODE() QDP_PUSH_PROFILE(QDP::getProfileLevel())
#define END_CODE()   QDP_POP_PROFILE()

#endif
*/

/********************** HACK ******************************/
// Primitive way for now to indicate the time direction

static int tDir() {return Nd-1;}
static Real xi_0() {return 1.0;}
static bool anisoP() {return false;}

/******************** END HACK ***************************/

#ifndef ONE_DEVIDE_SQRT3
#define ONE_DEVIDE_SQRT3 REAL(0.57735026918962576450914878050195744)
#endif

inline void fill_lambda(PColorMatrix< RComplex<REAL>, Nc>& lam , 
                                                    REAL r0,
                                                    REAL r1,
                                                    REAL r2,
                                                    REAL r3,
                                                    REAL r4,
                                                    REAL r5,
                                                    REAL r6,
                                                    REAL r7){
  lam.elem(0,0).real() = r2 + r7*ONE_DEVIDE_SQRT3;
  lam.elem(0,0).imag() = 0;

  lam.elem(1,1).real() = -r2 + r7*ONE_DEVIDE_SQRT3;
  lam.elem(1,1).imag() = 0;

  lam.elem(2,2).real() = -2 * r7 * ONE_DEVIDE_SQRT3;
  lam.elem(2,2).imag() = 0;

  lam.elem(1,0).real() = r0;
  lam.elem(1,0).imag() = r1;

  lam.elem(0,1).real() = r0;
  lam.elem(0,1).imag() = -r1;

  lam.elem(2,0).real() = r3;
  lam.elem(2,0).imag() = r4;

  lam.elem(0,2).real() = r3;
  lam.elem(0,2).imag() = -r4;

  lam.elem(1,2).real() = r5;
  lam.elem(1,2).imag() = -r6;

  lam.elem(2,1).real() = r5;
  lam.elem(2,1).imag() = r6;
}


//! dest	= gaussian	 under a subset
void  lambdafill(LatticeColorMatrix& d, const Subset& s, const multi1d<LatticeReal>& r)
{
	const int *tab = s.siteTable().slice();
	const int nodeSites = s.numSiteTable();

#pragma omp parallel for
	for(int j=0; j < nodeSites; ++j)
	{
		int i = tab[j];
		fill_lambda(d.elem(i).elem(), r[0].elem(i).elem().elem().elem(),
                           r[1].elem(i).elem().elem().elem(),
                           r[2].elem(i).elem().elem().elem(),
                           r[3].elem(i).elem().elem().elem(),
                           r[4].elem(i).elem().elem().elem(),
                           r[5].elem(i).elem().elem().elem(), 
                           r[6].elem(i).elem().elem().elem(),
                           r[7].elem(i).elem().elem().elem());

    /*
    fill_lambda(d.elem(tab[nodeSites/2-1-j]).elem(), 
                           -r[0].elem(i).elem().elem().elem(),
                           -r[1].elem(i).elem().elem().elem(),
                           -r[2].elem(i).elem().elem().elem(),
                           -r[3].elem(i).elem().elem().elem(),
                           -r[4].elem(i).elem().elem().elem(),
                           -r[5].elem(i).elem().elem().elem(), 
                           -r[6].elem(i).elem().elem().elem(),
                           -r[7].elem(i).elem().elem().elem());
                           */
	}
}


void getlambda(LatticeColorMatrix& lambda, const Real xi){
        multi1d<LatticeReal> randnum(8); 
        multi1d<LatticeReal> mean(8); 
        
  QDPIO::cout << "start generate Lambda ..." << std::endl;

  const Real xi_r = sqrt(xi);  
  LatticeReal sqrt_xi = xi_r;
  for(int i = 0; i < 8 ; ++i){
    gaussian(randnum[i]);
    mean[i] = sum(randnum[i])/Layout::vol();
    //using a linear transformation to keep sum 0 with the same exception and variance
    randnum[i] = sqrt(Layout::vol()/(Layout::vol()-1)) * (randnum[i] - mean[i]);
    randnum[i] *= sqrt_xi/2.0;
  }

  QDPIO::cout << "********** randnum[0].elem(0).elem().elem().elem() = " << randnum[0].elem(0).elem().elem().elem() << std::endl;
  QDPIO::cout << "********** randnum[1].elem(0).elem().elem().elem() = " << randnum[1].elem(0).elem().elem().elem() << std::endl;
  
  lambdafill(lambda, all, randnum);
  QDPIO::cout << "generate Lambda done..." << std::endl;

  //Debug:------------------------------- 
  /*     
        printf("Debug: lambda in lattice ilat = 7:\n");
        //LatticeColorMatrix debug = g * lambda_input; 
        for(int i = 0; i < Nc; i++){
          for (int j = 0; j < Nc; j++){
            std::cout << lambda.elem(7).elem().elem(i,j);
          }
          printf("\n");
        }
        printf("-----------------------------\n");
  */
  //end debug-------------------------------------
}




//! Perform a single gauge fixing iteration
/*!
 * \ingroup gfix
 *
 * Performs one gauge fixing 'iteration', one checkerboard and SU(2)
 * subgroup only, for gauge fixing to Coulomb gauge in slices perpendicular
 * to the direction "j_decay".
 *
 * \param g          Current (global) gauge transformation matrices ( Modify )
 * \param lambda_input	 gauge condition for gf (read)
 * \param u          original gauge field ( Read )
 * \param j_decay    direction perpendicular to slices to be gauge fixed ( Read )
 * \param su2_index  SU(2) subgroup index ( Read )
 * \param cb         checkerboard index ( Read )
 * \param ordo       use overrelaxation or not ( Read )
 * \param orpara     overrelaxation parameter ( Read ) 
 */



void Rxigrelax(LatticeColorMatrix& g,   //
      const LatticeColorMatrix& lambda_input,
	    const multi1d<LatticeColorMatrix>& u, 
	    int j_decay, int su2_index, int cb, bool ordo,
	    const Real& orpara)
{
  LatticeColorMatrix v;
  multi1d<LatticeReal> r(4);
  multi1d<LatticeReal> a(4);
  
  START_CODE();
      
  // Rotate the matrices using the current gauge rotation
  multi1d<LatticeColorMatrix> ug(Nd);

  for(int mu = 0; mu < Nd; ++mu)
    ug[mu] = g * u[mu] * shift(adj(g), FORWARD, mu) ;
    //ug[mu] = u[mu] ; //wrong !!!
  //XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  //ug[mu] = g * u[mu] * shift(adj(g), FORWARD, mu) - RxiGF:: muiltiimag(g * lambda_input) ;
  //NONONO!!! There's the problem. Cannot add RxiGF:: muiltiimag(g * lambda_input) here.

  /* Gather the Nd negative links attached to a site: */
  /* u_tmp(x,mu) = U(x-mu,mu) */
  multi1d<LatticeColorMatrix> u_neg(Nd);
  for(int mu=0; mu<Nd; ++mu)
    u_neg[mu][rb[cb]] = shift(ug[mu], BACKWARD, mu);
  
  /* Sum links to be gauge transformed on site x not in the direction */
  /* j_decay into matrix V: */
  v = 0;
  if (tDir() != j_decay)  //Landau gauge, for tDir 
  {
    v[rb[cb]] = ug[tDir()] + adj(u_neg[tDir()]);

    if (anisoP())
      v[rb[cb]] *= pow(xi_0(), 2);
  }

  for(int mu = 0; mu < Nd; ++mu)
  {
    if (mu != j_decay && mu != tDir())
      v[rb[cb]] += ug[mu] + adj(u_neg[mu]);
  }

  /*
  //Debug:-------------------------------
  printf("Debug: before -= RxiGF::extractlambda in lattice ilat = 7:\n");
  //LatticeColorMatrix debug = g * lambda_input; 
  for(int i = 0; i < Nc; i++){
    for (int j = 0; j < Nc; j++){
      std::cout << v.elem(7).elem().elem(i,j);
    }
    printf("\n");
  }
  printf("-----------------------------\n");
  //-------------------------------------
  */
  
  //LatticeComplex factor = cmplx( 0, );
  LatticeReal t1 = 0, t2 = 1.0;
  v[rb[cb]] -= cmplx(t1, t2) * lambda_input;
  //XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  //v[rb[cb]] -= cmplx(t1, t2) * (adj(g) * lambda_input); totally wrong!!
  //v[rb[cb]] -= cmplx(t1, t2) * (g * lambda_input);
  //NOTE:: borrow ug[0]

  /*
  //Debug:-------------------------------
  printf("Debug: after -= RxiGF::extractlambda in lattice ilat = 7 , cb = %d , rb[cb] = %d :\n", cb , rb[cb]);
  //LatticeColorMatrix debug = g * lambda_input; 
  for(int i = 0; i < Nc; i++){
    for (int j = 0; j < Nc; j++){
      std::cout << v.elem(7).elem().elem(i,j);
    }
    printf("\n");
  }
  printf("-----------------------------\n");
  //-------------------------------------
  */
  
  



  if (Nc > 1)
  {
                                    
    /* Extract components r_k proportional to SU(2) submatrix su2_index */
    /* from the SU(Nc) matrix V. The SU(2) matrix is parametrized in the */
    /* sigma matrix basis. */
    
    Chroma::su2Extract(r, v, su2_index, rb[cb]);

    //Chroma::su2Extract(r2, lambda, su2_index, rb[cb]);

    //r[0] = 
  
    /*
     * Now project onto SU(2)
     */
    LatticeReal r_l;
    r_l[rb[cb]] = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2] + r[3]*r[3]);
  
    // Normalize
    LatticeBoolean lbtmp;
    lbtmp[rb[cb]] = r_l > fuzz;
    //QDPIO::cout << "fuzz = " << fuzz <<std::endl;
    LatticeReal lftmp;
    lftmp[rb[cb]] = 1.0 / where(lbtmp, r_l, LatticeReal(1));

    // Fill   (r[0]/r_l, -r[1]/r_l, -r[2]/r_l, -r[3]/r_l) for r_l > fuzz
    //  and   (1,0,0,0)  for sites with r_l < fuzz
    multi1d<LatticeReal> a(4);
    a[0][rb[cb]] = where(lbtmp, r[0] * lftmp, LatticeReal(1));
    a[1][rb[cb]] = where(lbtmp, -(r[1] * lftmp), LatticeReal(0));
    a[2][rb[cb]] = where(lbtmp, -(r[2] * lftmp), LatticeReal(0));
    a[3][rb[cb]] = where(lbtmp, -(r[3] * lftmp), LatticeReal(0));
      
    /* Now do the overrelaxation, if desired */
    if( ordo )
    {
      /* get angle */
      LatticeReal theta_old;
      theta_old[rb[cb]] = acos(a[0]);

      /* old sin */
      LatticeReal oldsin;
      oldsin[rb[cb]] = sin(theta_old);

      /* overrelax, i.e. multiply by the angle */
      LatticeReal theta_new;
      theta_new[rb[cb]] = theta_old * orpara;

      /* compute sin(new)/sin(old) */
      /* set the ratio to 0, if sin(old) < FUZZ */
      lftmp[rb[cb]] = where(oldsin > fuzz, sin(theta_new) / oldsin, LatticeReal(0));

      /* get the new cos = a[0] */
      a[0][rb[cb]] = cos(theta_new);

      /* get the new a_k, k = 1, 2, 3 */
      a[1][rb[cb]] *= lftmp;
      a[2][rb[cb]] *= lftmp;
      a[3][rb[cb]] *= lftmp;
    }
  
          
    /* Now fill the SU(Nc) matrix V with the SU(2) submatrix 'su2_index' */
    /* paramtrized by a_k in the sigma matrix basis. */
    sunFill(v, a, su2_index, rb[cb]);

  }
  else		/* Nc = 1 */
  {
    r[0][rb[cb]] = real(peekColor(v, 0, 0));
    r[1][rb[cb]] = imag(peekColor(v, 0, 0));
    LatticeReal r_l;
    r_l[rb[cb]] = sqrt(r[0] * r[0] + r[1] * r[1]);

    // Normalize
    LatticeBoolean lbtmp;
    lbtmp[rb[cb]] = r_l > fuzz;
    LatticeReal lftmp;
    lftmp[rb[cb]] = 1.0 / where(lbtmp, r_l, LatticeReal(1));

    // Fill   (r[0]/r_l, -r[1]/r_l, -r[2]/r_l, -r[3]/r_l) for r_l > fuzz
    //  and   (1,0,0,0)  for sites with r_l < fuzz
    multi1d<LatticeReal> a(4);
    a[0][rb[cb]] = where(lbtmp, r[0] * lftmp, LatticeReal(1));
    a[1][rb[cb]] = where(lbtmp, -(r[1] * lftmp), LatticeReal(0));
      
    /* Now do the overrelaxation, if desired */
    if( ordo )
    {
      Real pi = 0.5 * twopi;

      /* get angle */
      LatticeReal theta;
      theta[rb[cb]] = acos(a[0]) + pi;     // Do I really want to add pi ??? This was in szin

      /* overrelax, i.e. multiply by the angle */
      theta[rb[cb]] *= orpara;

      /* get the new cos = a[0] */
      a[0][rb[cb]] = cos(theta);

      /* get the new sin */
      a[1][rb[cb]] = sin(theta);
    }

    pokeColor(v, cmplx(a[0],a[1]), 0, 0);
  }
    
  // Now do the gauge transformation with matrix V:
  // Multiply into the global "g" field
  LatticeColorMatrix u_tmp;
  u_tmp[rb[cb]] = v * g;
  g[rb[cb]] = u_tmp;
  
  END_CODE();
}

    void getdelta(LatticeColorMatrix& delta, const multi1d<LatticeColorMatrix>& u, const LatticeColorMatrix& lambda){
        delta = 0;
        LatticeColorMatrix delta_tmp; 
        for (int mu = 0; mu < Nd; mu++){
            delta_tmp = u[mu] - adj(u[mu]);
	    LatticeComplex tr=trace(delta_tmp)*1.0/3;
	    delta_tmp-=tr;
            delta = delta + delta_tmp - shift(delta_tmp, BACKWARD, mu);
        }

        //times -i/2
        LatticeReal t1 = 0, t2 = -0.5;
        delta *= cmplx(t1, t2); 
         
        /*
        //Debug:-------------------------------
            printf("Debug: no lambda delta in lattice ilat = 7:\n");
            //LatticeColorMatrix debug = g * lambda_input; 
            for(int i = 0; i < Nc; i++){
                for (int j = 0; j < Nc; j++){
                std::cout << delta.elem(7).elem().elem(i,j);
                }
                printf("\n");
            }
            printf("-----------------------------\n");
            //-------------------------------------
        */

        delta = delta - lambda;                  
    }


    Double convercriterion(const multi1d<LatticeColorMatrix>& u, const LatticeColorMatrix& lambda){
        LatticeColorMatrix delta;
        getdelta(delta , u , lambda);
        Double ret;
        ret = sum(real(trace(delta * adj(delta))) );
        ret /= Double(Layout::vol()*Nc);
        
        //debug: imag part == 0 , correct! So don't warry about real() here.
        /*
        Double impart;
        impart = sum(imag(trace(delta * adj(delta))) ) / Double(Layout::vol()*Nc);
        QDPIO:: cout << "imag part = " << impart << std::endl; 
        */
        /*
        if(toBool(ret < 1e-8)){
          //Debug:-------------------------------
            printf("Debug: the convergence delta in lattice ilat = 7:\n");
            //LatticeColorMatrix debug = g * lambda_input; 
            for(int i = 0; i < Nc; i++){
                for (int j = 0; j < Nc; j++){
                std::cout << delta.elem(7).elem().elem(i,j);
                }
                printf("\n");
            }
            printf("-----------------------------\n");
            //-------------------------------------
        }
        */
        return ret;
    }



//! a test code for R_xi gauge (and Landau) gauge fixing by ShiChunJiang 
/*!
 * \ingroup gfix
 *
 * Driver for gauge fixing to Coulomb gauge in slices perpendicular
 * to the direction "j_decay".
 * If j_decay >= Nd: fix to Landau gauge.
 * Note: as written this works only for SU(2) and SU(3)!

 * \param u        (gauge fixed) gauge field ( Modify )
 * \param lambda   the gauge condition Lambda(x) ( Read )
 * \param g        Gauge transformation matrices (Write)
 * \param n_gf     number of gauge fixing iterations ( Write )
 * \param j_decay  direction perpendicular to slices to be gauge fixed ( Read )    \\????
 * \param GFAccu   desired accuracy for gauge fixing ( Read )
 * \param GFMax    maximal number of gauge fixing iterations ( Read )
 * \param OrDo     use overrelaxation or not ( Read )
 * \param OrPara   overrelaxation parameter ( Read )
 */

void RxiGauge(multi1d<LatticeColorMatrix>& u, 
         const LatticeColorMatrix& lambda,
	       LatticeColorMatrix& g,
	       int& n_gf, 
	       int j_decay, const Real& GFAccu, int GFMax, 
	       bool OrDo, const Real& OrPara)
{
  Double tgfold;
  Double tgfnew;
  Double tgf_t;
  Double tgf_s;
  Double norm;
  int num_sdir;
  bool tdirp;

  //MY HACK
  bool MyOrDo = OrDo;
  Real MyOrPara = OrPara;
  

  //LatticeColorMatrix total_g = 1;

  Double xigftheta = 100;//xigf theta init = 100

  /*for modify u_tmp[mu]*/  
  multi1d<LatticeColorMatrix> u_iter(Nd);
  LatticeColorMatrix ltmp;
  REAL rtmp;

  START_CODE();

  //Set xi parameter, 
  Real xi_sq = pow(xi_0(),2);

  //judge time and space number, get norm
  if( j_decay >= 0 && j_decay < Nd )
  {
    if( tDir() >= 0 && tDir() != j_decay ) //tDir means time direction
    {
      num_sdir = Nd - 2;                   //num_sdir means number of space direction
      tdirp = true;                        //time direction perpendicular
      norm = Double(Layout::vol()*Nc) * (Double(num_sdir)+Double(xi_sq));
    }
    else
    {
      num_sdir = Nd - 1;                   // I gauss Couloum gauge
      tdirp = false;
      norm = Double(Layout::vol()*Nc*num_sdir);
    }
  }
  else
  {
    if( tDir() >= 0 && tDir() < Nd )
    {
      num_sdir = Nd - 1;
      tdirp = true;
      norm = Double(Layout::vol()*Nc) * (Double(num_sdir)+Double(xi_sq));
    }
    else
    {
      num_sdir = Nd;
      tdirp = false;
      norm = Double(Layout::vol()*Nc*num_sdir);
    }
  }

      
  /* Compute initial gauge fixing term: sum(trace(U_spacelike)); */
  tgf_t = 0;
  tgf_s = 0;
  for(int mu=0; mu<Nd; ++mu)
    //XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    if( mu != j_decay )
    {
      Double tgf_tmp = sum(real(trace(u[mu])));

      if( mu != tDir() )
	tgf_s += tgf_tmp;
      else
	tgf_t += tgf_tmp;
    }

  if( tdirp )
  {
    tgfold = (xi_sq*tgf_t+tgf_s )/norm;
    tgf_s = tgf_s/(Double(Layout::vol()*Nc*num_sdir));
    tgf_t = tgf_t/(Double(Layout::vol()*Nc));
  }
  else
  {
    tgf_s = tgf_s/(Double(Layout::vol()*Nc*num_sdir));
    tgfold = tgf_s;
  }

  // Rxi gauge fixing convergence criterion
  //xigftheta = convercriterion(u, lambda);
  
  // Gauge transf. matrices always start from identity
  g = 1; 

  /* Gauge fix until converged or too many iterations */
  n_gf = 0;
  bool wrswitch = true;    /* switch for writing of gauge fixing term */
  Double conver = 1;        /*old Coulomb gauge fixing convergence criterion */

  /*record for debug time */
  time_t t_gfstart, t_iteration;
  time(&t_gfstart);

  while( toBool(xigftheta > GFAccu)  &&  n_gf < GFMax )
  {
    n_gf = n_gf + 1;
    if( GFMax - n_gf < 11 ) 
      wrswitch = true;

    /*
    if(toBool(xigftheta < 1e-6) && OrDo){
      MyOrPara = 1.7;
      MyOrDo = true;
    }
    */
    


    /* Loop over checkerboards for gauge fixing */
    int N_cb = 2;
    for(int cb=0; cb<N_cb; ++cb)
    {
      if (Nc > 1)
      {
        /* Loop over SU(2) subgroup index */
        for(int su2_index=0; su2_index < Nc*(Nc-1)/2; ++su2_index)
        {
          /* Now do a gauge fixing relaxation step */
          //Rxigrelax(g, lambda, u, j_decay, su2_index, cb, OrDo, OrPara);
          Rxigrelax(g, lambda, u, j_decay, su2_index, cb, MyOrDo, MyOrPara); //HACK CHROMA
          //reunit(g);
        }   /* end su2_index loop */
      }
      else // Nc==1 case 
      {
        int su2_index = -1;
        /* Now do a gauge fixing relaxation step */
        //Rxigrelax(g, lambda, u, j_decay, su2_index, cb, OrDo, OrPara);
        Rxigrelax(g , lambda, u , j_decay, su2_index, cb, OrDo, OrPara); 
      }
    }     /* end cb loop */

    /* Reunitarize */
    reunit(g);
    /* now we get new iteration g(x) */


  //Debug:-------------------------------
  if (toBool(n_gf % 1000 ==0))
  {
    printf("Debug: the ColorMatrix (g) in lattice ilat = 7:\n");
    //LatticeColorMatrix debug = g * lambda_input; 
    for(int i = 0; i < Nc; i++){
      for (int j = 0; j < Nc; j++){
        std::cout << g.elem(7).elem().elem(i,j);
      }
      printf("\n");
    }
    printf("-----------------------------\n");
  }
  //-------------------------------------


  //???? question ????
  //total_g *= g;
  //total_g = g * total_g;

  //Debug:-------------------------------
  /*
  if (toBool(n_gf % 1000 == 0))
  {
    printf("Debug: the total rot total_g in lattice ilat = 7:\n");
    //LatticeColorMatrix debug = g * lambda_input; 
    for(int i = 0; i < Nc; i++){
      for (int j = 0; j < Nc; j++){
        std::cout << total_g.elem(7).elem().elem(i,j);
      }
      printf("\n");
    }
    printf("-----------------------------\n");
  }
  */
  //-------------------------------------

    /* Compute new gauge fixing term: sum(trace(U_spacelike)): */
    tgf_t = 0;
    tgf_s = 0;
    for(int mu=0; mu<Nd; ++mu)
      if( mu != j_decay )
      {
	Double tgf_tmp = sum(real(trace(g * u[mu] * shift(adj(g), FORWARD, mu))) );  

	if( mu != tDir() )
	  tgf_s += tgf_tmp;
	else
	  tgf_t += tgf_tmp;
      }

    if( tdirp )
    {
      tgfnew = (xi_sq*tgf_t+tgf_s + sum(imag(trace(g * lambda))) )/norm;
      //RxiGF:modified XXXXXXXXXXXXXXXXXXXXXXXXXX
      tgf_s = tgf_s/(Double(Layout::vol()*Nc*num_sdir));
      tgf_t = tgf_t/(Double(Layout::vol()*Nc));
    }
    else
    {
      tgf_s = tgf_s/(Double(Layout::vol()*Nc*num_sdir));
      /*Debug : Rxi gauge should not do this.*/
      printf("Error: tdirp == false !!! Rxi gauge should not do this \n");
      tgfnew = tgf_s;
    }

    //debug:
    //QDPIO::cout << "normed trace of(i*g*lambda)= " << (sum(imag(trace(g * lambda))) )/norm<<std::endl;
    
    /*Old convergence criterion, do not use it in RxiGf*/
    /* Normalized convergence criterion: */
    //conver = fabs((tgfnew - tgfold) / tgfnew);
    conver = fabs(tgfnew - tgfold);


  /*
  //Debug:-------------------------------
  printf("Debug: after rot, u[0] in lattice ilat = 7:\n");
  LatticeColorMatrix debug2 = u[0]; 
  for(int i = 0; i < Nc; i++){
    for (int j = 0; j < Nc; j++){
      std::cout << debug2.elem(7).elem().elem(i,j);
    }
    printf("\n");
  }
  printf("-----------------------------\n");
  //End debug-------------------------------------
  */

  // Rxi gauge fixing convergence criterion
  //all begin with xigf
  for(int mu = 0; mu < Nd; ++mu)
  {
    LatticeColorMatrix u_tmp1 = g * u[mu];
    u_iter[mu] = u_tmp1 * shift(adj(g), FORWARD, mu);
  }
  xigftheta = convercriterion(u_iter, lambda);

    time(&t_iteration);

    if( wrswitch && toBool(n_gf % 100 == 0)) 
    printf("\n");
      QDPIO::cout << "R_xiGAUGE: iter= " << n_gf
      << "  OrDo= " << MyOrDo 
      << "  OrPara= " << MyOrPara
      << "  theta= " << xigftheta 
      << "  delta= " << tgfnew - tgfold
      << "  time(sec)= " << (t_iteration-t_gfstart)
		  << "  tgfold= " << tgfold 
		  << "  tgfnew= " << tgfnew
		  << "  tgf_s= " << tgf_s 
		  << "  tgf_t= " << tgf_t << std::endl;

    
    tgfold = tgfnew;
  }       /* end while loop */

      
  if( wrswitch && toBool(n_gf % 500 == 0)) 
      QDPIO::cout << "R_xiGAUGE: iter= " << n_gf
      << "  theta= " << xigftheta  
      << "  OrDo= " << MyOrDo 
      << "  delta= " << tgfnew - tgfold
      << "  time(sec)= " << (t_iteration-t_gfstart)
		  << "  tgfold= " << tgfold 
		  << "  tgfnew= " << tgfnew
		  << "  tgf_s= " << tgf_s 
		  << "  tgf_t= " << tgf_t << std::endl;


  // Finally, gauge rotate the original matrices and overwrite them
    //！！！！！！！！！！！！！！！！！！！！！BUG！！
  for(int mu = 0; mu < Nd; ++mu)
  {
    LatticeColorMatrix u_tmp2 = g * u[mu];
    u[mu] = u_tmp2 * shift(adj(g), FORWARD, mu);
  }
    
#if 0
  /*+ debugging */
  XMLBufferWriter xml_out;
  push(xml_out,"Final_trace_max_in_R_xiGauge");
  write(xml_out, "j_decay", j_decay);
  write(xml_out, "t_dir", tDir());
  write(xml_out, "n_gf",n_gf);
  write(xml_out, "tgfold", tgfold);
  write(xml_out, "tgf_s", tgf_s);
  write(xml_out, "tgf_t", tgf_t);
  pop(xml_out);
#endif

  END_CODE();
}


} // Namespace Chroma

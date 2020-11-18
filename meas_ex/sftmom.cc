//  $Id: sftmom.cc,v 3.11 2015-05-24 11:11 sdcohen $
//  $Log: sftmom.cc,v $
//  Revision 3.11 2015-05-24 11:11 sdcohen
//  Added capability to specify a list of timeslices such that the source is
//  duplicated for all members of the list
//
//  Revision 3.10 2014-04-07 00:13 sdcohen
//  Added capability to specify a list of mom shells rather than merely the
//  radius of a ball in momentum space or a list of individual momenta; this
//  version includes the option to specify origin and shift the center-of-momenta
//
//  Revision 3.9  2009-02-17 16:33:38  edwards
//  Added a function that will return the multiplicity of a particular mom id.
//  This is number of equiv. canonical mom. Only non-zero when mom averaging
//  is enabled. Can use this to undo the averaging.
//
//  Revision 3.8  2008/12/21 21:03:43  edwards
//  Moved include of chromabase into the .h file.
//
//  Revision 3.7  2008/08/18 18:23:56  jbulava
//  Added an additional constructor to sftmom class that allows a list of momenta
//  to be calculated rather than all momenta less that a max value.
//
//  Revision 3.6  2007/08/23 21:28:31  edwards
//  Bug fix. Removed some invalid uses of resize within a multi2d.
//
//  Revision 3.5  2007/08/01 19:33:14  edwards
//  Removed check for origin_offset must be zero in case of momentum averaging.
//  The check is not needed. The phase of the origin is removed in the correct
//  way no matter whether the phases are added together or not.
//
//  Revision 3.4  2007/06/21 18:18:55  edwards
//  Added subset versions of "sft" function.
//
//  Revision 3.3  2007/06/12 16:10:01  edwards
//  Added a default constructor.
//
//  Revision 3.2  2006/08/30 02:10:19  edwards
//  Technically a bug fix. The test for a zero_offset should only be in directions
//  not in the fourier transform. E.g., there was a missing test of mu==decay_dir.
//
//  Revision 3.1  2006/08/19 19:29:33  flemingg
//  Changed the interface of the slow Fourier transform (SftMom) class to allow
//  for any lattice point to be chosen as the spatial origin.  Previously, this
//  meant that the origin was always implicitly (0,0,0,0), which lead to
//  various phase problems.  One example is 2pt correlation functions between a
//  point source not at the origin and a sink at non-zero momentum: the overall
//  phase depends upon the location of the source.  Using the new interface,
//  the phase can be made independent of the location of the source by choosing
//  the origin of the Fourier transform to be the location of the source.
//  Several things remain to be fixed:
//    (1) sequential sources at non-zero sink momentum,
//    (2) building blocks with link paths which cross boundaries,
//    (3) baryonic sequential sources where t_source > t_sink,
//        i.e. the boundary is between the source and sink.
//
//  Revision 3.0  2006/04/03 04:59:12  edwards
//  Major overhaul of fermion and gauge action interface. Basically,
//  all fermacts and gaugeacts now carry out  <T,P,Q>  template parameters. These are
//  the fermion type, the "P" - conjugate momenta, and "Q" - generalized coordinates
//  in the sense of Hamilton's equations. The fermbc's have been rationalized to never
//  be over multi1d<T>. The "createState" within the FermionAction is now fixed meaning
//  the "u" fields are now from the coordinate type. There are now "ConnectState" that
//  derive into FermState<T,P,Q> and GaugeState<P,Q>.
//
//  Revision 2.2  2006/02/26 21:17:34  edwards
//  Put TimeSliceFunc in an anonymous namespace.
//
//  Revision 2.1  2005/09/29 15:53:42  edwards
//  Updates to interface.
//
//  Revision 2.0  2005/09/25 21:04:44  edwards
//  Moved to version 2.0
//
//  Revision 1.10  2005/07/27 16:23:28  edwards
//  Added getter for momentum offset
//
//  Revision 1.9  2005/01/14 18:42:38  edwards
//  Converted all lib files to be in chroma namespace.
//
//  Revision 1.8  2004/04/28 18:55:38  edwards
//  Added numSites().
//
//  Revision 1.7  2004/02/03 20:04:05  edwards
//  Added sft for LatticeReal. Changed  getSubset() to the more accurate
//  getSet().
//
//  Revision 1.6  2003/12/17 04:50:42  edwards
//  Added function to return decay direction. Added some doxygen comments.
//
//  Revision 1.5  2003/04/02 22:28:22  edwards
//  Changed proto.h to qdp_util.h
//
//  Revision 1.4  2003/04/01 02:46:43  edwards
//  Added const qual.
//
//  Revision 1.3  2003/03/20 19:34:25  flemingg
//  Evolved formfac_w.cc to use SftMom class, which included some bug fixes
//  in features in SftMom which had been previously untested and evolution
//  of the corresponding test program.
//
//  Revision 1.2  2003/03/14 17:13:44  flemingg
//  SftMom::sft() now works.
//
//  Revision 1.1  2003/03/14 05:06:06  flemingg
//  Initial version of SftMom class
//

#include "meas_ex/sftmom.h"
#include "util/ft/single_phase.h"
#include "qdp_util.h"                 // part of QDP++, for crtesn()

namespace Chroma 
{
  //! input
  void read(XMLReader& xml, const std::string& path, SftMomSrcPos_t& input)
  {
    XMLReader inputtop(xml, path);

    // These bad names are left in for backward compatibility; remove when no longer used
    if (inputtop.count("src_pos") == 1)
      read(inputtop, "src_pos", input.src_pos);
    else
      read(inputtop, "src", input.src_pos);

    if (inputtop.count("t_min") == 1)
      read(inputtop, "t_min", input.t_min);
    else
      read(inputtop, "tA", input.t_min);

    if (inputtop.count("t_max") == 1)
      read(inputtop, "t_max", input.t_max);
    else
      read(inputtop, "tB", input.t_max);
  }

  //! output
  void write(XMLWriter& xml, const std::string& path, const SftMomSrcPos_t& input)
  {
    push(xml, path);

    write(xml, "src_pos", input.src_pos);
    write(xml, "t_min",  input.t_min);
    write(xml, "t_max",  input.t_max);

    pop(xml);
  }

  // Anonymous namespace
  namespace
  {
    //! Function object used for constructing the time-slice set
    class TimeSliceFunc : public SetFunc
    {
    public:
      TimeSliceFunc(int dir): dir_decay(dir) {}

      int operator() (const multi1d<int>& coordinate) const;

      int numSubsets() const;

    private:
      TimeSliceFunc() {}  // hide default constructor

      int dir_decay;
    };
  }

  int
  TimeSliceFunc::operator() (const multi1d<int>& coordinate) const
  {
    return ((dir_decay<0) || (dir_decay>=Nd)? 0 : coordinate[dir_decay]);
  }

  int
  TimeSliceFunc::numSubsets() const
  {
    if ((dir_decay<0) || (dir_decay>=Nd)) {
      return 1;
    } else {
      return Layout::lattSize()[dir_decay];
    }
  }

  // Construct ball about origin in momentum space
  SftMom_cohen::SftMom_cohen(int mom2_max, bool avg_mom, int j_decay)
  {
    multi1d<int> origin_off(Nd); origin_off = 0;
    multi1d<int> mom_off((j_decay<0) || (j_decay>=Nd)? Nd : Nd-1);
    mom_off = 0;

    init(mom2_max, origin_off, mom_off, avg_mom, j_decay);
  }

  // Construct only listed momenta about origin
  SftMom_cohen::SftMom_cohen(const multi2d<int> & moms, int j_decay)
  {
    multi1d<int> origin(Nd);
    origin = 0;
    
    init(0, origin, origin, false, j_decay);
    
    num_mom = moms.size2();
    mom_list = moms;
    phases.resize(num_mom);

    for (int m=0; m < num_mom; ++m)
    {
      phases[m] = singlePhase(origin, mom_list[m], j_decay);
    }
  }

  // Construct ball in momentum space with origin_offset in space
  SftMom_cohen::SftMom_cohen(int mom2_max, multi1d<int> origin_offset_, bool avg_mom,
                 int j_decay)
  {
    multi1d<int> mom_off(((j_decay<0) || (j_decay>=Nd)? Nd : Nd-1));
    mom_off = 0;

    init(mom2_max, origin_offset_, mom_off, avg_mom, j_decay);
  }

  SftMom_cohen::SftMom_cohen(int mom2_max, multi1d<SftMomSrcPos_t> origin_offsets_,
                 bool avg_mom, int j_decay)
  {
    multi1d<int> mom_off(((j_decay<0) || (j_decay>=Nd)? Nd : Nd-1));
    mom_off = 0;

    init(mom2_max, origin_offsets_, mom_off, avg_mom, j_decay);
  }

  int
  SftMom_cohen::numSites() const
  {
    int vol = 1;

    if ((decay_dir<0) || (decay_dir>=Nd))
      vol = Layout::vol();
    else 
    {
      for (int m=0; m < Nd; ++m)
        vol *= Layout::lattSize()[m];
    }

    return vol;
  }

  // Construct ball in momentum space around mom_offset with origin_offset
  void SftMom_cohen::init(int mom2_max, multi1d<int> origin_off, multi1d<int> mom_off,
                    bool avg_mom, int j_decay) {
    multi1d<SftMomSrcPos_t> multi_srcs(1);
    multi_srcs[0].src_pos = origin_off;
    multi_srcs[0].t_min = 0;
    multi_srcs[0].t_max = Layout::lattSize()[j_decay] - 1;

    // Invoke the more general initializer
    init(mom2_max, multi_srcs, mom_off, avg_mom, j_decay);
  }
  // Use specific momenta with origin_offset
  void SftMom_cohen::init(multi1d<int> mom2_list, multi1d<int> origin_off,
                    multi1d<int> mom_off, bool avg_mom, int j_decay) {
    multi1d<SftMomSrcPos_t> multi_srcs(1);
    multi_srcs[0].src_pos = origin_off;
    multi_srcs[0].t_min = 0;
    multi_srcs[0].t_max = Layout::lattSize()[j_decay] - 1;

    // Invoke the more general initializer
    init(mom2_list, multi_srcs, mom_off, avg_mom, j_decay);
  }

  // Construct ball in momentum space around mom_offset with origin_offsets
  void
  SftMom_cohen::init(int mom2_max, multi1d<SftMomSrcPos_t> origin_offs,
               multi1d<int> mom_off, bool avg_mom, int j_decay)
  {
    // Include every mom2 shell inside the ball
    multi1d<int> mom2_list(mom2_max+1);
    for (int q2=0; q2 <= mom2_max; q2++) mom2_list[q2] = q2;
    
    // Invoke the more general initializer
    init(mom2_list, origin_offs, mom_off, avg_mom, j_decay);
  }

  void
  SftMom_cohen::init(multi1d<int> mom2_list, multi1d<SftMomSrcPos_t> origin_offs,
               multi1d<int> mom_off, bool avg_mom, int j_decay)
  {
    // Copy input parameters into class members
    decay_dir      = j_decay;
    origin_offsets = origin_offs;
    mom_offset     = mom_off;
    avg_equiv_mom  = avg_mom;
    
    sft_set.make(TimeSliceFunc(j_decay));

    //--------------------------------------------------------------------------
    // determine the number of momenta with mom^2 in mom2_list
    // If avg_equiv_mom is true then only consider momenta with
    // mom[0] >= mom[1] >= ... >= mom[Nd] >= 0
    
    int mom_dim = ((j_decay<0) || (j_decay>=Nd)? Nd : Nd-1);
    multi1d<int> mom_size(mom_dim); // Dimensions of containing cube

    // Find size of cube containing momentum shells
    int L; for (L=1; L*L <= mom2_list[mom2_list.size()-1]; ++L);
    
    int mom_vol = 1;
    for (int mu=0; mu < mom_dim; ++mu) {
      mom_size[mu] = (avg_equiv_mom? L : 2*L+1);
      mom_vol *= mom_size[mu];
    }

    // Iterate over cube, counting momenta falling on shells
    num_mom = 0;
    for (int n=0; n < mom_vol; ++n) {
      multi1d<int> mom = crtesn(n, mom_size);

      int mom2 = 0;
      for (int mu=0; mu < mom_dim; ++mu) {
        if (!avg_equiv_mom) mom[mu] -= L;
        mom2 += mom[mu]*mom[mu];
      }

      for (int i=0; i<mom2_list.size() && mom2>=mom2_list[i]; ++i) {
        if (mom2 == mom2_list[i]) {
          bool skip = false;
          if (avg_equiv_mom) {
            // Ensure mom[0] >= mom[1] >= ... >= mom[Nd] >= 0
            for (int mu=0; mu < mom_dim-1; ++mu)
              for (int nu=mu+1; nu < mom_dim; ++nu)
                if (mom[nu] > mom[mu]) skip = true;
          }
          if (!skip) ++num_mom;
          break; // after finding the momentum in mom2_list
        } // end if in mom2_list
      } // end for loop over mom2_list
    }

    // After all those shenanigans to get num_mom, resize the momentum list
    mom_list.resize(num_mom, mom_dim);

    //--------------------------------------------------------------------------
    // Now we do exactly the same thing we did when counting, except this time
    // we can actually fill the list
    int num_filled = 0;
    for (int n=0; n < mom_vol; ++n) {
      multi1d<int> mom = crtesn(n, mom_size);

      int mom2 = 0;
      for (int mu=0; mu < mom_dim; ++mu) {
        if (!avg_equiv_mom) mom[mu] -= L;
        mom2 += mom[mu]*mom[mu];
      }

      for (int i=0; i<mom2_list.size() && mom2>=mom2_list[i]; ++i) {
        if (mom2 == mom2_list[i]) {
          bool skip = false;
          if (avg_equiv_mom) {
            // Ensure mom[0] >= mom[1] >= ... >= mom[Nd] >= 0
            for (int mu=0; mu < mom_dim-1; ++mu)
              for (int nu=mu+1; nu < mom_dim; ++nu)
                if (mom[nu] > mom[mu]) skip = true;
          }
          if (!skip) {
            for (int mu=0; mu < mom_dim; ++mu) {
              mom_list[num_filled][mu] = mom_offset[mu] + mom[mu];
            }
            ++num_filled;
          }
          break; // after finding the momentum in mom2_list
        } // end if in mom2_list
      } // end for loop over mom2_list
    }

    //--------------------------------------------------------------------------
    // Now resize and initialize the Fourier phase table.  Then, loop over
    // allowed momenta, optionally averaging over equivalent momenta.
    phases.resize(num_mom);
    phases = 0.;

    // Coordinates for sink momenta
    multi1d<LatticeInteger> my_coord(Nd);
    for (int mu=0; mu < Nd; ++mu)
      my_coord[mu] = Layout::latticeCoordinate(mu);

    // Keep track of |mom| degeneracy for averaging
    mom_degen.resize(num_mom);
    mom_degen = 0;

    // If averaging over equivalent momenta, we need redo mom_size and mom_vol
    // to allow both positive and negative momentum components
    if (avg_equiv_mom) {
      mom_vol = 1;
      for (int mu=0; mu < mom_dim; ++mu) {
        mom_size[mu] = 2*L+1;
        mom_vol *= mom_size[mu];
      }
    }

    // Index for iterating over the momentum list
    int mom_idx = 0;
    for (int n=0; n < mom_vol; ++n) {
      multi1d<int> mom = crtesn(n, mom_size);

      int mom2 = 0;
      for (int mu=0; mu < mom_dim; ++mu) {
        mom[mu] -= L;
        mom2 += mom[mu]*mom[mu];
      }

      // skip when (mom)^2 is not in the mom2_list
      for (int i=0; i<mom2_list.size() && mom2>=mom2_list[i]; ++i) {
        if (mom2 == mom2_list[i]) {

          // At this point, if (avg_equiv_mom == true) then we need to determine
          // mom_idx by a fairly time consuming process.
          // If (avg_equiv_mom == false) then (mom == mom_list[mom_idx])
          // (will double check this) and the momentum offset can be applied.

          if (avg_equiv_mom) {

            // determine mom_idx for entering table mom_list
            // put the momenta into canonical order
            multi1d<int> mom_tmp = canonicalOrder(mom);

            // mom_tmp should now contain a momentum that appears in mom_list.
            // scan through list until we find a match.
            mom_idx = -1;

            for (int k=0; k < num_mom; ++k) {
              bool match = true;
              for (int mu=0; mu < mom_tmp.size(); ++mu) {
                if (mom_list[k][mu] != mom_tmp[mu]) {
                  match = false;
                  break;
                }
              }
              if (match) {
                mom_idx = k;
                break;
              }
            }

            if (mom_idx < 0) {
              QDP_error_exit("SftMom_cohen: mom_idx < 0. Momentum not found on list.\n");
            }

            // increment degeneracy for this mom_idx
            ++(mom_degen[mom_idx]);
          } else /* (avg_equiv_mom == false) */ {

            // apply momentum offset
            for (int mu=0; mu < mom_dim; ++mu) {
              mom[mu] += mom_offset[mu];
            }

            // double check that (mom == mom_list[n])
            // this check should never fail and could be removed in the future
            for (int mu=0; mu < mom_dim; ++mu) {
              if (mom[mu] != mom_list[mom_idx][mu]) {
                // Should never get here !!!
                QDP_error_exit("SftMom_cohen: (mom != mom_list[mom_idx])\n");
              }
            }
          } // end if (avg_equiv_mom)

          // Different origin_offsets for different timeslices
          multi1d<LatticeInteger> origin_offsets_lat(Nd);
          Set timeslice;
          timeslice.make(TimeSliceFunc(j_decay));
          int tlen = Layout::lattSize()[j_decay];


          for (int i=0; i < origin_offsets.size(); ++i)
          {
            int tmin = origin_offsets[i].t_min;
            int tmax = origin_offsets[i].t_max;
            multi1d<int> origin = origin_offsets[i].src_pos;

            if (tmin <= tmax)
            {
              for (int tt=tmin; tt <= tmax; ++tt)
              {
                for (int mu=0; mu < Nd; ++mu) 
                  origin_offsets_lat[mu][timeslice[tt]] = origin[mu];
              }
            }
            else // tA > tB case: around the world,
                 // e.g. tlen=64, tmin=62, tmax=1; t={62, 63, 0, 1}
            {
              for (int tt=0; tt <= tmax; ++tt)
              {
                for (int mu=0; mu < Nd; ++mu) 
                  origin_offsets_lat[mu][timeslice[tt]] = origin[mu];
              }
              for (int tt=tmin; tt < tlen; ++tt)
              {
                for (int mu=0; mu < Nd; ++mu) 
                  origin_offsets_lat[mu][timeslice[tt]] = origin[mu];
              }
            } // end if tmax >= tmin
          } // end for loop over offsets

          //
          // Build the phase. 
          // RGE: the origin_offset works with or without momentum averaging
          //
          LatticeReal p_dot_x;
          p_dot_x = 0.;

          int j = 0;
          for (int mu = 0; mu < Nd; ++mu) {
            const Real twopi = 6.283185307179586476925286;

            if (mu == j_decay) continue;

            p_dot_x += LatticeReal(my_coord[mu] - origin_offsets_lat[mu]) * twopi *
                    Real(mom[j]) / Layout::lattSize()[mu];
            ++j;
          } // end for (mu)

          phases[mom_idx] += cmplx(cos(p_dot_x), sin(p_dot_x));

          // increment mom_idx for next valid momenta
          ++mom_idx;
          break; // after finding the momentum in mom2_list
        } // end if in mom2_list
      } // end for loop over mom2_list

    } // end for (int n=0; n < mom_vol; ++n)

    // Finish averaging
    // Momentum averaging works even in the presence of an origin_offset
    if (avg_equiv_mom) {
      for (int mom_idx=0; mom_idx < num_mom; ++mom_idx)
        phases[mom_idx] /= mom_degen[mom_idx];
    }
  }


  // Canonically order an array of momenta
  /* \return abs(mom[0]) >= abs(mom[1]) >= ... >= abs(mom[Nd]) >= 0 */
  multi1d<int> 
  SftMom_cohen::canonicalOrder(const multi1d<int>& mom) const
  {
    // first step: make all the components positive
    multi1d<int> mom_tmp = mom;
    for (int mu=0; mu < mom_tmp.size(); ++mu)
      if (mom_tmp[mu] < 0) mom_tmp[mu] = -mom_tmp[mu];

    // Initially, the first item is considered sorted.  mu divides mom
    // into a sorted region (<mu) and an unsorted one (>=mu)
    for (int mu=1; mu < mom_tmp.size(); ++mu) 
    {
      // Select the item at the beginning of the unsorted region
      int v = mom_tmp[mu];
      // Work backwards, finding where v should go
      int nu = mu;
      // If this element is less than v, move it up one
      while (mom_tmp[nu-1] < v) {
        mom_tmp[nu] = mom_tmp[nu-1];
        --nu;
        if (nu < 1) break;
      }
      // Stopped when mom_tmp[nu-1] >= v, so put v at postion nu
      mom_tmp[nu] = v;
    }

    return mom_tmp;
  }


  // Convert array of momenta to momenta id
  /* \return id in [0,numMom()-1] or -1 if not in list */
  int 
  SftMom_cohen::momToNum(const multi1d<int>& mom_in) const
  {
    multi1d<int> mom;

    // If mom avg is turned on, then canonicalize the input mom
    if (avg_equiv_mom)
      mom = canonicalOrder(mom_in);
    else
      mom = mom_in;

    // Search for the mom
    for (int mom_idx=0; mom_idx < num_mom; ++mom_idx) 
    {
      bool match = true;
      for (int mu=0; mu < mom.size(); ++mu)
      {
        if (mom_list[mom_idx][mu] != mom[mu]) 
        {
          match = false;
          break;
        }
      }
      if (match) return mom_idx;
    }
    return -1;
  }

  multi2d<DComplex>
  SftMom_cohen::sft(const LatticeComplex& cf) const
  {
    multi2d<DComplex> hsum(num_mom, sft_set.numSubsets());

    for (int mom_idx=0; mom_idx < num_mom; ++mom_idx)
      hsum[mom_idx] = sumMulti(phases[mom_idx]*cf, sft_set);

    return hsum;
  }

  multi2d<DComplex>
  SftMom_cohen::sft(const LatticeComplex& cf, int subset_color) const
  {
    int length = sft_set.numSubsets();
    multi2d<DComplex> hsum(num_mom, length);

    for (int mom_idx=0; mom_idx < num_mom; ++mom_idx)
    {
      hsum[mom_idx] = zero;
      hsum[mom_idx][subset_color] = sum(phases[mom_idx]*cf, sft_set[subset_color]);
    }

    return hsum;
  }

  multi2d<DComplex>
  SftMom_cohen::sft(const LatticeReal& cf) const
  {
    multi2d<DComplex> hsum(num_mom, sft_set.numSubsets());

    for (int mom_idx=0; mom_idx < num_mom; ++mom_idx)
      hsum[mom_idx] = sumMulti(phases[mom_idx]*cf, sft_set);

    return hsum;
  }

  multi2d<DComplex>
  SftMom_cohen::sft(const LatticeReal& cf, int subset_color) const
  {
    int length = sft_set.numSubsets();
    multi2d<DComplex> hsum(num_mom, length);

    for (int mom_idx=0; mom_idx < num_mom; ++mom_idx)
    {
      hsum[mom_idx] = zero;
      hsum[mom_idx][subset_color] = sum(phases[mom_idx]*cf, sft_set[subset_color]);
    }

    return hsum;
  }

#if BASE_PRECISION==32
  multi2d<DComplex>
  SftMom_cohen::sft(const LatticeComplexD& cf) const
  {
    multi2d<DComplex> hsum(num_mom, sft_set.numSubsets());

    for (int mom_idx=0; mom_idx < num_mom; ++mom_idx)
      hsum[mom_idx] = sumMulti(phases[mom_idx]*cf, sft_set);

    return hsum;
  }

  multi2d<DComplex>
  SftMom_cohen::sft(const LatticeComplexD& cf, int subset_color) const
  {
    int length = sft_set.numSubsets();
    multi2d<DComplex> hsum(num_mom, length);

    for (int mom_idx=0; mom_idx < num_mom; ++mom_idx)
    {
      hsum[mom_idx] = zero;
      hsum[mom_idx][subset_color] = sum(phases[mom_idx]*cf, sft_set[subset_color]);
    }

    return hsum;
  }
#endif

}  // end namespace Chroma

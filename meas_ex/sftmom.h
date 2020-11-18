// -*- C++ -*-
//  $Id: sftmom.h,v 3.11 2015-05-24 11:10 sdcohen $
/*! \file
 *  \brief Fourier transform phase factor support
 */

#ifndef __sftmom_cohen_h__
#define __sftmom_cohen_h__

#include "chromabase.h"
#include "util/ft/sftmom.h"

namespace Chroma {

//! Param struct for SftMom
/*!
 * \ingroup ft
 */
// Structure of source positions for FT origin
struct SftMomSrcPos_t {
  multi1d<int> src_pos;
  // This source will only have support on t_min <= t <= t_max
  int t_min, t_max;
};

//! input/output
void read(XMLReader& xml, const std::string& path, SftMomSrcPos_t& input);
void write(XMLWriter& xml, const std::string& path,
    const SftMomSrcPos_t& input);

//! Fourier transform phase factor support
/*!
 * \ingroup ft
 */
class SftMom_cohen {
public:
  //! Construct ball about origin in momentum space
  explicit SftMom_cohen(int mom2_max, bool avg_equiv_mom_ = false, int j_decay = -1);

  //! Construct only listed momenta about origin
  explicit SftMom_cohen(const multi2d<int>& moms, int j_decay = -1);

  //! Construct ball in momentum space with origin_offset in space
  SftMom_cohen(int mom2_max, multi1d<int> origin_offset_, bool avg_equiv_mom_ = false,
      int j_decay = -1);

  //! Construct ball in momentum space with origin offset in space for multiple times
  SftMom_cohen(int mom2_max, multi1d<SftMomSrcPos_t> origin_offsets_,
      bool avg_equiv_mom_ = false, int j_decay = -1);

  //! Construct listed momenta with origin offset
  template<typename T>
  SftMom_cohen(T mom_spec, multi1d<int> origin_offset_, bool avg_equiv_mom_ = false,
      int j_decay = -1) {
    multi1d<int> mom_offset(((j_decay < 0) || (j_decay >= Nd) ? Nd : Nd - 1));
    mom_offset = 0;

    init(mom_spec, origin_offset_, mom_offset, avg_equiv_mom_, j_decay);
  }

  //! Construct shells in momentum space around mom_offset with origin_offset
  SftMom_cohen(multi1d<int> mom2_list, multi1d<int> origin_offset_,
      multi1d<int> mom_offset_, bool avg_equiv_mom_ = false, int j_decay = -1) {
    init(mom2_list, origin_offset_, mom_offset_, avg_equiv_mom_, j_decay);
  }

  //! Construct ball in momentum space around mom_offset with origin_offset
  SftMom_cohen(int mom2_max, multi1d<int> origin_offset_, multi1d<int> mom_offset_,
      bool avg_equiv_mom_ = false, int j_decay = -1) {
    init(mom2_max, origin_offset_, mom_offset_, avg_equiv_mom_, j_decay);
  }

  //! Construct ball in momentum space with origin offset for multiple times
  template<typename T>
  SftMom_cohen(T mom_spec, multi1d<SftMomSrcPos_t> origin_offsets_,
      bool avg_equiv_mom_ = false, int j_decay = -1) {
    multi1d<int> mom_offset(((j_decay < 0) || (j_decay >= Nd) ? Nd : Nd - 1));
    mom_offset = 0;

    init(mom_spec, origin_offsets_, mom_offset, avg_equiv_mom_, j_decay);
  }

  //! Construct ball in momentum space with origin offset and mom_offset for multiple times
  template<typename T>
  SftMom_cohen(T mom_spec, multi1d<SftMomSrcPos_t> origin_offsets_,
      multi1d<int> mom_offset_, bool avg_equiv_mom_ = false, int j_decay = -1) {
    init(mom_spec, origin_offsets_, mom_offset_, avg_equiv_mom_, j_decay);
  }

  //! General constructor from parameter structure
  SftMom_cohen(const SftMomParams_t& p) {
    init(p.mom2_max, p.origin_offset, p.mom_offset, p.avg_equiv_mom,
      p.decay_dir);
  }

  //! The set to be used in sumMulti
  const Set& getSet() const {
    return sft_set;
  }

  //! Number of momenta
  int numMom() const {
    return num_mom;
  }

  //! Number of subsets - length in decay direction
  int numSubsets() const {
    return sft_set.numSubsets();
  }

  //! Number of sites in each subset
  int numSites() const;

  //! Decay direction
  int getDir() const {
    return decay_dir;
  }

  //! Are momenta averaged?
  bool getAvg() const {
    return avg_equiv_mom;
  }

  //! Momentum offset
  multi1d<int> getMomOffset() const {
    return mom_offset;
  }

  //! Convert momenta id to actual array of momenta
  multi1d<int> numToMom(int mom_num) const {
    return mom_list[mom_num];
  }

  //! Convert array of momenta to momenta id
  /*! \return id in [0,numMom()-1] or -1 if not in list */
  int momToNum(const multi1d<int>& mom_in) const;

  //! Canonically order an array of momenta
  /*! \return abs(mom[0]) >= abs(mom[1]) >= ... >= abs(mom[mu]) >= ... >= 0 */
  multi1d<int> canonicalOrder(const multi1d<int>& mom) const;

  //! Return the phase for this particular momenta id
  const LatticeComplex& operator[](int mom_num) const {
    return phases[mom_num];
  }

  //! Return the the multiplicity for this momenta id.
  /*! Only nonzero if momentum averaging is turned on */
  int multiplicity(int mom_num) const {
    return mom_degen[mom_num];
  }

  //! Do a sumMulti(cf*phases,getSet())
  multi2d<DComplex> sft(const LatticeComplex& cf) const;

  //! Do a sum(cf*phases,getSet()[my_subset])
  multi2d<DComplex> sft(const LatticeComplex& cf, int subset_color) const;

  //! Do a sumMulti(cf*phases,getSet())
  multi2d<DComplex> sft(const LatticeReal& cf) const;

  //! Do a sumMulti(cf*phases,getSet()[my_subset])
  multi2d<DComplex> sft(const LatticeReal& cf, int subset_color) const;

#if BASE_PRECISION==32
  multi2d<DComplex> sft(const LatticeComplexD& cf) const;
  //! Do a sum(cf*phases,getSet()[my_subset])
  multi2d<DComplex> sft(const LatticeComplexD& cf, int subset_color) const;
#endif

private:
  SftMom_cohen() {
  }  // hide default constructor

  void init(int mom2_max, multi1d<int> origin_offset, multi1d<int> mom_offset,
      bool avg_mom_ = false, int j_decay = -1);
  void init(multi1d<int> mom2_list, multi1d<int> origin_offset,
      multi1d<int> mom_offset, bool avg_mom_ = false, int j_decay = -1);
  void init(int mom2_max, multi1d<SftMomSrcPos_t> origin_offsets,
      multi1d<int> mom_offset, bool avg_mom_ = false, int j_decay = -1);
  void init(multi1d<int> mom2_list, multi1d<SftMomSrcPos_t> origin_offsets,
      multi1d<int> mom_offset, bool avg_mom_ = false, int j_decay = -1);

  multi2d<int> mom_list;
  bool avg_equiv_mom;
  int decay_dir;
  int num_mom;
  multi1d<SftMomSrcPos_t> origin_offsets;
  multi1d<int> mom_offset;
  multi1d<LatticeComplex> phases;
  multi1d<int> mom_degen;
  Set sft_set;
};

}  // end namespace Chroma

#endif

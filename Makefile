# define TOPDIR (for SCIDAC path) and CHROMA (for CHROMA path) in Makefile.local
include Makefile.local

CONFIG=$(CHROMA)/bin/chroma-config

MGLIBS= -lwilsonmg -lqopqdp -lqdp_common -lqdp_int -lqdp_f -lqdp_d -lqdp_df -lqdp_f2 -lqdp_d2 -lqdp_df2 -lqdp_f3 -lqdp_d3 -lqdp_df3 -lqdp_fn -lqdp_dn -lqdp_dfn -lqio -llime -lqla_c99 -lqla_cmath -lqla_d2 -lqla_d3 -lqla_d -lqla_df2 -lqla_df3 -lqla_df -lqla_dfn -lqla_dn -lqla_dq2 -lqla_dq3 -lqla_dq -lqla_dqn -lqla_f2 -lqla_f3 -lqla_f -lqla_fn -lqla_int -lqla_q2 -lqla_q3 -lqla_q -lqla_qn -lqla_random
MGLDFLAGS=-L$(TOPDIR)/qla/lib -L$(TOPDIR)/qopqdp/lib
MGCXXFLAGS=-L$(TOPDIR)/qla/lib -L$(TOPDIR)/qdp/lib -L$(TOPDIR)/qopqdp/lib

CXX=$(shell $(CONFIG) --cxx) $(CXXFLAGS_EX)
CXXFLAGS=$(shell $(CONFIG) --cxxflags) -I$(TOPDIR)/qla/include -I. $(MGCXXFLAGS) 
LDFLAGS=$(shell $(CONFIG) --ldflags) $(MGLDFLAGS)
LIBS=$(shell $(CONFIG) --libs) $(MGLIBS)

HDRS=\
     actions_ex/clover_utils_exp_w.h \
     actions_ex/eoprec_exp_clover_linop_w.h \
     actions_ex/eoprec_exp_clover_fermact_w.h \
     actions_ex/hwilson_eigenop.h \
     actions_ex/overlap_eigenop.h \
     actions_ex/target_quda/quda_multigrid_params.h \
     actions_ex/target_quda/quda_mg_utils.h \
     actions_ex/target_quda/syssolver_linop_clover_quda_multigrid_w.h \
     actions_ex/target_quda/syssolver_mdagm_clover_quda_multigrid_w.h \
     actions_ex/target_quda/syssolver_linop_clover_quda_w.h \
     meas_ex/momentum_volume_npr_source.h \
     meas_ex/grid_source.h \
     meas_ex/mom_quark_smearing.h \
     meas_ex/npr_vertex_pdf.h \
     meas_ex/sftmom.h\
     meas_ex/simple_baryon_seqsrc_w.h \
     meas_ex/Xigauge.h \
     inline_ex/inline_FixObc.h \
     inline_ex/inline_npr_vertex_pdf.h \
     inline_ex/inline_Xigauge.h \
     inline_ex/inline_qpropadd_w.h \
     inline_ex/inline_seqsource_fast_w.h \
     inline_ex/inline_erase_quda_multigrid_space.h \
     inline_ex/inline_loop_w.h \
     inline_ex/inline_hadspec_w.h \
     inline_ex/inline_building_blocks_w.h \
     inline_ex/inline_building_x_w.h \
     inline_ex/inline_pion_DA.h \
     inline_ex/inline_vev.h \
     inline_ex/inline_eigen_maker.h \
     inline_ex/inline_propagator_multi_eigen.h \

OBJS=chroma.o \
     io_ex/io_general.o\
     io_ex/io_general_class.o\
     io_ex/kyuqprop_io.o\
     actions_ex/eoprec_exp_clover_linop_w.o \
     actions_ex/eoprec_exp_clover_fermact_w.o \
     actions_ex/target_quda/quda_multigrid_params.o \
     actions_ex/target_quda/syssolver_linop_clover_quda_multigrid_w.o \
     actions_ex/target_quda/syssolver_mdagm_clover_quda_multigrid_w.o \
     actions_ex/target_quda/syssolver_linop_clover_quda_w.o \
     meas_ex/sftmom.o\
     meas_ex/simple_baryon_seqsrc_w.o \
     meas_ex/Xigauge.o \
     meas_ex/momentum_volume_npr_source.o \
     meas_ex/grid_source.o \
     meas_ex/mom_quark_smearing.o \
     meas_ex/npr_vertex_pdf.o \
     inline_ex/inline_FixObc.o \
     inline_ex/inline_npr_vertex_pdf.o \
     inline_ex/inline_Xigauge.o \
     inline_ex/inline_qpropadd_w.o \
     inline_ex/inline_seqsource_fast_w.o \
     inline_ex/inline_erase_quda_multigrid_space.o \
     inline_ex/inline_loop_w.o \
     inline_ex/inline_hadspec_w.o \
     inline_ex/inline_building_blocks_w.o \
     inline_ex/inline_building_x_w.o \
     inline_ex/inline_pion_DA.o \
     inline_ex/inline_vev.o \
     inline_ex/inline_eigen_maker.o \
     inline_ex/inline_propagator_multi_eigen.o \
     

OBJS2=hmc.o\
      actions_ex/target_quda/syssolver_mdagm_clover_quda_multigrid_w.o\
      actions_ex/target_quda/quda_multigrid_params.o \
      actions_ex/target_quda/syssolver_linop_clover_quda_multigrid_w.o \
      io_ex/io_general.o\
      io_ex/io_general_class.o\
      meas_ex/grid_source.o \
      inline_ex/inline_hadspec_w.o \

chroma: $(OBJS)
	$(CXX) -o $@ $(CXXFLAGS) $(OBJS) $(LDFLAGS) $(LIBS)

hmc: $(OBJS2)
	$(CXX) -o $@ $(CXXFLAGS) $(OBJS2) $(LDFLAGS) $(LIBS)


%.o: %.cc $(HDRS)
	$(CXX) $(CXXFLAGS) -c $< 

clean:
	rm -rf chroma hmc $(OBJS) *~ io_ex/*~ actions_ex/*~ actions_ex/target_quda/*~ meas_ex/*~ inline_ex/*~

#                                                                           c
#
include make.inc
#
#==========================================#
# List object files needed for compilation #
#  OBJS0 is used for make a                #
#  OBJS1 is used for make b                #
#==========================================#
OBJS0   = testing.o makefilename.o inputdata.o getfilename.o openfile.o numbernodeseqs.o writesetup.o dmacheps.o legendrep.o quadtest.o condnumber.o rotg.o quatmult.o qnorm.o vertex.o spharm.o  sobolev.o eqquad.o eqquadpw.o neqest.o toarray.o coefjac.o coefjacRF.o eqcount.o stepnewt_testing.o fcenter.o scenter.o toCart.o output.o xyzpair.o  
OBJS1   = testing.o makefilename.o inputfromfile.o getfilename.o openfile.o numbernodeseqs.o writesetup.o dmacheps.o legendrep.o quadtest.o condnumber.o rotg.o quatmult.o qnorm.o vertex.o spharm.o spharmN.o sobolev.o eqquad.o eqquadpw.o neqest.o toarray.o coefjac.o coefjacRF.o eqcount.o stepnewt_testing.o fcenter.o scenter.o toCart.o output.o xyzpair.o  
#
#  
#==========================#
# Custum libraries to link #
#==========================#
#LIBRARY  = $(DIRLIB)/libio.a $(DIRLIB)/libqrmin.a 
LIBRARY  = $(DIRLIB)/libio.a 
#
#================================================
#
#
#=====================================#
# "DPdriver-a" is the executable name #
#=====================================#
a: $(OBJS0) $(LIBRARY)
	$(FORTRAN)  $(FFLAGS)  $(OTHERLIB) -o DPdriver-a  $(OBJS0) $(LIBRARY)
#================================================
#
#================================================
#
#
b: $(OBJS1) $(LIBRARY)
	$(FORTRAN)  $(FFLAGS) -o DPdriver-b  $(OBJS1) $(LIBRARY)  $(OTHERLIB)
#================================================
#
#
.f.o:
	$(FORTRAN) -c $(FFLAGS) $*.f 

































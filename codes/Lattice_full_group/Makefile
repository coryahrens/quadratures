#vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvc
#                                                                           c
#        This file is a part of                                             c
#        Fast  Transform Library                                            c
#        Contains proprietary information supplied by GB Consulting.        c
#        Copyright (C), 1993-2008 GB Consulting. All rights reserved        c
#                                                                           c
#vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvc
#                                                                           c
#
include make.inc
#
#==========================================#
# List object files needed for compilation #
#  OBJS0 is used for make a                #
#  OBJS1 is used for make b                #
#==========================================#
OBJS0 = main.o fcenter.o bccoords.o bctoxyz.o makefilename.o
#  
#==========================#
# Custum libraries to link #
#==========================#
LIBRARY  = $(DIRLIB)/libio.a 
#
#================================================
#
#
#=====================================#
# "DPdriver-a" is the executable name #
#=====================================#
a: $(OBJS0) $(LIBRARY)
	$(FORTRAN)  $(FFLAGS) -o main  $(OBJS0) $(LIBRARY) $(OTHERLIB)
#================================================
#
clean:
	rm *.o
#
.f.o:
	$(FORTRAN) -c $(FFLAGS) $*.f 

































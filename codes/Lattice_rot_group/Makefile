
#	$Log: Makefile,v $
# 
#vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvc
#                                                                           c
#        This file is a part of                                             c
#        Fast  Transform Library                                            c
#        Contains proprietary information supplied by GB Consulting.        c
#        Copyright (C), 1993-2008 GB Consulting. All rights reserved        c
#                                                                           c
#vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvc
#                                                                           c
# Makefile for compiling programs.
#
include make.inc
#
#================
# 
OBJS0   = main-fg.o fcenter.o bccoords.o bctoxyz.o makefilename.o
#
LIBRARY  = $(DIRLIB)/libio.a  
#
#================================================
#
#
t: $(OBJS0) $(LIBRARY)
	$(FORTRAN)  $(FFLAGS) -o main  $(OBJS0) $(LIBRARY)  $(otherlib)
	rm -f *.o
#================================================
#
#================================================
#
#
#s: $(OBJS1) $(LIBRARY)
#	$(FORTRAN)  $(FFLAGS) -o testing  $(OBJS0) $(LIBRARY)  $(otherlib)
#================================================
#
#
.f.o:
	$(FORTRAN) -c $(FFLAGS) $*.f 

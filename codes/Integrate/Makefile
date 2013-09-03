
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
OBJS0    = driver.o cart2sphere.o inputdata.o getfilename.o openfile.o lpn.o spharm.o 
#
LIBRARY  = $(DIRLIB)/libio.a 
#
#================================================
#
#
t: $(OBJS0) $(LIBRARY)
	$(FORTRAN)  $(FFLAGS) -o driver  $(OBJS0) $(LIBRARY)  
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

































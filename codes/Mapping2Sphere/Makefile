
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
OBJS0   = main.o gen2D3D.o ico2plane.o findbcc3d.o sphere2ico.o  scalerotate.o findbcc2d.o bctocart.o spDist.o dmacheps.o findneigh.o coarse2sphere.o cross.o
#
LIBRARY  = $(DIRLIB)/libio.a  
#
#================================================
#
#
t: $(OBJS0) $(LIBRARY)
	$(FORTRAN)  $(FFLAGS) -o main  $(OBJS0) $(LIBRARY)  $(otherlib)
#================================================
#
.f.o:
	$(FORTRAN) -c $(FFLAGS) $*.f 

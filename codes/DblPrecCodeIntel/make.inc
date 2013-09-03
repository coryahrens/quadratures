##############################
#  make include file         #
#  to change parameters      #
#  for different platforms   #
##############################
#
#
#=================================#
# Name of Fortran and C compilers #
#=================================#
#FORTRAN = gfortran
#FORTRAN = lf95
FORTRAN = ifort
CC      = gcc 
ENDIAN  = LITTLE
#
#
#================#
# Compiler flags #
#================#
##FFLAGS = -pg 
##FFLAGS = -L$(MKLPATH) -I$(MKLINCLUDE)
##FFLAGS=-O3 -parallel -L/usr/local/composerxe-2011.5.220/mkl/lib/intel64/ -L$(HOME)/usr/lib/intel
FFLAGS = -O3 
##
##
#
#====================#
# Personal libraries #
#====================#
DIRLIB  = $(HOME)/usr/lib/intel
#
#=====================#
# 3rd party libraries #
#=====================#
OTHERLIB = -lnsl -lblas -llapack
#OTHERLIB = -lnsl  -lblas -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
#
#




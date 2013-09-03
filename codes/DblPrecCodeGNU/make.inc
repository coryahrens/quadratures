##############################
#  make include file         #
#  to change parameters      #
#  for different platforms   #
##############################
#
#
#====================#
# Compiler and flags #
#====================#
F90 := gfortran
F90LINK := $(F90)
F90FLAGS := -fdefault-integer-8 
F90LINKFLAGS :=  -fopenmp
F90LINKLIBS := -lrt -ldl

#=======Location of ACML installation=======#
ACMLDIR := /opt/acml5.3.1/gfortran64_mp_int64
LIBACML := $(ACMLDIR)/lib/libacml_mp.a

OMP_NUM_THREADS := 16

EXPORT_PATH := export LD_LIBRARY_PATH='$(ACMLDIR)/lib:$(LD_LIBRARY_PATH)

 
#================#
# Compiler flags #
#================#
#
#====================#
# Personal libraries #
#====================#
DIRLIB  = $(HOME)/usr/lib/gnu
#
#=====================#
# 3rd party libraries #
#=====================#
OTHERLIB = -lnsl
#
#




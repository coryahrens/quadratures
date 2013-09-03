##############################
#  make include file         #
#  to change parameters      #
#  for different platforms   #
##############################
#
#
FORTRAN  = ifort
CC       = gcc 
#
ENDIAN   = LITTLE
#
# intel chips
#
#FFLAGS     =  -march=pentium4 -msse2  -malign-double -O3 -funroll-loops
#FFLAGS   = --chk  -g  --info  
FFLAGS    = 
####
####
CFLAGS     = $(FFLAGS) 
#
DIRLIBRARY =../libs
DIRLIB     = $(HOME)/usr/lib/intel
#DIRBIN     =../bin
#
otherlib=
#-lnsl -lblas -llapack
#
#otherlib=-lxnet
#




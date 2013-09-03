##############################
#  make include file         #
#  to change parameters      #
#  for different platforms   #
##############################
#
#
FORTRAN = lf95
CC       = gcc 
ENDIAN   = LITTLE
#
# intel chips
#
FFLAGS = 
#FFLAGS     =  -march=pentium4 -msse2  -malign-double -O3 -funroll-loops
#FFLAGS     =  -extend_source 80 -check all -O3
####
####
CFLAGS     = $(FFLAGS) 
#
#
DIRLIBRARY=../libs
DIRLIB  = $(HOME)/usr/lib/lf
DIRBIN  =../bin
#
otherlib=-lnsl -lblas -llapack
#
#otherlib=-lxnet
#




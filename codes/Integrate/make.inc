##############################
#  make include file         #
#  to change parameters      #
#  for different platforms   #
##############################
#
#
FORTRAN  = lf95
#CC       = gcc 
#CC       = egcs -DLINUX
##AR       = ar
ENDIAN   = LITTLE
#
# intel chips
#
FFLAGS     = --o2 --fast 
####
####
CFLAGS     = $(FFLAGS) 
##ARFLAGS    = r
#
DIRLIB  = $(HOME)/usr/lib/lf
#
otherlib=-lnsl -lblas -llapack
#




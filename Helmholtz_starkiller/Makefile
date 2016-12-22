#	Makefile for the Helmholtz equation of state

Eos += eos_helmConstData.o  eos_vecData.o eos_helmData.o\
       eos_readHfet.o eos_writeHfet.o \
       eos_helm.o eos_vecAlloc.o eos_vecDealloc.o \
       eos_helmInterface.o

eos_vecAlloc.o : eos_vecData.o  
eos_vecDealloc.o : eos_vecData.o  
Eos.o : eos_vecData.o  Eos_data.o 

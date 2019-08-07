#Makefile for stereotomography
include ./Makefile.config
B=../bin
CFLAGS=                 \
	-I../include	\
	-I$(CWPROOT)/include        \
        -I./include     \
	-I$(SUITROOT)/include	\
	-L../lib       \
	-L/usr/local/gfortran/lib	\
	-L$(SUITROOT)/lib	\
	-L$(CWPROOT)/lib       \
	-L./lib

LFLAGS=                 \
	-lspqr -lsuitesparseconfig -lcholmod -lamd -lcolamd	\
	-lm -llapack -lblas -ltbb	\
	-lsub -lsu -lpar -lcwp -lstdc++

PROGS=			\
	$B/kinematic	\
	$B/dynamic
ALL: GCC

r:
	cd ./src; $(MAKE)

GCC: $(PROGS)
$(PROGS): $B/%: %.c
	$(CC) $(CFLAGS) $(@F).c $(LFLAGS) -o $@ -g


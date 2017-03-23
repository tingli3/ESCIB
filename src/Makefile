GCC	:= g++


TARGETS := io countPoints clusters
OBJS    := $(TARGETS:=.o)
SRCS    := $(TARGETS:=.c)
HDRS    := $(TARGETS:=.h)



all: ESCIB_Bernoulli ESCIB_Poisson DBSCAN

$(OBJS): %.o: %.c %.h
	$(GCC) -o $@ -c $<

ESCIB_Bernoulli.o: ESCIB_Bernoulli.c
	$(GCC) -o $@ -c $<

ESCIB_Poisson.o: ESCIB_Poisson.c
	$(GCC) -o $@ -c $<

DBSCAN.o: DBSCAN.c
	$(GCC) -o $@ -c $<

ESCIB_Bernoulli: ESCIB_Bernoulli.o $(OBJS)
	$(GCC) -o ../$@ $+

ESCIB_Poisson: ESCIB_Poisson.o $(OBJS)
	$(GCC) -o ../$@ $+

DBSCAN: DBSCAN.o $(OBJS)
	$(GCC) -o ../$@ $+

clean: 
	rm -f ../ESCIB_Bernoulli ../ESCIB_Poisson ../DBSCAN *.o 

#Makefile of s2kit10 interface for Determination library
CC = gcc

FFTWINCLUDE ?= /usr/local/Cellar/fftw/3.3.5/include
FFTWDIR ?= /usr/local/Cellar/fftw/3.3.5/lib

FFTWINC = -I$(FFTWINCLUDE)
FFTWLIB = -L$(FFTWDIR) -lfftw3

s2kitdir = s2kit10

S2KITSRC = $(s2kitdir)/FST_semi_memo.c $(s2kitdir)/cospmls.c $(s2kitdir)/csecond.c \
	$(s2kitdir)/makeweights.c $(s2kitdir)/naive_synthesis.c $(s2kitdir)/pmls.c $(s2kitdir)/primitive.c \
	$(s2kitdir)/seminaive.c

CFLAGS = -O3 -fPIC $(FFTWINC)
LDFLAGS = -lm

all: s2kit_interface.c
	$(CC) -shared $(CFLAGS) $(S2KITSRC) s2kit_interface.c \
	$(FFTWLIB) $(LDFLAGS) -o libsh.so

clean:
	rm *.so

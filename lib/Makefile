##  Linux versions
# -lrt needed for clock_gettime(), Linux only
#      but only used in check.c
#  assumes GSL stored in regular paths
CFLAGS = -Wall -DNDEBUG
#CFLAGS = -Wall -DNDEBUG -O5
# CFLAGS = -g

SRC = sxtable.c stable.c digamma.c arms.c sapprox.c sslice.c \
	sampleb.c samplea.c yaps.c lgamma.c sympoly.c
HFILES = sxtable.h stable.h arms.h sapprox.h psample.h digamma.h \
	rng.h yaps.h lgamma.h sympoly.h
OBJ = $(SRC:.c=.o)
LIBRARY = libstb.a

library:  $(OBJ)
	ar rcs $(LIBRARY) $(OBJ)

clean: 
	rm -f *.o $(LIBRARY)

etags:  $(SRC) $(HFILES)
	etags  $(SRC) $(HFILES)

depend:
	makedepend $(FILES) 

# DO NOT DELETE THIS LINE -- make depend depends on it. 

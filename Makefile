MF=	Makefile

CC=	mpicc
CFLAGS= -fastsse

LFLAGS=	-lm

EXE=	main

SRC= \
	main.c \
	cio.c 

#
# No need to edit below this line
#

.SUFFIXES:
.SUFFIXES: .c .o

OBJ=	$(SRC:.c=.o)

.c.o:
	$(CC) $(CFLAGS) -c $<

all:	$(EXE)

$(EXE):	$(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ) $(LFLAGS)

$(OBJ):	$(MF)

clean:
	rm -f $(OBJ) $(EXE) core

#depend:
#       makedepend -f $(MAKEFILE) $(SRC)

# DO NOT DELETE


main.o: cio.h


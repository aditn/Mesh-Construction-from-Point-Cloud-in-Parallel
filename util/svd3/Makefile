PROG=svd3demo 
SRC=main.c svd3.c
OBJ=$(SRC:.c=.o)
#CFLAGS=-std=c99 -DMATLAB_FMT
CFLAGS=-std=c99
LDFLAGS=-lm

.SUFFIXES : .c .o

$(PROG) : $(OBJ)
	$(CC) -o $(PROG) $(OBJ) $(LDFLAGS)

.c.o :
	$(CC) -c $(CFLAGS) $<

clean :
	rm -f $(OBJ)


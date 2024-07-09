CC = mpicc
CFLAGS = -Wall
LIBS = -lm -lpthread
MAIN = main.c
NODES = nodes.c
BASE = base.c
all: run

run: $(MAIN) $(NODES) $(BASE)
	$(CC) $(CFLAGS) $^ $(LIBS) -o $@


clean:
	rm -f run

.PHONY: all clean
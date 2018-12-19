CC = gcc # c- compiler
CFLAGS = # compiler flags such as -O3 -W
LIBS = -lm -lgsl -lgslcblas # libraries such as -lm -lgsl -lgslcblas
OBJECTS =  main1.o # object files to produce
PROGRAM = main1 # name of executable to be built
# build .o files
%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS)
# compile program
all : $(PROGRAM)
$(PROGRAM): $(OBJECTS)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

clean:
		rm -f *.o
		touch *.c
		rm -rf *.dat

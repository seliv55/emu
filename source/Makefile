CFLAGS = -static -O3 -c -I../include

SRC = source.cpp

OBJECTS = $(SRC:.cpp=.o)

CEL = source.a
all: $(CEL)
$(CEL): $(OBJECTS)
	rm -f $(CEL)
	ar r $(CEL) *.o

.cpp.o:
	g++ $(CFLAGS) $< -o $@


clean:
	rm -f *.o *~ $(CEL)


CFLAGS = -static -O3 -c -I../include -I../integrT

SRC = nv.cpp solvers.cpp readexp.cpp

OBJECTS = $(SRC:.cpp=.o)

CEL = code.a

all: $(CEL)

$(CEL): $(OBJECTS)
	rm -f $(CEL)
	ar r $(CEL) *.o

.cpp.o:
	g++ $(CFLAGS) $< -o $@


clean:
	rm -f *.o *~ $(CEL)


CFLG = -c -O3 -Iinclude

LFLG = -O3 -lgfortran

SRC = rmod.cpp

LIBS = source/source.a

COBJ = $(SRC:.cpp=.o)

CEL = emu.out

all: $(CEL)

$(CEL): $(COBJ) $(LIBS)
	g++ -std=c++11 $(COBJ) $(LIBS) $(LFLG) -o $(CEL)

.cpp.o:
	g++ -std=c++11 $(CFLG) $< -o $@

source/source.a:
	make -C source

clean:
	rm -f emu.out *.o *~
	make clean -C source


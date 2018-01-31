CFLG = -c -O3 -Iinclude -IintegrT

LFLG = -O3 -lgfortran

SRC = main.cpp asdf.cpp

LIBS = code/code.a integrT/integT.a fort/fort.a nrused/nr.a

COBJ = $(SRC:.cpp=.o)

CEL = emu.out

all: $(CEL)

$(CEL): $(COBJ) $(LIBS)
	g++ $(COBJ) $(LIBS) $(LFLG) -o $(CEL)

.cpp.o:
	g++ $(CFLG) $< -o $@

code/code.a:
	make -C code

integrT/integT.a:
	make -C integrT

clean:
	rm -f emu.out *.o *~
	make clean -C code
#	make clean -C nrused
#	make clean -C integrT
cleanc:
	make clean -C code





INCLUDE   = -Iinclude
LIBS      = 
CXXFLAGS  = -c
	


COMPILER_CXX = mpicxx -g -O2
COMPILER_CC  = mpicc -g -O2
# PETSC library
INCLUDE      += -I/usr/local/lib/petsc/3.5.2/include
LIBS         += -L/usr/local/lib/petsc/3.5.2/lib -lpetsc
# Mesquite library
INCLUDE      += -I/home/camata/local/include
LIBS         += -L/home/camata/local/lib -lmesquite
# GTS
INCLUDE      += -I/usr/include/glib-2.0 -I/usr/lib/x86_64-linux-gnu/glib-2.0/include
LIBS += -lgts -glib-2.0



# CFLAGS  Must have INCLUDES to get  Implicit Rules working
CXXFLAGS    += $(INCLUDE)

headerfiles 	:= $(wildcard *.h)
srccfiles       := $(wildcard src/*.cpp)
objects         := $(patsubst %.cpp,  %.o, $(srccfiles))
testfiles       := $(wildcard tests/*.cpp)
testobjects     := $(patsubst %.cpp, %.o, $(testfiles))



#
# How to compile C
#
%.o : %.cpp
	@echo "Compiling C "$<"..."
	$(COMPILER_CXX) $(CXXFLAGS) $< -o $@

all: sedmesher

sedmesher: main.cpp $(objects) 
	$(COMPILER_CXX) $(CXXFLAGS) $< -o main.o $(DEFS)
	$(COMPILER_CXX) -o sedmesher main.o $(objects) ${LIBS}
clean:
	rm $(objects)
	rm test_mesh.o
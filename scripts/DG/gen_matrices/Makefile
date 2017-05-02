.DEFAULT_GOAL=all

files=gl_integration.f90 lu.f90 dubiner.f90 elementary.f90 gen_matrices.f90 main.f90

out=gen_matrices

TOP_DIR := $(shell pwd)
SRC_DIR=$(TOP_DIR)/src
TEST_DIR=$(TOP_DIR)/tests

sources=$(addprefix $(SRC_DIR)/,$(files))
objects=$(sources:.f90=.o)

#LIB=-llapack -lblas

comp = gfortran
cflags = -c $(LIB) -J$(SRC_DIR) -ffree-line-length-none -fimplicit-none
ldflags = -g $(LIB) -ffree-line-length-none -fimplicit-none

all : $(out)


$(out):	$(objects) test
	$(comp) $(ldflags) $(objects) -o $(out)


%.o: 	%.f90
	$(comp) $(cflags) $^ -o $@

test: 	
	make -C $(TEST_DIR)

clean: 
	rm -f $(out)
	rm -f $(SRC_DIR)/*.o
	rm -f $(SRC_DIR)/*.mod
	rm -f $(TEST_DIR)/*.o
	rm -f $(TEST_DIR)/*.mod
	make -C $(TEST_DIR) clean


export SRC_DIR
export TEST_DIR
export LIB
export files




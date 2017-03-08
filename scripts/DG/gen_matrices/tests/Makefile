.DEFAULT_GOAL=test

ifndef TOP_DIR
TOP_DIR := ..
endif

ifndef SRC_DIR
SRC_DIR=$(TOP_DIR)/src
endif

ifndef TEST_DIR
TEST_DIR=.
endif

$(info $(TEST_DIR))

source_files=$(filter-out main.f90,$(files))

test_files=unit.f90 test_defs.f90
test_main=test_main.f90

source_objects_temp=$(source_files:.f90=.o)
source_objects=$(addprefix $(SRC_DIR)/,$(source_objects_temp))

test_objects_temp=$(test_files:.f90=.o)
test_objects=$(addprefix $(TEST_DIR)/,$(test_objects_temp))


comp = gfortran
cflags = -c -fimplicit-none -ffree-line-length-none
ldflags = -g $(LIB) -I$(SRC_DIR) -fimplicit-none -ffree-line-length-none

TEST_EXE=test.x


%.o: 	%.f90
	$(comp) $(cflags) $^ -o $@

test: 	$(TEST_EXE)
	./$(TEST_EXE)

$(TEST_EXE): $(test_objects) $(objects)
	$(comp) -o $@ $(test_main) $(ldflags) $(test_objects) $(source_objects)


clean: 
	rm -f *.o
	rm -f *.mod
	rm -f $(TEST_EXE)






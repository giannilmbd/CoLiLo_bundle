.PHONY: all clean

# Define the default target
all: test_fit_ls.out

# Define the targets and their dependencies
test_fit_ls.out: fit_ls.o my_kinds_mod.mod test_fit_ls.f90
	ifort -O3 -r8 fit_ls.o test_fit_ls.f90 -o test_fit_ls.out -llapack

fit_ls.o: fit_ls.f90 my_kinds_mod.mod
	ifort -c fit_ls.f90

my_kinds_mod.mod: ../my_kinds_mod.f90
	ifort -c ../my_kinds_mod.f90

# Define a clean target to remove generated files
clean:
	rm -f *.out *.mod *.o

#####################################
## libraries and compilers
## update to your installation
####################################
# NETCDF INC and LIB dirs
NETCDF_INC=/usr/local/include
NETCDF_LIB=/usr/local/lib

# fortrangis INC and LIB dir
PROJ_INC=/usr/local/include
PROJ_LIB=/usr/local/lib

# fortran compiler
F90=gfortran
# fortran flag
FFLAGS=-fcheck=all

# source directory
SRC=../elmersrc

#######################################
LIBS= NetcdfUGRIDOutputSolver USFIcyMask 

all: $(LIBS)

%:$(SRC)/%.F90
	elmerf90  $(FFLAGS)  $^ -o $@

NetcdfUGRIDOutputSolver: $(SRC)/NetcdfUGRIDOutputSolver.F90
	elmerf90 $(FFLAGS) -DHAVE_PROJ -I$(NETCDF_INC)  -I$(PROJ_INC) $^ -o $@ -L$(NETCDF_LIB) -lnetcdf -lnetcdff -L$(PROJ_LIB) -lproj -lfortranc -lfortrangis

clean:
	rm -f $(LIBS)


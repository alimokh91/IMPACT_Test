ifeq ($(machine),artorg)
  ALL = local
else 
ifeq ($(machine),daint)
  ALL = xc30
else
  ALL = portable
endif
endif


execute:	$(ALL)

# target for local system (PC)
local:	COMP = mpifort
local:	INCL = -I/usr/include -I/opt/hdf5_gcc_native_shared/include -I/opt/openmpi_gcc_native_shared/include
local:	LIBS = -L/usr/lib/x86_64-linux-gnu -lm -ldl -llapack -L/opt/openmpi_gcc_native_shared/lib -lmpi_mpifh -lmpi_usempif08 -L/opt/hdf5_gcc_native_shared/lib -lhdf5_fortran -lhdf5hl_fortran -lhdf5_hl -lhdf5 -lz -lm -ldl
local:	OPT1 = -ffree-form -ffree-line-length-none -cpp -fdefault-real-8 -fdefault-double-8 \
-O3 -ffast-math -funroll-all-loops -fpeel-loops -ftracer -funswitch-loops \
#-ftree-vectorize -march=native -funsafe-math-optimizations -fomit-frame-pointer -fno-math-errno -fno-stack-limit 
	#-fopenmp # -J$(DST) 
local:	OPT2 = $(OPT1) -xf95-cpp-input
local:	$(BUILD_TARGET)
#local:  ftopy

# target for laptop
portable: COMP = mpifort
portable: INCL = -I$(HDF5_DIR)/include -I$(MPI_HOME)/include
portable: LIBS = -L/usr/lib -lm -ldl -lz -lblas -llapack -L$(MPI_HOME)/lib -L$(HDF5_DIR)/lib -lmpifort -lhdf5_fortran -lhdf5_hl -lhdf5
portable:       OPT1 = -ffree-form -ffree-line-length-none -cpp -fdefault-real-8 -O2
#-O3 -ffast-math -funroll-all-loops -fpeel-loops -ftracer -funswitch-loops \
#-ftree-vectorize -funsafe-math-optimizations -fomit-frame-pointer -fno-math-errno
portable:       OPT2 = $(OPT1) #-xf95-cpp-input
portable:       $(BUILD_TARGET)
#portable: f2py


# target for laptop
# target for Cray XC30 (Piz Daint, 5272 8-core SandyBridge 64-bit CPU compute nodes, 5272 NVidia Tesla K20X w/ 6GB memory)
xc30:	COMP = ftn
xc30:	INCL = -I$(HDF5_DIR)/include 
xc30:	LIBS = -L$(HDF5_DIR)/lib -lhdf5_fortran -lhdf5_hl -lhdf5 -lz -lm -ldl
xc30:	OPT1 = -ffree-form -ffree-line-length-none -cpp -fdefault-real-8 \
	-g -O2
#	-O3 -ffast-math -funroll-all-loops -fpeel-loops -ftracer -funswitch-loops \
#	-ftree-vectorize -funsafe-math-optimizations -fomit-frame-pointer -fno-math-errno
xc30:	OPT2 = $(OPT1) -xf95-cpp-input
xc30:	$(BUILD_TARGET)
#xc30:	f2py

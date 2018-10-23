
set -x

# Loading ncep environment
module load ncep/1.0

# Loading Intel Compiler Suite
module load PrgEnv-intel/5.2.56

module swap craype-sandybridge craype-haswell

module use /gpfs/hps/nco/ops/nwprod/lib/modulefiles

# Loading IOBUF
module load iobuf

module load cce/8.3.10
module load craype/2.3.0
module load cray-libsci/13.0.3

# Loading nceplibs modules
module load jasper-gnu-haswell/1.900.1
module load png-intel-haswell/1.2.49
module load zlib-gnu-haswell/1.2.7

module load bacio-cray-haswell/2.0.1 
module load w3emc-cray-haswell/2.2.0
module load w3nco-cray-haswell/2.0.6

module list

make -f Makefile_iobuf_cray

: "${PREFIX:=$HOME/gkylsoft}"
: "${MPI_HOME:=$HOME/gkylsoft/openmpi}"
./configure CC=cc --prefix=$PREFIX --use-lua=yes --use-mpi=yes

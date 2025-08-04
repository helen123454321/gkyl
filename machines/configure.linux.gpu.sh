: "${PREFIX:=$HOME/gkylsoft}"
./configure CC=nvcc --prefix=$PREFIX --use-lua=yes CUDA_ARCH=90 --cudamath-libdir=$LD_LIBRARY_PATH

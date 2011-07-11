export GAUSS_DIR=/home/tgraen/Software/g03
export GAUSS_EXEDIR=$GAUSS_DIR
export GAUSS_EXE=g03
export g03root=/home/tgraen/Software
export GAUSS_SCRDIR=$TMPDIR
export DEVEL_DIR=/home/tgraen/Software/g03/devel
export TMPDIR=$TMPDIR
export LD_LIBRARY_PATH=$GAUSS_DIR:/usr/lib64
export LD_LIBRARY_PATH=/home/tgraen/Software/pgi/linux86-64/7.2-1/libso:$GAUSS_DIR:/usr/lib64
sh /usr/local/pgi/pgi.sh

/usr/local/g03/g03 <mol.com >mol.out


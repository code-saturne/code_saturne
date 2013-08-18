#/bin/sh

ARCH=/home/B43051/saturne/arch

SALOME=/mingw/opt/salome-7.2/MODULES

if test -d "$SALOME" ; then

  OMNIORB=/mingw/opt/salome-7.2/PRODUCTS/omniORB-4.1.6
  export PATH=$OMNIORB/bin/x86_win32:$PATH
  export PYTHONPATH=$OMNIORB/bin/x86_win32:$PATH
  export PYTHONPATH=$OMNIORB/bin/x86_win32:$OMNIORB/lib/python/omniidl_be:$PYTHONPATH

  WITH_SALOME="--with-salome-kernel=$SALOME/KERNEL/RELEASE/KERNEL_INSTALL --with-salome-gui=$SALOME/GUI/RELEASE/GUI_INSTALL OMNIIDL=$OMNIORB/bin/x86_win32/omniidl.exe OMNIIDLPYBE=$OMNIORB/lib/python/omniidl_be SALOMEENVCMD=:"

fi


../configure \
  --prefix=$ARCH \
  --without-blas \
  --without-adf --without-ccm \
  --without-scotch \
  --disable-mpi-io \
  --enable-silent-rules \
  --enable-debug \
  --enable-relocatable $WITH_SALOME

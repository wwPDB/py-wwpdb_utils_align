#!/bin/bash
#
LIBDIR=$(python -c "import distutils.sysconfig as sysconfig; import os; print(os.path.join(sysconfig.get_config_var('LIBDIR'), sysconfig.get_config_var('LDLIBRARY')))")
INCDIR=$(python -c "from distutils.sysconfig import get_python_inc; print(get_python_inc())")
#
echo "Using INCDIR=${INCDIR}"
echo "Using LIBDIR=${LIBDIR}"
#
rm -rf build
mkdir build
cd build
#cmake ..
cmake -DPYTHON_LIBRARY=${LIBDIR}  -DPYTHON_INCLUDE_DIR=${INCDIR} ..
#
make
#
cd lib
ls -la
python -c "import alignlib" | grep "Symbol not found"
#
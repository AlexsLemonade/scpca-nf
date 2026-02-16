#!/bin/bash
set -euo pipefail

# Install bcl2fastq per illumina docs
BASE=${PWD}
SOURCE=${BASE}/bcl2fastq
BUILD=${BASE}/bcl2fastq2-build 
INSTALL_DIR=/usr/local

# cmake won't find libraries without this
export C_INCLUDE_PATH=/usr/include/x86_64-linux-gnu

tar -xzf bcl2fastq2.tar.gz

mkdir -p ${BUILD}
cd ${BUILD}
chmod ugo+x ${SOURCE}/src/configure
chmod ugo+x ${SOURCE}/src/cmake/bootstrap/installCmake.sh 
${SOURCE}/src/configure --prefix=${INSTALL_DIR}

cd ${BUILD}
make
make install


# clean up
cd ${BASE}
rm -rf ${SOURCE}
rm -rf ${BUILD}

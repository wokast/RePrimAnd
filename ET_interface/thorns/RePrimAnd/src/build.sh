#! /bin/bash

################################################################################
# Build
################################################################################

# Set up shell
if [ "$(echo ${VERBOSE} | tr '[:upper:]' '[:lower:]')" = 'yes' ]; then
    set -x                      # Output commands
fi
set -e                          # Abort on errors



# Set locations
THORN=RePrimAnd
NAME=RePrimAnd
SRCDIR="$(dirname $0)"
BUILD_DIR=${SCRATCH_BUILD}/build/${THORN}
if [ -z "${REPRIMAND_INSTALL_DIR}" ]; then
    INSTALL_DIR=${SCRATCH_BUILD}/external/${THORN}
else
    echo "BEGIN MESSAGE"
    echo "Installing REPRIMAND into ${REPRIMAND_INSTALL_DIR} "
    echo "END MESSAGE"
    INSTALL_DIR=${REPRIMAND_INSTALL_DIR}
fi
DONE_FILE=${SCRATCH_BUILD}/done/${THORN}
REPRIMAND_DIR=${INSTALL_DIR}

# Set up environment
#unset LIBS
#if echo '' ${ARFLAGS} | grep 64 > /dev/null 2>&1; then
#    export OBJECT_MODE=64
#fi

echo "REPRIMAND: Preparing directory structure..."
cd ${SCRATCH_BUILD}
mkdir build external done 2> /dev/null || true
rm -rf ${BUILD_DIR} ${INSTALL_DIR}
mkdir ${BUILD_DIR} ${INSTALL_DIR}

echo "REPRIMAND: Unpacking archive..."
pushd ${BUILD_DIR} >/dev/null
${TAR?} xf ${SRCDIR}/../dist/${NAME}.tar

echo "REPRIMAND: Configuring..."
cd ${NAME}/library

echo "REPRIMAND: Building..."
${MAKE} -f ${SRCDIR}/../dist/Makefile DESTDIR=${REPRIMAND_DIR}

echo "REPRIMAND: Installing..."
${MAKE} -f ${SRCDIR}/../dist/Makefile DESTDIR=${REPRIMAND_DIR} install
popd >/dev/null

echo "REPRIMAND: Cleaning up..."
rm -rf ${BUILD_DIR}

date > ${DONE_FILE}
echo "REPRIMAND: Done."

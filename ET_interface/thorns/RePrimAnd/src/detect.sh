#! /bin/bash

################################################################################
# Prepare
################################################################################

# Set up shell
if [ "$(echo ${VERBOSE} | tr '[:upper:]' '[:lower:]')" = 'yes' ]; then
    set -x                      # Output commands
fi
set -e                          # Abort on errors

. $CCTK_HOME/lib/make/bash_utils.sh

# Take care of requests to build the library in any case
REPRIMAND_DIR_INPUT=$REPRIMAND_DIR
if [ "$(echo "${REPRIMAND_DIR}" | tr '[a-z]' '[A-Z]')" = 'BUILD' ]; then
    REPRIMAND_BUILD=yes
    REPRIMAND_DIR=
else
    REPRIMAND_BUILD=
fi


THORN=RePrimAnd

################################################################################
# Build
################################################################################

if [ -n "$REPRIMAND_BUILD" -o -z "${REPRIMAND_DIR}" ]; then
    echo "BEGIN MESSAGE"
    echo "Using bundled REPRIMAND..."
    echo "END MESSAGE"
    
    check_tools "tar"

    # Set locations
    BUILD_DIR=${SCRATCH_BUILD}/build/${THORN}
    if [ -z "${REPRIMAND_INSTALL_DIR}" ]; then
        INSTALL_DIR=${SCRATCH_BUILD}/external/${THORN}
    else
        echo "BEGIN MESSAGE"
        echo "Installing REPRIMAND into ${REPRIMAND_INSTALL_DIR} "
        echo "END MESSAGE"
        INSTALL_DIR=${REPRIMAND_INSTALL_DIR}
    fi
    REPRIMAND_BUILD=1
    REPRIMAND_DIR=${INSTALL_DIR}
    REPRIMAND_INC_DIRS="$REPRIMAND_DIR/include"
    REPRIMAND_LIB_DIRS="$REPRIMAND_DIR/lib"
    REPRIMAND_LIBS="RePrimAnd"
else
    REPRIMAND_BUILD=
    DONE_FILE=${SCRATCH_BUILD}/done/${THORN}
    if [ ! -e ${DONE_FILE} ]; then
        mkdir ${SCRATCH_BUILD}/done 2> /dev/null || true
        date > ${DONE_FILE}
    fi
fi

################################################################################
# Configure Cactus
################################################################################

# Pass configuration options to build script
echo "BEGIN MAKE_DEFINITION"
echo "REPRIMAND_BUILD       = ${REPRIMAND_BUILD}"
echo "REPRIMAND_INSTALL_DIR = ${REPRIMAND_INSTALL_DIR}"
echo "END MAKE_DEFINITION"

# Set default options if not all variables were set by user input
if [ "${REPRIMAND_DIR}" != 'NO_BUILD' ]; then
    : ${REPRIMAND_INC_DIRS="${REPRIMAND_DIR}/include"}
    : ${REPRIMAND_LIB_DIRS="${REPRIMAND_DIR}/lib"}
fi
: ${REPRIMAND_LIBS='RePrimAnd'}

REPRIMAND_LIB_DIRS="$(${CCTK_HOME}/lib/sbin/strip-libdirs.sh $REPRIMAND_LIB_DIRS)"
REPRIMAND_INC_DIRS="$(${CCTK_HOME}/lib/sbin/strip-incdirs.sh $REPRIMAND_INC_DIRS)"

# Pass options to Cactus
echo "BEGIN MAKE_DEFINITION"
echo "REPRIMAND_DIR      = ${REPRIMAND_DIR}"
echo "REPRIMAND_INC_DIRS = ${REPRIMAND_INC_DIRS}"
echo "REPRIMAND_LIB_DIRS = ${REPRIMAND_LIB_DIRS}"
echo "REPRIMAND_LIBS     = ${REPRIMAND_LIBS}"
echo "END MAKE_DEFINITION"

echo 'INCLUDE_DIRECTORY $(REPRIMAND_INC_DIRS)'
echo 'LIBRARY_DIRECTORY $(REPRIMAND_LIB_DIRS)'
echo 'LIBRARY           $(REPRIMAND_LIBS)'

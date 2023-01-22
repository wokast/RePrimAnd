#!/bin/sh

set -e

/opt/python/cp39-cp39/bin/python3 -m venv /opt/env_meson
.  /opt/env_meson/bin/activate
pip install meson ninja

export BOOST_ROOT=/usr/local
BUILD_DIR=/opt/build/reprimand
INST_DIR=/usr

mkdir -p $BUILD_DIR
cd /project
meson setup --prefix=$INST_DIR -Dbuild_tests=false -Dbuild_benchmarks=false -Dbuild_documentation=false $BUILD_DIR
cd $BUILD_DIR
ninja 
ninja install
deactivate

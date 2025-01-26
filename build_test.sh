#!/bin/bash

PROJECT_ROOT=$(dirname "$(realpath "$0")")
BUILD_DIR="$PROJECT_ROOT/build"
GTEST_EXEC="$BUILD_DIR/openbnsl_test"

RED='\033[31m'
GREEN='\033[32m'
YELLOW='\033[33m'
BLUE='\033[34m'
RESET='\033[0m'

if [ "$1" == "clean" ] && [ -d "$BUILD_DIR" ]; then
    echo -e "${YELLOW}Cleaning build directory...${RESET}"
    rm -rf "$BUILD_DIR"
fi

if [ ! -d "$BUILD_DIR" ]; then
    mkdir -p "$BUILD_DIR"
fi

cd "$BUILD_DIR"

echo -e "${BLUE}Configuring build with tests...${RESET}"
cmake -DBUILD_TESTS=ON ..

echo -e "${BLUE}Building...${RESET}"
make -j$(nproc)
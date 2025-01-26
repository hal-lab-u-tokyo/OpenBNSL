#!/bin/bash

PROJECT_ROOT=$(dirname "$(realpath "$0")")
BUILD_DIR="$PROJECT_ROOT/build"
GTEST_EXEC="$BUILD_DIR/openbnsl_test"

RED='\033[31m'
GREEN='\033[32m'
YELLOW='\033[33m'
BLUE='\033[34m'
RESET='\033[0m'

cd "$BUILD_DIR"

if [ -x "$GTEST_EXEC" ]; then
    echo -e "${BLUE}Running tests...${RESET}"
    "$GTEST_EXEC"
else
    echo -e "${RED}Error: Test executable not found at $GTEST_EXEC${RESET}"
    exit 1
fi
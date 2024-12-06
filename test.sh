#!/bin/bash

PROJECT_ROOT=$(dirname "$(realpath "$0")")
BUILD_DIR="$PROJECT_ROOT/build"
GTEST_EXEC="$BUILD_DIR/openbnsl_test"

# Backend test (GoogleTest)
if [ -x "$GTEST_EXEC" ]; then
    echo "Running GoogleTest..."
    "$GTEST_EXEC"
else
    echo "Error: GoogleTest executable not found at $GTEST_EXEC"
    exit 1
fi

# Frontend test (Pytest)
echo "Running Pytest..."
pytest

echo "All tests completed."

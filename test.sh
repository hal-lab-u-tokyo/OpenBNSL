#!/bin/bash

PROJECT_ROOT=$(dirname "$(realpath "$0")")
BUILD_DIR="$PROJECT_ROOT/build"
TEST_EXEC="$BUILD_DIR/openbn_test"

# Backend test (GoogleTest)
if [ -x "$TEST_EXEC" ]; then
    echo "Running GoogleTest..."
    "$TEST_EXEC"
else
    echo "Error: GoogleTest executable not found at $TEST_EXEC"
    exit 1
fi

# Frontend test (Pytest)
echo "Running Pytest..."
pytest

echo "All tests completed."

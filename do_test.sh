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
if [ -d "$PROJECT_ROOT/tests" ]; then
    echo "Running Pytest..."
    pytest --rootdir="$PROJECT_ROOT/tests"
else
    echo "No tests directory found for Pytest. Skipping..."
fi

echo "All tests completed."

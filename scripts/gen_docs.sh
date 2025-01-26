#!/bin/bash

SCRIPT_DIR=$(dirname "$(realpath "$0")")
PROJECT_ROOT=$(realpath "$SCRIPT_DIR/..")

DOXYGEN_DIR="$PROJECT_ROOT/doxygen"
if [ -d "$DOXYGEN_DIR" ]; then
    echo "Moving to Doxygen directory: $DOXYGEN_DIR"
    cd "$DOXYGEN_DIR" || { echo "Error: Failed to enter $DOXYGEN_DIR"; exit 1; }
    if [ -f "Doxyfile" ]; then
        echo "Generating Doxygen documentation..."
        doxygen Doxyfile
    else
        echo "Error: Doxyfile not found in $DOXYGEN_DIR"
        exit 1
    fi
else
    echo "Error: Doxygen directory not found at $DOXYGEN_DIR"
    exit 1
fi

cd "$PROJECT_ROOT" || { echo "Error: Failed to return to $PROJECT_ROOT"; exit 1; }

SPHINX_DIR="$PROJECT_ROOT/sphinx"
if [ -d "$SPHINX_DIR" ]; then
    echo "Moving to Sphinx directory: $SPHINX_DIR"
    cd "$SPHINX_DIR" || { echo "Error: Failed to enter $SPHINX_DIR"; exit 1; }
    echo "Generating Sphinx HTML documentation..."
    make html
else
    echo "Error: Sphinx directory not found at $SPHINX_DIR"
    exit 1
fi

echo "Documentation generation completed."

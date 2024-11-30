#!/bin/bash

PROJECT_ROOT=$(dirname "$(realpath "$0")")

BUILD_DIR="$PROJECT_ROOT/build"

if [ -d "$BUILD_DIR" ]; then
    echo "Cleaning existing build directory..."
    rm -rf "$BUILD_DIR"
fi
mkdir -p "$BUILD_DIR"

echo "Configuring the project with CMake..."
cd "$BUILD_DIR"
cmake ..

echo "Building the project..."
make -j$(nproc)

echo "Exporting PYTHONPATH..."
export PYTHONPATH="$BUILD_DIR:$PYTHONPATH"

echo "Build completed. PYTHONPATH has been set to: $BUILD_DIR"
echo "You can now run Python scripts that use the 'openbnsl' module."

echo "Testing Python import..."
python3 -c "import openbnsl; print('openbnsl module loaded successfully.')"

echo "Done."

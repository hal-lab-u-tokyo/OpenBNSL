#!/bin/bash

TARGET_DIRS=("src" "include" "gtest" "bindings")
FILE_EXTENSIONS=("*.c" "*.cpp" "*.cu" "*.h" "*.hpp")

PYTHON_TARGET_DIRS=("scripts" "pytest")
PYTHON_EXTENSION="*.py"

if [ ! -f .clang-format ]; then
  echo "Error: .clang-format file not found in the project root."
  exit 1
fi

for dir in "${TARGET_DIRS[@]}"; do
  if [ -d "$dir" ]; then
    echo "Formatting C/C++ files in directory: $dir"
    for ext in "${FILE_EXTENSIONS[@]}"; do
      find "$dir" -type f -name "$ext" -exec clang-format -i {} +
    done
  else
    echo "Warning: Directory '$dir' does not exist. Skipping."
  fi
done

if ! command -v black &> /dev/null; then
  echo "Error: black is not installed. Install it using 'pip install black'."
  exit 1
fi

if [ ! -f pyproject.toml ]; then
  echo "Error: pyproject.toml file not found in the project root for Black."
  exit 1
fi

for dir in "${PYTHON_TARGET_DIRS[@]}"; do
  if [ -d "$dir" ]; then
    echo "Formatting Python files in directory: $dir"
    find "$dir" -type f -name "$PYTHON_EXTENSION" -exec black {} +
  else
    echo "Warning: Directory '$dir' does not exist. Skipping."
  fi
done

echo "Formatting complete."

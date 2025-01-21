#!/bin/bash

TARGET_DIRS=("src" "include" "gtest" "bindings")
FILE_EXTENSIONS=("*.c" "*.cpp" "*.cu" "*.h" "*.hpp")

PYTHON_TARGET_DIRS=("scripts" "pytest" "modules")
PYTHON_EXTENSION="*.py"

CHECK_MODE=false
if [[ "$1" == "--check" ]]; then
  CHECK_MODE=true
fi

if [ ! -f .clang-format ]; then
  echo "Error: .clang-format file not found in the project root."
  exit 1
fi

for dir in "${TARGET_DIRS[@]}"; do
  if [ -d "$dir" ]; then
    echo "Processing C/C++ files in directory: $dir"
    for ext in "${FILE_EXTENSIONS[@]}"; do
      if $CHECK_MODE; then
        find "$dir" -type f -name "$ext" -exec clang-format --dry-run --Werror {} +
      else
        find "$dir" -type f -name "$ext" -exec clang-format -i {} +
      fi
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
    echo "Processing Python files in directory: $dir"
    if $CHECK_MODE; then
      find "$dir" -type f -name "$PYTHON_EXTENSION" -exec black --check {} +
    else
      find "$dir" -type f -name "$PYTHON_EXTENSION" -exec black {} +
    fi
  else
    echo "Warning: Directory '$dir' does not exist. Skipping."
  fi
done

if $CHECK_MODE; then
  echo "Format check complete. No issues found."
else
  echo "Formatting complete."
fi

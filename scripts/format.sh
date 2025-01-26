#!/bin/bash

SCRIPT_DIR=$(dirname "$(realpath "$0")")
PROJECT_ROOT=$(realpath "$SCRIPT_DIR/..")

TARGET_DIRS=("src" "include" "gtest" "bindings")
FILE_EXTENSIONS=("*.c" "*.cpp" "*.cu" "*.h" "*.hpp")

PYTHON_TARGET_DIRS=("benchmarks" "pytest" "modules")
PYTHON_EXTENSION="*.py"

CHECK_MODE=false
if [[ "$1" == "--check" ]]; then
  CHECK_MODE=true
fi

if [ ! -f "$PROJECT_ROOT/.clang-format" ]; then
  echo "Error: .clang-format file not found in the project root."
  exit 1
fi

for dir in "${TARGET_DIRS[@]}"; do
  TARGET_PATH="$PROJECT_ROOT/$dir"
  if [ -d "$TARGET_PATH" ]; then
    echo "Processing C/C++ files in directory: $TARGET_PATH"
    for ext in "${FILE_EXTENSIONS[@]}"; do
      if $CHECK_MODE; then
        find "$TARGET_PATH" -type f -name "$ext" -exec clang-format --dry-run --Werror {} +
      else
        find "$TARGET_PATH" -type f -name "$ext" -exec clang-format -i {} +
      fi
    done
  else
    echo "Warning: Directory '$TARGET_PATH' does not exist. Skipping."
  fi
done

if ! command -v black &> /dev/null; then
  echo "Error: black is not installed. Install it using 'pip install black'."
  exit 1
fi

if [ ! -f "$PROJECT_ROOT/pyproject.toml" ]; then
  echo "Error: pyproject.toml file not found in the project root for Black: $PROJECT_ROOT"
  exit 1
fi

for dir in "${PYTHON_TARGET_DIRS[@]}"; do
  TARGET_PATH="$PROJECT_ROOT/$dir"
  if [ -d "$TARGET_PATH" ]; then
    echo "Processing Python files in directory: $TARGET_PATH"
    if $CHECK_MODE; then
      find "$TARGET_PATH" -type f -name "$PYTHON_EXTENSION" -exec black --check {} +
    else
      find "$TARGET_PATH" -type f -name "$PYTHON_EXTENSION" -exec black {} +
    fi
  else
    echo "Warning: Directory '$TARGET_PATH' does not exist. Skipping."
  fi
done

if $CHECK_MODE; then
  echo "Format check complete. No issues found."
else
  echo "Formatting complete."
fi

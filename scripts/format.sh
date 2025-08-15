#!/bin/bash
set -euo pipefail

RED='\033[31m'
GREEN='\033[32m'
YELLOW='\033[33m'
BLUE='\033[34m'
RESET='\033[0m'

SCRIPT_DIR=$(dirname "$(realpath "$0")")
PROJECT_ROOT=$(realpath "$SCRIPT_DIR/..")

TARGET_DIRS=("src" "include" "gtest" "bindings")
FILE_EXTENSIONS=("*.c" "*.cpp" "*.cu" "*.h" "*.hpp")

PYTHON_TARGET_DIRS=("benchmarks" "pytest" "helpers")
PYTHON_EXTENSION="*.py"

CHECK_MODE=false
if [[ "${1:-}" == "--check" ]]; then
  CHECK_MODE=true
fi

if [ ! -f "$PROJECT_ROOT/.clang-format" ]; then
  echo -e "${RED}Error: .clang-format file not found in the project root.${RESET}"
  exit 1
fi

for dir in "${TARGET_DIRS[@]}"; do
  TARGET_PATH="$PROJECT_ROOT/$dir"
  if [ -d "$TARGET_PATH" ]; then
    echo -e "${BLUE}Processing C/C++ files in: $TARGET_PATH${RESET}"
    for ext in "${FILE_EXTENSIONS[@]}"; do
      if $CHECK_MODE; then
        find "$TARGET_PATH" -type f -name "$ext" -exec clang-format --dry-run --Werror {} + 2>/dev/null
      else
        find "$TARGET_PATH" -type f -name "$ext" -exec clang-format -i {} + 2>/dev/null
      fi
    done
  else
    echo -e "${YELLOW}  Warning: Directory '$TARGET_PATH' does not exist. Skipping.${RESET}"
  fi
done

if ! command -v black &> /dev/null; then
  echo -e "${RED}Error: black is not installed. Install it using 'pip install black'.${RESET}"
  exit 1
fi

if [ ! -f "$PROJECT_ROOT/pyproject.toml" ]; then
  echo -e "${RED}Error: pyproject.toml file not found in the project root.${RESET}"
  exit 1
fi

for dir in "${PYTHON_TARGET_DIRS[@]}"; do
  TARGET_PATH="$PROJECT_ROOT/$dir"
  if [ -d "$TARGET_PATH" ]; then
    echo -e "${BLUE}Processing Python files in: $TARGET_PATH${RESET}"
    if $CHECK_MODE; then
      black --check "$TARGET_PATH"
    else
      black "$TARGET_PATH"
    fi
  else
    echo -e "${YELLOW}  Warning: Directory '$TARGET_PATH' does not exist. Skipping.${RESET}"
  fi
done

if $CHECK_MODE; then
  echo -e "${GREEN}Format check complete.${RESET}"
else
  echo -e "${GREEN}Formatting complete.${RESET}"
fi
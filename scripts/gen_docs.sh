#!/bin/bash

RED='\033[31m'
GREEN='\033[32m'
YELLOW='\033[33m'
BLUE='\033[34m'
RESET='\033[0m'

SCRIPT_DIR=$(dirname "$(realpath "$0")")
PROJECT_ROOT=$(realpath "$SCRIPT_DIR/..")
DOCS_DIR="$PROJECT_ROOT/docs"

echo -e "${BLUE}Starting documentation generation...${RESET}"

DOXYGEN_DIR="$DOCS_DIR/doxygen"
if [ -d "$DOXYGEN_DIR" ] && [ -f "$DOXYGEN_DIR/Doxyfile" ]; then
    echo -e "${BLUE}Generating Doxygen documentation...${RESET}"
    (cd "$DOXYGEN_DIR" && doxygen Doxyfile) \
        && echo -e "${GREEN}Doxygen documentation generated successfully!${RESET}" \
        || echo -e "${RED}Failed to generate Doxygen documentation.${RESET}"
else
    echo -e "${RED}Error: Doxygen directory or Doxyfile not found.${RESET}"
fi

SPHINX_DIR="$DOCS_DIR/sphinx"
if [ -d "$SPHINX_DIR" ]; then
    echo -e "${BLUE}Generating Sphinx HTML documentation...${RESET}"
    (cd "$SPHINX_DIR" && make html) \
        && echo -e "${GREEN}Sphinx HTML documentation generated successfully!${RESET}" \
        || echo -e "${RED}Failed to generate Sphinx documentation.${RESET}"
else
    echo -e "${RED}Error: Sphinx directory not found.${RESET}"
fi

echo -e "${BLUE}Documentation generation process completed.${RESET}"


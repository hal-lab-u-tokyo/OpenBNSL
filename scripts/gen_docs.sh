#!/bin/bash
set -euo pipefail

# ===== Colors =====
RED='\033[31m'
GREEN='\033[32m'
YELLOW='\033[33m'
BLUE='\033[34m'
RESET='\033[0m'

# ===== Paths =====
SCRIPT_DIR=$(dirname "$(realpath "$0")")
PROJECT_ROOT=$(realpath "$SCRIPT_DIR/..")
DOCS_DIR="$PROJECT_ROOT/docs"
DOXYGEN_DIR="$DOCS_DIR/doxygen"
SPHINX_DIR="$DOCS_DIR/sphinx"

# ===== Helpers =====
die() { echo -e "${RED}Error: $*${RESET}"; exit 1; }
info() { echo -e "${BLUE}$*${RESET}"; }
ok() { echo -e "${GREEN}$*${RESET}"; }

trap 'echo -e "'"${RED}"'Failed at line $LINENO. Aborting.'"${RESET}"'" >&2' ERR

info "Starting documentation generation..."

# ===== Doxygen =====
command -v doxygen >/dev/null 2>&1 || die "doxygen is not installed."
[ -d "$DOXYGEN_DIR" ] || die "Doxygen directory not found: $DOXYGEN_DIR"
[ -f "$DOXYGEN_DIR/Doxyfile" ] || die "Doxyfile not found: $DOXYGEN_DIR/Doxyfile"

info "Generating Doxygen documentation..."
(
  cd "$DOXYGEN_DIR"
  doxygen Doxyfile
)
ok "Doxygen documentation generated successfully!"

# ===== Sphinx =====
command -v make >/dev/null 2>&1 || die "make is not installed."
[ -d "$SPHINX_DIR" ] || die "Sphinx directory not found: $SPHINX_DIR"

[ -f "$SPHINX_DIR/Makefile" ] || die "Sphinx Makefile not found: $SPHINX_DIR/Makefile"

info "Generating Sphinx HTML documentation..."
(
  cd "$SPHINX_DIR"
  make html
)
ok "Sphinx HTML documentation generated successfully!"

ok "Documentation generation process completed."

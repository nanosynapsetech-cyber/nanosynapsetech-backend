#!/usr/bin/env bash
# setup_turso.sh — Install libsql Python client in the correct environment
# Run from WSL in the backend directory

set -e

echo "=== NanoSynapse Turso Setup ==="
echo ""

# Find the Python being used to run uvicorn/fastapi
PYTHON_BIN=$(which python3 2>/dev/null || which python3.14 2>/dev/null || which python 2>/dev/null)
echo "Found Python: $PYTHON_BIN"
$PYTHON_BIN --version

echo ""
echo "=== Installing libsql-experimental ==="
$PYTHON_BIN -m pip install libsql-experimental --break-system-packages 2>/dev/null \
  || $PYTHON_BIN -m pip install libsql-experimental 2>/dev/null \
  || pip3 install libsql-experimental 2>/dev/null \
  || pip install libsql-experimental 2>/dev/null

echo ""
echo "=== Verifying installation ==="
$PYTHON_BIN -c "import libsql_experimental as libsql; print('libsql OK:', libsql.__version__ if hasattr(libsql, '__version__') else 'installed')"

echo ""
echo "=== Done! ==="
echo "Now update .env.local with your Turso URL and token, then restart the backend."

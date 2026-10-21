#!/usr/bin/env bash
# Build the SageMath environment used by the four Sage scripts cited in Paper II.
#
#   ./setup_sage_env.sh [ENV_NAME]        (default ENV_NAME = sage-nsb)
#
# Creates a conda environment with SageMath 10.9 from conda-forge (native build,
# no Rosetta needed on Apple silicon), then installs abelfunctions at a fixed
# public revision with the small compatibility patch in this folder.  Nothing is
# installed outside the new environment and a temporary directory.
#
# Requires: conda (miniconda, miniforge or mambaforge), git, a C compiler
# (Xcode command-line tools on macOS, gcc on Linux), network access to
# conda-forge, PyPI and GitHub.  Disk: about 4-5 GB.
set -euo pipefail

ENV_NAME="${1:-sage-nsb}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ABEL_REPO="https://github.com/abelfunctions/abelfunctions.git"
ABEL_REV="d98843e4d25cebe9dc2993c86f469776142c5820"   # master, 24 May 2026

command -v conda >/dev/null || { echo "conda not found"; exit 1; }
START=$SECONDS
step () { echo; echo "== [$1/5] $2   (elapsed $(( (SECONDS-START)/60 )) min $(( (SECONDS-START)%60 )) s)"; }
CONDA_BASE="$(conda info --base)"
ENV_PREFIX="$CONDA_BASE/envs/$ENV_NAME"

if [ -d "$ENV_PREFIX" ]; then
  step 1 "environment $ENV_NAME already exists: reusing it"
else
  step 1 "creating $ENV_NAME with SageMath 10.9 (conda-forge; the long step, about 4 GB to download)"
  conda create -y -n "$ENV_NAME" --override-channels -c conda-forge "sage=10.9"
fi
PY="$ENV_PREFIX/bin/python"
"$PY" -c "import sage.version as v; assert v.version.startswith('10.9'), v.version; print('Sage', v.version)"

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

step 2 "cloning abelfunctions at $ABEL_REV and applying the patch"
git clone -q "$ABEL_REPO" "$WORK/abelfunctions"
git -C "$WORK/abelfunctions" checkout -q "$ABEL_REV"
# Sage 10.x removed the is_LaurentSeries / is_PowerSeries(Ring) helpers.
git -C "$WORK/abelfunctions" apply "$HERE/abelfunctions_sage10.patch"

# Cython >= 3.1 generates C for this package that does not compile
# (undeclared __pyx_tp_new_vectorcall_* functions).  Build with Cython 3.0 taken
# from a throw-away directory, so the environment keeps its own Cython.
step 3 "building abelfunctions with Cython 3.0 (temporary, not installed in the environment)"
"$PY" -m pip install -q --target "$WORK/cy30" "cython>=3.0,<3.1"
( cd "$WORK/abelfunctions" && \
  PYTHONPATH="$WORK/cy30" PATH="$WORK/cy30/bin:$PATH" \
  "$PY" -m pip install --no-cache-dir --no-deps --no-build-isolation --progress-bar on . )
step 4 "abelfunctions installed"

# Import from a neutral directory, inside Sage (abelfunctions needs sage.all).
step 5 "checking the import inside Sage"
( cd "$WORK" && "$ENV_PREFIX/bin/sage" -c "import abelfunctions; from abelfunctions import RiemannTheta; print('abelfunctions OK')" )
echo "== done.  Run the checks with:"
echo "   $PY $HERE/run_sage_checks.py --env $ENV_NAME"

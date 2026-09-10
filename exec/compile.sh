#!/bin/bash
set -e

# Verify Basilisk environment
if [ -z "$BASILISK" ]; then
    echo "Error: BASILISK environment variable is not set."
    exit 1
fi

module load gcc

echo "Translating Basilisk C to C99 with MPI..."
qcc -source -disable-dimensions -grid=octree -D_MPI=1 bubble_in_wedge.c

echo "Compiling binary..."
mpicc -Wall -std=c99 -O2 -D_FORTIFY_SOURCE=2 -D_XOPEN_SOURCE=700 _bubble_in_wedge.c \
      -o bubble_in_wedge -L"$BASILISK/gl" -lglutils -lfb_tiny -lm

echo "Build complete: ./bubble_in_wedge"

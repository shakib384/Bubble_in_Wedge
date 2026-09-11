# Bubble_in_Wedge
# Confined Bubble Dynamics in a Wedge Channel Under Microgravity

This repository contains the computational fluid dynamics (CFD) solver and automation scripts used to simulate three-dimensional, capillary-driven bubble migration in a wedge geometry. The code is built on the [Basilisk](http://basilisk.fr/) framework and utilizes C and MPI for high-performance parallel execution.

This work supports research on phase-change phenomena, microgravity fluid systems, and multiphase flow transport.

## 🚀 Features
* **Custom Solver:** Extends the Navier-Stokes equations with Volume-of-Fluid (VOF) interface tracking.
* **Immersed Boundary Method (IBM):** Utilizes an in-house IBM-GCM module to accurately map the complex wedge geometry.
* **Adaptive Mesh Refinement (AMR):** Dynamically refines the octree grid based on gas volume fraction and velocity gradients to minimize computational cost.
* **HPC Ready:** Configured for parallel execution on Linux-based clusters using OpenMPI.
* **Automated Post-Processing:** Includes shell scripts for rapid FFmpeg video rendering of simulation dumps.

## 📂 Repository Structure
```text
├── exec/
│   ├── bubble_in_wedge.c    # Main Basilisk C solver
│   ├── compile.sh           # Build script (GCC, MPI)
│   ├── job.sh               # Slurm batch submission script
│   └── ibm/                 # Lab-developed Immersed Boundary headers
├── out/
│   ├── makemovie.sh         # FFmpeg script to render .mp4 animations
│   └── clean.sh             # Utility to clear simulation output dumps
└── README.md
🛠️ Prerequisites
Basilisk: (Configured for C99 and OpenMPI)

Compiler: GCC (via module load gcc on HPC clusters)

MPI: OpenMPI (mpicc and mpirun)

Media: FFmpeg (for animation rendering)

Note: The required in-house Immersed Boundary Method (IBM-GCM) headers are bundled in the exec/ibm/ directory to ensure standalone compilation without modifying the upstream Basilisk tree.

💻 Compilation and Execution
1. Build the Executable
Navigate to the execution directory and run the compile script:

Bash
cd exec
./compile.sh
2. Run the Simulation
For local testing, run via standard MPI:

Bash
mpirun -np 8 ./bubble_in_wedge [parameters]
For HPC cluster execution, submit the batch script:

Bash
sbatch job.sh
📊 Post-Processing
Outputs, including interface geometry .dat files and cross-sectional .png renders, are generated in the out/ directory.

To compile the image sequences into high-quality MP4 animations (H.264 codec):

Bash
cd ../out
./makemovie.sh
To safely wipe the output data before a fresh run:

Bash
cd ../out
./clean.sh
📚 Associated Publications
The physics, boundary conditions, and numerical methods implemented in this repository are detailed in the following works:

Ahmed, S., Tryggvason, G., and Ling, Y., "Migration and breakup of confined bubbles in a wedge under microgravity," Physical Review of Fluids, (Submitted) 2026.

Ahmed, S., Tryggvason, G., and Ling, Y., "Rise of a confined bubble in a wedge," International Journal of Multiphase Flow, 2026.

Ahmed, S., Tryggvason, G., and Ling, Y., "Bubble Motion in a Wedge Channel," APS Division of Fluid Dynamics Annual Meeting, 2025.

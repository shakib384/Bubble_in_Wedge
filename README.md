# Confined Bubble Dynamics in a Wedge Channel Under Microgravity

This repository contains the computational fluid dynamics (CFD) solver and automation scripts used to simulate three-dimensional, capillary-driven bubble migration in a wedge geometry. The code is built on the [Basilisk](http://basilisk.fr/) framework and utilizes C and MPI for high-performance parallel execution.

This work supports research on phase-change phenomena, microgravity fluid systems, and multiphase flow transport.

## 🚀 Features
* **Custom Solver:** Extends the Navier-Stokes equations with Volume-of-Fluid (VOF) interface tracking.
* **Immersed Boundary Method (IBM):** Utilizes an in-house IBM-GCM module to accurately map the complex wedge geometry.
* **Adaptive Mesh Refinement (AMR):** Dynamically refines the octree grid based on gas volume fraction and velocity gradients to minimize computational cost.
* **HPC Ready:** Configured for parallel execution on Linux-based clusters using OpenMPI.
* **Automated Post-Processing:** Includes shell scripts for rapid FFmpeg video rendering of simulation dumps.

## 🔬 Simulation Phases
The solver is designed to handle multiphase transitions, specifically for the study *"Migration and breakup of confined bubbles in a wedge under microgravity"*:
* **Phase 1 (Initialization & Equilibrium):** The bubble rises and reaches an equilibrium state within the confined wedge channel under baseline gravitational conditions.
* **Phase 2 (Microgravity Migration):** The simulation is restarted from the Phase 1 equilibrium dump file. Gravity is turned off (microgravity conditions), and capillary forces exclusively drive the bubble's migration and subsequent breakup out of the wedge.

## 📂 Repository Structure
```text
├── exec/
│   ├── bubble_in_wedge.c           # Main Basilisk C solver
│   ├── compile.sh                  # Build script (GCC, MPI)
│   ├── job_phase1_equilibrium.sh   # Slurm script for baseline gravity
│   ├── job_phase2_microgravity.sh  # Slurm script for capillary-driven restart
│   └── ibm/                        # Lab-developed Immersed Boundary headers
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
(To run Phase 2, ensure the Phase 1 dump file is present in the output directory and adjust your parameter flags accordingly before submission).

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

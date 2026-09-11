#!/bin/sh
#!/bin/bash
#SBATCH --job-name=bubble_wedge
#SBATCH --output=run_%j.out
#SBATCH --error=run_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --time=04:00:00
#SBATCH --partition=standard

module list
pwd
date

echo "Job working directory: $SLURM_SUBMIT_DIR"
cd $SLURM_SUBMIT_DIR


max_level=9
L=29.5e-3
t_out=0.0001
t_end=7.0001
rhod=1.24;
rhof=1614.;
mud=0.018e-3;
muf=1.25e-3;
sigma=16.2e-3;
Uf=0.0;
H=4.0828e-3;
femax=1e-4;
uemax=5e-5;
R0=2.05779e-3;
z0=14.35e-3;
slope=0.1361;
x0=0e-3;
film=0.1e-4;
vmax=0.01;
pmax=5;
init_grav=-9.81;
grav_plus=-9.81;
t_grav=7.0;
min_level=3;

mpirun -np 48 ./bubble_in_wedge $max_level $L $t_out $t_end $rhod $rhof $mud $muf $sigma $Uf $H $femax $uemax $R0 $z0 $slope $x0 $film $vmax $pmax $init_grav $grav_plus $t_grav $min_level 1>out 2>log


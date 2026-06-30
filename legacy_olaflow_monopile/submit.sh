#!/bin/bash
#SBATCH --job-name=monopile_wave
#SBATCH --output=output.log
#SBATCH --account=def-yongshen
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000M

module load openfoam/v2412
source $WM_PROJECT_DIR/etc/bashrc

blockMesh > blockMesh.log
snappyHexMesh -overwrite > snappyHexMesh.log
rm -fr 0
cp -r 0.org 0
setFields > setFields.log
decomposePar -force > decomposePar.log
mpirun -np 8 olaFlow -parallel > olaFlow.log

echo "Simulation complete."

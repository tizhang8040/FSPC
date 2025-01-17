# Experimental Gmsh

# export GMSH_DIR=/mnt/Data/Devs/Gmsh/build/bin
# export PYTHONPATH=${GMSH_DIR}/lib:${PYTHONPATH}
# export INCLUDE=${GMSH_DIR}/include:${INCLUDE}
# export PATH=${GMSH_DIR}/bin:${PATH}
# export LIB=${GMSH_DIR}/lib:${LIB}

# Path to the external library

export PYTHONPATH=${PWD}:${PYTHONPATH}
export PYTHONPATH=${PWD}/../Metafor_dev/oo_meta:${PYTHONPATH}
export PYTHONPATH=${PWD}/../Metafor_dev/MetaBin/bin:${PYTHONPATH}
export PYTHONPATH=${PWD}/../PFEM3D/build/bin:${PYTHONPATH}

# Path to the Python script

# export SCRIPT=${PWD}/examples/2D/carsherWall/main.py
# export SCRIPT=${PWD}/examples/2D/coolingDisk/main.py
# export SCRIPT=${PWD}/examples/2D/damBreak/main.py
# export SCRIPT=${PWD}/examples/2D/elasticFunnel/main.py
# export SCRIPT=${PWD}/examples/2D/flowContact/main.py
# export SCRIPT=${PWD}/examples/2D/hydroStatic/main.py
# export SCRIPT=${PWD}/examples/2D/pureConduction/main.py
export SCRIPT=${PWD}/examples/2D/rubberGate/main.py
# export SCRIPT=${PWD}/examples/2D/staticAxisym/main.py
# export SCRIPT=${PWD}/examples/2D/thermoSquare/main.py
# export SCRIPT=${PWD}/examples/2D/vonKarman/main.py

# export SCRIPT=${PWD}/examples/3D/coolingDisk/main.py
# export SCRIPT=${PWD}/examples/3D/crossFlow/main.py
# export SCRIPT=${PWD}/examples/3D/damBreak/main.py
# export SCRIPT=${PWD}/examples/3D/hydroStatic/main.py

#export SCRIPT=${PWD}/examples/3D/pureCond_3Plates/main.py
export SCRIPT=${PWD}/examples/3D/alfeweld_1Plate/main.py

# Clean output folder

rm -rf workspace
mkdir workspace
cd workspace

# Runs the code

export CPU_PER_PROC=4
export MKL_NUM_THREADS=${CPU_PER_PROC}
export OMP_NUM_THREADS=${CPU_PER_PROC}
export OPTION="-map-by node:PE=${CPU_PER_PROC}"
mpiexec ${OPTION} -n 2 python3 ${SCRIPT} -k ${CPU_PER_PROC} 2>&1 | tee workspace.txt

:: Set input file

set CPU_PER_PROC=4
set INPUT=examples\damNcompRBF\input.py

:: Clean output folder

rd /s /q workspace
mkdir workspace

:: Runs the code

set OPTION=-cores %CPU_PER_PROC%
set MKL_NUM_THREADS=%CPU_PER_PROC%
set OMP_NUM_THREADS=%CPU_PER_PROC%
mpiexec %OPTION% -n 2 python main.py -k %CPU_PER_PROC% %INPUT%
::mpiexec %OPTION% -n 2 python main.py -k 1 %INPUT%
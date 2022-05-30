rm -rf sol_file*txt
rm -rf poissonfortran2D_Neumann
make poissonfortran2D_Neumann 
mpirun -np 8 ./poissonfortran2D_Neumann  -ksp_type fgmres -ksp_monitor_short
cat sol_file*.txt > pressure.txt
rm -rf sol_file*.txt

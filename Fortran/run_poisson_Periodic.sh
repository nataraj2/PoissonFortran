rm -rf poissonfortran2D_Periodic
make poissonfortran2D_Periodic 
mpirun -np 8 ./poissonfortran2D_Periodic  -ksp_type fgmres -ksp_monitor_short
cat sol_file*.txt > pressure.txt
rm -rf sol_file*.txt

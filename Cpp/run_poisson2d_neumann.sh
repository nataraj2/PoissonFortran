rm -rf poisson2d_neumann
make poisson2d_neumann
mpirun -np 8 ./Poisson2D_Neumann -ksp_monitor_short
#-pc_type mg -pc_mg_type full -ksp_type fgmres -ksp_monitor_short -pc_mg_levels 3 -mg_coarse_pc_factor_shift_type nonzero
cat sol_file*.txt > pressure.txt
rm -rf sol_file*.txt

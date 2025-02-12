rm -rf test_matzerorows
make test_matzerorows
mpirun -np 4 ./test_matzerorows #-ksp_monitor_short
#-pc_type mg -pc_mg_type full -ksp_type fgmres -ksp_monitor_short -pc_mg_levels 3 -mg_coarse_pc_factor_shift_type nonzero

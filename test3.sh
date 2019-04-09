mkdir ./sim/MARR_trunc$1/unfolding


echo "Modifying cfg file"
sed "s+/n/home03/amirbitran/ShakhnovichLab/MCPU/sim/MARR/files/marr_0.100_Emin.pdb+/n/home03/amirbitran/ShakhnovichLab/MCPU/sim/MARR_trunc$1/files/marr_truncated$1.pdb+g" ./src_mpi/cfg_MARR_unfolding > ./src_mpi/cfg_MARR_t$1_unfolding
sed -i "s+/n/home03/amirbitran/ShakhnovichLab/MCPU/sim/MARR/unfolding/marr+/n/home03/amirbitran/ShakhnovichLab/MCPU/sim/MARR_trunc$1/unfolding/marr_t$1+g" ./src_mpi/cfg_MARR_t$1_unfolding
sed -i "s+/n/home03/amirbitran/ShakhnovichLab/MCPU/sim/MARR/files/MARR.triple+/n/home03/amirbitran/ShakhnovichLab/MCPU/sim/MARR_trunc$1/files/MARR.triple+g" ./src_mpi/cfg_MARR_t$1_unfolding
sed -i "s+/n/home03/amirbitran/ShakhnovichLab/MCPU/sim/MARR/files/MARR.sctorsion+/n/home03/amirbitran/ShakhnovichLab/MCPU/sim/MARR_trunc$1/files/MARR.sctorsion+g" ./src_mpi/cfg_MARR_t$1_unfolding
sed -i "s+/n/home03/amirbitran/ShakhnovichLab/MCPU/sim/MARR/files/MARR.sec_str+/n/home03/amirbitran/ShakhnovichLab/MCPU/sim/MARR_trunc$1/files/MARR.sec_str+g" ./src_mpi/cfg_MARR_t$1_unfolding


#The extension -i tells you to replace in the current file, rather than create a new output file
#In OsX, you need blank quotes after -i. In Odyssey, we do not
#We don't use -i the first time because we haven't created te new file yet




echo "Modifying run.sh"
sed -i "s+./cfg_MARR_unfolding+./cfg_MARR_t$1_unfolding+g" ./src_mpi/runMARR_unfolding.sh > ./src_mpi/runMARR_t$1_unfolding.sh
 


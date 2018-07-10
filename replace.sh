#! /bin/bash
echo "Modifying backbone.c"
sed "s+/n/home03/amirbitran/ShakhnovichLab/MCPU/sim/DHFR/replica/dhfr+/n/home03/amirbitran/ShakhnovichLab/MCPU/sim/DHFR_trunc$1/replica/dhfr_t$1+g" ./src_mpi/backbone.c > ./src_mpi/backbone_t$1.c

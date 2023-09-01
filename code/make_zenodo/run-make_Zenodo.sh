#!/bin/bash

suite=$1

#mkdir -p ../../data/Zenodo_archive/$suite/atm/
#python subset_atm.py $suite

ocean_at='/nesi/nobackup/nesi00442/thoma97p/cylc-run/'$suite'/share/data/History_Data/NEMOhist/archive_ready'

mkdir -p ../../data/Zenodo_archive/$suite/ocean_T/
for filename in $ocean_at/*_1m_*T.nc
do
echo ${filename:98}
python subset_Tgrid.py $filename ../../data/Zenodo_archive/$suite/ocean_T/${filename:98}
done

mkdir -p ../../data/Zenodo_archive/$suite/ocean_U/
for filename in $ocean_at/*_1m_*U.nc
do
echo ${filename:98}
python subset_Ugrid.py $filename ../../data/Zenodo_archive/$suite/ocean_U/${filename:98}
done

mkdir -p ../../data/Zenodo_archive/$suite/ocean_V/
for filename in $ocean_at/*_1m_*V.nc
do
echo ${filename:98}
python subset_Vgrid.py $filename ../../data/Zenodo_archive/$suite/ocean_V/${filename:98}
done

ice_at='/nesi/nobackup/nesi00442/thoma97p/cylc-run/'$suite'/share/data/History_Data/CICEhist/archive_ready'
mkdir -p ../../data/Zenodo_archive/$suite/ice
for filename in $ice_at/*_1m*.nc
do
echo ${filename:98}
python subset_CICE.py $filename ../../data/Zenodo_archive/$suite/ice/${filename:98}
done


mkdir -p ../../data/Zenodo_archive/$suite/atm/
python subset_atm.py $suite



#!/bin/bash


now="$(date +'%Y%m%d')"
results="../../results/"$now"/ssp585_ISMIPfw.py"
`mkdir -p $results`
`cp ssp585_ISMIPfw.py $results`
`cp run-ssp585_ISMIPfw-Apr23_resub.sh $results`
rm ISMIPfw*.png

# fw timeseries
python ssp585_ISMIPfw.py --analysis=fw_timeseries --suite=ci501,cm483,cn043,cn077
mv fw_timeseries.png $results 
mv ISMIPfw-data_for_fw_timeseries*.nc $results

# thetao 400m
python ssp585_ISMIPfw.py --variable=thetao --depth=400 --analysis=map --suite=ci501,cm483,cn043,cn077 --lats=-90,0 --abs_lims=-2.5,7.5 --diff_lims=-1.5,1.5 --plot_sia_contour=False
mv ISMIPfw-map-thetao.png $results/ISMIPfw-thetao-z400.png && echo "done thetao z400"
mv ISMIPfw-map-thetao-colorbars.png $results
mv ISMIPfw-data_for_map*.nc $results/ 

# uo 400m
python ssp585_ISMIPfw.py --variable=uo --depth=400 --analysis=map --suite=ci501,cm483,cn043,cn077 --lats=-90,0 --abs_lims=-0.1,0.1 --diff_lims=-0.1,0.1 --plot_sia_contour=False
mv ISMIPfw-map-uo.png $results/ISMIPfw-uo-z400.png && echo "done uo z400"
mv ISMIPfw-data_for_map*.nc $results/ 

# uo 200m
python ssp585_ISMIPfw.py --variable=uo --depth=200 --analysis=map --suite=ci501,cm483,cn043,cn077 --lats=-90,0 --abs_lims=-0.1,0.1 --diff_lims=-0.1,0.1 --plot_sia_contour=False
mv ISMIPfw-map-uo.png $results/ISMIPfw-uo-z200.png && echo "done uo z200"
mv ISMIPfw-data_for_map*.nc $results/

# uo 100m
python ssp585_ISMIPfw.py --variable=uo --depth=100 --analysis=map --suite=ci501,cm483,cn043,cn077 --lats=-90,0 --abs_lims=-0.1,0.1 --diff_lims=-0.1,0.1 --plot_sia_contour=False
mv ISMIPfw-map-uo.png $results/ISMIPfw-uo-z100.png && echo "done uo z100"
mv ISMIPfw-data_for_map*.nc $results/

# uo 0m
python ssp585_ISMIPfw.py --variable=uo --depth=0 --analysis=map --suite=ci501,cm483,cn043,cn077 --lats=-90,0 --abs_lims=-0.1,0.1 --diff_lims=-0.1,0.1 --plot_sia_contour=False
mv ISMIPfw-map-uo.png $results/ISMIPfw-uo-z0.png && echo "done uo z0"
mv ISMIPfw-map-uo-colorbars.png $results
mv ISMIPfw-data_for_map*.nc $results/

# u wind
python ssp585_ISMIPfw.py --variable=10m_wind_speed_u --depth=0 --analysis=map --suite=ci501,cm483,cn043,cn077 --lats=-90,0 --abs_lims=-10,10 --diff_lims=-2.5,2.5 --plot_0_wind=True --land_mask=True 
mv ISMIPfw-map-10m_wind_speed_u.png $results && echo "done u wind"
mv ISMIPfw-map-10m_wind_speed_u-colorbars.png $results
mv ISMIPfw-data_for_map*.nc $results/

python ssp585_ISMIPfw.py --analysis=mackie 
mv mackie_maps.png $results

python ssp585_ISMIPfw.py  --suite=ci501,cm483,cn043,cn077 --depths=0,1000 --variable=so --analysis=section --lats=-78.2,-70 --lons=176,178
var_dir=$results"/section/177E" && `mkdir -p $var_dir` && echo $var_dir
`cp ISMIPfw-section-so.png $var_dir`
`cp ISMIPfw-section-thetao.png $var_dir`
mv ISMIPfw-data_for_section*.nc $var_dir/

python ssp585_ISMIPfw.py  --suite=ci501,cm483,cn043,cn077 --depths=0,1000 --variable=so --analysis=section --lats=-76,-68 --lons=251,253        
var_dir=$results"/section/Thwaites" && `mkdir -p $var_dir` && echo $var_dir
`cp ISMIPfw-section-so.png $var_dir`
`cp ISMIPfw-section-thetao.png $var_dir`
mv ISMIPfw-data_for_section*.nc $var_dir/

python ssp585_ISMIPfw.py  --suite=ci501,cm483,cn043,cn077 --depths=0,1000 --variable=so --analysis=section --lats=-69,-65 --lons=70,72
var_dir=$results"/section/Amery" && `mkdir -p $var_dir` && echo $var_dir
`cp ISMIPfw-section-so.png $var_dir`
`cp ISMIPfw-section-thetao.png $var_dir`
mv ISMIPfw-data_for_section*.nc $var_dir/

python ssp585_ISMIPfw.py --analysis=ekman --lats=-90,90 --abs_lims=-10e-6,10e-6 --diff_lims=-1.5e-6,1.5e-6 --suite=ci501,cm483,cn043,cn077 > ekman_results.txt
var_dir=$results"/ekman/" && `mkdir -p $var_dir` && echo $var_dir
mv ekman_results.txt $var_dir/
mv ISMIPfw*ekman*.png $var_dir/
mv ISMIPfw-data_for_map-ekman*.nc $var_dir/

# total_precipitation
python ssp585_ISMIPfw.py --variable=total_precipitation --depth=0 --analysis=map --suite=ci501,cm483,cn043,cn077 --lats=-90,0 --abs_lims=0,0.5e-4 --diff_lims=-0.1e-4,0.1e-4 
mv ISMIPfw-map-total_precipitation.png $results/ISMIPfw-total_precipitation.png && echo "done precip"
mv ISMIPfw-map-total_precipitation-colorbars.png $results
mv ISMIPfw-data_for_map*.nc $results/

python ssp585_ISMIPfw.py --analysis=sidmassth --lats=-90,90 --abs_lims=-1.5e-4,1.5e-4 --diff_lims=-1.5e-4,1.5e-4 --suite=ci501,cm483,cn043,cn077 > sidmassth_results.txt
var_dir=$results"/sidmassth/" && `mkdir -p $var_dir` && echo $var_dir
mv sidmassth_results.txt $var_dir/
mv ISMIPfw*sidmassth*.png $var_dir/
mv ISMIPfw-data_for_map-sidmassth*.nc $var_dir/

python ssp585_ISMIPfw.py --analysis=sea_ice_area --lats=-90,90 --abs_lims=0,1 --diff_lims=-1,1 --suite=ci501,cm483,cn043,cn077 > sia_results.txt
var_dir=$results"/sia/" && `mkdir -p $var_dir` && echo $var_dir
mv sia_results.txt $var_dir/
mv ISMIPfw*soicecov*.png $var_dir/
mv ISMIPfw-data_for_map-sia*.nc $var_dir/

results_mackie=$results"/mackie_sections/" && `mkdir -p $var_dir` && echo $var_dir
rm ISMIPfw-section*.png

for suite in ay187 az576 bb819
    do
        eval "python ssp585_ISMIPfw.py  --suite=mackie,$suite --depths=0,1000 --variable=so --analysis=section --lats=-78.2,-70 --lons=176,178"
        var_dir=$results_mackie"/section/177E" && `mkdir -p $var_dir` && echo $var_dir
        cp ISMIPfw-section-thetao.png $var_dir/thetao-$suite".png"
        cp ISMIPfw-section-so.png $var_dir/so-$suite".png"

        eval "python ssp585_ISMIPfw.py  --suite=mackie,$suite --depths=0,1000 --variable=so --analysis=section --lats=-76,-68 --lons=251,253"
        var_dir=$results_mackie"/section/Thwaites" && `mkdir -p $var_dir` && echo $var_dir
        cp ISMIPfw-section-thetao.png $var_dir/thetao-$suite".png"
        cp ISMIPfw-section-so.png $var_dir/so-$suite".png"

        eval "python ssp585_ISMIPfw.py  --suite=mackie,$suite --depths=0,1000 --variable=so --analysis=section --lats=-69,-65 --lons=70,72"
        var_dir=$results_mackie"/section/Amery" && `mkdir -p $var_dir` && echo $var_dir
        cp ISMIPfw-section-thetao.png $var_dir/thetao-$suite".png"
        cp ISMIPfw-section-so.png $var_dir/so-$suite".png"
    done

python ssp585_ISMIPfw.py --analysis=latent_heat --lons=210,290 --lats=-90,-60 --depths=0,1000
python ssp585_ISMIPfw.py --analysis=opottemprmadvect --lons=210,290 --lats=-90,-60 --depths=0,1000
python ssp585_ISMIPfw.py --analysis=opt_lh_plot
var_dir=$results"/lh_opottemprmadvect/" && `mkdir -p $var_dir` && echo $var_dir
mv opottemprmadvect_latentHeat.png $var_dir
mv cmip-opot-spatial.nc $var_dir
mv maui-opot-spatial.nc $var_dir
mv lh_spatial-sum1000.0_lons210.0to290.0_shelf.nc $var_dir/lh_spatial-AmBe-sum1000-shelf.nc

python ssp585_ISMIPfw.py --analysis=input_regions
mv input_regions.png $results

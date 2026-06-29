#!/bin/bash
/bin/python3 plot_ttr_vs_time.py \
	--case "hB_case1_TTR" .  \
	--field TTR \
	--title "N2 V-T Relaxation" \
	--output n2_vt_ttr.png \
	--export-csv

/bin/python3 plot_ttr_vs_time.py \
	--case "hB_case1_TVib" .  \
	--field TVib \
	--title "N2 V-T Relaxation" \
	--output n2_vt_tvib.png \
	--export-csv

/bin/python3 plot_relaxation.py \
    --case "$ T_{tr} $ (nexaFoam)" hB_case1_TTR.csv \
    --case "$ T_{vib} $ (nexaFoam)" hB_case1_TVib.csv \
    --field "Temp" \
    --title "N2 V-T Relaxation Comparison" \
    --output n2_vt_relaxation.png \
    --export-csv

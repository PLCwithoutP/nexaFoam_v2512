#!/bin/bash
/bin/python3 plot_ttr_vs_time.py \
	--case "hB_case4_TTR" .  \
	--field TTR \
	--title "N2-N V-T & C-V Relaxation" \
	--output n2n_vt_cv_ttr.png \
	--export-csv

/bin/python3 plot_ttr_vs_time.py \
	--case "hB_case4_TVib" .  \
	--field TVib \
	--title "N2-N V-T & C-V Relaxation" \
	--output n2n_vt_cv_tvib.png \
	--export-csv

/bin/python3 plot_relaxation.py \
    --case "$ T_{tr} $ (nexaFoam)" hB_case4_TTR.csv \
    --case "$ T_{vib, N_2} $ (nexaFoam)" hB_case4_TVib.csv \
    --field "Temp" \
    --title "N2-N V-T & C-V Relaxation Comparison" \
    --output n2n_vt_cv_relaxation.png \
    --export-csv

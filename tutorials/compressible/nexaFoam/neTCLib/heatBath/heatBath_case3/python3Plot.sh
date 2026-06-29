#!/bin/bash
/bin/python3 plot_ttr_vs_time.py \
	--case "hB_case3_TTR" .  \
	--field TTR \
	--title "N2-O2 V-T & V-V Relaxation" \
	--output n2o2_vt_ttr.png \
	--export-csv

/bin/python3 plot_ttr_vs_time.py \
	--case "hB_case3_TVib_N2" .  \
	--field TVib_Y_N2 \
	--title "N2-O2 V-T & V-V Relaxation" \
	--output n2o2_vt_tvibn2.png \
	--export-csv

/bin/python3 plot_ttr_vs_time.py \
	--case "hB_case3_TVib_O2" .  \
	--field TVib_Y_O2 \
	--title "N2-O2 V-T & V-V Relaxation" \
	--output n2o2_vt_tvibo2.png \
	--export-csv

/bin/python3 plot_relaxation.py \
    --case "$ T_{tr} $ (nexaFoam)" hB_case3_TTR.csv \
    --case "$ T_{vib, N_2} $ (nexaFoam)" hB_case3_TVib_N2.csv \
    --case "$ T_{vib, O_2} $ (nexaFoam)" hB_case3_TVib_O2.csv \
    --field "Temp" \
    --title "N2-O2 V-T & V-V Relaxation Comparison" \
    --output n2o2_vt_relaxation.png \
    --export-csv

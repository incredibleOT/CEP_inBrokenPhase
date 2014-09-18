##expects only one file per kappa
fileName=""
function fileName_func { kappa_name=$( echo "${kappa}" | sed "s:-:m:g" | sed "s:\.:p:g" );  fileName=F_${kappa_name}; }
declare -a colornames=( "red"  "blue"  "orange"  "dark-green"  "magenta"  "green"  "grey"  "cyan"  "black"  "dark-violet"  )
##general numbers
lambda_6=0.01
LS=64
LT=128


# lambda=-0.0472
lambda=-0.0476

allDataFiles=$( ls ../potential*_l_${lambda}_k_*_scan*.txt )

plot_file=tmp_plot_pot.plt

declare -a all_kappas=( ` echo ${allDataFiles} | grep -oE "_k_[[:digit:].-]+" | sed "s:_k_::g" | sort -g `)


# # echo "${all_kappas[@]}"
# # for k_index in ${!all_kappas[@]}; do
# # 	echo "index: ${k_index}"
# # done
# # 
# # exit 0

if [[ ( -e ${plot_file} ) ]]; then
	rm ${plot_file}
fi
touch ${plot_file}

echo "reset" >> ${plot_file}
echo ""      >> ${plot_file}

for k_index in ${!all_kappas[@]}; do kappa=${all_kappas[$k_index]}; fileName_func ${kappa}  
	file_loc=$( ls ../potential*_l_${lambda}_k_${kappa}_scan*.txt )
	echo "${fileName}=\"${file_loc}\"" >> ${plot_file}
done
echo "" >> ${plot_file}

echo "set pointsize 0.5" >> ${plot_file}
echo "set key noautotitle" >> ${plot_file}
echo "set xlabel \"\$\\\\varphi\$\"" >> ${plot_file}
echo "set ylabel \"\$U_{CEP}(\\\\varphi)\$\"" >> ${plot_file}
echo "set title sprintf(\"CEP \$\\\\lambda_6=${lambda_6}\$, \$\\\\lambda=${lambda}\$, Lattice: \$${LS}^3 \\\\times ${LT}\$\")" >> ${plot_file}

echo "" >>${plot_file}


##stat all files
for k_index in ${!all_kappas[@]}; do kappa=${all_kappas[$k_index]}; fileName_func ${kappa}  
	echo "stats     ${fileName}   using (\$1):(\$2)   name \"info_${fileName}\"" >> ${plot_file}
done
##largest kappa determines bounds of plot
fileName_func ${all_kappas[${#all_kappas[@]}-1]}
echo "" >>${plot_file}
echo "os_y = info_${fileName}_min_y" >>${plot_file}
##stats for max y
echo "" >>${plot_file}
echo "stats     ${fileName}   using (\$2):(\$1)   name \"inv_${fileName}\"" >> ${plot_file}
echo "y_max = (inv_${fileName}_pos_min_y - os_y)" >>${plot_file}

echo "" >>${plot_file}
echo "set yrange [0.0 : y_max]" >>${plot_file}
echo "" >>${plot_file}

echo "set format y \"$%.e$\"" >>${plot_file}
echo "" >>${plot_file}

##mark lines at minima
for k_index in ${!all_kappas[@]}; do kappa=${all_kappas[$k_index]}; fileName_func ${kappa}  
	echo "set arrow from info_${fileName}_pos_min_y, 0.0    to info_${fileName}_pos_min_y, y_max nohead lt 2 lw 2 lc rgb \"${colornames[\k_index]}\"" >>${plot_file}
done


##write the plotting commands
echo -n "plot " >>  ${plot_file}
for k_index in ${!all_kappas[@]}; do kappa=${all_kappas[$k_index]}; fileName_func ${kappa}  
	if [[ k_index -eq 0 ]]; then echo -e "\\" >> ${plot_file}; else echo -e ",\\" >> ${plot_file}; fi
	echo -n "     ${fileName}   using (\$1):(\$2-os_y)   with points lc rgbcolor \"${colornames[\k_index]}\" title \"\$\\\\kappa=${kappa}\$\"" >> ${plot_file}
done
echo "" >> ${plot_file}

# echo "size of array: ${#all_kappas[@]}"
# 
# echo "last plus 1:   ${all_kappas[${#all_kappas[@]}-1]}"
# for kappa in ${all_kappas}; do
# 	echo "kappa:   ${kappa}"
# # 	kappa_name=$( echo "${kappa}" | sed "s:-:m:g" | sed "s:.:p:g" )
# 	fileName_func ${kappa}
# 	echo "fileName:   ${fileName}"
# done

#Script for launching multiple runs (a loop) of a Geant4 application.
#Runs can differ for random seed and for settings.


######################### INPUT ############################
t_array=($(seq 8.0 1.0 16.0)) #mm
#t_array=(12.0) #mm

#mis_array=(0.000 $(seq 0.001 0.001 0.010)) #rad
mis_array=(0.000) #rad

#pot_array=("0.050A", "0.065A")
pot_array=("0.050A")

binsize=0.5 #mm

energy='2.86GeV'

macro='macros/run_FS_crystalline.mac'
############################################################


#definition of functions
ceil() {
    local num=$1
    result=$(echo '(' $num '* 2 + 1 - 0.0) / 2' | bc)
    return $result
}

isEven() { 
	n=$1;
	if [ `expr $n % 2` == 0 ]
	then
		return 1;
	else
		return 0;
	fi
}


#elaboration
n_t=${#t_array[@]}
n_mis=${#mis_array[@]}
n_pot=${#pot_array[@]}

n_runs=$(($n_t * $n_mis * $n_pot))
echo "A loop of $n_runs runs will be executed"

Ncores=$(grep -c ^processor /proc/cpuinfo)
echo "number of available cores on this computer: $Ncores"

if [ ! -d output ] ; then
	mkdir output 
fi

count=0
for (( i=0; i<$n_t; i++ ))
do
	for (( j=0; j<$n_mis; j++ ))
	do
		for (( k=0; k<$n_pot; k++ ))
	    do
			angX=$(echo ${mis_array[$j]} '* 0.149438' | bc) #mm (angleX=mis*cos(atan(slopeYX=20/3)))
			angY=$(echo ${mis_array[$j]} '* 0.988770' | bc) #mm (angleY=mis*sin(atan(slopeYX=20/3)))
			echo Running simulation $(($count + 1)) at $energy with: crystal_th "=" ${t_array[$i]} mm, mis "=" 0${mis_array[$j]} rad [angleX "=" 0$angX rad, angleY "=" 0$angY rad], pot "=" ${pot_array[$k]}...
		
			sed -i "s+/crystal/setCrystalSize .*+/crystal/setCrystalSize 20. 20. ${t_array[$i]} mm+" $macro
			sed -i "s+/crystal/setCrystalAngleX .*+/crystal/setCrystalAngleX 0$angX+" $macro
			sed -i "s+/crystal/setCrystalAngleY .*+/crystal/setCrystalAngleY 0$angY+" $macro
			sed -i "s+/crystal/setPotentialPath .*+/crystal/setPotentialPath ../Potentials/W111/u1_${pot_array[$k]}/+" $macro
			sed -i "s+/run/setfilenamesave output/output_${energy}_W.*+/run/setfilenamesave output/output_${energy}_W${t_array[$i]}mm_crystalline_mis${mis_array[$j]}rad_pot${pot_array[$k]}.root+" $macro
			
			halfsize=$(echo ${t_array[$i]} '* 0.5' | bc) #mm
			sed -i "s+/score/mesh/boxSize .*+/score/mesh/boxSize 10. 10. $halfsize mm+" $macro
			var=$(echo ${t_array[$i]} '/' $binsize | bc)
			ceil $var
			int_var=$?
			isEven $int_var
			is_int_var_even=$?
			if [ $is_int_var_even == 1 ] ; then
				int_var=$(($int_var + 1))
			fi
			NbinZ=$int_var			
			sed -i "s+/score/mesh/nBin .*+/score/mesh/nBin 41 41 $NbinZ+" $macro
			sed -i "s+/score/mesh/translate/xyz 0. 0. .*+/score/mesh/translate/xyz 0. 0. -$halfsize mm+" $macro
				
			sed -i "s+/score/dumpQuantityToFile boxMesh1 Edep output/Edep_${energy}_W.*+/score/dumpQuantityToFile boxMesh1 Edep output/Edep_${energy}_W${t_array[$i]}mm_crystalline_mis${mis_array[$j]}rad_pot${pot_array[$k]}.txt+" $macro
			sleep 1
			
			#RUN YOUR MACRO HERE 
			./positronSource $macro $Ncores > output/output_${energy}_W${t_array[$i]}mm_crystalline_mis${mis_array[$j]}rad_pot${pot_array[$k]}.out

			count=$(($count + 1))
			echo -e "Run $count done!"
		done
	done
done

#Script for launching multiple runs (a loop) of a Geant4 application.
#Runs can differ for random seed and for settings.


######################### INPUT ############################
t_array=($(seq 5.0 1.0 20.0)) #mm
#t_array=(15.0) #mm

binsize=0.5 #mm

energy='2.86GeV'

macro='macros/run_conventional.mac'
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

n_runs=$(($n_t))
echo "A loop of $n_runs runs will be executed"

Ncores=$(grep -c ^processor /proc/cpuinfo)
echo "number of available cores on this computer: $Ncores"

if [ ! -d output ] ; then
	mkdir output 
fi

count=0
for (( i=0; i<$n_t; i++ ))
do
	echo Running simulation $(($count + 1)) at $energy with: volume_th "=" ${t_array[$i]} mm ...

	sed -i "s+/crystal/setCrystalSize .*+/crystal/setCrystalSize 20. 20. ${t_array[$i]} mm+" $macro
	sed -i "s+/run/setfilenamesave output/output_${energy}_W.*+/run/setfilenamesave output/output_${energy}_W${t_array[$i]}mm_conventional.root+" $macro
	
	halfsize=$(echo ${t_array[$i]} '* 0.5' | bc) #mm
	sed -i "s+/score/mesh/boxSize .*+/score/mesh/boxSize 10. 10. $halfsize mm+" $macro
	var=$(echo "scale=2; ${t_array[$i]} / $binsize" | bc -l) #mm
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
	
	sed -i "s+/det/setSlices .*+/det/setSlices $NbinZ+" $macro
		
	sed -i "s+/score/dumpQuantityToFile boxMesh1 Edep output/Edep_${energy}_W.*+/score/dumpQuantityToFile boxMesh1 Edep output/Edep_${energy}_W${t_array[$i]}mm_conventional.txt+" $macro
	sleep 1
	
	#RUN YOUR MACRO HERE 
	./positronSource $macro $Ncores > output/output_${energy}_W${t_array[$i]}mm_conventional.out

	count=$(($count + 1))
	echo -e "Run $count done!"
done

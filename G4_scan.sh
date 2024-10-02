#Script for launching multiple runs (a loop) of a Geant4 application.
#Runs can differ for random seed and for settings.


######################### INPUT ############################
#t_array=($(seq 1.0 0.1 2.0)) #mm
#t_array=(1.50 2.00 2.50) #mm
t_array=(2.0) #mm

t_conv_array=($(seq 6.0 1.0 12.0)) #mm
#t_conv_array=(9.00 9.50 10.00) #mm
#t_conv_array=(10.0) #mm

D_array=(0) #cm

binsize=0.5 #mm

energy='2.86GeV'

macro='macros/run_FS.mac'
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
n_t_conv=${#t_conv_array[@]}
n_D=${#D_array[@]}

n_runs=$(($n_t * $n_t_conv * $n_D))
echo "A loop of $n_runs runs will be executed"

Ncores=$(grep -c ^processor /proc/cpuinfo)
echo "number of available cores on this computer: $Ncores"

if [ ! -d output ] ; then
	mkdir output 
fi

count=0
for (( i=0; i<$n_t; i++ ))
do
	for (( j=0; j<$n_t_conv; j++ ))
	do
		for (( k=0; k<$n_D; k++ ))
		do
			echo Running simulation $(($count + 1)) at $energy with: rad_th "=" ${t_array[$i]} mm, D "=" ${D_array[$k]} cm, conv_th "=" ${t_conv_array[$j]} mm ...
			
			sed -i "s+/crystal/setCrystalSize .*+/crystal/setCrystalSize 20. 20. ${t_array[$i]} mm+" $macro
			sed -i "s+/det/setRadiatorConverterSepDistance .*+/det/setRadiatorConverterSepDistance ${D_array[$k]} cm+" $macro
			sed -i "s+/det/setConverterSize .*+/det/setConverterSize 100. 100. ${t_conv_array[$j]} mm+" $macro
				
			sed -i "s+/run/setfilenamesave output/output_6GeV_W.*+/run/setfilenamesave output/output_${energy}_W${t_array[$i]}mm_D${D_array[$k]}cm_target${t_conv_array[$j]}mm.root+" $macro		
			
			halfsize_conv=$(echo ${t_conv_array[$j]} '* 0.5' | bc) #mm
			sed -i "s+/score/mesh/boxSize .* #converter+/score/mesh/boxSize 50. 50. $halfsize_conv mm #converter+" $macro
			var=$(echo ${t_conv_array[$j]} '/' $binsize | bc) #mm
			ceil $var
			int_var=$?
			isEven $int_var
			is_int_var_even=$?
			if [ $is_int_var_even == 1 ] ; then
				int_var=$(($int_var + 1))
			fi
			NbinZ_conv=$int_var			
			sed -i "s+/score/mesh/nBin .* #converter+/score/mesh/nBin 201 201 $NbinZ_conv #converter+" $macro	
			
			defZ0=$(echo ${t_conv_array[$j]} ' * 0.5 + 0.02' | bc) #mm
			#echo defZ0: $defZ0 mm		
			if [ ${D_array[$k]} > 0 ] ; then
				defZ=$(echo $defZ0 '+ 0.01' | bc) #mm
			else
				defZ=$defZ0
			fi
			convZ=$(echo ${D_array[$k]} '* 10 +' $defZ | bc) #mm
			sed -i "s+/score/mesh/translate/xyz .* #converter+/score/mesh/translate/xyz 0. 0. $convZ mm #converter+" $macro
			
			halfsize=$(echo ${t_array[$i]} '* 0.5' | bc) #mm
			sed -i "s+/score/mesh/boxSize .* #radiator+/score/mesh/boxSize 10. 10. $halfsize mm #radiator+" $macro
			var=$(echo ${t_array[$i]} '/' $binsize | bc) #mm
			ceil $var
			int_var=$?
			isEven $int_var
			is_int_var_even=$?
			if [ $is_int_var_even == 1 ] ; then
				int_var=$(($int_var + 1))
			fi
			NbinZ=$int_var			
			sed -i "s+/score/mesh/nBin .* #radiator+/score/mesh/nBin 41 41 $NbinZ #radiator+" $macro
			sed -i "s+/score/mesh/translate/xyz .* #radiator+/score/mesh/translate/xyz 0. 0. -$halfsize mm #radiator+" $macro	
			
			sed -i "s+/score/dumpQuantityToFile boxMesh1 Edep output/Edep_${energy}_W.*+/score/dumpQuantityToFile boxMesh1 Edep output/Edep_${energy}_W${t_array[$i]}mm_D${D_array[$k]}cm_target${t_conv_array[$j]}mm.txt+" $macro		
			sed -i "s+/score/dumpQuantityToFile boxMesh2 Edep output/Edep_${energy}_W.*+/score/dumpQuantityToFile boxMesh2 Edep output/Edep_${energy}_W${t_array[$i]}mm_D${D_array[$k]}cm_target${t_conv_array[$j]}mm_rad.txt+" $macro		
			sleep 1
			
			#RUN YOUR MACRO HERE 
			./positronSource $macro $Ncores > output/output_${energy}_W${t_array[$i]}mm_D${D_array[$k]}cm_target${t_conv_array[$j]}mm.out
			
			count=$(($count + 1))
			echo -e "Run $count done!"
		done
	done
done

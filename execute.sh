#Ncores=$1

Ncores=$(grep -c ^processor /proc/cpuinfo)
echo "number of available cores on this computer: $Ncores"

./positronSource macros/run.mac $Ncores > output/output.out

#rm output/*_t*.root

#nohup sh execute.sh > output/nohup &
#nohup sh G4_scan_conventional.sh > output/nohup &
#nohup sh G4_scan.sh > output/nohup &
#nohup sh G4_scan_crystalline.sh > output/nohup &
#nohup sh G4_scan_GT.sh > output/nohup &

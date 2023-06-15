#! /bin/bash

N=20
K=1
Q=1
Ham=7 #3 simplified, 4 quadratic
Alpha=0.025
Cm=0.025
Cgap=1
maxT=50
samplingStep=0.01
complexM=40
capacity=1
initState=1 #1 everything in 0-mode, 2 N-subState in 0-mode - subState in 1-mode 
subState=0 #only important for initState = 2
complexity=false 
alphaN = 0.9

for N2 in $(seq 10 20 200)
do
Alpha2=$(echo "scale=6', $alphaN/$N2" | bc)
./quantumBreaker --N $N2 --K $K --Q $Q --Ham $Ham --Alpha $Alpha2 --Cm $Cm --Cgap $Cgap --maxT $maxT --samplingStep $samplingStep --complexM $complexM --capacity $capacity --initState $initState --subState $subState --complexity $complexity
done


for N2 in $(seq 250 50 1000)
do
Alpha2=$(echo "scale=6', $alphaN/$N2" | bc)
./quantumBreaker --N $N2 --K $K --Q $Q --Ham $Ham --Alpha $Alpha2 --Cm $Cm --Cgap $Cgap --maxT $maxT --samplingStep $samplingStep --complexM $complexM --capacity $capacity --initState $initState --subState $subState --complexity $complexity
done
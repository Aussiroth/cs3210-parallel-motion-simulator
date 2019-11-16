echo > output.txt
for i in {1..4}
do
	echo $i >> output.txt
	for j in {1..2}
	do
		mpirun -machinefile machinefile.${i} -np 64 --map-by node ./a.out < input/random5000.txt | grep "Time taken" >> output.txt
	done
done

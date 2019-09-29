echo > output.txt
threads='1 2 4 8 10 12 16 20'
for i in $threads
do
	export OMP_NUM_THREADS=$i
	echo Thread count = $i >> output.txt
	for j in {1..5}
	do
		./a.out < random4000.txt | grep "Time taken" >> output.txt
	done
done


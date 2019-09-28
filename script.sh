echo > output.txt
for i in {1..20}
do
	export OMP_NUM_THREADS=$i
	echo Thread count = $i >> output.txt
	for j in {1..5}
	do
		./a.out < random.txt | grep "Time taken" >> output.txt
	done
done


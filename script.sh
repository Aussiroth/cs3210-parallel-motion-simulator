echo > output.txt
for i in {1..5}
do
	echo $i >> output.txt
	for j in {1..5}
	do
		./a.out < input/random${i}000.txt | grep "Time taken" >> output.txt
	done
done

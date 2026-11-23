# for r in `seq 10 10 100`
# do
	# ./mcp -i ~/Documents/data/REAL/all.txt -o ANS_DEC/REAL_$r -k $r ga -s
	# ./mcp -i ~/Documents/data/REAL/all.txt -o ANS_DEC/REAL_$r  -d ANS_DEC/REAL_$r  -s -b -f none
# done
for r in `seq 400 100 1000`
do
	./mcp -i ~/Documents/data/REAL/all.txt -o ANS_DEC/REAL_$r -k $r ga -s
	# ./mcp -i ~/Documents/data/REAL/all.txt -o ANS_DEC/REAL_$r  -d ANS_DEC/REAL_$r  -s -b -f none
done

.PHONY: clean test bpv2 ampv4

clean:
	rm -rf puma-out

test:
	pytest -xv test.py

bpv2: clean
	./run_puma.py -i bpv2.fa -d ../data_dir -D debug

ampv4: clean
	./run_puma.py -i ampv4.fa -f fasta -d ../data_dir -D debug	

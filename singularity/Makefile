SINGULARITY = /usr/local/bin/singularity
IMG = puma-1.2.0.img
DEF = image.def

run:
	sudo $(SINGULARITY) run $(IMG)

clean:
	rm -f $(IMG)

img: clean
	sudo $(SINGULARITY) build $(IMG) $(DEF)

sandbox: clean
	sudo $(SINGULARITY) build --sandbox $(IMG) $(DEF)

shell:
	sudo $(SINGULARITY) shell --writable -B /data:/work -B $(HOME):/home1 -B $(shell pwd):/tmp $(IMG)

test:
	sudo $(SINGULARITY) exec -B /data:/work -B $(HOME):/home1 $(IMG) run_puma.py \
    -o /work/kyclark/hpv/puma-out \
    -i /work/kyclark/hpv/BPV2.fa \
    -d /work/kyclark/hpv/data_dir -D debug

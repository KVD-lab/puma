SINGULARITY = /usr/local/bin/singularity
IMG = puma-0.3.0.img
DEF = image.def
SIZE = 2048

run:
	sudo $(SINGULARITY) run $(IMG)

clean:
	rm -f $(IMG)

img: clean
	#sudo $(SINGULARITY) build $(IMG) $(DEF)
	sudo $(SINGULARITY) create --size $(SIZE) $(IMG)
	sudo $(SINGULARITY) bootstrap $(IMG) $(DEF)

shell:
	sudo $(SINGULARITY) shell --writable -B /data:/work -B $(HOME):/home1 -B $(shell pwd):/tmp $(IMG)

test:
	sudo $(SINGULARITY) exec -B /data:/work -B $(HOME):/home1 $(IMG) puma.py \
    -o /work/kyclark/hpv/puma-out \
    -i /work/kyclark/hpv/BPV2.fa \
    -d /work/kyclark/hpv/data_dir

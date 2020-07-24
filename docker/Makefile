TAG = kvdlab/puma:1.2.0
DATA_DIR = $(shell pwd)/../data_dir
DEV_DIR = $(shell pwd)/../development
OUT_DIR = $(shell pwd)

img:
	docker build --tag=$(TAG) .

shell:
	docker run --rm -it $(TAG) /bin/bash

push: 
	docker push $(TAG)

# -v $(OUT_DIR):/out \

test:
	docker run --rm \
		-v $(DEV_DIR):/input \
		-v $(DATA_DIR):/data \
		$(TAG) \
		run_puma.py -i /input/BPV2_new.fa -d /data 

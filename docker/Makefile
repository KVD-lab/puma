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

test:
	docker run --rm \
		-v $(DATA_DIR):/data \
		-v $(DEV_DIR):/development \
		-v $(OUT_DIR):/out \
		$(TAG) \
		run_puma.py -o /out -i /development/BPV2_new.fa -d /data 

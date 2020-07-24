TAG = kvdlab/puma:1.2.0

img:
	docker build --tag=$(TAG) .

shell:
	docker run --rm -it $(TAG) /bin/bash

push: 
	docker push $(TAG)

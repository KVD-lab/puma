TAG = kvdlab/puma:1.2.0

img:
	docker build --tag=$(TAG) .

shell:
	docker run -it $(TAG) bash

push: 
	docker push $(TAG)

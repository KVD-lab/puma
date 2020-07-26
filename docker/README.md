# How to use the Docker app

Be sure you have Docker installed on your machine.
Do the following to pull the image locally:

```bash
docker pull kvdlab/puma:1.2.0
```

Now you should be able to run the image and see the following output:

```bash
$ docker run --rm -it kvdlab/puma:1.2.0
usage: run_puma.py [-h] -i FILE [-f FORMAT] [-d DIR] [-o DIR] [-e FLOAT]
                   [-m NUM] [-D STR] [-L FILE]
run_puma.py: error: the following arguments are required: -i/--input
```

To run Puma on your data, you will need to [mount](https://docs.docker.com/storage/bind-mounts/) your local input and output directories.
For instance, from within this "docker" directory, we can run the program like so:

```bash
$ docker run --rm \
    -v `pwd`/../data_dir:/data \           <1>
    -v `pwd`/../input_and_output:/in_out \ <2>
    kvdlab/puma:1.2.0 \                    <3>
    run_puma.py -i /in_out/BPV2_new.fa -o /in_out/puma_out -d /data <4>
```

1. The "data_dir" dir will be mounted as "/data"
2. The "development" dir will be mounted as "/development"
3. This is the tag of the Docker image to run
4. This is the command to execute Puma with the input and output arguments

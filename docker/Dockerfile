FROM python:3.8.5
RUN apt-get -y update && \
        apt-get install -y locales git python3 python3-pip ncbi-blast+ \
                wget zlib1g-dev muscle && \
        locale-gen en_US.UTF-8


WORKDIR /app
RUN git clone https://github.com/KVD-lab/puma.git && \
    python3 -m pip install -r /app/puma/requirements.txt

ENV PATH=/app/puma/scripts:$PATH

RUN wget http://meme-suite.org/meme-software/5.0.2/meme-5.0.2.tar.gz && \
        tar xvf meme-5.0.2.tar.gz && \
        cd meme-5.0.2 && \
        ./configure --prefix=/usr/local && make && make install

# Run when the container launches
CMD ["/app/puma/scripts/run_puma.py"]

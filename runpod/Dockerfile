FROM nvidia/cuda:12.2.2-runtime-ubuntu22.04

RUN apt-get update && apt-get install -y \
    python3 python3-pip python3-venv wget curl git \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install --no-cache-dir \
    openmm==8.6.0 \
    pdbfixer \
    mdtraj \
    numpy \
    scipy \
    matplotlib \
    requests

WORKDIR /workspace
COPY . /workspace/

CMD ["python3", "runpod_membrane_sim.py"]

#!/bin/bash 

CUDA=$1

echo "General DESC-jl installation..." 

if [ "$#" -eq 0 ]; then
    echo "No arguments provided. One of {cpu, cuda11, cuda12} required."
    exit 1
fi

if [ "$1" == "cpu" ]; then 
    echo "Installing JAX with CPU-only. This may substantially slow \
        down performance. " 
    pip install --upgrade pip
    pip install --upgrade "jax[cpu]"
fi 

if [ "$1" != "cuda11" ] && [ "$1" != "cuda12" ]; then 
    echo "Provided installation version of is not supported. Currently supported versions are one of {cpu, cuda11, cuda12}" 
        exit 1 
fi 

conda install pip

echo "Beginning with installing jax with CUDA=$CUDA"

pip install --upgrade "jax[${CUDA}_pip]" -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html

pip install numpy 
pip install scipy 
pip install pandas 
pip install desc-opt
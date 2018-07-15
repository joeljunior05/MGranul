FROM jupyter/tensorflow-notebook:2c80cf3537ca

USER root

RUN apt-get update && \
    apt-get install -y clang && \
    apt-get install -y libopencv-dev python-opencv && \
    apt-get install -y cmake

USER $NB_UID
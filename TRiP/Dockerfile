FROM ubuntu:18.04

LABEL version="0.6"

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y wget git-all \
    && rm -rf /var/lib/apt/lists/*

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm Miniconda3-latest-Linux-x86_64.sh

RUN git clone https://github.com/HugoAi2bc/TRiP.git

ENV BASH_ENV ~/.bashrc

ENV PATH=/root/miniconda3/bin:$PATH
ENV PATH=/root/miniconda3/envs/all_TRiP/bin:$PATH

SHELL ["/bin/bash", "-c"]

RUN conda env create -f /TRiP/TRiP/all_TRiP.yml
RUN conda init bash
RUN echo "conda activate all_TRiP" >> /root/.bashrc \
    && conda info --envs \
    && conda list --name all_TRiP

CMD ["bash","/TRiP/TRiP/TRiP.sh"]

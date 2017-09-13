FROM ubuntu:latest

USER root

RUN apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv EA312927
RUN apt-key adv --keyserver pgp.skewed.de --recv-key 612DEFB798507F25
RUN echo "deb http://repo.mongodb.org/apt/ubuntu xenial/mongodb-org/3.2 multiverse" | tee /etc/apt/sources.list.d/mongodb-3.2.list
RUN echo "deb http://downloads.skewed.de/apt/xenial xenial universe" | tee /etc/apt/sources.list.d/graph-tool.list
RUN apt-get update

RUN apt-get install -y build-essential python-pip python virtualenv python-dev zlib1g-dev \
    liblapack-dev gfortran libfreetype6-dev libpng-dev fontconfig pkg-config mongodb-org \
    python-graph-tool python-matplotlib python-matplotlib-venn python-seaborn

RUN adduser --disabled-password --gecos "" mgkit

USER mgkit

WORKDIR /home/mgkit/

RUN pip install --user mgkit jupyter ipython statsmodels patsy scikit-learn numexpr ipywidgets pandas pysam
RUN echo 'export PYTHONPATH=$HOME/.local/:$PTYHONPATH' >> .bashrc
RUN echo 'export PATH=$HOME/.local/bin:$PATH' >> .bashrc


EXPOSE 8888

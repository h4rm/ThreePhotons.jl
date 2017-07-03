#To run this docker, use:
#docker run --rm --ipc="host" harmonics/julia_cluster_performance

FROM ubuntu:14.04
MAINTAINER Benjamin von Ardenne

#Basic julia installation 0.4.5
RUN apt-get update
RUN apt-get install -y software-properties-common build-essential make libfftw3-dev libgsl0-dev curl wget python2.7
RUN apt-add-repository ppa:staticfloat/juliareleases
RUN apt-add-repository ppa:staticfloat/julia-deps
RUN apt-get update
RUN apt-get install -y julia

#Copy over all important files
RUN mkdir home/standalone
COPY .* home/standalone/

#Let's recompile the s2kit in the machine
RUN cd home/standalone/src && make clean && make

#Customize julia for my needs
RUN julia home/standalone/start_fresh.jl

#Let julia run one time to prepare all the packages
RUN cd home/standalone && julia -e 'include("determination.jl")'

#ENTRYPOINT ["./home/standalone/performance_run.sh"]

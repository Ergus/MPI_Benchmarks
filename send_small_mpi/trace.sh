#!/usr/bin/env bash

export EXTRAE_ON=1
export EXTRAE_CONFIG_FILE=extrae.xml
export LD_PRELOAD=${EXTRAE_HOME}/lib/libnanosmpitrace.so

taskset -c 0-8 ./$@

#!/usr/bin/env bash

# Copyright (C) 2021  Jimmy Aguilar Mena

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


export EXTRAE_ON=1
export EXTRAE_CONFIG_FILE=extrae.xml

export ASAN_OPTIONS=@ASAN_OPTIONS@

# Remove previous traces
if  [ -z ${MPI_LOCALRANKID} ] || [ $MPI_LOCALRANKID -eq 0 ]; then
	echo "Deleting old traces"
	rm -rf TRACE.* set-0/ *.prv *.row *.pcf
fi

export LD_PRELOAD="@LIBASAN@:${EXTRAE_HOME}/lib/libompitrace.so"
$@
unset LD_PRELOAD

# Create the mpits
# if  [ -z ${MPI_LOCALRANKID} ] || [ $MPI_LOCALRANKID -eq 0 ]; then
# 	echo "Creating mpits file manually"
# 	ls ${PWD}/set-0/*.mpit | sed 's/$/ named/' > TRACE.mpits
# fi

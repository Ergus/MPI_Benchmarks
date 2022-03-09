# Copyright (C) 2022  Jimmy Aguilar Mena

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

# Compilation
ADD_EXECUTABLE(Node Node.cxx ${SOURCES})

# Submitter scripts and filters
CONFIGURE_FILE(test_spawn.sh test_spawn.sh COPYONLY)
CONFIGURE_FILE(test_shrink.sh test_shrink.sh COPYONLY)

CONFIGURE_FILE(gnuplot.plt gnuplot.plt COPYONLY)

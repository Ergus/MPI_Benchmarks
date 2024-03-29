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

# Add the executables and link


FIND_PACKAGE(PkgConfig REQUIRED)
PKG_CHECK_MODULES(SLURM REQUIRED slurm)

set(executable "malleability")

message("- Add benchmark: ${executable}")

add_executable(${executable} Malleability.cxx)
target_include_directories(${executable} PUBLIC ${SLURM_INCLUDE_DIRS})
target_link_libraries(${executable} ${SLURM_LINK_LIBRARIES})


# submitter
configure_file("${PROJECT_SOURCE_DIR}/utils/submit_dim.sh" submit_dim.sh @ONLY)
configure_file("${PROJECT_SOURCE_DIR}/utils/submitter_dim.sh" submitter_dim.sh @ONLY)

# Extrae
configure_file("${PROJECT_SOURCE_DIR}/utils/submit_extrae.sh" submit_extrae.sh @ONLY)
configure_file("${PROJECT_SOURCE_DIR}/utils/trace.sh" trace.sh @ONLY)
configure_file("${PROJECT_SOURCE_DIR}/utils/extrae_template.xml" extrae.xml @ONLY)

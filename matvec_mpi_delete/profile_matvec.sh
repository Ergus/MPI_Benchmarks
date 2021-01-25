#!/bin/bash

l_nodes=( 1 2 4 8 16 )
l_sizes=( 1024 2048 4096 8192 16384 32768 )
l_threads=( 1 2 4 8 16 )

id=$(date +%a%d%m%Y_%H%M)
iterations=1

res_dir="@CMAKE_BINARY_DIR@/results/matvec_${id}"
exec_dir="@CMAKE_BINARY_DIR@/matvec_mpi"

with_extrae=false   # do not use extrae by default

print_help()
{
    echo -e "Available swithces:\n
	 -i|--iterations
		\tExtrae profiling
	 -e|--extrae
		\tExtrae profiling
	 -h|--help
		\tPrint this message"
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    key=$1
    case $key in
	-i|--iterations)
	    iterations=$2
	    shift
	    ;;
	-h|--help)
	    print_help
	    exit 1
	    ;;
	-e|--extrae)
	    with_extrae=true
	    res_dir="@CMAKE_BINARY_DIR@/results/matvec-extrae_${id}"
	    ;;
	*)
	    echo "Unknown option"
	    print_help
	    exit 1
	    ;;
    esac
    shift   # Parse next
done

echo "===== SUBMITING EXPERIMENT ====="
echo "USE_EXTRAE: ${with_extrae}"
echo "Iterations: ${iterations}"
echo "nodes: ${l_nodes[@]}"
echo "sizes: ${l_sizes[@]}"
echo "threads: ${l_threads[@]}"
echo "output_dir: ${res_dir}"
echo "================================"

mkdir -p $res_dir

for size in ${l_sizes[@]}; do
    for nodes in ${l_nodes[@]}; do
	for threads in ${l_threads[@]}; do
	    for ((it=0; it<iterations; ++it)); do
		prefix="sz_${size}_nr_${nodes}_tr_${threads}_it_${it}"
		ofile="${res_dir}/${prefix}.out"
		efile="${res_dir}/${prefix}.err"

		echo "Submiting: sz: ${size} nr: ${nodes} tr: ${threads} it: ${it}"
		
		if ( $with_extrae ); then	   
		    sed -e "s|PREFIX|${prefix}|"                  \
			-e "s|TMPDIR|${res_dir}/${prefix}|"       \
			-e "s|FINALDIR|${res_dir}/${prefix}|"     \
			@CMAKE_BINARY_DIR@/extrae_template.xml > ${res_dir}/${prefix}.xml
		fi

		sed -e "s|SIZE|${size}|" 		       	\
		    -e "s|NODES|${nodes}|" 	               	\
		    -e "s|THREADS|${threads}|"		       	\
		    -e "s|ITS|100|" 		               	\
		    -e "s|EXTRAEFILE|${res_dir}/${prefix}.xml|"	\
		    -e "s|WITHEXTRAE|${with_extrae}|"  	        \
		    ${exec_dir}/run_matvec.sh | bsub -W 01:00 -n $nodes    \
						     -J MV${size}_${nodes} \
						     -e $efile -o $ofile
	    done # it
	done # threads
    done # nodes
done # size


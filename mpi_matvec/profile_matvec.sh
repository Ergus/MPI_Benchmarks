#!/bin/bash

conf_procs=( 1 2 4 8 16 )
conf_size=( 1024 2048 4096 8192 16384 32768 )

res_dir="@CMAKE_BINARY_DIR@/results/matvec"
exec_dir="@CMAKE_BINARY_DIR@/mpi_matvec"

tasks_per_node=1    # default setup
with_extrae=false   # do not use extrae by default

print_help()
{
    echo -e "Available swithces:\n
    	 -t|--tasks-per-node
		\tNumber of tasks to create per node. By default 1
	 -e|--extrae
		\tExtrae profiling
	 -h|--help
		\tPrint this message"

}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    key="$1"

    case $key in
	-t|--tasks-per-node)
	    tasks_per_node=$2
	    shift
	    ;;
	-h|--help)
	    print_help
	    exit 1
	    ;;
	-e|--extrae)
	    with_extrae=true
	    res_dir="@CMAKE_BINARY_DIR@/results/matvec-extrae"
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
echo "TASKS_PER_NODE: ${tasks_per_node}"
echo "USE_EXTRAE: ${with_extrae}"
echo "================================"

mkdir -p $res_dir

for size in "${conf_size[@]}"; do
    for proc_nr in "${conf_procs[@]}"; do
	task_size=$(( size / (proc_nr * tasks_per_node) ))
	
	prefix="sz_${size}_nr_${proc_nr}_ts_${tasks_per_node}"
	
	ofile="${res_dir}/${prefix}.out"
	efile="${res_dir}/${prefix}.err"
	rm -f ${ofile} ${efile}

	if ( $with_extrae ); then	   
	    sed -e "s|PREFIX|${prefix}|"                  \
		-e "s|TMPDIR|${res_dir}/${prefix}|"       \
		-e "s|FINALDIR|${res_dir}/${prefix}|"     \
		@CMAKE_BINARY_DIR@/extrae_template.xml > ${res_dir}/${prefix}.xml
	fi
	
	sed -e "s/NR_PROCS/${proc_nr}/" 	           \
	    -e "s/ROWS/${size}/" 		           \
	    -e "s/TS/${task_size}/"		           \
	    -e "s/ITS/100/" 		 	           \
	    -e "s|EXTRAEFILE|${res_dir}/${prefix}.xml|"         \
	    ${exec_dir}/run_matvec.sh | bsub -W 01:00 -n $proc_nr         \
					    -J matvec_${size}_${proc_nr} \
					    -e $efile -o $ofile
    done #proc_nr
done #size

echo "Wait for all jobs to finish"
while : ; do
    sleep 10
    [[ "$((bjobs) 2>&1)" != "No unfinished job found" ]] || break
done

echo "Parsing results..."
awk -v odir=${res_dir}            \
    '/algorithm_time/ {           \
        split(FILENAME,a,"[_.]"); \
        print a[4], a[2]/(a[4]*a[6]), $5 >> odir"/stats_matvec_"(a[2])"_"(a[6])".csv"; \
	}'                        \
    ${res_dir}/*.out

# Sort the files content (optional)
for file in ${res_dir}/*.csv; do
    sort -n $file -o $file &
done
wait



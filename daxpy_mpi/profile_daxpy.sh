#!/bin/bash

conf_procs=( 1 2 4 8 16 )
conf_size=( 100000 200000 400000 800000 )

res_dir="@CMAKE_BINARY_DIR@/results/daxpy"
exec_dir="@CMAKE_BINARY_DIR@/daxpy_mpi"

tasks_per_node=1	#default setup
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

#parse arguments
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
	    res_dir="@CMAKE_BINARY_DIR@/results/daxpy-extrae"
	    ;;		
	*)
	    echo "Unknown option"
	    print_help
	    exit 1
	    ;;
    esac
    shift
done

echo "===== SUBMITING EXPERIMENT ====="
echo "TASKS_PER_NODE:${tasks_per_node}"
echo "USE_EXTRAE: ${with_extrae}"
echo "================================"

alpha=0.7
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
	
	sed -e "s/NR_PROCS/${proc_nr}/"		        \
	    -e "s/SIZE/${size}/"			\
	    -e "s/TS/${task_size}/"			\
	    -e "s/ALPHA/${alpha}/" 			\
	    -e "s|EXTRAEFILE|${res_dir}/${prefix}.xml|" \
	    -e "s|WITHEXTRAE|${with_extrae}|"               \
	    ${exec_dir}/run_daxpy.sh | bsub -W 30 -n $proc_nr 	         \
					    -J daxpy_${size}_${proc_nr}  \
			                    -e $efile -o $ofile
    done #size
done #proc_nr

echo "Wait for all jobs to finish"
while : ; do
    sleep 10
    [[ "$((bjobs) 2>&1)" != "No unfinished job found" ]] || break
done

echo "Parsing results..."
awk -v odir=${res_dir}            \
    '/algorithm_time/ {           \
    	split(FILENAME,a,"[_.]"); \
	print a[4], a[2]/(a[4]*a[6]), $5 >> odir"/stats_daxpy_"(a[2])"_"(a[6])".csv"; \
	}'                        \
    ${res_dir}/*.out

# Sort the files (optional)
for file in ${res_dir}/*.csv; do
    sort -n $file -o $file &
done
wait


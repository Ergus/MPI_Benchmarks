#!/bin/bash

conf_procs=( 1 2 4 8 16 )
conf_size=( 1024 2048 4096 8192 16384 32768 )
res_dir="@CMAKE_BINARY_DIR@/results/matvec"
exec_dir="@CMAKE_BINARY_DIR@/mpi_matvec/"

tasks_per_node=1    #default setup

print_help()
{
    echo -e "Available swithces:\n
    	 -t|--tasks-per-node
			\tNumber of tasks to create per node. By default 1
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
	*)
	    echo "Unknown option"
	    print_help
	    exit 1
	    ;;
    esac
    shift   # Parse next
done

echo "===== SUBMITING EXPERIMENT ====="
echo "TASKS_PER_NODE:${tasks_per_node}"
echo "================================"

mkdir -p $res_dir

for size in "${conf_size[@]}"; do
    for proc_nr in "${conf_procs[@]}"; do
		task_size=$(( size / (proc_nr * tasks_per_node) ))
		ofile="${res_dir}/sz_${size}_nr_${proc_nr}_ts_${tasks_per_node}.out"
		efile="${res_dir}/sz_${size}_nr_${proc_nr}_ts_${tasks_per_node}.err"
		rm -f ${ofile} ${efile}
		sed -e "s/NR_PROCS/${proc_nr}/" 			    \
			-e "s/ROWS/${size}/" 				    	\
			-e "s/TS/${task_size}/"				    	\
			-e "s/ITS/100/" 				            \
			${exec_dir}/run_matvec.sh | bsub -W 01:00 -n $proc_nr         \
											 -J matvec_${size}_${proc_nr} \
											 -e $efile -o $ofile
    done #proc_nr
done #size

#Wait for all jobs to finish
sleep 10;
while [[ "$((bjobs) 2>&1)" != "No unfinished job found" ]]; do
    sleep 10
done

echo "Parsing results..."
awk -v odir=${res_dir}            \
	'/algorithm_time/ {           \
		split(FILENAME,a,"[_.]"); \
		print a[4], a[2]/(a[4]*a[6]), $5 >> odir"/stats_matvec_"(a[2])"_"(a[6])".csv" \
		}'                        \
	${res_dir}/*.out

# Sort the files content (optional)
for file in ${res_dir}/*.csv; do
	sort -n $file -o $file &
done
wait

# for size in "${conf_size[@]}"; do
#     stats_file="${res_dir}/stats_matvec_${size}_${tasks_per_node}.csv"
#     echo "#procs task_size algorithm_time(nsec)" >> ${stats_file}
#     echo "#TASKS_PER_NODE: ${tasks_per_node}" > ${stats_file}
#     for proc_nr in "${conf_procs[@]}"; do
# 		task_size=$(( size / (proc_nr * tasks_per_node) ))
# 		ofile="${res_dir}/sz_${size}_nr_${proc_nr}_ts_${tasks_per_node}.out"
# 		algorithm_time=$(awk '/algorithm_time/{print $5}' $ofile)
# 		echo "${proc_nr} ${task_size} ${algorithm_time}" >> ${stats_file}
#     done #proc_nr
# done #size

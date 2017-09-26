#!/bin/bash

conf_procs=( 1 2 4 8 16 )
conf_size=( 100000 200000 400000 800000 )
res_dir="@CMAKE_BINARY_DIR@/results/daxpy"
exec_dir="@CMAKE_BINARY_DIR@/mpi_daxpy/"

tasks_per_node=1	#default setup

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
	shift
done

echo "===== SUBMITING EXPERIMENT ====="
echo "TASKS_PER_NODE:${tasks_per_node}"
echo "================================"

alpha=0.7

mkdir -p $res_dir
for proc_nr in "${conf_procs[@]}"; do
	for size in "${conf_size[@]}"; do
		task_size=$(( size / (proc_nr * tasks_per_node) ))
		ofile="${res_dir}/sz_${size}_nr_${proc_nr}_ts_${tasks_per_node}.out"
		efile="${res_dir}/sz_${size}_nr_${proc_nr}_ts_${tasks_per_node}.err"
		rm -f ${ofile} ${efile}
		sed -e "s/NR_PROCS/$proc_nr/" 		\
		    -e "s/SIZE/$size/" 				\
		    -e "s/TS/$task_size/"			\
		    -e "s/ALPHA/$alpha/" 			\
		    -e "s/PRINT/0/" 				\
		    ${exec_dir}/run_daxpy.sh | bsub -W 30 -n $proc_nr 	         \
											-J daxpy_${size}_${proc_nr}  \
			                                -e $efile -o $ofile
	done #size
done #proc_nr

#Wait for all jobs to finish
sleep 10
while [[ "$((bjobs) 2>&1)" != "No unfinished job found" ]]; do
	sleep 10
done

awk -v odir=${res_dir}            \
	'/algorithm_time/ {               \
		split(FILENAME,a,"[_.]"); \
		print a[4], a[2]/(a[4]*a[6]), $5 >> odir"/stats_daxpy_"(a[2])"_"(a[6])".csv" \
		}'                        \
	${res_dir}/*.out

# Sort the files (optional)
for file in *.csv; do
	sort -n $file -o $file &
done
wait

# echo "Parsing results..."
# for size in "${conf_size[@]}"; do
# 	stats_file="${res_dir}/stats_daxpy_${size}_${tasks_per_node}.csv"
# 	echo "#procs task_size algorithm_time(nsec)" >> ${stats_file}
#     echo "#TASKS_PER_NODE: ${tasks_per_node}" > ${stats_file}
# 	for proc_nr in "${conf_procs[@]}"; do
# 		task_size=$(( size / (proc_nr * tasks_per_node) ))
# 		ofile="${res_dir}/sz_${size}_nr_${proc_nr}_ts_${tasks_per_node}.out"
# 		algorithm_time=$(awk '/algorithm_time/{print $5}' $ofile)
# 		echo "${$proc_nr} ${task_size} ${algorithm_time}" >> ${stats_file}
# 	done #proc_nr
# done #size

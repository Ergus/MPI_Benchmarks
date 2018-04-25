

BEGIN {
	first=1
}
/WSIZE/||/^SIZE/||/^BSIZE/ {
	sizes[$1]=$3
}
/^TIME/ {
	array[$1]=$3
}

ENDFILE {
	if (length(array) > 0) {
		if (first) {
			printf "NR_PROCS ROWS TASK_SIZE "
			for (x in array)
				printf "%s ", x
			printf "PERFORMANCE"
			printf "\n";
		first = 0
		}
		size=sizes["SIZE"]
		printf "%d %d %d ", sizes["WSIZE"], sizes["SIZE"], sizes["BSIZE"]
		for (x in array)
			printf "%g ", array[x];
		printf "%g\n", (size^3) / (array["TIME(cholesky)"] * 3.0e+3);

		delete array;
		delete sizes;

	} else {
		failed[cont++] = FILENAME
	}
}
END{
	print "# Failed:"
	for (x in failed) {
		print "#", failed[x]
	}
}

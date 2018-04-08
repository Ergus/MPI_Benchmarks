

BEGIN {
	print "# Here start"
}
BEGINFILE {
	#split(FILENAME,arr,"[|,/,_]")
}
/WSIZE/||/^SIZE/||/^BSIZE/||/^PERFORMANCE/||/^TIME/ {
	array[$1]=$3
}
ENDFILE {
	if (length(array) > 0) {
		for (x in array) {
			printf "%s\t%s\t", x, array[x]};
		printf "\n";
		delete array;
	} else {
		failed[cont++] = FILENAME
	}
}
END{
	print "Failed:"
	for (x in failed) {
		print failed[x]
	}
}

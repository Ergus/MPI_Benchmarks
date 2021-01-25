#!/usr/bin/awk -f

# mean
function mean(myarray)
{
	ret = 0
	cont = 0
	for (x in myarray) {
		ret += myarray[x]
		cont++
	}
	return ret / cont
}

# standard deviation
function devest(myarray, mean)
{
	dev = 0
	cont = 0
	for (x in myarray) {
		dev += (myarray[x] - mean)^2
		cont++
	}
	return sqrt(dev / (cont - 1) )
}

# output for babis format
function getoutput(file)
{
	n = split(file,tmp,"/")
	filename = tmp[n]
	sub("/"filename, "", file)
	return file
}

# Filter to temporal array
/WSIZE/||/^SIZE/||/^BSIZE/||/^TIME/ { array[$1]=$3 }

ENDFILE { # save to global array at the end of every file
	if (length(array) > 0) {
		iti = ti[array["WSIZE"]][array["SIZE"]][array["BSIZE"]]++
		size = array["SIZE"]
		time = array["TIME(cholesky)"]
		perf = (size^3) / (time * 3.0e+3) # For Babis format
		times[array["SIZE"]][array["WSIZE"]][array["BSIZE"]][iti] = time
		perfo[array["SIZE"]][array["WSIZE"]][array["BSIZE"]][iti] = perf

		delete array;
	} else { # The file failed, save the name to check latter
		failed[cont++] = FILENAME
	}
}
END { # Print at the very end in my format or babis' format
	if (!VERBOSE) { # Babis output (everything to a file(size))
		header = "NR_PROCS ROWS TASK_SIZE TIME PERFORMANCE"
		for (sz in times) { # size
			output = getoutput(FILENAME) "_" sz ".out"
			print header > output
			for (np in times[sz]) { # Both arrays have the same indices, wsize
				for (bsz in times[sz][np]) { #bsize
					for (iti in times[sz][np][bsz]) {
						printf "%d %d %d %g %g\n", np, sz, bsz,			\
							times[sz][np][bsz][iti], perfo[sz][np][bsz][iti] > output
					}
				}
			}
		}
	} else if (VERBOSE == 1) { # My output everything to stdout averaged
		print "NR_PROCS SIZE BSIZE TIME ERROR PERFORMANCE ERROR"
		for (sz in times) { # size
			for (np in times[sz]) { # Both arrays have the same indices, wsize
				for (bsz in times[sz][np]) { #bsize
					meant = mean(times[sz][np][bsz])
					meanp = mean(perfo[sz][np][bsz])
					printf "%d %d %d %g %g %g %g\n", np, sz, bsz,			\
						meant, devest(times[sz][np][bsz], meant),			\
						meanp, devest(perfo[sz][np][bsz], meanp)
				}
			}
		} # Also print the failed benchmarks
		if (length(failed) > 0) {
			print "#", "Failed"
			for (x in failed) {
				print "#", failed[x]
			}
		}
	}
}

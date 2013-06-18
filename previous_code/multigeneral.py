import subprocess
import sys

for i in range(12):
	process = subprocess.Popen("./generalstatistics " + str(i+int(sys.argv[1])) + " " + str(0), shell=True, stdout=subprocess.PIPE)
	for line in process.stdout:
	    print line
	process.wait()

import os
import sys

GSE = sys.argv[1]
stringtie_dir = sys.argv[2]
files = os.listdir(stringtie_dir)
with open(f"lst_{GSE}.txt", "w") as f:
	for file in files:
        	line = f"{file.split('.')[0]} {stringtie_dir}/{file}"
        	f.write(line+"\n")

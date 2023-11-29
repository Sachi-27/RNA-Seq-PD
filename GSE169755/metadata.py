import os
import sys
import pandas as pd

files = os.listdir("stringtie")
files = [file[:-4] for file in files]
query = " OR ".join(files)
os.system(f"esearch -db sra -query '{query}' | efetch -format runinfo > runinfo.csv")

file_path = sys.argv[1]
sample_dict = {}

with open(file_path, 'r') as file:
    for line in file:
        if line.startswith('!Sample_title'):
            titles = line.strip().split()[1:]
            titles = [title[1:-1] for title in titles]
        elif line.startswith('!Sample_geo_accession'):
            geo_accessions = line.strip().split()[1:]
            geo_accessions = [geo_accession[1:-1] for geo_accession in geo_accessions]

for i in range(len(titles)):
    sample_dict[geo_accessions[i]] = titles[i]
    
    
data = pd.read_csv("runinfo.csv")
runs = data["Run"].values
samples = data["SampleName"].values

conditions = []
for i in range(len(samples)):
    conditions.append(sample_dict[samples[i]])
    
conditions = ["Control" if "CTRL" in condition else "Disease" for condition in conditions]

os.system("rm runinfo.csv")

df = pd.DataFrame({'ID': runs, 'SampleName': samples, 'Condition': conditions})
output_file = 'metadata.csv'
df.to_csv(output_file, index=False)
print(f'Data has been written to {output_file}')

import os
import sys
import pandas as pd

#files = os.listdir("stringtie")
#files = [file[:-4] for file in files]


lst_file_path = sys.argv[1]
files = []
with open(lst_file_path, "r") as f:
    for line in f:
        files.append(f.strip().split()[0])

query = " OR ".join(files)
os.system(f"esearch -db sra -query '{query}' | efetch -format runinfo > runinfo.csv")


series_file_path = sys.argv[2]
sample_dict = {}

with open(series_file_path, 'r') as file:
    for line in file:
        if line.startswith('!Sample_title'):
            titles = line.strip().split()[1:]
            titles = [title[1:-1] for title in titles]
        elif line.startswith('!Sample_geo_accession'):
            geo_accessions = line.strip().split()[1:]
            geo_accessions = [geo_accession[1:-1] for geo_accession in geo_accessions]
        elif line.startswith('!Sample_characteristics_ch11'):
            diseases = line.strip().split()[1:]
            diseases = [disease[10:-1] for disease in diseases]
        elif line.startswith('!Sample_characteristics_ch12'):
            regions = line.strip().split()[1:]
            regions = [region[9:-1] for region in regions]
        elif line.startswith('!Sample_characteristics_ch13'):
            genders = line.strip().split()[1:]
            genders = [gender[-2:-1] for gender in genders]
        elif line.startswith('!Sample_characteristics_ch14'):
            rins = line.strip().split()[1:]
            rins = [float(rin[-4:-1]) for rin in rins]
        elif line.startswith('!Sample_characteristics_ch15'):
            dementias = line.strip().split()[1:]
            dementias = [dementia[11:-1] for dementia in dementias]
        elif line.startswith('!Sample_characteristics_ch16'):
            dyskinesias = line.strip().split()[1:]
            dyskinesias = [dyskinesia[13:-1] for dyskinesia in dyskinesias]
        elif line.startswith('!Sample_characteristics_ch16'):
            cholinesterases = line.strip().split()[1:]
            cholinesterases = [cholinesterase[28:-1] for cholinesterase in cholinesterases] 

for i in range(len(titles)):
    sample_dict[geo_accessions[i]] = i
    
    
data = pd.read_csv("runinfo.csv")
runs = data["Run"].values
samples = data["SampleName"].values

ids = []
for i in range(len(samples)):
    ids.append(sample_dict[samples[i]])
diseases = [diseases[i] for i in ids]
regions = [regions[i] for i in ids]
genders = [genders[i] for i in ids]
rins = [rins[i] for i in ids]
dementias = ["Control" for i in ids if dementias[i] == "CTRL" else "Disease"]
dyskinesias = ["Control" for i in ids if dyskinesias[i] == "CTRL" else "Disease"]
cholinesterases = ["Control" for i in ids if cholinesterases[i] == "CTRL" else "Disease"

os.system("rm runinfo.csv")

df = pd.DataFrame({'ID': runs, 'SampleName': samples, 'Condition': diseases, "Region":regions, "Gender":genders, "RIN":rins, "Dementia":dementias, "Dyskinesia":dyskinesias, "Chlonisterase":chlonisterases})
output_file = 'metadata.csv'
df.to_csv(output_file, index=False)
print(f'Data has been written to {output_file}')

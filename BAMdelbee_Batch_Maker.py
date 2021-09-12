from glob import glob
from os import system
import pandas as pd

try:
    system("mkdir temp")
except:
    "temp folder already present"
bams = glob("*.bam")
print("Located BAM files:\n" +"\n".join(bams) + "\n\nProcessing...")

users = [bam[:-4] for bam in bams]
open("temp/samples.txt", 'w').write('\n'.join(users))
def depthCounter(chromosome):
    try:
        for bam in bams:
            system("samtools depth -a -r chr" + chromosome + ' ' + bam + " | cut -f2,3 | awk '$1%50==0' | head -n -1 | sort -n -s -k1,1 > temp/" + bam[:-4] + "_chr" + chromosome + ".counted")
        df = pd.read_csv("temp/" + bams[0][:-4] + "_chr" + chromosome + ".counted", sep = '\t', header = None)
        df.columns = ['Coordinate', bams[0][:-4]]
        for bam in bams[1:]:
            df2 = pd.read_csv("temp/" + bam[:-4] + "_chr" + chromosome + ".counted", sep = '\t', header = None)
            df2.columns = ['Coordinate', bam[:-4]]
            df = df.merge(df2, on='Coordinate')
        df['coverages'] = df.apply(lambda x : x.values[1:], axis = 1) 
        df = df[df['coverages'].apply(lambda x : (sum(x) != 0) and (0 in x) and (max(x) >= 5))]
        del df['coverages']
        df.to_csv("temp/chr" + chromosome + ".coverage", sep = '\t', index = False)
        system("rm -f temp/*hr" + chromosome + ".*ounted")
    except:
        print("Error occured while processing chromosome " + chromosome)
    print("Finished processing chromosome " + chromosome)

x = list(map(depthCounter, [str(chromosome) for chromosome in range(1,23)] + ["X","Y"]))
system("mv temp results")
system("zip rename_me.zip results/*")
system("rm -rf results")



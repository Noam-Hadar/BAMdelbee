from glob import glob
import os
import pandas as pd
import subprocess
import shutil
from time import sleep

print('Please ignore messages such as "samtools view: writing to standard output failed: Broken pipe" or "samtools view: error closing standard output: -1", they do not affect the analysis.')
sleep(1)
print('\n')

try:
    os.mkdir('temp')
except:
    'temp folder already present'
bams = glob('**/*.bam', recursive = True)
print('Located BAM files:\n' + '\n'.join(bams) + '\n')


def getSampleName(name):
    sample = name[:-4]
    if '/' in sample:
        sample = sample.split('/')[-1]
    if '\\' in sample:
        sample = sample.split('\\')[-1]
    return sample

samples = [getSampleName(bam) for bam in bams]
open('temp/samples.txt', 'w').write('\n'.join(samples))

def depthCounter(chromosome):
    print('Processing chromosome ' + chromosome + '...')
    try:
        for bam in bams:
            sample = getSampleName(bam)
            if '\\' in bam:
                bam = bam.replace('\\','/')
            chrAnnotation = 'samtools view ' + bam + ' | head | cut -f3 | uniq'
            if os.name == 'nt':
                head = subprocess.check_output(['wsl'] + chrAnnotation.split(' '))
                if 'chr' in str(head).lower():
                    command = 'samtools depth -a -r chr' + chromosome + ' ' + bam + " | cut -f2,3 | awk '$1%50==0' | head -n -1 | sort -n -s -k1,1 > temp/" + sample + '_chr' + chromosome + '.counted'
                else:
                    command = 'samtools depth -a -r ' + chromosome + ' ' + bam + " | cut -f2,3 | awk '$1%50==0' | head -n -1 | sort -n -s -k1,1 > temp/" + sample + '_chr' + chromosome + '.counted'
                subprocess.call(['wsl'] + command.split(' '), stderr = subprocess.STDOUT)
            else:
                os.system(chrAnnotation + ' > temp.txt')
                head = open('temp.txt', 'r').read()
                os.remove('temp.txt')
                if 'chr' in str(head).lower():       
                    command = 'samtools depth -a -r chr' + chromosome + ' ' + bam + " | cut -f2,3 | awk '$1%50==0' | head -n -1 | sort -n -s -k1,1 > temp/" + sample + '_chr' + chromosome + '.counted'
                else:
                    command = 'samtools depth -a -r ' + chromosome + ' ' + bam + " | cut -f2,3 | awk '$1%50==0' | head -n -1 | sort -n -s -k1,1 > temp/" + sample + '_chr' + chromosome + '.counted'
                subprocess.run(command, shell = True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
            sample = getSampleName(bams[0])
            df = pd.read_csv('temp/' + sample + '_chr' + chromosome + '.counted', sep = '\t', header = None)
        df.columns = ['Coordinate', sample ]
        for bam in bams[1:]:
            sample = getSampleName(bam)
            df2 = pd.read_csv('temp/' + sample  + '_chr' + chromosome + '.counted', sep = '\t', header = None)
            df2.columns = ['Coordinate', sample ]
            df = df.merge(df2, on = 'Coordinate')
        df['coverages'] = df.apply(lambda x : x.values[1:], axis = 1) 
        df = df[df['coverages'].apply(lambda x : (sum(x) != 0) and (0 in x) and (max(x) >= 5))]
        del df['coverages']
        df.to_csv('temp/chr' + chromosome + '.coverage', sep = '\t', index = False)
        for coverage_file in glob('temp/*hr' + chromosome + '.*ounted'):
            os.remove(coverage_file)
    except:
        print('Error occured while processing chromosome ' + chromosome)
    print('Finished processing chromosome ' + chromosome)

x = list(map(depthCounter, [str(chromosome) for chromosome in range(1,23)] + ["X"]))

shutil.make_archive('rename_me', 'zip', 'temp')
shutil.rmtree('temp', ignore_errors=True)
input('Done, press the enter key to exit')

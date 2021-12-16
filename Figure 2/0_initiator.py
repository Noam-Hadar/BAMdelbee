from os import system
from time import sleep

for i in range(1,150):
    system('python3 run1_fig2_fastq_generator.py ' + str(i))
sleep(180)
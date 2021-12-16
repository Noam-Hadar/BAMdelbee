import pandas as pd, matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = [10, 5]
plt.rcParams['figure.autolayout'] = True
with plt.xkcd():
    font = 'Times New Roman'
    for n, i in enumerate(['del', 'ins']):
        ax = plt.subplot(1, 2, n + 1)
        bwa_mem_df = pd.read_csv(i + '_bwa_mem.tsv', sep = '\t')
        bt2_df = pd.read_csv(i + '_bt2.tsv', sep = '\t')
        plt.plot(bt2_df['Size'], bt2_df['Percent aligned'], color = 'red')
        plt.plot(bwa_mem_df['Size'], bwa_mem_df['Percent aligned'], color = 'navy')
        plt.ylabel('Aligned reads (%)', fontname = font, fontsize = 14)
        plt.legend(['Bowtie2', 'BWA-MEM'], prop = font)
        xlimRight = len([i for i in bwa_mem_df['Percent aligned'] + bt2_df['Percent aligned'] if i == 0])
        plt.xlim(0,150 - xlimRight + 5)
        plt.tight_layout()
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        if i == 'del':
            plt.xlabel('Deletion size (base pairs)', fontname = font, fontsize = 14)
        else:
            plt.xlabel('Insertion size (base pairs)', fontname = font, fontsize = 14)
    plt.savefig('Fig2.tif', dpi = 320)
    plt.savefig('Fig2.svg')


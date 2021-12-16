for f in ('del_bwa_mem.tsv', 'ins_bwa_mem.tsv', 'del_bt2.tsv', 'ins_bt2.tsv'):
    open(f, 'w').write('\t'.join(['Size','Reads', 'Aligned reads', 'Percent aligned']))
for i in range(1,149):
    print(i)
    i = str(i)
    for s in ('del_', 'ins_'):
        for e in ('bwa_mem', 'bt2'):
            lines = open(s + i + '.' + e, 'r').read().split('\n')
            if s == 'del_':
                d = 'D'
            else:
                d = 'I'
            hits = len([line for line in lines if i + d in line])
            lines = len(lines)
            percentage = str(round(hits/lines*100))
            open(s + e + '.tsv', 'a').write('\n' + '\t'.join([i, str(lines), str(hits), percentage]))

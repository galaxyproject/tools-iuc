import sys


if __name__ == '__main__':
    with open(sys.argv[1]) as i:
        getmasked_output = i.readline().strip()

    if not getmasked_output:
        print()
        print('No affected primer binding sites found!')
    else:
        masked_primers = getmasked_output.split('\t')
        with open(sys.argv[2]) as i:
            amplicon_data = [line.strip().split('\t') for line in i]

        masked_complete = []
        for primer in masked_primers:
            for amplicon in amplicon_data:
                if primer in amplicon:
                    masked_complete += amplicon
        result = '\t'.join(sorted(set(masked_complete)))
        print()
        print('Removing reads primed with any of:')
        print(result)
        with open(sys.argv[1], 'w') as o:
            o.write(result + '\n')

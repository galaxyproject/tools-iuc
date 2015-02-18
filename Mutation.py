import random
def insertion(mut_seq,changes_list):
    i = random.randint(0,len(mut_seq)-1)
    base = ['A','T','C','G'][random.randint(0,3)]
    new_seq = mut_seq[:i] + str(base) + mut_seq[i:]
    changes_list += ['Insertion of '+str(base)+' at position'+str(i)]
    return new_seq,changes_list

def deletion(mut_seq,changes_list):
    i = random.randint(0,len(mut_seq)-1)
    new_seq = mut_seq[:i] + mut_seq[i+1:]
    changes_list += ['Deletion at position'+str(i)]
    return new_seq,changes_list

def inversion(mut_seq,changes_list):
    i = random.randint(0,len(mut_seq)-1)
    j = random.randint(0,len(mut_seq)-1)
    if i=j:
        break
    elif i<j:
        endpoints = (i,j)
    else:
        endpoints = (j,i)
    new_seq = mut_seq[:endpoints[0]] + mut_seq[endpoints[0]:endpoints[1]:-1] + mut_seq[endpoints[1]:]
    changes_list += ['Inversion between positions ' + str(endpoints[0]) + ' and ' + str(endpoints[1])]
    return new_seq,changes_list

def translocation(mut_seq,changes_list):
    i = random.randint(0,len(mut_seq)-1)
    j = random.randint(0,len(mut_seq)-1)
    k = random.randint(0,len(mut_seq)-1)
    l = random.randint(0,len(mut_seq)-1)
    if i=j:
        break
    elif i<j:
        endpoints = (i,j)
    else:
        endpoints = (j,i)
    if i=j:
        break
    elif i<j:
        new_pos = (i,j)
    else:
        new_pos = (j,i)
    transposon = mut_seq[endpoints[0]:endpoints[1]]
    
    changes_list += ['Inversion between positions ' + str(endpoints[0]) + ' and ' + str(endpoints[1])]
    return new_seq,changes_list




def main():
    bases = ['A','T','C','G']
    base_seq= ''
    for iteration in range(500):
        base_seq += bases[random.randint(0,3)]
    mut_seq = base_seq
    print base_seq
    print mut_seq
    changes_list = []
    mut_seq,changes_list = insertion(mut_seq,changes_list)
    print 'Sequence:',mut_seq
    for change in changes_list:
        print change
    mut_seq,changes_list = deletion(mut_seq,changes_list)
    print 'Sequence:',mut_seq
    for change in changes_list:
        print change
    mut_seq,changes_list = translocation(mut_seq,changes_list)
    print 'Sequence:',mut_seq
    for change in changes_list:
        print change
    mut_seq,changes_list = inversion(mut_seq,changes_list)
    print 'Sequence:',mut_seq
    for change in changes_list:
        print change

main()

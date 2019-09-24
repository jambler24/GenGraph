l =[['Aligned_Read_1', 'Aligned_Read_2'], ['Aligned_Read_2', 'Aligned_Read_3'], ['Aligned_Read_3', 'Aligned_Read_4']]

out = []
while len(l)>0:
    first, *rest = l
    first = set(first)

    lf = -1
    while len(first)>lf:
        lf = len(first)

        rest2 = []
        for r in rest:
            if len(first.intersection(set(r)))>0:
                first |= set(r)
            else:
                rest2.append(r)
        rest = rest2

    out.append(first)
    l = rest

print(list(out[0]))
t = list(out[0])
t.sort()
print(t)

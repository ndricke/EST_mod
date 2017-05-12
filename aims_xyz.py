import sys



infile = sys.argv[1]
outfile = sys.argv[2]

coord_list = []
with open(infile,'r') as f:
    for line in f:
        if line[0] == '#':
            continue
        elif line.split()[0] == "trust_radius":
            break
        else:
            coord = line.split()[1:]
            coord.insert(0,coord.pop())
            coord_list.append(coord)

line_num = len(coord_list)

print coord_list
with open(outfile,'w') as f:
    f.write(str(line_num)+'\n')
    f.write('100\n')
    for item in coord_list:
        f.write(' '.join(item)+'\n')


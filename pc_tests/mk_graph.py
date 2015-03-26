def build_MRF():
    edges = {}
    with open(sys.argv[1]) as fp:
        for line in fp:
            data = line.split()
            edges[int(data[0])] = set(int(x) for x in data[2:])
    new_edges = []
    for i, conns in edges.items():
        for j in conns:
            if i < j:
            #if i in edges[j]:
                print i, j
        
    return

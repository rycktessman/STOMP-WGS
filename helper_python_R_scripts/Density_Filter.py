def ParseSNP(snp_fh, field, sep):
    '''GENERATOR: Parse snp_file and yield tuples of lines
    as strings and positions'''

    for line in snp_fh:
        sline = line.rstrip().split(sep)
        pos = int(sline[field])
        yield (line, pos)

def GetPrefix(infile, prefix="not_provided", sep="."):
    '''Try to get prefix from a file name'''
    if prefix == "not_provided":
        prefix = infile.split(sep)
        # prefix = prefix[0:-1] # remove trailing character
        return prefix
    else:
        return prefix

def GetDensityPos(snp_pos, window, max_density):
    candidates = []
    density_pos = set()
    density = 1
    for i in range(len(snp_pos) - 1): # Voy recorriendo la lista de posiciones
        posA = snp_pos[i] # Cogiendo cada posicion
        posB = snp_pos[i+1] # La que le sigue
        distAB = posB - posA # Y calculando su distancia
        if distAB <= window: # Si es menor a la permitida
            density += 1 # Se aumenta el valor de densidad
            candidates.append(posA) # y la posicion primera se anyade como candidata
        else: # En caso contrario se checkea si hasta ese punto se habia alcanzado
              # la densidad maxima permitida para filtrar los SNPs hasta ese punto
            if density >= max_density:
                for pos in candidates:
                    density_pos.add(pos) # Se anyaden los candidatos hasta el momento
                density_pos.add(posA) # Y el ultimo posA, que era el posB de la
                                     # anterior iteracion
            candidates = [] # Se resetean los candidatos
            density = 1 # y la densidad

    return density_pos


def FilterByDensity(args):
    '''Load deleted positions and filter snp_file
    according to those deletions'''
    import sys

    field = args.field
    field = int(field) - 1
    try:
        prefix = GetPrefix(args.snp_file)
    except IOError:
        sys.exit("{} does not exist. Filtering aborted".format(args.snp_file))

    with open("{}.densF".format(".".join(prefix)), "w") as outfile:
        with open(args.snp_file) as infile:
            header = infile.readline()
            if not header.startswith("Chrom"):
                header = "Chrom\tPosition\tRef\tCons\tReads1\tReads2\tVarFreq\tStrands1\tStrands2\tQual1\tQual2\tPvalue\tMapQual1\tMapQual2\tReads1Plus\tReads1Minus\tReads2Plus\tReads2Minus\tVarAllele\n"
            outfile.write(header)

            snp_pos = []
            for line, pos in ParseSNP(infile, field, args.snp_sep):
                snp_pos.append(pos)

        density_pos = GetDensityPos(snp_pos, args.window, args.max_density)
        with open(args.snp_file) as infile:
            header = infile.readline()
            for line, pos in ParseSNP(infile, field, args.snp_sep):
                if pos not in density_pos:
                    outfile.write(line)



    return 0


import argparse

parser=argparse.ArgumentParser(description="testing")

parser.add_argument("--snp", "--snpfile", dest="snp_file", required=True)
parser.add_argument("-f", "--field", dest="field", default="2", metavar="Field of positions in file. DEFAULT=2")
parser.add_argument("--snpsep", dest="snp_sep", metavar="field delimiter for .snp file. DEFAULT=Tab", default="\t")
parser.add_argument("-m", "--maxdensity", dest="max_density", default=3, type=int, metavar="Max distance between SNPs to be considered as candidate false positives. DEFAULT=3")
parser.add_argument("-w", "--window", dest="window", default=10, type=int, metavar="Max distance between SNPs to be considered as candidate false positives. DEFAULT=10")

args=parser.parse_args()
print(args)

FilterByDensity(args)

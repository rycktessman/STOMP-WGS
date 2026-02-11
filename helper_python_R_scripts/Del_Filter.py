def length_after_minus(input_string):
    # Split the string by the minus sign
    parts = input_string.split('-', 1)
    
    # Return the length of the segment after the minus sign
    if len(parts) > 1:
        return len(parts[1])
    else:
        return 0 

def ParseSNP(snp_fh, field, sep):
    '''GENERATOR: Parse snp_file and yield tuples of lines
    as strings and positions'''

    for line in snp_fh:
        sline = line.rstrip().split(sep)
        pos = int(sline[field])
        yield (line, pos)

def ParseDel(del_file, field, sep, distance, reads_field):
    '''Return a list with deleted positions parsed from
    .del files'''

    with open(del_file) as del_fh:
        deletions = set() # deleted positions
        # Skip header
        del_fh.readline()
        for line in del_fh:
            line = line.rstrip().split(sep)
            #print(line)
            pos = int(line[field])
            # This is just the start of the deletion
            # Also include subsequent deleted positions
            # For deletions only - NOT INSERTIONS
            del_length=length_after_minus(line[reads_field])
            del_set=set(range(0,del_length))
            del_set = {x + pos for x in del_set}
            #print(del_set)
            if(len(del_set)==0):
                #if there is only 1 position deleted
                deletions.add(pos)
                # Now also include all positions that are "distance" nt
                # away upstream and downstream
                for i in range(1, distance+1):
                    deletions.add(pos-i)
                    deletions.add(pos+i)
            else: #if there are multiple positions deleted
                deletions.update(del_set)
                # Now also include all positions that are "distance" nt
                # away upstream and downstream
                for i in range(1, distance+1):
                    deletions.add(min(del_set)-i)
                    deletions.add(max(del_set)+i)
            
    return deletions

def GetPrefix(infile, prefix="not_provided", sep="."):
    '''Try to get prefix from a file name'''
    if prefix == "not_provided":
        prefix = infile.split(sep)
        # prefix = prefix[0:-1] # remove trailing character
        return prefix
    else:
        return prefix

def FilterByDel(args):
    '''Load deleted positions and filter snp_file
    according to those deletions'''
    import sys
    from shutil import copyfile

    field = args.field
    field = int(field) - 1
    try:
        prefix = GetPrefix(args.snp_file)
    except IOError:
        sys.exit("{} does not exist. Filtering aborted".format(args.snp_file))
    try:
        deletions = ParseDel(args.del_file, field, args.del_sep, args.distance, args.reads_field)
    except IOError:
        # Create anyway a delF file that
        copyfile(args.snp_file, "{}.NoDel.delF".format(".".join(prefix)))
        sys.exit("{} does not exist. Filtering aborted".format(args.del_file))

    with open("{}.delF".format(".".join(prefix)), "w") as outfile:
        with open(args.snp_file) as infile:
            header = infile.readline()
            if not header.startswith("Chrom"):
                header = "Chrom\tPosition\tRef\tCons\tReads1\tReads2\tVarFreq\tStrands1\tStrands2\tQual1\tQual2\tPvalue\tMapQual1\tMapQual2\tReads1Plus\tReads1Minus\tReads2Plus\tReads2Minus\tVarAllele\n"
            outfile.write(header)

            for line, pos in ParseSNP(infile, field, args.snp_sep):
                if pos not in deletions:
                    outfile.write(line)

    

    return 0



import argparse

parser=argparse.ArgumentParser(description="testing")

parser.add_argument("--snp", "--snpfile", dest="snp_file", required=True)
parser.add_argument("--delfile", dest="del_file", required=True)
parser.add_argument("-f", "--field", dest="field", default="2", metavar="Field of positions in file. DEFAULT=2")
parser.add_argument("--reads_field", dest="reads_field", default=18)
parser.add_argument("-d", "--distance", dest="distance", default=2, type=int, metavar="Distance of a SNP from a deleted position to be filtered out. DEFAULT=2") 
parser.add_argument("--delsep", dest="del_sep", metavar="field delimiter for .del file. DEFAULT=Space", default="\t")
parser.add_argument("--snpsep", dest="snp_sep", metavar="field delimiter for .snp file. DEFAULT=TAB", default="\t")

args=parser.parse_args()
print(args)

FilterByDel(args)
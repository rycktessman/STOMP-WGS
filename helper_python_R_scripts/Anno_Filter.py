def LoadAnnotation(anno_file):
    '''Load the annotation file in .ppt format
    to store each position with all its
    info in a dictionary

    This function is meant to load by default
    the H37Rv annotation under
    /data/ThePipeline/PipeModules/PipeScripts/H37Rv_annotation2sytems.ptt
    '''
    
    #Note from TR: despite the above comment, this code only works for the .tab file currently used in the pipeline, not the .ptt file
    annotation = {}
    with open(anno_file) as infile:
        # First skip lines until the header line that
        # starts with << Location >> string
        header = infile.readline()
        for line in infile:
            tokens = line.rstrip().split("\t")
            start = int(tokens[1])
            end = int(tokens[2])
            for pos in range(start, end+1):
                annotation[pos] = tokens

    return annotation

def GetPrefix(infile, prefix="not_provided", sep="."):
    '''Try to get prefix from a file name'''
    if prefix == "not_provided":
        prefix = infile.split(sep)
        # prefix = prefix[0:-1] # remove trailing character
        return prefix
    else:
        return prefix

def ParseSNP(snp_fh, field, sep):
    '''GENERATOR: Parse snp_file and yield tuples of lines
    as strings and positions'''

    for line in snp_fh:
        if not line.startswith("#"):
            sline = line.rstrip().split(sep)
            pos = int(sline[field])
            yield (line, pos)

def FilterSnps(args):
    import sys

    field = args.field
    field = int(field) - 1
    annotation = LoadAnnotation(args.annofile)
    try:
        prefix = GetPrefix(args.snp_file)
    except IOError:
        sys.exit("{} does not exist. Filtering aborted".format(args.snp_file))

    with open("{}.annoF".format(".".join(prefix)), "w") as outfile:
        with open(args.snp_file) as infile:
            if args.file_type=="mpileup":
                if not args.noheader:
                    header = infile.readline()
                    if not header.startswith("Chrom"):
                        header = "Chrom\tPosition\tRef\tCons\tReads1\tReads2\tVarFreq\tStrands1\tStrands2\tQual1\tQual2\tPvalue\tMapQual1\tMapQual2\tReads1Plus\tReads1Minus\tReads2Plus\tReads2Minus\tVarAllele\n"
                    outfile.write(header)
            elif args.file_type=="vcf":
                header="Chrom\tPosition\tID\tRef\tCons\tQual\tFilter\tInfo\tFormat\tPL"
                outfile.write(header)

            for line, pos in ParseSNP(infile, field, args.sep):
            	# It could be that position is not in annotation file, in that case
            	# keep anyway
            	if pos not in annotation:
            		outfile.write(line)
            	else:
    	            tokens = annotation[pos]
    	            if tokens[-1] == "KEEP":
    	                outfile.write(line)


    return 0


import argparse

parser=argparse.ArgumentParser(description="testing")

parser.add_argument("--anno", "--annofile", dest="annofile")
parser.add_argument("--snp", "--snpfile", dest="snp_file", required=True)
parser.add_argument("--no-header", dest="noheader", action="store_true", default=0)
parser.add_argument("-t", "--filetype", dest="file_type", default="mpileup")
parser.add_argument("--sep", dest="sep", metavar="field delimiter. DEFAULT=TAB", default="\t")
parser.add_argument("-f", "--field", dest="field", default="2", metavar="Field of positions in file. DEFAULT=2")


args=parser.parse_args()
print(args)

FilterSnps(args)

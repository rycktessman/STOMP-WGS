#CREDIT: https://gitlab.com/tbgenomicsunit/ThePipeline; Iñaki Comas's lab. 

def MeanCoverage(prefix, depth4cov):
    ''' Calculate the mean coverage per pb in a sample from
    prefix.coverage file '''

    positions = 0.0
    coverage = 0
    covered_positions = 0.0
    # Add cov list to sort and calc median
    covlist = []
    with open("{}.coverage".format(prefix)) as infile:
        for line in infile:
            if line == "":
                positions = 0.0
                break
            ref, pos, cov = line.rstrip().split()
            cov  = int(cov)
            positions += 1.0
            coverage += cov
            if cov >= depth4cov:
                covered_positions += 1.0
            covlist.append(int(cov))

    if positions > 0:
        # Calc mean
        mean_cov = coverage / positions
        # Calc median
        # Recall that position X in a list is accesed with X-1
        # so for example the 6th element is accessed as covlist[5]
        covlist.sort()
        if len(covlist) % 2 == 0:
            m1 = int(len(covlist) / 2)
            cov1 = covlist[m1 - 1]
            cov2 = covlist[m1]
            median_cov = (cov1 + cov2) / 2.0
        else:
            m1 = int(len(covlist) / 2)
            # not neccessary tu sum 1 as the m1 +1 element is
            # accessed in fact with just m1
            median_cov = covlist[m1]

        genome_coverage = covered_positions / positions
        return (mean_cov, median_cov, genome_coverage)
    else:
        return (0,0,0)

def MeanSampleCoverage(prefix, keepcoverage, depth4cov):
    import os

    meancov, mediancov, genome_coverage = MeanCoverage(prefix, depth4cov)
    with open("{}.meancov".format(prefix), "w") as outfh:
        outfh.write("{}_mean_depth\t{}\n".format(prefix, meancov))
        outfh.write("{}_median_depth\t{}\n".format(prefix, mediancov))
        outfh.write("{}_genome_coverage\t{}\n".format(prefix, genome_coverage))

    # Remove .coverage files
    if not keepcoverage:
        os.remove("{}.coverage".format(prefix))

    return meancov, mediancov, genome_coverage

def FiltByCov(prefix, mean, median, cov, minmean, minmedian, mincov):
    '''If a sample do not pass coverage/depth thresholds, move all its files
    to NoPassCov folder'''
    import os, shutil
    from glob import glob

    if mean < minmean or median < minmedian or cov < mincov:
        # Create NoPassCov if it does not exist
        try:
            os.mkdir("NoPassCov")
        except OSError:
            pass

        files = glob("{}*".format(prefix))
        for file in files:
            shutil.move(file, "./NoPassCov/{}".format(file))

    return 0


def CalcCoverage(args):
    #from Repository import Programs, Data

    #data = Data()

    prefix = args.prefix
    reference = args.reference
    keepcoverage = args.keepcoverage
    depth4cov = args.depth4cov

    mean, median, cov = MeanSampleCoverage(prefix, keepcoverage, depth4cov)

    # Filter by Coverage in case it is specified
    if args.filter:
        FiltByCov(args.prefix, mean, median, cov,
                  args.minmean, args.minmedian, args.mincov)

    return 0


import argparse

parser=argparse.ArgumentParser(description="testing")
parser.add_argument("-p", "--prefix", dest="prefix")
parser.add_argument("-r", "--reference", dest="reference")
parser.add_argument("-d", "--depth4cov", dest="depth4cov", default=10)
parser.add_argument("-k", "--keepcoverage", dest="keepcoverage", default=1)
parser.add_argument("-f", "--filter", dest="filter")
parser.add_argument("-e", "--extension", dest="extension", default="bam")
parser.add_argument("-mmn", "--minmean", dest="minmean", default=0)
parser.add_argument("-mmd", "--minmedian", dest="minmedian", default=20)
parser.add_argument("-mc", "--mincov", dest="mincov", default=0.95)

args=parser.parse_args()
print(args)
print(args.prefix)
print(args.reference)
print(args.depth4cov)
print(args.keepcoverage)
print(args.filter)
print(args.extension)
CalcCoverage(args)
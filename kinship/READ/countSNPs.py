import pandas as pd
from argparse import ArgumentParser

parser = ArgumentParser(description="Count intersecting SNPs between READ individual pairs.")
parser.add_argument("-f", dest="filename", required=True,
                    help="READ_output_ordered file", metavar="FILE")
args = parser.parse_args()

df = pd.read_csv(args.filename, sep="\t")
result = pd.DataFrame({'SNPcount' : df.groupby("PairIndividuals")["SNVperWindow"].sum()}).reset_index()

result.to_csv('~/SNPCount.{}'.format(args.filename),sep = '\t', index=False)

print('End of run for {}'.format(args.filename))

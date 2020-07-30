import re
import csv

proteomepattern = re.compile(r'UP\d+(?=:)')
inputtree = snakemake.input.tree
proteomefile = snakemake.input.proteomes
newname = snakemake.config["proteome_name"]
outputtree = snakemake.output[0]

def removeIllegalChars(string):
    for char in "()[]-":
        string = string.replace(char, "")
    return string

with open(inputtree, 'rt') as treestream:
    tree = treestream.read(-1)
    treeproteomes = set(re.findall(proteomepattern, tree))

proteomedict = {}
with open(proteomefile, 'rt') as proteomestream:
    proteomereader = csv.DictReader(proteomestream, delimiter='\t', quotechar='"')
    proteomeinfo = [line for line in proteomereader if line['Proteome'] in treeproteomes]

for proteome in proteomeinfo:
    tree = tree.replace(proteome['Proteome'], removeIllegalChars(newname.format(**proteome)))

with open(outputtree, 'wt') as treewriter:
    treewriter.write(tree)
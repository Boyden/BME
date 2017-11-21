import re
regx = r"\|CY1\_(\w+)"
regx = re.compile(regx)
path = "C:\\Users\\acer\\Desktop\\uniprot-Cytochrome+C1.fasta\\CY1.fasta"
with open(path, "rt") as f:
    data = f.read()
newFilePath = "C:\\Users\\acer\\Desktop\\uniprot-Cytochrome+C1.fasta\\newCY1.fasta"
li = data.split(">")
for elem in li:
    match = regx.findall(elem)
    if len(match) != 0:
        fileName = match[0]
        templi = elem.split("\n")[1:]
        s = ">" + fileName + "\n" + "".join(templi) + "\n"
        with open(newFilePath, "at") as f:
            f.write(s)

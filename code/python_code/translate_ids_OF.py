# read the contents of dict file into a dictionary
with open("assets/species.dct", "r") as f:
    b_dict = dict(line.strip().split() for line in f)

# create an output file to write the translated values
with open("outputs/functional_annotation/TF_annotation/TFannot_orthofinder/orthogroups_longformat_translated.tsv", "w") as out:
    # iterate over the lines in gfams file
    with open("outputs/functional_annotation/TF_annotation/TFannot_orthofinder/orthogroups_longformat.tsv", "r") as f:
        for line in f:
            # split the line into two columns
            columns = line.strip().split("\t")
            # get the value from file "B" based on the key in column 2 of file "A"
            b_value = b_dict.get(columns[1], columns[1])
            # write the translated line to the output file
            out.write(columns[0] + "\t" + b_value + "\n")

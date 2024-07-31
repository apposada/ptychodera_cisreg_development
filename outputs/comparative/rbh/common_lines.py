import sys

if len(sys.argv) != 4:
    print("Usage: python common_lines.py <file1> <file2> <outfile>")
    sys.exit(1)

file1 = sys.argv[1]
file2 = sys.argv[2]
outfile = sys.argv[3]

with open(file1) as f1:
    with open(file2) as f2:
        with open(outfile, "w") as output:
            # read lines from both files
            lines1 = set(f1.readlines())
            lines2 = set(f2.readlines())
            # find the common lines
            common_lines = lines1.intersection(lines2)
            # write the common lines to the output file
            output.writelines(common_lines)


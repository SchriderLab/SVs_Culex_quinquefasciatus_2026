import pandas as pd


df = pd.read_csv("/proj/dschridelab/remi/ANNA_RELAB/02_2026_FINAL2_AT/08_SVs/05_FINAL_FILES/01_INT_FILES/cg_chr3_99.txt", sep="\t")


locations = [int(x) for x in df["location"]]

if not locations:
    exit()

ranges = []
start = locations[0]
end = locations[0]

for loc in locations[1:]:
    if loc - end == 1000:
        end = loc
    else:
	ranges.append((start, end))
        start = end = loc

ranges.append((start, end))


for start, end in ranges:
    print(f"{start}\t{end}")

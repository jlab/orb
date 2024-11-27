import re
import pandas as pd
import sys

def parse_bowtie2_log(filepath):
    stats = {
        "total_reads": None,
        "paired_reads": None,
        "concordant_0": None,
        "concordant_1": None,
        "concordant_gt1": None,
        "discordant_1": None,
        "unaligned_pairs": None,
        "mates_unaligned_0": None,
        "mates_unaligned_1": None,
        "mates_unaligned_gt1": None,
        "overall_alignment_rate": None,
    }

    with open(filepath, 'r') as file:
        for line in file:
            line = line.strip()
            if match := re.match(r"^(\d+) reads; of these:", line):
                stats["total_reads"] = int(match.group(1))
            elif match := re.match(r"^(\d+) \([\d.]+%\) were paired", line):
                stats["paired_reads"] = int(match.group(1))
            elif match := re.match(r"^(\d+) \([\d.]+%\) aligned concordantly 0 times", line):
                stats["concordant_0"] = int(match.group(1))
            elif match := re.match(r"^(\d+) \([\d.]+%\) aligned concordantly exactly 1 time", line):
                stats["concordant_1"] = int(match.group(1))
            elif match := re.match(r"^(\d+) \([\d.]+%\) aligned concordantly >1 times", line):
                stats["concordant_gt1"] = int(match.group(1))
            elif match := re.match(r"^(\d+) pairs aligned concordantly 0 times; of these:", line):
                stats["unaligned_pairs"] = int(match.group(1))
            elif match := re.match(r"^(\d+) \([\d.]+%\) aligned discordantly 1 time", line):
                stats["discordant_1"] = int(match.group(1))
            elif match := re.match(r"^(\d+) mates make up the pairs; of these:", line):
                pass
            elif match := re.match(r"^(\d+) \([\d.]+%\) aligned 0 times", line):
                stats["mates_unaligned_0"] = int(match.group(1))
            elif match := re.match(r"^(\d+) \([\d.]+%\) aligned exactly 1 time", line):
                stats["mates_unaligned_1"] = int(match.group(1))
            elif match := re.match(r"^(\d+) \([\d.]+%\) aligned >1 times", line):
                stats["mates_unaligned_gt1"] = int(match.group(1))
            elif match := re.match(r"([\d.]+)% overall alignment rate", line):
                stats["overall_alignment_rate"] = float(match.group(1))

    return stats

log_file = sys.argv[1]
colname = sys.argv[2]

logs_dict = parse_bowtie2_log(log_file)
logs_df = pd.DataFrame.from_dict(logs_dict, orient='index')
logs_df.columns = [colname]
logs_df.to_csv(sys.stdout, sep='\t')
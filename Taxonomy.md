### **Extracting Taxonomy Information Using Diamond Blast Report (IMG/VR Database)**

**Note**: Run this command in the root folder where the SOVAP pipeline is executed.

This command processes all directories in the current folder. For each directory, it extracts data from `output.diamond.tsv`, joins it with taxonomy information from the `IMGVR_all_Sequence_information.tsv` file, and outputs a `.taxo` file with the extracted taxonomy data.

```bash
for folder in */; do \
    parent_dir=${folder%/}; \
    cut -f1,2 ${folder}6_Diamond-Taxonomy/output.diamond.tsv | cut -d "|" -f1 | \
        awk '{print $2 "\t" $1}' | \
        sort -k1,1 | \
        join -a1 -t $'\t' - /path/to/IMGVR_all_Sequence_information_sorted.tsv > \
        ${folder}6_Diamond-Taxonomy/$parent_dir.taxo; \
done
```

**Details**:

1. The script loops through each subdirectory in the root folder.
2. For each directory:
   - **Extract relevant columns** from `output.diamond.tsv` (column 1 and 2).
   - **Clean the taxonomic information**, removing everything after the first `|` symbol.
   - **Format** the output by swapping the order of columns, and sort the data by the first column.
   - **Join** the extracted data with the `IMGVR_all_Sequence_information.tsv` file, which contains sequence metadata. The `-a1` option ensures that all lines from the first file (the blast output) are included, even if no match is found in the second file.
3. **Save the output** to a file called `foldername.taxo` within the `6_Diamond-Taxonomy/` subdirectory of each folder.

The resulting `.taxo` file will contain taxonomic details for each sequence in the corresponding `output.diamond.tsv` file, ready for downstream analysis.

  
**Important**: 
- **Ensure the `IMGVR_all_Sequence_information.tsv` file is pre-sorted by the first column** (sequence ID) for the join operation to work correctly. Using an unsorted file may result in missing matches or errors. If the file is not sorted, you can pre-sort it as follows:

```bash
sort -k1,1 /path/to/IMGVR_all_Sequence_information.tsv -o /path/to/IMGVR_all_Sequence_information_sorted.tsv
```

Alternatively, **if the file is not pre-sorted and sorting it is not feasible due to its size**, you can sort the file during execution using the following command. However, this is **not recommended** when working with many directories, as it significantly increases runtime:

```bash
for folder in */; do \
    parent_dir=${folder%/}; \
    cut -f1,2 ${folder}6_Diamond-Taxonomy/output.diamond.tsv | cut -d "|" -f1 | \
        awk '{print $2 "\t" $1}' | \
        sort -k1,1 | \
        join -a1 -t $'\t' - <(sort -k1,1 /path/to/IMGVR_all_Sequence_information.tsv) > \
        ${folder}6_Diamond-Taxonomy/$parent_dir.taxo; \
done
```

**Explanation**: 
- In this command, the `IMGVR_all_Sequence_information.tsv` file is sorted on-the-fly during the join operation. This adds overhead, particularly when processing many directories, as the file is sorted every time the command is run.

**Due to the large size of IMGVR_all_Sequence_information.tsv file, this step is provided as an optional side script. However, it is recommended to perform this step for reproducing graphs in the manuscript and analyzing data in R for taxonomy, diversity indices, and etc.**  

---------------------
IMGVR_all_Sequence_information.tsv:
a table listing the characteristics of each viral sequence such as its origin, affiliation, and predicted host (tsv format).  

More information: 
https://genome.jgi.doe.gov/portal/IMG_VR/IMG_VR.home.html

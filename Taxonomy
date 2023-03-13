**To extract taxonomy information using Diamond blast report (when using IMG/VR db) and using IMGVR_all_Sequence_information.tsv file:**

**Note: It should be run in the root folder where the SOVAP pipeline is executed.**

This is a shell script that loops over all directories (indicated by */) in the current directory, extracts information from a file called output.diamond.tsv, and joins it with information from a file called IMGVR_all_Sequence_information.tsv. 

```
for folder in */; do \
    parent_dir=${folder%/}; \
    cut -f1,2 ${folder%}6_Diamond-Taxonomy/output.diamond.tsv | cut -d "|" -f1 | \
        paste <(cut -f2) <(cut -f1) | \
        sort -k1,1 | \
        join -t $'\t' -a1 - /path/to/IMGVR_all_Sequence_information.tsv > ${folder%}6_Diamond-Taxonomy/$parent_dir.taxo ; \
done
```

_Details:_
The script loops through each subdirectory in the root folder and performs the following steps for each subdirectory:

- Define the `parent_dir` variable as the subdirectory name without the trailing slash.
- Extract the first and second columns from the output.diamond.tsv file using `cut` command and separate taxonomic information up to the first "|" character using cut again.
- Merge the taxonomic information into a single column using the paste command.
- Sort the merged column using `sort` command.
- Use the `join` command to join the sorted merged column with IMGVR_all_Sequence_information.tsv file. The -t option specifies the tab as the delimiter, and the `-a1` option tells join to print unpairable lines from the first file. The `-` character specifies to use standard input as the first file for `join`.
- Write the output to a file named after the subdirectory and with the extension .taxo.
- The resulting output files will contain taxonomic information for each sequence in the subdirectory's output.diamond.tsv file.

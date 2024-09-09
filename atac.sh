#!/bin/bash

. CONFIG_ATAC

create_dirs() {
    local dirs=("$@")
    for dir in "${dirs[@]}"; do
        mkdir -p "${dir}"
    done
}

TEMP="$ProjDir/test"
Data="$ProjDir/test"
ana="$ProjDir/Results"
MACS2_PEAKS="$ProjDir/macs2_peaks"
FILTERED_PEAKS="$ProjDir/macs2_peaks/filtered_p10"
MOTIFS="$Data/motifs"
FOOTPRINTING="$Data/footprinting"
BINDETECTOUT="$Data/binDetect"

create_dirs "${TEMP}" "${Data}"  "${ana}" \
            "${MACS2_PEAKS}" "${FILTERED_PEAKS}" "${LOSS_MOTIFS}" \
            "${GAIN_MOTIFS}" "${FOOTPRINTING}" "${BINDETECT}"

generate_process_script() {
    cat > runProcess.sh <<EOI
#!/bin/bash
#SBATCH --job-name=ATAC_seq_analysis
#SBATCH --chdir=${ProjDir}
#SBATCH --output=${ProjDir}/ATAC.out.log
#SBATCH --error=${ProjDir}/ATAC.err.log
#SBATCH --mail-type=FAIL
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=200GB
srun preprocessing.sh
srun peakCalling.sh
srun heatmap.sh
srun transcriptionFactors.sh
srun footprinting.sh
EOI
}

declare -A samples_by_condition

# while IFS=$'\t' read -r sampleID condition; do
#     samples_by_condition["$condition"]+="$sampleID "
# done < sampleID.txt

generate_process_script

cat > preprocessing.sh <<EOI
#!/bin/bash

while IFS=$'\t' read -r sampleID condition; do
    samples_by_condition["\$condition"]+="\$sampleID "
done < sampleID.txt

# Iterate over conditions and process each set of samples
for condition in "\${!samples_by_condition[@]}"; do
    sampleIDs="\${samples_by_condition[\${condition}]}"
    replicate_num=1

    # FastQC for samples under the current condition
    for sampleID in \$sampleIDs; do
        ${FASTQC} ${rawData}/\${sampleID}_R1.fastq.gz -o ${Data}
        ${FASTQC} ${rawData}/\${sampleID}_R2.fastq.gz -o ${Data}
    # done

    # Trimming for samples under the current condition
    # for sampleID in \$sampleIDs; do
        ${TRIM} ${rawData}/\${sampleID}_R1.fastq.gz ${rawData}/\${sampleID}_R2.fastq.gz -o ${TEMP}
    # done

    # FastQC after trimming
    # for sampleID in \$sampleIDs; do
        ${FASTQC} ${TEMP}/\${sampleID}_R1_val_1.fq.gz -o ${Data}
        ${FASTQC} ${TEMP}/\${sampleID}_R2_val_2.fq.gz -o ${Data}
    # done

    # Mapping for samples under the current condition
    # for sampleID in \$sampleIDs; do
        bam_file="${Data}/\${condition}_rep\${replicate_num}.bam"
        ${MAP} -x ${genomeDir}/mm10_bt2 -1 ${TEMP}/\${sampleID}_R1_val_1.fq.gz -2 ${TEMP}/\${sampleID}_R2_val_2.fq.gz \
            | ${SAMVIEW} -bS - > "\$bam_file"
    # done

    # BAM file stats
    stats_file="${Data}/\${condition}_rep\${replicate_num}_stats.txt"
    ${STATS} "\${bam_file}" > "\${stats_file}"

    # Sorting the BAM files
    sort_prefix="${TEMP}/\${condition}_rep\${replicate_num}.sort"
    sort_bam="${TEMP}/\${condition}_rep\${replicate_num}.sort.bam"
    ${SAMSORT} "\${bam_file}" -T "\${sort_prefix}" -o "${sort_bam}"

    # Remove duplicates
    rmdup="${TEMP}/\${condition}_rep\${replicate_num}.sort.rmdup.bam"
    dup_metrics="${Data}/\${condition}_rep\${replicate_num}.rmdup_metrics.txt"
    ${PICARD} INPUT="\${sort_bam}" OUTPUT="\${rmdup}" METRICS_FILE="\${dup_metrics}"
    
    # Remove mitochondrial reads
    rmmito="${TEMP}/\${condition}_rep\${replicate_num}.sort.rmdup.rmmito.bam"
    ${SAMVIEW} -o "\${rmmito}" -e 'rname != "chrM" && rname !~ "_random"' "\${rmdup}"

    # Remove blacklist
    rmblacklist="${TEMP}/\${condition}_rep\${replicate_num}.sort.rmdup.rmmito.rmblacklist.bam"
    ${BEDTOOLS} -abam "\${rmmito}" -b "\${blacklist}" -v > "\${rmblacklist}"

    # Filter reads by quality
    filtered="${Data}/\${condition}_rep\${replicate_num}.filtered.bam"
    ${SAMVIEW} -b -f 2 -F 4 -q 30 "\${rmblacklist}" -o "\${filtered}"

    # Make index file
    index="${Data}/\${condition}_rep\${replicate_num}.sort.rmdup.rmmito.rmblacklist.bai"
    ${SAMINDEX} "\${filtered}" "\${index}"
    replicate_num=\$((replicate_num + 1))
done

# Cleanup
#rm ${TEMP}/*.fq.gz

EOI

cat > peakCalling.sh <<EOI
#!/bin/bash

while IFS=$'\t' read -r sampleID condition; do
    samples_by_condition["\$condition"]+="\$sampleID "
done < sampleID.txt

for condition in "\${!samples_by_condition[@]}"; do
    # sampleIDs="\${samples_by_condition[\${condition}]}"
    replicate_num=1
    
    # MACS2 sample-wise Peak calling
    filtered="${Data}/\${condition}_rep\${replicate_num}.filtered.bam"
    ${MACS2} -t "\${filtered}" -f BAM -g 2.6e9 -q 0.01 --outdir ${MACS2_PEAKS} -n "\${condition}_rep\${replicate_num}"

    # Filter peaks by fold enrichment and p-value
    peak_file="${FILTERED_PEAKS}/\${condition}_rep\${replicate_num}_peaks.bed"
    ${TAIL} "${MACS2_PEAKS}/\${condition}_rep\${replicate_num}_peaks.xls" | LC_ALL=C awk '${AWK}' > "\${peak_file}"
    replicate_num=\$((replicate_num + 1))
done

# Merge the replicates
for condition in "\${!samples_by_condition[@]}"; do
    ${BEDTOOLS} -a "${FILTERED_PEAKS}/\${condition}_rep1_peaks.bed" -b "${FILTERED_PEAKS}/\${condition}_rep2_peaks.bed" -f 0.5 -wa -wb > "${FILTERED_PEAKS}/\${condition}_merged.bed"
done

${BEDTOOLS} -a "condition1_merged" -b "condition2_merged" -f 0.5 -wa -wb > ${commonPeaks}

for condition in "\${!samples_by_condition[@]}"; do
    ${BEDTOOLS} -a "${FILTERED_PEAKS}/\${condition}_rep1_peaks.bed" -b "${FILTERED_PEAKS}/\${condition}_rep2_peaks.bed" -f 0.5 -wa -wb > "${FILTERED_PEAKS}/\${condition}_merged.bed"
done

for condition1 in "\${!samples_by_condition[@]}"; do
    for condition2 in "\${!samples_by_condition[@]}"; do
        if [ "\${condition1}" != "\${condition2}" ]; then
            condition1_merged="${FILTERED_PEAKS}/\${condition1}_merged.bed"
            condition2_merged="${FILTERED_PEAKS}/\${condition2}_merged.bed"
            ${BEDTOOLS} -a "\${condition1_merged}" -b "\${condition2_merged}" -f 0.5 -wa -wb > "${FILTERED_PEAKS}/\${condition1}_\${condition2}_common_peaks.bed"
            ${MACS2} -c "${Data}/\${condition1}_rep1.filtered.bam" "${Data}/\${condition1}_rep2.filtered.bam" -t "${Data}/\${condition2}_rep1.filtered.bam" "${Data}/\${condition2}_rep2.filtered.bam" -f BAM -g mm -q 0.01 --outdir ${MACS2_PEAKS} -n "Loss"
            ${MACS2} -c "${Data}/\${condition2}_rep1.filtered.bam" "${Data}/\${condition2}_rep2.filtered.bam" -t "${Data}/\${condition1}_rep1.filtered.bam" "${Data}/\${condition1}_rep2.filtered.bam" -f BAM -g mm -q 0.01 --outdir ${MACS2_PEAKS} -n "Gain"
        fi
    done
done

# Filter peaks by fold enrichment and p-value
${TAIL} ${MACS2_PEAKS}/Loss_peaks.xls | LC_ALL=C awk '${AWK}' > "${FILTERED_PEAKS}/Loss_peaks.bed"
${TAIL} ${MACS2_PEAKS}/Gain_peaks.xls | LC_ALL=C awk '${AWK}' > "${FILTERED_PEAKS}/Gain_peaks.bed"

# Filter peaks overlapping with common peaks
${BEDTOOLS} -a ${FILTERED_PEAKS}/Loss_peaks.bed -b "${FILTERED_PEAKS}/*_common_peaks.bed" -f 0.5 -v > ${FILTERED_PEAKS}/Loss_macs2.bed
${BEDTOOLS} -a ${FILTERED_PEAKS}/Gain_peaks.bed -b "${FILTERED_PEAKS}/*_common_peaks.bed" -f 0.5 -v > ${FILTERED_PEAKS}/Gain_macs2.bed

# DiffBind
Rscript diffBindAnalysis.R

# Merge differential peaks from MACS2 and DiffBind
cat ${FILTERED_PEAKS}/Loss_macs2.bed ${FILTERED_PEAKS}/Loss_diffBind.bed | cut -f1-3 | sort -k1,1 -k2,2n | ${MERGE} -i – > ${FILTERED_PEAKS}/Loss.bed
cat ${FILTERED_PEAKS}/Gain_macs2.bed ${FILTERED_PEAKS}/Gain_diffBind.bed | cut -f1-3 | sort -k1,1 -k2,2n | ${MERGE} -i – > ${FILTERED_PEAKS}/Gain.bed

EOI

cat > heatmap.sh<<EOI

for condition in "\${!samples_by_condition[@]}"; do
    replicate_num=1
    bw_file="${Data}/\${condition}_rep\${replicate_num}.bw"
    ${COVERAGE} -b "${Data}/\${condition1}_rep1.filtered.bam" -o "\${bw_file}" --normalizeUsing RPKM --effectiveGenomeSize 2652783500 -p 12
    replicate_num=\$((replicate_num + 1))
done

for condition1 in "\${!samples_by_condition[@]}"; do
    for condition2 in "\${!samples_by_condition[@]}"; do
        if [ "\${condition1}" != "\${condition2}" ]; then
            matrix="${Data}/\${condition1}_\${condition2}_matrix.tab.gz"
            heatmapPdf="${Data}/\${condition1}_\${condition2}_heatmap.pdf"
            ${COMPUTEMATRIX} -S "${Data}/\${condition1}_rep1.bw" "${Data}/\${condition1}_rep2.bw" "${Data}/\${condition2}_rep1.bw" "${Data}/\${condition2}_rep2.bw"  \
                -R ${FILTERED_PEAKS}/Loss.bed ${FILTERED_PEAKS}/Gain.bed -a 3000 -b 3000 -out "\${matrix}"
            ${PLOTHEATMAP} -m "\${matrix}" -out "\${heatmapPdf}" --regionsLabel "Loss" "Gain"
        fi
    done
done

EOI

cat > transcriptionFactors.sh<<EOI

for status in Loss Gain; do
    # HOMER motif prediction
    ${FINDMOTIFS} ${FILTERED_PEAKS}/\${status}.bed mm10 "\${MOTIFS}/${status}"

    # TF co-occurence
    ${BEDTOOLS} -wb -wa -a ${chipAtlas} -b "${FILTERED_PEAKS}/\${status}.bed" -F 0.5 | cut -f 1-4,10-12 | awk -v OFS="\t" '{ print $5,$6,$7,$1,$2,$3,$4 }' | cut -d';' -f 1,2,4 | sed 's/%20/ /g;s/;/\t/g' | awk '!seen[$4,$6,$7]++' > "${Data}/\${status}_TF.bed"

    # Plot for TF co-occurence
    ${RScript} ChIP-Atlas_Viz.R ${Data}/\${status}_TF.bed "${FILTERED_PEAKS}/\${status}.bed"

    # WhichTF
    ${WHICHTF} "${FILTERED_PEAKS}/\${status}.bed" mm10 --outfile "${Data}/\${status}_whichTF.tsv"

done

EOI

cat > footprinting.sh<<EOI

# Merge peaks replicate-wise
for condition in "\${!samples_by_condition[@]}"; do
    cat "${FILTERED_PEAKS}/\${condition}_rep1_peaks.bed" "${FILTERED_PEAKS}/\${condition}_rep1_peaks.bed" | cut -f1-3 | ${SORT} | ${MERGE} -d 5 | awk '{ print $0,"\${condition}" }' OFS='\t' > "${FILTERED_PEAKS}/\${condition}_union.bed"
    ${SAMMERGE} -o "${Data}/\${condition}_merged.bam" "${Data}/\${condition}_rep1.filtered.bam" "${Data}/\${condition}_rep2.filtered.bam"
done

# Merge peaks sample-wise
for condition1 in "\${!samples_by_condition[@]}"; do
    for condition2 in "\${!samples_by_condition[@]}"; do
        if [ "\${condition1}" != "\${condition2}" ]; then
            cat ${FILTERED_PEAKS}/\${condition1}_union.bed" "${FILTERED_PEAKS}/\${condition2}_union.bed" | sort -k1,1 -k2,2n | ${MERGE} -d 5 -c 4 -o distinct > "${FILTERED_PEAKS}/peaks_union.bed"
        fi
    done
done

for condition in "\${!samples_by_condition[@]}"; do
    corrected="${FOOTPRINTING}/\${condition}_corrected.bw"
    ${CORRECT} --bam "${Data}/\${condition}_merged.bam" --genome ${genome} --peaks "${FILTERED_PEAKS}/\${condition}_union.bed" --outdir ${FOOTPRINTING}
    ${FOOTPRINT} --signal "\${corrected}" --regions "${FILTERED_PEAKS}/peaks_union.bed" --output "${FOOTPRINTING}/\${condition}_footprints.bw"
done
TOBIAS BINDetect --signals ${ctrlfootprint} ${treatfootprint} --genome ${genome} --peaks ${unionPeak} --motifs ${motifs} --outdir ${BINDETECTOUT} --cond_names ${ctrlsampleID} ${treatsampleID} --cores 12

for condition1 in "\${!samples_by_condition[@]}"; do
    for condition2 in "\${!samples_by_condition[@]}"; do
        if [ "\${condition1}" != "\${condition2}" ]; then
            ${BINDETECT}  --signals "${FOOTPRINTING}/\${condition1}_footprints.bw" "${FOOTPRINTING}/\${condition2}_footprints.bw" \
                --genome ${genome} --peaks "${FILTERED_PEAKS}/peaks_union.bed" --motifs ${motifs} --outdir ${BINDETECTOUT} --cond_names "\${condition1}" "\${condition2}"
        fi
    done
done

EOI

chmod a+x *.sh
sbatch runProcess.sh
#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


date=new Date().format( 'yyMMdd' )
user="$USER"
runID="${date}.${user}"


//////////////////////////// SWITCHES ///////////////////////////////// 

switch (params.gatk) {

    case 'danak':
    gatk_image="gatk419.sif";
    break;
    case 'new':
    gatk_image="gatk4400.sif";
    break;
    default:
    gatk_image="gatk4400.sif";
    break;
}




switch (params.server) {
    case 'lnx01':
        s_bind="/data/:/data/,/lnx01_data2/:/lnx01_data2/";
        simgpath="/data/shared/programmer/simg";
        params.intervals_list="/data/shared/genomes/hg38/interval.files/WGS_splitIntervals/wgs_splitinterval_BWI_subdivision3/*.interval_list";
        tmpDIR="/data/TMP/TMP.${user}/";
        gatk_exec="singularity run -B ${s_bind} ${simgpath}/${gatk_image} gatk";
        refFilesDir="/data/shared/genomes";
    break;
    default:
        s_bind="/data/:/data/,/lnx01_data2/:/lnx01_data2/,/fast/:/fast/,/lnx01_data3/:/lnx01_data3/,/lnx01_data4/:/lnx01_data4/";
        simgpath="/data/shared/programmer/simg";
        params.intervals_list="/data/shared/genomes/hg38/interval.files/WGS_splitIntervals/wgs_splitinterval_BWI_subdivision3/*.interval_list";
        tmpDIR="/fast/TMP/TMP.${user}/";
        gatk_exec="singularity run -B ${s_bind} ${simgpath}/${gatk_image} gatk";
        refFilesDir="/fast/shared/genomes";
    break;
}
/*
switch ($user) {
    case 'mmaj':
        permissions="full";
    break;
    
    case 'raspau':
        permissions="full";
    break;
    
    default:
        permissions="reduced";
    break;

}
*/
switch (params.genome) {
    case 'hg19':
        assembly="hg19"
        // Genome assembly files:
        genome_fasta = "/data/shared/genomes/hg19/human_g1k_v37.fasta"
        genome_fasta_fai = "/data/shared/genomes/hg19/human_g1k_v37.fasta.fai"
        genome_fasta_dict = "/data/shared/genomes/hg19/human_g1k_v37.dict"
        genome_version="V1"
        break;


    case 'hg38':
        assembly="hg38"
        // Genome assembly files:
        if (params.hg38v1) {
        genome_fasta = "/data/shared/genomes/hg38/GRCh38.primary.fa"
        genome_fasta_fai = "/data/shared/genomes/hg38/GRCh38.primary.fa.fai"
        genome_fasta_dict = "/data/shared/genomes/hg38/GRCh38.primary.dict"
        genome_version="hg38v1"
        cnvkit_germline_reference_PON="/data/shared/genomes/hg38/inhouse_DBs/hg38v1_primary/cnvkit/wgs_germline_PON/jgmr_45samples.reference.cnn"
        cnvkit_inhouse_cnn_dir="/data/shared/genomes/hg38/inhouse_DBs/hg38v1_primary/cnvkit/wgs_persample_cnn/"
        inhouse_SV="/data/shared/genomes/hg38/inhouse_DBs/hg38v1_primary/"
        }
        
        if (params.hg38v2){
        genome_fasta = "/data/shared/genomes/hg38/ucsc.hg38.NGS.analysisSet.fa"
        genome_fasta_fai = "/data/shared/genomes/hg38/ucsc.hg38.NGS.analysisSet.fa.fai"
        genome_fasta_dict = "/data/shared/genomes/hg38/ucsc.hg38.NGS.analysisSet.dict"
        genome_version="hg38v2"
        }
        // Current hg38 version (v3): NGC with masks and decoys.
        if (!params.hg38v2 && !params.hg38v1){
        genome_fasta = "/data/shared/genomes/hg38/GRCh38_masked_v2_decoy_exclude.fa"
        genome_fasta_fai = "/data/shared/genomes/hg38/GRCh38_masked_v2_decoy_exclude.fa.fai"
        genome_fasta_dict = "/data/shared/genomes/hg38/GRCh38_masked_v2_decoy_exclude.dict"
        genome_version="hg38v3"
        cnvkit_germline_reference_PON="/data/shared/genomes/hg38/inhouse_DBs/hg38v3_primary/cnvkit/hg38v3_109samples.cnvkit.reference.cnn"
        cnvkit_inhouse_cnn_dir="/data/shared/genomes/hg38/inhouse_DBs/hg38v3_primary/cnvkit/wgs_persample_cnn/"
        inhouse_SV="/data/shared//genomes/hg38/inhouse_DBs/hg38v3_primary/"
        }


        AV1_ROI="/data/shared/genomes/${params.genome}/interval.files/panels/av1.hg38.ROI.v2.bed"
        CV1_ROI="/data/shared/genomes/${params.genome}/interval.files/panels/cv3.hg38.ROI.bed"
        CV2_ROI="/data/shared/genomes/${params.genome}/interval.files/panels/cv3.hg38.ROI.bed"
        CV3_ROI="/data/shared/genomes/${params.genome}/interval.files/panels/cv3.hg38.ROI.bed"
        CV4_ROI="/data/shared/genomes/${params.genome}/interval.files/panels/cv4.hg38.ROI.bed"
        CV5_ROI="/data/shared/genomes/${params.genome}/interval.files/panels/cv5.hg38.ROI.bed"
        GV3_ROI="/data/shared/genomes/${params.genome}/interval.files/panels/gv3.hg38.ROI.v2.bed"
        NV1_ROI="/data/shared/genomes/${params.genome}/interval.files/panels/nv1.hg38.ROI.bed"
        WES_ROI="/data/shared/genomes/hg38/interval.files/exome.ROIs/211130.hg38.refseq.gencode.fullexons.50bp.SM.bed"
        MV1_ROI="/data/shared/genomes/${params.genome}/interval.files/panels/muc1.hg38.coordinates.bed"
 
        break;
}

dataStorage="/lnx01_data3/storage/"
multiqc_config="/data/shared/programmer/configfiles/multiqc_config_tumorBoard.yaml"
/*
if (!params.useBasesMask && !params.RNA) {
  dnaMask="Y*,I8nnnnnnnnn,I8,Y*"
}
if (!params.useBasesMask && params.RNA) {
  dnaMask="Y*,I8nnnnnnnnnnn,I8nn,Y*"
  rnaMask="Y*,I10nnnnnnnnn,I10,Y*"
}

if (params.useBasesMask) {
  dnaMask=params.useBasesMask
  params.DNA=true
}
*/

if (params.miniseq) {
    umiConvertDNA = "Y151;I8;I8;Y151"
    umiConvertRNA = "Y151;I8;I8;Y151"
}
else if (params.useBasesMask) {
    umiConvertDNA = params.useBasesMask
    umiConvertRNA = params.useBasesMask
}
else if (params.RNA) {
    umiConvertDNA = "Y151;I8N2U9;I8N2;Y151"
    umiConvertRNA = "Y151;I10U9;I10;Y151"
}
else {
    umiConvertDNA = "Y151;I8U9;I8;Y151"
}

if (!params.DNA && !params.miniseq && !params.useBasesMask) {
    umiConvertDNA = "Y151;I10U9;I10;Y151"
}


if (params.localStorage ) {
aln_output_dir="${params.outdir}/"
fastq_dir="${params.outdir}/"
qc_dir="${params.outdir}/QC/"
}

if (!params.localStorage) {
aln_output_dir="${dataStorage}/alignedData/${params.genome}/novaRuns/2026/"
fastq_dir="${dataStorage}/fastqStorage/novaRuns/"
qc_dir="${dataStorage}/fastqStorage/demultiQC/2026"
}


/////////////////////////////////////////////////////////////////////////////
////////////////////////// DEMULTI PROCESSES: ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////

process prepare_DNA_samplesheet {

    input:
    tuple val(samplesheet_basename), path(samplesheet)
    path(runinfo)

    output:
    path("*.DNA_SAMPLES.csv"), emit: std
    path("*.UMI.csv"), emit: umi

    script:
    """
    # ---------------------------------------------------------
    # Keep DNA samples only
    # ---------------------------------------------------------

    grep -v "RV1" ${samplesheet} > ${samplesheet_basename}.DNA_SAMPLES.csv


    # ---------------------------------------------------------
    # Read run structure from RunInfo.xml
    #
    # DNA layout:
    #   Index 1 = 8 bp sample index + optional N + 9 bp UMI
    #   Index 2 = 8 bp sample index + optional N
    #
    # Examples:
    #   17 / 8  -> I8U9 / I8
    #   19 / 10 -> I8N2U9 / I8N2
    # ---------------------------------------------------------

    python3 - <<'PY'
import csv
import xml.etree.ElementTree as ET
from pathlib import Path

samplesheet = Path("${samplesheet_basename}.DNA_SAMPLES.csv")
output      = Path("${samplesheet_basename}.DNA_SAMPLES.UMI.csv")
runinfo     = Path("${runinfo}")


# ---------------------------------------------------------
# Parse RunInfo.xml
# ---------------------------------------------------------

root = ET.parse(runinfo).getroot()

reads = root.findall(".//Read")

if len(reads) < 4:
    raise RuntimeError(
        f"Expected at least 4 reads in RunInfo.xml, found {len(reads)}"
    )

r1 = int(reads[0].attrib["NumCycles"])
i1 = int(reads[1].attrib["NumCycles"])
i2 = int(reads[2].attrib["NumCycles"])
r2 = int(reads[3].attrib["NumCycles"])


# Determine DNA mask from the actual run structure.
#
# Supported layouts:
#   8 / 8   -> I8 / I8               (no UMI; e.g. MiniSeq/non-UMI runs)
#   17 / 8  -> I8U9 / I8             (DNA UMI)
#   19 / 10 -> I8N2U9 / I8N2         (DNA UMI with 2 spacer cycles)
#
# A manually supplied --useBasesMask still takes precedence.

forced_mask = "${params.useBasesMask ?: ''}".strip()

if forced_mask:
    override_cycles = forced_mask

elif i1 == 8 and i2 == 8:
    override_cycles = f"Y{r1};I8;I8;Y{r2}"

elif i1 >= 17 and i2 >= 8:
    pad1 = i1 - 17   # 8 sample-index + pad + 9 UMI
    pad2 = i2 - 8    # 8 sample-index + optional pad

    i1_mask = "I8"
    if pad1:
        i1_mask += f"N{pad1}"
    i1_mask += "U9"

    i2_mask = "I8"
    if pad2:
        i2_mask += f"N{pad2}"

    override_cycles = f"Y{r1};{i1_mask};{i2_mask};Y{r2}"

else:
    raise RuntimeError(
        f"Unsupported DNA index lengths from RunInfo.xml: "
        f"I1={i1}, I2={i2}"
    )

print(
    f"DNA RunInfo: R1={r1}, I1={i1}, I2={i2}, R2={r2}"
)
print(
    f"DNA OverrideCycles: {override_cycles}"
)


# ---------------------------------------------------------
# Instrument
# ---------------------------------------------------------

instrument_node = root.find(".//Instrument")
instrument = (
    instrument_node.text.strip()
    if instrument_node is not None and instrument_node.text
    else ""
)

# Your MN02212 data show that the index pair reaching
# BCL Convert is:
#
#   RC(index2), RC(index)
#
# Therefore correct the DNA samplesheet only for the MN 19/10 layout
# that was verified from Top_Unknown_Barcodes.csv.  Do not apply this
# transformation blindly to 8/8 or 17/8 runs.

swap_revcomp = instrument.upper().startswith("MN")

print(f"Instrument: {instrument or 'UNKNOWN'}")
print(f"Swap/reverse-complement DNA indexes: {swap_revcomp}")


def revcomp(seq):
    table = str.maketrans(
        "ACGTNacgtn",
        "TGCANtgcan"
    )
    return seq.translate(table)[::-1]


# ---------------------------------------------------------
# Read samplesheet as raw rows
# ---------------------------------------------------------

with samplesheet.open(newline="") as fh:
    rows = list(csv.reader(fh))


# Find Data header containing Sample_ID/index/index2
header_idx = None

for n, row in enumerate(rows):
    stripped = [x.strip() for x in row]

    if (
        "Sample_ID" in stripped
        and "index" in stripped
        and "index2" in stripped
    ):
        header_idx = n
        break

if header_idx is None:
    raise RuntimeError(
        "Could not locate Sample_ID/index/index2 header in samplesheet"
    )


header = rows[header_idx]

sample_col = header.index("Sample_ID")
index_col  = header.index("index")
index2_col = header.index("index2")


# ---------------------------------------------------------
# Transform MN index orientation
# ---------------------------------------------------------

if swap_revcomp:

    print("Converting DNA indexes:")
    print("  new index  = RC(old index2)")
    print("  new index2 = RC(old index)")

    for row in rows[header_idx + 1:]:

        if not row or len(row) <= max(index_col, index2_col):
            continue

        # Skip section headers / empty data
        if not row[sample_col].strip():
            continue

        old_i1 = row[index_col].strip()
        old_i2 = row[index2_col].strip()

        if not old_i1 or not old_i2:
            continue

        new_i1 = revcomp(old_i2)
        new_i2 = revcomp(old_i1)

        print(
            f"  {row[sample_col]}: "
            f"{old_i1}/{old_i2} -> {new_i1}/{new_i2}"
        )

        row[index_col]  = new_i1
        row[index2_col] = new_i2


# ---------------------------------------------------------
# Insert / replace BCLConvert settings
# ---------------------------------------------------------

settings = [
    ["OverrideCycles", override_cycles],
    ["NoLaneSplitting", "true"],
]

# TrimUMI is only valid when OverrideCycles actually defines UMI bases.
# BCL Convert 4.3.6 aborts if TrimUMI is present for a non-UMI mask
# such as Y151;I8;I8;Y151.
if "U" in override_cycles.upper():
    settings.append(["TrimUMI", "0"])


# Remove old versions of settings we control.
remove_keys = {
    "OverrideCycles",
    "NoLaneSplitting",
    "TrimUMI",
    "CreateFastqIndexForReads",
    "ReverseComplement",
}

cleaned = []

for row in rows:
    if row and row[0].strip() in remove_keys:
        continue

    cleaned.append(row)

rows = cleaned


# Find Settings section
settings_idx = None

for n, row in enumerate(rows):
    if row and row[0].strip() in {
        "[Settings]",
        "[BCLConvert_Settings]",
    }:
        settings_idx = n
        break

if settings_idx is None:
    raise RuntimeError(
        "Could not locate [Settings] or [BCLConvert_Settings] section"
    )


for offset, setting in enumerate(settings, start=1):
    rows.insert(settings_idx + offset, setting)


with output.open("w", newline="") as fh:
    writer = csv.writer(fh, lineterminator="\\n")
    writer.writerows(rows)

print(f"Wrote corrected samplesheet: {output}")

PY
    """
}

process prepare_RNA_samplesheet {

    input:
    tuple val(samplesheet_basename), path(samplesheet)// from original_samplesheet2
    output:
    path("*.RNA_SAMPLES.csv"), emit: std
    path("*.UMI.csv"),emit: umi// into rnaSS1    

    script:
    """

    cat ${samplesheet} | grep "RV1" > ${samplesheet_basename}.RNAsamples.intermediate.txt
    sed -n '1,/Sample_ID/p' ${samplesheet} > ${samplesheet_basename}.HEADER.txt
    cat ${samplesheet_basename}.HEADER.txt ${samplesheet_basename}.RNAsamples.intermediate.txt > ${samplesheet_basename}.RNA_SAMPLES.csv 

    if [[ "${umiConvertRNA}" == *U* ]]; then
        sed 's/Settings]/&\\nOverrideCycles,${umiConvertRNA}\\nNoLaneSplitting,true\\nTrimUMI,0/' ${samplesheet_basename}.RNA_SAMPLES.csv > ${samplesheet_basename}.RNA_SAMPLES.UMI.csv
    else
        sed 's/Settings]/&\\nOverrideCycles,${umiConvertRNA}\\nNoLaneSplitting,true/' ${samplesheet_basename}.RNA_SAMPLES.csv > ${samplesheet_basename}.RNA_SAMPLES.UMI.csv
    fi
    """
}

process bclConvert_DNA {
    tag "$runfolder_simplename"
    errorStrategy 'ignore'

    publishDir "${fastq_dir}/${runfolder_simplename}_umi", mode: 'copy', pattern:"*.fastq.gz"
    publishDir "${qc_dir}/", mode: 'copy', pattern:"*.DNA.html"

    input:
    tuple val(runfolder_simplename), path(runfolder)// from runfolder_ch2
    path(dnaSS) // from dnaSS1
    path(runinfo) // from xml_ch

    output:
    path("*.fastq.gz"), emit: dna_fastq// into (dna_fq_out,dna_fq_out2)
    path("*.Multiqc.DNA.html")
    script:
    """
    bcl-convert \
    --sample-sheet ${dnaSS} \
    --bcl-input-directory ${runfolder} \
    --output-directory ${runfolder_simplename}_umi/
    
    singularity run -B ${s_bind} ${simgpath}/multiqc.sif \
    -c ${multiqc_config} \
    -f -q ${runfolder_simplename}_umi/ \
    -n ${runfolder_simplename}.DemultiplexRunStats.Multiqc.DNA.html
    
    mv ${runfolder_simplename}_umi/*.fastq.gz .
    rm -rf Undetermined*
    """
}

process bclConvert_RNA {
    tag "$runfolder_simplename"
    errorStrategy 'ignore'

    publishDir "${fastq_dir}/${runfolder_simplename}_umi/", mode: 'copy', pattern:"*.fastq.gz"
    publishDir "${qc_dir}/", mode: 'copy', pattern:"*.RNA.html"

    input:
    tuple val(runfolder_simplename), path(runfolder)// from runfolder_ch2
    path(rnaSS) // from dnaSS1
    path(runinfo) // from xml_ch

    output:

    path("*.fastq.gz"), emit: rna_fastq// into (dna_fq_out,dna_fq_out2)
    path("*.Multiqc.RNA.html")

    script:
    """
    bcl-convert \
    --sample-sheet ${rnaSS} \
    --bcl-input-directory ${runfolder} \
    --output-directory ${runfolder_simplename}_umi/
    
    singularity run -B ${s_bind} ${simgpath}/multiqc.sif \
    -c ${multiqc_config} \
    -f -q ${runfolder_simplename}_umi/ \
    -n ${runfolder_simplename}.DemultiplexRunStats.Multiqc.RNA.html

    mv ${runfolder_simplename}_umi/*.fastq.gz .
    rm -rf Undetermined*
    """
}


///////////////////////////////// PREPROCESS MODULES - STANDARD //////////////////////// 


process fastq_to_ubam {
    errorStrategy 'ignore'
    tag "$meta.id"
    //publishDir "${params.outdir}/unmappedBAM/", mode: 'copy',pattern: '*.{bam,bai}'
    //publishDir "${params.outdir}/${runfolder_basename}/fastq_symlinks/", mode: 'link', pattern:'*.{fastq,fq}.gz'
    cpus 2
    maxForks 30

    input:
    tuple val(meta), path(data) // Meta [npn,id(npn_superpanel),superpanel,panel,subpanel,runfolder], data: [r1,r2]

    output:
    tuple val(meta), path("${meta.id}.unmapped.from.fq.bam"),emit:ubam// into (ubam_out1, ubam_out2)
    
    script:
    """
    ${gatk_exec} FastqToSam \
    -F1 ${data[0]} \
    -F2 ${data[1]} \
    -SM ${meta.id} \
    -PL illumina \
    -PU KGA_PU \
    -RG KGA_RG \
    -O ${meta.id}.unmapped.from.fq.bam
    """
}

process markAdapters {
    tag "$meta.id"
    errorStrategy 'ignore'

    input:
    tuple val(meta), path(uBAM)

    output:
    tuple val(meta), path("${meta.id}.ubamXT.bam"), path("${meta.id}.markAdapterMetrics.txt")// into (ubamXT_out,ubamXT_out2)

    script:
    """
    ${gatk_exec} MarkIlluminaAdapters \
    -I ${uBAM} \
    -O ${meta.id}.ubamXT.bam \
    --TMP_DIR ${tmpDIR} \
    -M ${meta.id}.markAdapterMetrics.txt
    """
}


process align {
    tag "$meta.id"

    maxForks 10
    errorStrategy 'ignore'
    cpus 60

    input:
    tuple val(meta), path(uBAMXT), path(metrics)//  from ubamXT_out

    output:
    tuple val(meta), path("${meta.id}.${genome_version}.QNsort.BWA.clean.bam")
    
    script:
    """
    ${gatk_exec} SamToFastq \
    -I ${uBAMXT} \
    -INTER \
    -CLIP_ATTR XT \
    -CLIP_ACT 2 \
    -NON_PF \
    -F /dev/stdout \
    |  singularity run -B ${s_bind} ${simgpath}/bwa0717.sif bwa mem \
    -t ${task.cpus} \
    -p \
    ${genome_fasta} \
    /dev/stdin \
    | ${gatk_exec} MergeBamAlignment \
    -R ${genome_fasta} \
    -UNMAPPED ${uBAMXT} \
    -ALIGNED /dev/stdin \
    -MAX_GAPS -1 \
    -ORIENTATIONS FR \
    -SO queryname \
    -O ${meta.id}.${genome_version}.QNsort.BWA.clean.bam
    """
}


process markDup_cram {
    errorStrategy 'ignore'
    maxForks 16
    tag "$meta.id"
    publishDir "${aln_output_dir}/${meta.runfolder}_umi/", mode: 'copy', pattern: "*.cra*"
    conda '/lnx01_data3/shared/programmer/miniconda3/envs/samblasterSambamba/'
    input:
    tuple val(meta), path(bam)
    
    output:
    tuple val(meta), path("${meta.id}.${genome_version}.cram"), path("${meta.id}.${genome_version}*crai")
    
    script:
    """
    samtools view -h ${bam} \
    | samblaster | sambamba view -t 8 -S -f bam /dev/stdin | sambamba sort -t 8 --tmpdir=${tmpDIR} -o /dev/stdout /dev/stdin \
    |  samtools view \
    -T ${genome_fasta} \
    -C \
    -o ${meta.id}.${genome_version}.cram -

    samtools index ${meta.id}.${genome_version}.cram
    """
}


///////////////////////////////// PREPROCESS MODULES - UMI samples //////////////////////// 

process fastq_to_ubam_umi {
    errorStrategy 'ignore'
    tag "$sampleID"
    //publishDir "${params.outdir}/unmappedBAM/", mode: 'copy',pattern: '*.{bam,bai}'
    //publishDir "${params.outdir}/${runfolder_basename}/fastq_symlinks/", mode: 'link', pattern:'*.{fastq,fq}.gz'
    cpus 2
    maxForks 30
    conda '/lnx01_data3/shared/programmer/miniconda3/envs/fgbio/'

    input:
    tuple val(meta), path(data)
    
    output:
    tuple val(meta), path("${meta.id}.unmapped.umi.from.fq.bam"),   emit: unmappedBam
    tuple val(meta), path(data),                                    emit: fastq       

    
    script:
    """
    fgbio FastqToBam \
    -i ${data} \
    -n \
    --sample ${meta.id} \
    --library ${meta.npn} \
    --read-group-id KGA_RG \
    -O ${meta.id}.unmapped.umi.from.fq.bam
    """
}



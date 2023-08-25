# SEGLH ExomeDepth pipeline

## input files

Exomdepth runs on a collection of intervals that must be provided in BED format.
The target intervals must be annotated with gene symbol and exon number as follows:

```
3	128199861	128200161	GATA2_5
3	128200661	128200787	GATA2_4
3	128202702	128202848	GATA2_3
3	128204569	128205211	GATA2_2
3	128205645	128205874	GATA2_1
```
If intronic or intergenic regions are included they can be named accordingly using the same naming principle (eg. GATA2-INTRON_1). Exomedepth uses the exon/intron numbers only to label the CNVs in the report output and they are not part of the applied mathematical model.

## Output files

While only single output files are to be defined in the read count and exomedepth step, both scripts produce a set of output files. Be aware that only the PDF report masks CNVs which are outside of the BED regions supplied to the ExomeDepth step. ``VCF`` and ``BED`` files contain the full output and should be filtered accordingly to mask incidental findings.

### Read count

- `readCount.RData` - Read count data and selected references per sample
- `readCount.csv` - Model parameters and QC metrics output (can be used to build a QC classifier, see below)

### Exomedepth

- `output.pdf` - Exomedepth CNV report with all QC information
- `output.bed` - CNVs in BED format (whole panel)
- `output.vcf` - CNVs in VCF format (whole panel, with out-of-scope filter tags, see VCF file header)

## Usage

Build the docker container which contains all scripts with `make`.
The built docker image is tagged as `seglh/exomedepth:latest`.

The docker image contains all required packages and scripts but no reference genome.

### Step 1 - Create read count data

```
docker run -it \
	-v /path_to_data:/data \
	-v /path_to_genome:/resources \
	seglh/exomedepth:latest \
	readCount.R \
	/data/readCount.RData \
	/genome/human_g1k_v37_decoy.fasta \
	/data/targets.bed \
	/data/sample1.bam \
	/data/sample2.bam \
	/data/sample3.bam \
	/data/sample4.bam
```



### Step 2 - Run Exomedepth

This steps will pick the most appropriate reference samples and run exomedepth for a given sample.
Reference samples are picked from within the same batch. If any BAM file is prefixed with _NORMAL_ it will be considered part of a panel of normals which will be used as reference instead of other batched samples (requires at least 2 normal samples).

```
docker run -it \
	-v /path_to_data:/data \
	-v /path_to_genome:/resources \
	seglh/exomedepth:latest \
	exomeDepth.R \
	VERSIONSTRING \
	/data/output.pdf \
	/data/targets.bed:PANELNAME \
	/data/readCount.RData \
	/data/sample1.bam:SampleLabel \
	/data/qcmodel.RData
```
NB: The last argument is optional (see below).

The script will reuse the same genome reference file for VCF generation. Make sure its location is available (mounted into the container). Additional output files are produced as described above.

## Reference sample modes

Exomedepth, as default, runs with intra-batch normalisation (normals are picked from batch). Alternatively a panel of normals (PoN) can be used. This is highly recommended if the batch sizes are small or known negatives are available to ensure sensitivity to clinically relevant variants.

### Intra-batch normalisation

### Panel of Normals (PoN)

To run exomedepth in PoN mode, supply normal BAMs to the readCount step as any other BAM file, but ensure their name is prefixed with _NORMAL_. This will automatically switch to PoN mode. The used normalisation mode is also indicated on the PDF reports, alongside the picked normal samples.

Any normals that have low coverage (below 100X) in _any_ of the supplied target regions will automatically be dropped as a safety measure.

If a Panel of Normals is supplied, the `readCount.R` script will pick reference sets for both, the PoN and the supplied batch, under the condition that there are at least 3 reference samples to choose from.

To avoid recounting reads for the normal samples, a `readCount.RData` file from a previous run can be supplied and the normals contained in this run will be merged into the counts table. Normal BAMs supplied will not be recounted if they are in the RData file.

```
docker run -it \
	-v /path_to_data:/data \
	-v /path_to_genome:/resources \
	seglh/exomedepth:latest \
	readCount.R \
	/data/readCount.RData \
	/genome/human_g1k_v37_decoy.fasta \
	/data/targets.bed \
	/data/sample1.bam \
	/data/sample2.bam \
	/data/NORMAL1.bam \ <- will not be recounted if in RData file below
	/data/previous_readCount.RData  <- output from a previous run containing the PoN
```

Ideally a panel of normals is created as a standalone file without any test samples included, and then supplied to the ``readCount`` step as an `.RData` file as follows:

```
docker run -it \
	-v /path_to_data:/data \
	-v /path_to_genome:/resources \
	seglh/exomedepth:latest \
	readCount.R \
	/data/normals.RData \
	/genome/human_g1k_v37_decoy.fasta \
	/data/targets.bed \
	/data/NORMAL1.bam \
	/data/NORMAL2.bam \
	/data/NORMAL3.bam
```


```
docker run -it \
	-v /path_to_data:/data \
	-v /path_to_genome:/resources \
	seglh/exomedepth:latest \
	readCount.R \
	/data/readCount.RData \
	/genome/human_g1k_v37_decoy.fasta \
	/data/targets.bed \
	/data/sample1.bam \
	/data/sample2.bam \
	/data/sample3.bam \
	/data/normals.RData
```

#### Minimal reference count for PoN
If the maximum Bayes Factor (BF) of the selected reference set is for a single reference sample, the next higher BF will be selected. This forces the reference to be composed of at least 2 reference samples and sacrifices some statistical power for enhanced specificity (lower FNR and FPR).
If this rule is applied it will be highlighted on the QC report.

#### Implicit sex matching for PoN

If a PoN is used and target regions are located on sex chromosomes, the readcount step will automatically try to match test samples with normals of the same sex.
Sex is extracted from input file names (BAM) using the regex `_[MF]_`.
No sex matching will occur of not all samples in the reference (PoN) and the test sample do not match this regular expression.

## Integrated Quality Control

The report contains pertinent information for Quality control. Two plots visualise the RPKM correlation with other samples in the same batch, and the coefficient of variation.

A table summarises if QC thresholds have been met (PASS), failed a soft threshold (CAUTION) or a hard threshold (FAIL). If a threshold is not met, an interpretation help is provided in the report.

Optionally a random forest classifier can be provided which will indicate the PASS/FAIL classification in the report. This classification is only a recommendation and the scientist has ultimate authority to classify the sample as failed.
### QC Thresholds

Exomdedepth implements default QC thresholds as follows:

| Metric                              | Soft Threshold | Hard Threshold |
| :---------------------------------- | -------------: | -------------: |
| Median Correlation in Batch         |             NA |           0.90 |
| Maximum Correlation in Batch        |           0.95 |           0.90 |
| Coefficent of Variation             |             30 |             35 |
| Correlation with selected reference |           0.95 |           0.90 |
| Selected reference samples          |              3 |              1 |

These can be overridden by providing an RData file containing thresholds as follows:

```
limits<-list(
  medcor=c(NA, 0.90),    # median correlation within batch
  maxcor=c(0.95, 0.90),  # max correlation within batch
  refcor=c(0.95, 0.90),  # reference set correlation
  refcount=c(5, 3),      # refernce set size (selected reference samples)
  coeffvar=c(30, 35),    # coefficient of variation
  coverage=c(100),       # Coverage warning threshold
  expectedbf=c(5.0, NA), # expected BF of reference choice (power)
  minrefs=c(2,Inf)       # minimum reference set size
)
```

The `configure.R` script contains an example how these can be specified.

### Sex check on reference set

The report will check if test sample and chosen reference samples are sex-matched using the regex `_[MF]_` on the BAM file names.

### Additional annotations for CNVs

Additional annotations of CNVs can be added by generating a configuration file. Any BED4 file supplied to the `configure.R` script will be used to annotation any found and reported CNVs (PDF only).
CNVs will be annotated if they overlap at least 10 percent with any additional annotation supplied.

### Random Forest QC classifier

SEGLH-Exomedepth can optionally provide a QC classification based on previous labeled data. The classfier is automatically optimised for the number of random variables considered and runs 200 trees by default. All labeled samples are part of the training set.

The model can be built as follows:

1. Amend a column with header `status` to the CSV metrics output from the `readCount.R` step.
2. Label the samples with a known classification (eg. FAIL or PASS) in this column. Providing a label is optional. Missing values will not be included in the model generation.
3. Build the model with `Rscript configure.R qcmodel.RData file1.csv [file2.csv ...]`
4. Supply the `qcmodel.RData` file to the exomedepth command line as last parameter. The report will automaticcly be amended with the RF QC classification and an explanation.

It is recommended to inspect the generated Random Forest classifier in terms of error rate convergence and out-of-bag error. Also consider adjusting of the number of trees!

## Bootstrapping a set of normal samples

Often the complete negative status of a potential normal sample is not known. the `bootstrapCnv.R` script allows to bootstrap CNVs from a pool of likely-normals which then can be used to exclude samples which do not meet user set criteria for being a normal sample (eg presence of no CNVs in a certain region).

### Process

The following steps are recommended to further reduce a set of likely-normals to a more stringent set of normals, by excluding samples that do have reproducible CNVs in a within-batch normalisation approach.

The `utils/getNormals.sh` script implements an example of the steps below.

#### Read counting for pool of normals

Run the readcount step as normal supplying the PoN candidates.

#### Bootstrap CNVs

Run the resampling script. As default, 30 samples are picked from the pool 100 times, for all samples in sequence. For large sets of normals, the runs can be split into multiple parallel runs by specifying start and end samples.

`Rscript bootstrapCnv.R readcounts.RData [BOOTSTRAP:100] [SAMPLESIZE:30] ([FROM:1] [TO:samplecount])`

A BED file with known highly variable regions (common CNVs) can be supplied as the last argument.

The output files are named according to the input file and time stamped:

- `readcounts.210912T130216.scores` - The CNVs from all runs (aggregated)
- `readcounts.210912T130216.RData` - Saved state 

#### Filtering

The user can either manually assess the returned bootstrapped CNV calls or rely on the integrated filtering.
The filtering aims to identify samples that have highly reproducible CNV calls and therefore do not meet the normal expectation.

To run the filtering, supply a previous _scores_ file to the command line:

`Rscript bootstrapCnv.R readcounts.210912T130216.RData 100 30 readcounts.210912T130216.scores`

NB: If previous runs have been partitioned, ensure the scores files are concatenated before initiating the filtering step.

The following filter steps are applied in order:

1. Remove any common CNVs (according to supplied BED file)
2. Remove CNVs with high entropy (>2) and prevalence (>0.4) in the bootstrap (highly variable CNVs)
3. Remove CNVs with less than 0.5 bootstrap frequency (unstable CNVs or noise).

The excluded samples are printed on screen.

#### Re-test normals (iterate)
If the excluded samples are considered appropriate, one can rerun the bootrapping processing while excluding these samples by appending `RUN` as an additional argument.

`Rscript bootstrapCnv.R readcounts.210912T130216.RData 100 30 readcounts.210912T130216.scores RUN`

This process can be iterated multiple times. However, this could lead to overfitting the normal set.


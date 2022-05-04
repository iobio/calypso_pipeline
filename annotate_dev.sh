#! /bin/bash
set -eou pipefail

# Inputs are a compressed, indexed vcf file and a ped file
VCF=$1
PED=$2

# Generate names of output files
VCFBASE=`echo $VCF | rev | cut -d '/' -f 1 | rev`
NORMVCF=${VCFBASE/vcf.gz/norm.vcf.gz}
CLEANVCF=${VCFBASE/vcf.gz/clean.vcf.gz}
ANNOTATEDVCF=${VCFBASE/vcf.gz/annotated.vcf.gz}
SLIVAR1=${VCFBASE/vcf.gz/slivar1.vcf.gz}
SLIVAR2=${VCFBASE/vcf.gz/slivar2.vcf.gz}
FINALVCF=${VCFBASE/vcf.gz/final.vcf.gz}

# Get the path where this script resides
ROOT_PATH=`echo ${BASH_SOURCE[0]} | rev | cut -d '/' -f 2- | rev`

# Required files
REF=$ROOT_PATH/data/GRCh38/human_g1k_v38_decoy_phix.fasta
GFF=$ROOT_PATH/data/Homo_sapiens.GRCh38.105.chr.gff3.gz
# GFF must be reprocessed from Ensembl site file to add chr prefixes to file, code below would likely do it
# zgrep ^# Homo_sapiens.GRCh38.105.gff3.gz > Homo_sapiens.GRCh38.105.chr.gff3 && zgrep -v ^# Homo_sapiens.GRCh38.105.gff3.gz | sed 's/^/chr/g' >> Homo_sapiens.GRCh38.105.chr.gff3
# bgzip Homo_sapiens.GRCh38.105.chr.gff3
GNOMAD=$ROOT_PATH/data/gnomad.genomes.v3.1.sites.slivar.zip
JS=$ROOT_PATH/data/slivar-functions.js
# slivar_static must be 0.2.4 or greater

# Create the toml file for vcfanno
TOML=$ROOT_PATH/data/annotations.toml
# Add REVEL annotations to TOML
echo "[[annotation]]" > $TOML
echo "file=\"$ROOT_PATH/data/revel_scores_grch38.vcf.gz\"" >> $TOML
echo "fields=[\"REVEL\"]" >> $TOML
echo "ops=[\"self\"]" >> $TOML
echo "" >> $TOML
# Add ClinVar annotations to TOML
echo "[[annotation]]" >> $TOML
echo "file=\"$ROOT_PATH/data/ClinVar_GRCh38_Mar8_2022.vcf.gz\"" >> $TOML
echo "fields=[\"CLNDN\", \"CLNSIG\"]" >> $TOML
echo "ops=[\"self\", \"self\"]" >> $TOML
echo "" >> $TOML
# Add pLI annotations to TOML
echo "[[annotation]]" >> $TOML
echo "file=\"$ROOT_PATH/data/gnomad.v2.1.1.lof_metrics.by_gene.txt.gz.pli.grch38.sort.bed.gz\"" >> $TOML
echo "columns=[4]" >> $TOML
echo "names=[\"pLI\"]" >> $TOML
echo "ops=[\"self\"]" >> $TOML
echo "" >> $TOML
# Add CCR annotations to TOML
echo "[[annotation]]" >> $TOML
echo "file=\"$ROOT_PATH/data/ccrs.autosomes.v2.20180420.grch38.sort.bed.gz\"" >> $TOML
echo "columns=[4]" >> $TOML
echo "names=[\"CCR\"]" >> $TOML
echo "ops=[\"self\"]" >> $TOML
echo "" >> $TOML
# Add EVE annotations to TOML
echo "[[annotation]]" >> $TOML
echo "file=\"$ROOT_PATH/data/EVE.all.variants.chr.vcf.gz\"" >> $TOML
echo "fields=[\"EVE\"]" >> $TOML
echo "ops=[\"self\"]" >> $TOML
echo "" >> $TOML
# Add MutScore annotations to TOML
echo "[[annotation]]" >> $TOML
echo "file=\"$ROOT_PATH/data/mutscore-v1.0-hg38.vcf.gz\"" >> $TOML
echo "fields=[\"MutScore\"]" >> $TOML
echo "ops=[\"self\"]" >> $TOML
echo "" >> $TOML

# normalize and subset original VCF
echo "Starting annotation..."
echo "Subsetting and normalizing input VCF..."

tail -n+2 $PED | cut -f 2 | sort -u > samples.txt
bcftools norm -m - -w 10000 -f $REF $VCF \
	| bcftools view -a -c 1 -S samples.txt -O z -o $NORMVCF

bcftools index -t $NORMVCF

# strip existing VCF annotations
echo "Removing any existing annotations from input VCF..."

bcftools annotate -x INFO $NORMVCF -O z -o $CLEANVCF

bcftools index -t $CLEANVCF

# annotate with bcftools csq and add ClinVar and REVEL annotations with vcfanno
echo "Annotating cleaned VCF..."

bcftools csq -f $REF \
	--ncsq 40 \
	-l \
	-g $GFF \
	$CLEANVCF \
	| vcfanno -p 16 $TOML /dev/stdin \
	| bcftools view -O z -o $ANNOTATEDVCF -
#	| bgzip -c > $ANNOTATEDVCF

slivar_static expr \
	--vcf $ANNOTATEDVCF \
	--ped $PED \
	--js $JS \
	-g $GNOMAD \
	--info 'INFO.gnomad_popmax_af < 0.01 && variant.FILTER == "PASS" && variant.ALT[0] != "*"' \
	--family-expr 'denovo:fam.every(segregating_denovo) && INFO.gnomad_popmax_af < 0.001' \
	--family-expr 'x_denovo:(variant.CHROM == "X" || variant.CHROM == "chrX") && fam.every(segregating_denovo_x) && INFO.gnomad_popmax_af < 0.001' \
	--family-expr 'recessive:fam.every(segregating_recessive)' \
	--family-expr 'dominant:fam.every(segregating_dominant)' \
	-o $SLIVAR1

slivar_static expr \
	--pass-only \
	--vcf $ANNOTATEDVCF \
	--ped $PED \
	--js $JS \
	-g $GNOMAD \
	--family-expr 'denovo:fam.every(segregating_denovo) && INFO.gnomad_popmax_af < 0.001' \
	--trio 'comphet_side:comphet_side(kid, mom, dad) && INFO.gnomad_popmax_af < 0.005' \
	| slivar_static compound-hets \
		--vcf /dev/stdin \
		--skip NONE \
		-s comphet_side \
		-s denovo \
		-p $PED \
		-o $SLIVAR2

bcftools index -t $SLIVAR1
bcftools index -t $SLIVAR2

echo "Combining VCFs and cleaning things up..."

bcftools concat -a -d none $SLIVAR1 $SLIVAR2 -O z -o $FINALVCF

bcftools index -t $FINALVCF

rm $SLIVAR1* $SLIVAR2* $ANNOTATEDVCF* $CLEANVCF* $NORMVCF*

echo "Everything completed! Annotated VCF written to $FINALVCF"

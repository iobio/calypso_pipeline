import os

# Use Slivar to filter the annotated vcf file based on a trio family structure
def annotateTrio(resourceInfo, bashFile):

  # Generate the list of annotations to be extracted from the CSQ INFO field for use in slivar
  annotationString = ''
  for annotation in resourceInfo['resources']['vep']['slivar']: annotationString += annotation + ','
  annotationString = annotationString.rstrip(',')

  # Print out status messages
  print('# Perform filtering using Slivar. This is a trio', file = bashFile)
  print('echo -n "Filtering trio variants..."', file = bashFile)

  # The slivar filter needs some of the VEP annotations to appear as standalone annotations in the INFO field. These are defined in the
  # resource file. bcftools +split-vep is used to extract the required annotations from the CSQ field, and the resulting vcf is passed to
  # slivar
  print('$BCFTOOLS +split-vep -O z -a CSQ \\', file = bashFile)
  print('  -c ', annotationString, ' \\', sep = '', file = bashFile)
  print('  $ANNOTATEDVCF \\', file = bashFile)
  print('  | $SLIVAR expr \\', file = bashFile)
  print('  --vcf - \\', file = bashFile)
  print('  --ped $PED \\', file = bashFile)
  print('  --pass-only \\', file = bashFile)
#  print('  -g $SLIVAR_GNOMAD_V2 \\', file = bashFile)
#  print('  -g $SLIVAR_GNOMAD_V3 \\', file = bashFile)
  print('  --js $JS \\', file = bashFile)
  print('  -o $SLIVAR1VCF \\', file = bashFile)
#  print('  --info "(INFO.impactful || (\'SpliceAI_pred_DS_AG\' in INFO && INFO.SpliceAI_pred_DS_AG >= 0.1) || (\'SpliceAI_pred_DS_AL\' in INFO && INFO.SpliceAI_pred_DS_AL >= 0.1) || \\', file = bashFile)
#  print('    (\'SpliceAI_pred_DS_DG\' in INFO && INFO.SpliceAI_pred_DS_DG >= 0.1) || (\'SpliceAI_pred_DS_DL\' in INFO && INFO.SpliceAI_pred_DS_DL >= 0.1) || \\', file = bashFile)
#  print('    ((\'MaxEntScan_diff\' in INFO && INFO.MaxEntScan_diff >= 1.15) && (\'MaxEntScan_alt\' in INFO && INFO.MaxEntScan_alt < 6.2)) || \\', file = bashFile)
#  print('    ((\'MaxEntScan_diff\' in INFO && INFO.MaxEntScan_diff < 0) && (\'MaxEntScan_alt\' in INFO && INFO.MaxEntScan_alt > 8.5)) || \\', file = bashFile)
#  print('    (\'SpliceRegion\' in INFO && INFO.SpliceRegion.includes(\'splice\')) || \\', file = bashFile)
#  print('    (\'ada_score\' in INFO && INFO.ada_score >= 0.6) || (\'rf_score\' in INFO && INFO.rf_score >= 0.6) || \\', file = bashFile)
#  print('    (\'GeneSplicer\' in INFO && INFO.GeneSplicer.includes(\'/High\')) || \\', file = bashFile)
#  print('    (\'five_prime_UTR_variant_consequence\' in INFO && INFO.five_prime_UTR_variant_consequence.includes(\'uAUG\')) || \\', file = bashFile)
#  print('    (\'five_prime_UTR_variant_consequence\' in INFO && INFO.five_prime_UTR_variant_consequence.includes(\'uSTOP\')) || \\', file = bashFile)
#  print('    (\'five_prime_UTR_variant_consequence\' in INFO && INFO.five_prime_UTR_variant_consequence.includes(\'uFrame\'))) && \\', file = bashFile)
#  print('    INFO.gg2_1_1_controls_AF_popmax < 0.02 && INFO.gg2_1_1_AF_popmax < 0.02 && INFO.gg3_1_2_AF_controls < 0.02 && INFO.gg3_1_2_AF < 0.02 && \\', file = bashFile)
#  print('    INFO.gg3_1_2_AF_nfe < 0.02 && INFO.gg3_1_2_AF_amr < 0.02 && INFO.gg3_1_2_AF_afr < 0.02 && INFO.gg3_1_2_AF_eas < 0.02 && INFO.gg3_1_2_AF_sas < 0.02" \\', file = bashFile)
#  print('  --trio "denovo:denovo(kid, mom, dad) && INFO.gg2_1_1_controls_AF_popmax < 0.001 && INFO.gg2_1_1_AF_popmax < 0.001 && INFO.gg3_1_2_AF_controls < 0.001 && INFO.gg3_1_2_AF < 0.001 && \\', file = bashFile)
#  print('    INFO.gg3_1_2_AF_nfe < 0.001 && INFO.gg3_1_2_AF_amr < 0.001 && INFO.gg3_1_2_AF_afr < 0.001 && INFO.gg3_1_2_AF_eas < 0.001 && INFO.gg3_1_2_AF_sas < 0.001 && \\', file = bashFile)
#  print('    (kid.sex != \'male\' || variant.CHROM != \'chrX\' || variant.POS <= 2781479 || variant.POS >= 155701383)" \\', file = bashFile)
#  print('  --trio "x_denovo:x_denovo(kid, mom, dad) && kid.hom_alt && INFO.gg2_1_1_controls_AF_popmax < 0.001 && INFO.gg2_1_1_AF_popmax < 0.001 && INFO.gg3_1_2_AF_controls < 0.001 && INFO.gg3_1_2_AF < 0.001 && \\', file = bashFile)
#  print('    INFO.gg3_1_2_AF_nfe < 0.001 && INFO.gg3_1_2_AF_amr < 0.001 && INFO.gg3_1_2_AF_afr < 0.001 && INFO.gg3_1_2_AF_eas < 0.001 && INFO.gg3_1_2_AF_sas < 0.001 && \\', file = bashFile)
#  print('    kid.sex == \'male\' && variant.CHROM == \'chrX\' && variant.POS > 2781479 && variant.POS < 155701383" \\', file = bashFile)
#  print('  --trio "hom:hom(kid, mom, dad) && INFO.gg2_1_1_controls_nhomalt < 2 && INFO.gg2_1_1_nhomalt < 3 && INFO.gg3_1_2_nhomalt_controls < 2 && INFO.gg3_1_2_nhomalt < 3 && \\', file = bashFile)
#  print('    (kid.sex != \'male\' || variant.CHROM != \'chrX\' || variant.POS <= 2781479 || variant.POS >= 155701383)" \\', file = bashFile)
#  print('  --trio "x_hemi:x_hemi(kid, mom, dad) && kid.hom_alt && INFO.gg2_1_1_controls_AF_popmax < 0.001 && INFO.gg2_1_1_AF_popmax < 0.001 && INFO.gg3_1_2_AF_controls < 0.001 && INFO.gg3_1_2_AF < 0.001 && \\', file = bashFile)
#  print('    INFO.gg3_1_2_AF_nfe < 0.001 && INFO.gg3_1_2_AF_amr < 0.001 && INFO.gg3_1_2_AF_afr < 0.001 && INFO.gg3_1_2_AF_eas < 0.001 && INFO.gg3_1_2_AF_sas < 0.001 && \\', file = bashFile)
#  print('    INFO.gg2_1_1_AC_male_controls < 2 && INFO.gg2_1_1_AC_male < 3 && INFO.gg3_1_2_AC_controls_XY < 2 && INFO.gg3_1_2_AC_XY < 3 && \\', file = bashFile)
#  print('    INFO.gg2_1_1_controls_nhomalt < 2 && INFO.gg2_1_1_nhomalt < 3 && INFO.gg3_1_2_nhomalt_controls < 2 && INFO.gg3_1_2_nhomalt < 3 && \\', file = bashFile)
#  print('    kid.sex == \'male\' && variant.CHROM == \'chrX\' && variant.POS > 2781479 && variant.POS < 155701383" \\', file = bashFile)
#  print('  --trio "comphet_side:comphet_side(kid, mom, dad) && INFO.gg2_1_1_controls_nhomalt < 2 && INFO.gg2_1_1_nhomalt < 3 && INFO.gg3_1_2_nhomalt_controls < 2 && INFO.gg3_1_2_nhomalt < 3" \\', file = bashFile)
  print('  --info \'(INFO.impactful || ("SpliceAI_pred_DS_AG" in INFO && INFO.SpliceAI_pred_DS_AG >= 0.1) || ("SpliceAI_pred_DS_AL" in INFO && INFO.SpliceAI_pred_DS_AL >= 0.1) ||', file = bashFile)
  print('    ("SpliceAI_pred_DS_DG" in INFO && INFO.SpliceAI_pred_DS_DG >= 0.1) || ("SpliceAI_pred_DS_DL" in INFO && INFO.SpliceAI_pred_DS_DL >= 0.1) ||', file = bashFile)
  print('    (("MaxEntScan_diff" in INFO && INFO.MaxEntScan_diff >= 1.15) && ("MaxEntScan_alt" in INFO && INFO.MaxEntScan_alt < 6.2)) ||', file = bashFile)
  print('    (("MaxEntScan_diff" in INFO && INFO.MaxEntScan_diff < 0) && ("MaxEntScan_alt" in INFO && INFO.MaxEntScan_alt > 8.5)) ||', file = bashFile)
  print('    ("SpliceRegion" in INFO && INFO.SpliceRegion.includes("splice")) ||', file = bashFile)
  print('    ("ada_score" in INFO && INFO.ada_score >= 0.6) || ("rf_score" in INFO && INFO.rf_score >= 0.6) ||', file = bashFile)
  print('    ("GeneSplicer" in INFO && INFO.GeneSplicer.includes("High")) ||', file = bashFile)
  print('    ("five_prime_UTR_variant_consequence" in INFO && INFO.five_prime_UTR_variant_consequence.includes("uAUG")) ||', file = bashFile)
  print('    ("five_prime_UTR_variant_consequence" in INFO && INFO.five_prime_UTR_variant_consequence.includes("uSTOP")) ||', file = bashFile)
  print('    ("five_prime_UTR_variant_consequence" in INFO && INFO.five_prime_UTR_variant_consequence.includes("uFrame"))) &&', file = bashFile)
  print('    (!("gg2_1_1_controls_AF_popmax" in INFO) || ("gg2_1_1_controls_AF_popmax" in INFO && INFO.gg2_1_1_controls_AF_popmax < 0.02)) &&', file = bashFile)
  print('    (!("gg2_1_1_AF_popmax" in INFO) || ("gg2_1_1_AF_popmax" in INFO && INFO.gg2_1_1_AF_popmax < 0.02)) &&', file = bashFile)
  print('    (!("gg3_1_2_AF_controls" in INFO) || ("gg3_1_2_AF_controls" in INFO && INFO.gg3_1_2_AF_controls < 0.02)) &&', file = bashFile)
  print('    (!("gg3_1_2_AF" in INFO) || ("gg3_1_2_AF" in INFO && INFO.gg3_1_2_AF < 0.02)) &&', file = bashFile)
  print('    (!("gg3_1_2_AF_nfe" in INFO) || ("gg3_1_2_AF_nfe" in INFO && INFO.gg3_1_2_AF_nfe < 0.02)) &&', file = bashFile)
  print('    (!("gg3_1_2_AF_amr" in INFO) || ("gg3_1_2_AF_amr" in INFO && INFO.gg3_1_2_AF_amr < 0.02)) &&', file = bashFile)
  print('    (!("gg3_1_2_AF_afr" in INFO) || ("gg3_1_2_AF_afr" in INFO && INFO.gg3_1_2_AF_afr < 0.02)) &&', file = bashFile)
  print('    (!("gg3_1_2_AF_eas" in INFO) || ("gg3_1_2_AF_eas" in INFO && INFO.gg3_1_2_AF_eas < 0.02)) &&', file = bashFile)
  print('    (!("gg3_1_2_AF_sas" in INFO) || ("gg3_1_2_AF_sas" in INFO && INFO.gg3_1_2_AF_sas < 0.02))\' \\', file = bashFile)
  print('  --trio \'denovo:denovo(kid, mom, dad) &&', file = bashFile)
  print('    (!("gg2_1_1_controls_AF_popmax" in INFO) || ("gg2_1_1_controls_AF_popmax" in INFO && INFO.gg2_1_1_controls_AF_popmax < 0.001)) &&', file = bashFile)
  print('    (!("gg2_1_1_AF_popmax" in INFO) || ("gg2_1_1_AF_popmax" in INFO && INFO.gg2_1_1_AF_popmax < 0.001)) &&', file = bashFile)
  print('    (!("gg3_1_2_AF_controls" in INFO) || ("gg3_1_2_AF_controls" in INFO && INFO.gg3_1_2_AF_controls < 0.001)) &&', file = bashFile)
  print('    (!("gg3_1_2_AF" in INFO) || ("gg3_1_2_AF" in INFO && INFO.gg3_1_2_AF < 0.001)) &&', file = bashFile)
  print('    (!("gg3_1_2_AF_nfe" in INFO) || ("gg3_1_2_AF_nfe" in INFO && INFO.gg3_1_2_AF_nfe < 0.001)) &&', file = bashFile)
  print('    (!("gg3_1_2_AF_amr" in INFO) || ("gg3_1_2_AF_amr" in INFO && INFO.gg3_1_2_AF_amr < 0.001)) &&', file = bashFile)
  print('    (!("gg3_1_2_AF_afr" in INFO) || ("gg3_1_2_AF_afr" in INFO && INFO.gg3_1_2_AF_afr < 0.001)) &&', file = bashFile)
  print('    (!("gg3_1_2_AF_eas" in INFO) || ("gg3_1_2_AF_eas" in INFO && INFO.gg3_1_2_AF_eas < 0.001)) &&', file = bashFile)
  print('    (!("gg3_1_2_AF_sas" in INFO) || ("gg3_1_2_AF_sas" in INFO && INFO.gg3_1_2_AF_sas < 0.001)) &&', file = bashFile)
  print('    (kid.sex != "male" || variant.CHROM != "chrX" || variant.POS <= 2781479 || variant.POS >= 155701383)\' \\', file = bashFile)
  print('  --trio \'x_denovo:x_denovo(kid, mom, dad) && kid.hom_alt &&', file = bashFile)
  print('    (!("gg2_1_1_controls_AF_popmax" in INFO) || ("gg2_1_1_controls_AF_popmax" in INFO && INFO.gg2_1_1_controls_AF_popmax < 0.001)) &&', file = bashFile)
  print('    (!("gg2_1_1_AF_popmax" in INFO) || ("gg2_1_1_AF_popmax" in INFO && INFO.gg2_1_1_AF_popmax < 0.001)) &&', file = bashFile)
  print('    (!("gg3_1_2_AF_controls" in INFO) || ("gg3_1_2_AF_controls" in INFO && INFO.gg3_1_2_AF_controls < 0.001)) &&', file = bashFile)
  print('    (!("gg3_1_2_AF" in INFO) || ("gg3_1_2_AF" in INFO && INFO.gg3_1_2_AF < 0.001)) &&', file = bashFile)
  print('    (!("gg3_1_2_AF_nfe" in INFO) || ("gg3_1_2_AF_nfe" in INFO && INFO.gg3_1_2_AF_nfe < 0.001)) &&', file = bashFile)
  print('    (!("gg3_1_2_AF_amr" in INFO) || ("gg3_1_2_AF_amr" in INFO && INFO.gg3_1_2_AF_amr < 0.001)) &&', file = bashFile)
  print('    (!("gg3_1_2_AF_afr" in INFO) || ("gg3_1_2_AF_afr" in INFO && INFO.gg3_1_2_AF_afr < 0.001)) &&', file = bashFile)
  print('    (!("gg3_1_2_AF_eas" in INFO) || ("gg3_1_2_AF_eas" in INFO && INFO.gg3_1_2_AF_eas < 0.001)) &&', file = bashFile)
  print('    (!("gg3_1_2_AF_sas" in INFO) || ("gg3_1_2_AF_sas" in INFO && INFO.gg3_1_2_AF_sas < 0.001)) &&', file = bashFile)
  print('    kid.sex == "male" && variant.CHROM == "chrX" && variant.POS > 2781479 && variant.POS < 155701383\' \\', file = bashFile)
  print('  --trio \'hom:hom(kid, mom, dad) &&', file = bashFile)
  print('    (!("gg2_1_1_controls_nhomalt" in INFO) || ("gg2_1_1_controls_nhomalt" in INFO && INFO.gg2_1_1_controls_nhomalt < 2)) &&', file = bashFile)
  print('    (!("gg2_1_1_nhomalt" in INFO) || ("gg2_1_1_nhomalt" in INFO && INFO.gg2_1_1_nhomalt < 3)) &&', file = bashFile)
  print('    (!("gg3_1_2_nhomalt_controls" in INFO) || ("gg3_1_2_nhomalt_controls" in INFO && INFO.gg3_1_2_nhomalt_controls < 2)) &&', file = bashFile)
  print('    (!("gg3_1_2_nhomalt" in INFO) || ("gg3_1_2_nhomalt" in INFO && INFO.gg3_1_2_nhomalt < 3)) &&', file = bashFile)
  print('    (kid.sex != "male" || variant.CHROM != "chrX" || variant.POS <= 2781479 || variant.POS >= 155701383)\' \\', file = bashFile)
  print('  --trio \'x_hemi:x_hemi(kid, mom, dad) && kid.hom_alt &&', file = bashFile)
  print('    (!("gg2_1_1_controls_AF_popmax" in INFO) || ("gg2_1_1_controls_AF_popmax" in INFO && INFO.gg2_1_1_controls_AF_popmax < 0.001)) &&', file = bashFile)
  print('    (!("gg2_1_1_AF_popmax" in INFO) || ("gg2_1_1_AF_popmax" in INFO && INFO.gg2_1_1_AF_popmax < 0.001)) &&', file = bashFile)
  print('    (!("gg3_1_2_AF_controls" in INFO) || ("gg3_1_2_AF_controls" in INFO && INFO.gg3_1_2_AF_controls < 0.001)) &&', file = bashFile)
  print('    (!("gg3_1_2_AF" in INFO) || ("gg3_1_2_AF" in INFO && INFO.gg3_1_2_AF < 0.001)) &&', file = bashFile)
  print('    (!("gg3_1_2_AF_nfe" in INFO) || ("gg3_1_2_AF_nfe" in INFO && INFO.gg3_1_2_AF_nfe < 0.001)) &&', file = bashFile)
  print('    (!("gg3_1_2_AF_amr" in INFO) || ("gg3_1_2_AF_amr" in INFO && INFO.gg3_1_2_AF_amr < 0.001)) &&', file = bashFile)
  print('    (!("gg3_1_2_AF_afr" in INFO) || ("gg3_1_2_AF_afr" in INFO && INFO.gg3_1_2_AF_afr < 0.001)) &&', file = bashFile)
  print('    (!("gg3_1_2_AF_eas" in INFO) || ("gg3_1_2_AF_eas" in INFO && INFO.gg3_1_2_AF_eas < 0.001)) &&', file = bashFile)
  print('    (!("gg3_1_2_AF_sas" in INFO) || ("gg3_1_2_AF_sas" in INFO && INFO.gg3_1_2_AF_sas < 0.001)) &&', file = bashFile)
  print('    (!("gg2_1_1_AC_male_controls" in INFO) || ("gg2_1_1_AC_male_controls" in INFO && INFO.gg2_1_1_AC_male_controls < 2)) &&', file = bashFile)
  print('    (!("gg2_1_1_AC_male" in INFO) || ("gg2_1_1_AC_male" in INFO && INFO.gg2_1_1_AC_male < 3)) &&', file = bashFile)
  print('    (!("gg3_1_2_AC_controls_XY" in INFO) || ("gg3_1_2_AC_controls_XY" in INFO && INFO.gg3_1_2_AC_controls_XY < 2)) &&', file = bashFile)
  print('    (!("gg3_1_2_AC_XY" in INFO) || ("gg3_1_2_AC_XY" in INFO && INFO.gg3_1_2_AC_XY < 3)) &&', file = bashFile)
  print('    (!("gg2_1_1_controls_nhomalt" in INFO) || ("gg2_1_1_controls_nhomalt" in INFO && INFO.gg2_1_1_controls_nhomalt < 2)) &&', file = bashFile)
  print('    (!("gg2_1_1_nhomalt" in INFO) || ("gg2_1_1_nhomalt" in INFO && INFO.gg2_1_1_nhomalt < 3)) &&', file = bashFile)
  print('    (!("gg3_1_2_nhomalt_controls" in INFO) || ("gg3_1_2_nhomalt_controls" in INFO && INFO.gg3_1_2_nhomalt_controls < 2)) &&', file = bashFile)
  print('    (!("gg3_1_2_nhomalt" in INFO) || ("gg3_1_2_nhomalt" in INFO && INFO.gg3_1_2_nhomalt < 3)) &&', file = bashFile)
  print('    kid.sex == "male" && variant.CHROM == "chrX" && variant.POS > 2781479 && variant.POS < 155701383\' \\', file = bashFile)
  print('  --trio \'comphet_side:comphet_side(kid, mom, dad) &&', file = bashFile)
  print('    (!("gg2_1_1_controls_nhomalt" in INFO ) || ("gg2_1_1_controls_nhomalt" in INFO && INFO.gg2_1_1_controls_nhomalt < 2)) &&', file = bashFile)
  print('    (!("gg2_1_1_nhomalt" in INFO) || ("gg2_1_1_nhomalt" in INFO && INFO.gg2_1_1_nhomalt < 3)) &&', file = bashFile)
  print('    (!("gg3_1_2_nhomalt_controls" in INFO) || ("gg3_1_2_nhomalt_controls" in INFO && INFO.gg3_1_2_nhomalt_controls < 2)) &&', file = bashFile)
  print('    (!("gg3_1_2_nhomalt" in INFO) || ("gg3_1_2_nhomalt" in INFO && INFO.gg3_1_2_nhomalt < 3))\' \\', file = bashFile)
  print('  >> $STDOUT 2>> $STDERR', file = bashFile)
  print('echo "complete"', file = bashFile)
  print(file = bashFile)

  # Remove comphet_sides from the Slivar filtered vcf
  print('# Removing comphet_sides from the filtered vcf', file = bashFile)
  print('echo -n "Removing comphet_sides from filtered vcf..."', file = bashFile)
  print('$SLIVAR expr \\', file = bashFile)
  print('  --vcf $SLIVAR1VCF \\', file = bashFile)
  print('  --ped $PED \\', file = bashFile)
  print('  --js $JS \\', file = bashFile)
  print('  --pass-only \\', file = bashFile)
  print('  --trio "denovo:denovo(kid, mom, dad) && (kid.sex != \'male\' || variant.CHROM != \'chrX\' || variant.POS <= 2781479 || variant.POS >= 155701383)" \\', file = bashFile)
  print('  --trio "x_denovo:x_denovo(kid, mom, dad) && kid.sex == \'male\' && variant.CHROM == \'chrX\' && variant.POS > 2781479 && variant.POS < 155701383" \\', file = bashFile)
  print('  --trio "hom:hom(kid, mom, dad) && (kid.sex != \'male\' || variant.CHROM != \'chrX\' || variant.POS <= 2781479 || variant.POS >= 155701383)" \\', file = bashFile)
  print('  --trio "x_hemi:x_hemi(kid, mom, dad) && kid.sex == \'male\' && variant.CHROM == \'chrX\' && variant.POS > 2781479 && variant.POS < 155701383" \\', file = bashFile)
  print('  2>> $STDERR \\', file = bashFile)
  print('| $BCFTOOLS view -O z -o $SLIVAR2VCF - \\', file = bashFile)
  print('  >> $STDOUT 2>> $STDERR', file = bashFile)
  print('echo "complete"', file = bashFile)
  print(file = bashFile)

  # Now filter the slivar filtered vcf file for paired compound gets
  print('# Searching for paired compound hets in the slivar filtered vcf', file = bashFile)
  print('echo -n "Searching for compound heterozygotes..."', file = bashFile)
  print('$SLIVAR compound-hets \\', file = bashFile)
  print('  --vcf $SLIVAR1VCF \\', file = bashFile)
  print('  --ped $PED \\', file = bashFile)
  print('  --sample-field denovo \\', file = bashFile)
  print('  --sample-field x_denovo \\', file = bashFile)
  print('  --sample-field comphet_side \\', file = bashFile)
  print('  --skip intergenic_region,upstream_gene,downstream_gene,non_coding,non_coding_transcript,non_coding_transcript_exon,NMD_transcript,3_prime_UTR \\', file = bashFile)
  print('  2>> $STDERR \\', file = bashFile)
  print('| $BCFTOOLS view -O z -o $COMPHETSVCF - \\', file = bashFile)
  print('  >> $STDOUT 2>> $STDERR', file = bashFile)
  print('echo "complete"', file = bashFile)
  print(file = bashFile)

  # Index the files before concatenating
  print('# Index the slivar files', file = bashFile)
  print('echo -n "Indexing slivar vcf files..."', file = bashFile)
  print('$BCFTOOLS index -t $SLIVAR2VCF', file = bashFile)
  print('$BCFTOOLS index -t $COMPHETSVCF', file = bashFile)
  print('echo "complete"', file = bashFile)
  print(file = bashFile)

  # Determine the number of variants in the final slivar files. Concatting files with zero variants will fail, so
  # this step is dependent on the results for this case
  print('# Determine the number of variants in the slivar files', file = bashFile)
  print('echo -n "Determining number of variants in slivar output files..."', file = bashFile)
  print('NO_SLIVAR1=`$BCFTOOLS view -H $SLIVAR1VCF | wc -l`', file = bashFile)
  print('NO_SLIVAR2=`$BCFTOOLS view -H $SLIVAR2VCF | wc -l`', file = bashFile)
  print('NO_COMPHETS=`$BCFTOOLS view -H $COMPHETSVCF | wc -l`', file = bashFile)
  print('if [[ $NO_SLIVAR2 == 0 ]]; then', file = bashFile)
  print('  if [[ $NO_COMPHETS == 0 ]]; then', file = bashFile)

  # If there are no variants in either of the files to be concatted, and there are less than 30 variants in the
  # original filtered slivar file, upload these variants. If there are more than 30 variants in the original
  # filtered slivar vcf and none in either of the final 2 files, no prioritized variats will be uploaded
  print('    if [[ $NO_SLIVAR1 < 10 ]]; then', file = bashFile)
  print('      echo "mv $SLIVAR1 $SLIVARFILTEREDVCF"', file = bashFile)
  print('    fi', file = bashFile)

  # If only the comp hets file has variants, use this as the final slivar output
  print('  else', file = bashFile)
  print('    echo "mv $COMPHETSVCF $SLIVARFILTEREDVCF"', file = bashFile)
  print('  fi', file = bashFile)

  # If only the first filtered slivar file has variants, e.g. no comp hets, use this file
  print('else', file = bashFile)
  print('  if [[ $NO_COMPHETS == 0 ]]; then', file = bashFile)
  print('    echo "mv $SLIVAR2VCF $SLIVARFILTEREDVCF"', file = bashFile)

  # If bothe final slivar files have results, concat them
  print('  else', file = bashFile)
  print('    echo -n "Concatenating files..."', file = bashFile)
  print('    $BCFTOOLS concat -O u -a -D $SLIVAR2VCF $COMPHETSVCF 2>> $STDERR \\', file = bashFile)
  print('    | $BCFTOOLS sort -O z -o $SLIVARFILTEREDVCF - >> $STDOUT 2>> $STDERR', file = bashFile)
  print('  fi', file = bashFile)
  print('fi', file = bashFile)
  print('echo "complete"', file = bashFile)
  print(file = bashFile)

# If the script fails, provide an error message and exit
def fail(message):
  print(message, sep = "")
  exit(1)
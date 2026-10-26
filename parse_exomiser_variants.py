from os.path import exists

import argparse

# Parse an Exomiser variants.tsv file into a tsv that can be uploaded to Mosaic. The Exomiser annotations are
# semi-public annotations that already exist in Mosaic: their values are specific to the project they are
# uploaded to, so no project specific (private) annotations need to be created. Exomiser runs four ways (SNVs /
# Indels, or SNVs / Indels + SVs, each with or without HPO terms), and each run has its own set of the same 10
# annotations, distinguished by a uid suffix
def main():

  # Parse the command line
  args = parse_command_line()

  # Check that the exomiser output file exists
  if not exists(args.input):
    fail('The file ' + str(args.input) + ' was not found')

  # Read the exomiser variants. A variant can appear multiple times (e.g. under different modes of inheritance).
  # Mosaic allows a single value per variant, so keep the values from the best (lowest) rank, but list all of
  # the modes of inheritance
  variants = {}
  columns = False
  with open(args.input, 'r') as variants_file:
    for line in variants_file:
      fields = line.rstrip('\n').split('\t')

      # Get the column positions from the header line
      if line.startswith('#'):
        header = [field.lstrip('#') for field in fields]
        missing_columns = [column for column in required_columns if column not in header]
        if missing_columns:
          fail('The exomiser file (' + str(args.input) + ') is missing the columns: ' + ', '.join(missing_columns))
        columns = {column: header.index(column) for column in required_columns}
        continue
      if not columns:
        fail('The exomiser file (' + str(args.input) + ') does not have a header line')

      variant = get_variant(fields, columns)
      values = [fields[columns[column]] for column in exomiser_fields.values()]
      moi = fields[columns['MOI']]
      rank = float(fields[columns['RANK']])
      if variant not in variants:
        variants[variant] = {'rank': rank, 'values': values, 'moi': [moi]}
      else:
        if rank < variants[variant]['rank']:
          variants[variant]['rank'] = rank
          variants[variant]['values'] = values
        if moi not in variants[variant]['moi']:
          variants[variant]['moi'].append(moi)

  # Write the variants to the output tsv. The header uses the uid for this exomiser run and the version
  moi_index = list(exomiser_fields.keys()).index('moi')
  with open(args.output, 'w') as output_file:
    print('CHROM', 'START', 'END', 'REF', 'ALT', '\t'.join(str(uid) + str(args.uid_suffix) + '@' + str(args.version) for uid in exomiser_fields), sep = '\t', file = output_file)
    for variant in variants:
      values = variants[variant]['values']
      values[moi_index] = ','.join(variants[variant]['moi'])
      print('\t'.join(variant), '\t'.join(values), sep = '\t', file = output_file)

# Input options
def parse_command_line():
  parser = argparse.ArgumentParser(description='Process the command line')

  # Required arguments
  parser.add_argument('--input', '-i', required = True, metavar = 'string', help = 'The exomiser variants.tsv file')
  parser.add_argument('--output', '-o', required = True, metavar = 'string', help = 'The output tsv file to upload to Mosaic')
  parser.add_argument('--uid_suffix', '-s', required = True, metavar = 'string', choices = allowed_suffixes, help = 'The suffix of the Mosaic annotation uids for this exomiser run. Allowed values: ' + ', '.join(allowed_suffixes))

  # Optional arguments
  parser.add_argument('--version', '-v', required = False, metavar = 'string', default = 'v15.1.0_data-2602_v1', help = 'The Mosaic annotation version to upload to')

  return parser.parse_args()

# Get the Mosaic variant identifier (CHROM, START, END, REF, ALT) for an exomiser variant
def get_variant(fields, columns):

  # Mosaic does not use the "chr" prefix and the end position is exclusive
  chrom = fields[columns['CONTIG']]
  chrom = chrom[3:] if chrom.startswith('chr') else chrom
  start = fields[columns['START']]
  end = str(int(fields[columns['END']]) + 1)
  ref = fields[columns['REF']]
  alt = fields[columns['ALT']]

  # Structural variants are identified in Mosaic by the SV type in the ALT column, e.g. <DEL> or <DUP:TANDEM> are
  # DEL and DUP, and breakends are BND
  if alt.startswith('<'):
    alt = alt.strip('<>').split(':')[0]
  elif '[' in alt or ']' in alt:
    alt = 'BND'

  # For small variants, if the ref or alt allele are 'N', change them to '*'
  else:
    ref = '*' if ref == 'N' else ref
    alt = '*' if alt == 'N' else alt

  return (chrom, start, end, ref, alt)

# If the script fails, provide an error message and exit
def fail(message):
  print('ERROR: ', message, sep = '')
  exit(1)

# The suffixes of the Mosaic uids for each exomiser run
allowed_suffixes = ['_snv_indel', '_snv_indel_no_hpo', '_snv_indel_sv', '_snv_indel_sv_no_hpo']

# The Mosaic annotation uids (without the suffix) and the exomiser column they are populated from
exomiser_fields = {'rank': 'RANK',
                   'p_value': 'P-VALUE',
                   'moi': 'MOI',
                   'variant_score': 'EXOMISER_VARIANT_SCORE',
                   'gene_variant_score': 'EXOMISER_GENE_VARIANT_SCORE',
                   'gene_phenotype_score': 'EXOMISER_GENE_PHENO_SCORE',
                   'gene_combined_score': 'EXOMISER_GENE_COMBINED_SCORE',
                   'contributing_variant': 'CONTRIBUTING_VARIANT',
                   'acmg_classification': 'EXOMISER_ACMG_CLASSIFICATION',
                   'acmg_evidence': 'EXOMISER_ACMG_EVIDENCE'}

# All the columns required from the exomiser variants file
required_columns = list(exomiser_fields.values()) + ['CONTIG', 'START', 'END', 'REF', 'ALT']

if __name__ == "__main__":
  main()

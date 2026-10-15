import argparse
import math
import os
import subprocess
import sys

import read_resource_jsons as read_resources

# The SV type is written to the ALT column and is also output as an annotation
SVTYPE = 'SVTYPE'

# CNV records have an SVTYPE of CNV, so the SV type is taken from the ALT allele instead
CNV_ALTS = {'<DEL>': 'DEL', '<DUP>': 'DUP'}

# Mosaic cannot store numeric values outside of these bounds
MIN_VALUE = 1e-37
MAX_VALUE = 1e37

def main():

  # Parse the command line
  args = parseCommandLine()
  if not args.sv_vcf and not args.cnv_vcf:
    fail('At least one of --sv_vcf (-v) or --cnv_vcf (-c) must be supplied')

  # Define the executable bcftools command
  bcftools = os.path.join(args.tools_directory, 'bcftools', 'bcftools') if args.tools_directory else 'bcftools'

  # Read the mosaicJson file and get the name, uid and type of every annotation to be uploaded to Mosaic. SVTYPE is
  # handled separately as it is always extracted (to populate the ALT column), and is only output if the record has
  # other annotations
  mosaic_info = read_resources.read_mosaic_json(args.mosaic_json, args.reference)
  annotations = {}
  svtype_uid = False
  for resource in mosaic_info['resources']:
    for annotation, info in mosaic_info['resources'][resource]['annotations'].items():
      if annotation == SVTYPE:
        if svtype_uid:
          fail('Annotation ' + SVTYPE + ' is defined for multiple resources in ' + str(args.mosaic_json))
        svtype_uid = info['uid']
        continue
      if annotation in annotations:
        fail('Annotation ' + str(annotation) + ' is defined for multiple resources in ' + str(args.mosaic_json))
      if str(info['type']) not in ['integer', 'float', 'string']:
        fail('Annotation type (' + str(info['type']) + ') for annotation ' + str(annotation) + ', is unknown')
      annotations[annotation] = {'uid': info['uid'], 'type': str(info['type'])}

  if not svtype_uid:
    fail('The Mosaic json (' + str(args.mosaic_json) + ') does not define the ' + SVTYPE + ' annotation')

  # Remove any output from a previous run so a stale tsv is never uploaded
  if os.path.exists(args.output_tsv):
    os.remove(args.output_tsv)

  # If there are no annotations to upload, there is nothing to do
  if not annotations:
    print('No annotations (other than ' + SVTYPE + ') are defined in ' + str(args.mosaic_json) + '. No tsv will be generated')
    return

  # Parse the vcf files and build the records for the tsv. Records with the same coordinates and SV type are only
  # output once, including across the SV and CNV vcfs
  records = {}
  for vcf in [args.sv_vcf, args.cnv_vcf]:
    if vcf:
      process_vcf(bcftools, vcf, annotations, records)

  # Only generate the tsv if there are records with annotations
  if not records:
    print('No variants with annotations found. No tsv will be generated')
    return
  with open(args.output_tsv, 'w') as output_file:
    print('CHROM\tSTART\tEND\tREF\tALT\t', svtype_uid, '\t', '\t'.join(str(x['uid']) for x in annotations.values()), sep = '', file = output_file)

    # The SV type is the last field of the identifier (the ALT column)
    for identifier, values in records.items():
      print('\t'.join(identifier + identifier[-1:] + values), file = output_file)
  print('Wrote ', len(records), ' records to ', args.output_tsv, sep = '')

# Input options
def parseCommandLine():
  parser = argparse.ArgumentParser(description = 'Generate a tsv of SV / CNV annotations for upload to Mosaic')
  parser.add_argument('--sv_vcf', '-v', required = False, metavar = 'string', help = 'The input SV vcf file')
  parser.add_argument('--cnv_vcf', '-c', required = False, metavar = 'string', help = 'The input CNV vcf file')
  parser.add_argument('--output_tsv', '-o', required = True, metavar = 'string', help = 'The output tsv file')
  parser.add_argument('--reference', '-r', required = True, metavar = 'string', help = 'The genome reference file used')
  parser.add_argument('--tools_directory', '-s', required = False, metavar = 'string', help = 'The path to the directory where the tools live. If not set, bcftools must be in the PATH')
  parser.add_argument('--mosaic_json', '-m', required = True, metavar = 'string', help = 'The json file describing the Mosaic parameters')

  return parser.parse_args()

# Extract the annotations from the vcf with bcftools query, check the values and add the record to the records dictionary
# if it has at least one annotation
def process_vcf(bcftools, vcf, annotations, records):
  if not os.path.exists(vcf):
    fail('Input vcf file (' + str(vcf) + ') does not exist')

  # bcftools query fails if a requested INFO field is not in the header, so only query annotations that are present. Any
  # missing annotations are output as empty values
  info_fields = get_info_fields(bcftools, vcf)
  if SVTYPE not in info_fields:
    fail('The vcf file (' + str(vcf) + ') does not contain the ' + SVTYPE + ' INFO field')
  present = [annotation for annotation in annotations if annotation in info_fields]
  for annotation in annotations:
    if annotation not in info_fields:
      print('WARNING: Annotation ', annotation, ' is not present in ', vcf, ' and will be left empty', sep = '', file = sys.stderr)
  if not present:
    print('WARNING: None of the annotations are present in ', vcf, sep = '', file = sys.stderr)
    return

  # Stream the records rather than reading the whole output into memory
  query_format = '%CHROM\\t%POS\\t%END\\t%REF\\t%ALT\\t%INFO/' + SVTYPE + ''.join('\\t%INFO/' + x for x in present) + '\\n'
  process = subprocess.Popen([bcftools, 'query', '-f', query_format, vcf], stdout = subprocess.PIPE, text = True)
  duplicates = 0
  for record in process.stdout:
    fields = record.rstrip('\n').split('\t')
    chrom, start, end, ref, alt, svtype = fields[0:6]
    queried = dict(zip(present, fields[6:]))

    # Check the values for each annotation and determine if the record has any annotations. Records with no annotations
    # (all values are missing, zero or empty) are skipped
    values = []
    has_annotation = False
    for annotation, info in annotations.items():
      value = check_value(queried.get(annotation, '.'), info['type'], annotation, record)
      if value_present(value, info['type']):
        has_annotation = True
      values.append(value)
    if not has_annotation:
      continue

    # Replace the CNV SV type with DEL or DUP based on the ALT allele. If the ALT is anything else (e.g. multiple
    # alleles), the CNV SV type is kept
    if svtype == 'CNV':
      svtype = CNV_ALTS.get(alt, svtype)

    # Mosaic does not use the "chr" prefix and the end position is exclusive
    identifier = (chrom.removeprefix('chr'), start, str(int(end) + 1), ref, svtype)
    if identifier in records:
      duplicates += 1
    else:
      records[identifier] = tuple(values)

  process.stdout.close()
  if process.wait() != 0:
    fail('bcftools query failed for ' + str(vcf))
  if duplicates:
    print('WARNING: ', duplicates, ' duplicate records in ', vcf, ' were skipped', sep = '', file = sys.stderr)

# Get the ids of all INFO fields in the vcf header
def get_info_fields(bcftools, vcf):
  try:
    result = subprocess.run([bcftools, 'view', '-h', vcf], capture_output = True, text = True)
  except FileNotFoundError:
    fail('bcftools could not be found at ' + str(bcftools))
  if result.returncode != 0:
    fail('Could not read the header of ' + str(vcf) + ': ' + result.stderr.strip())
  info_fields = set()
  for line in result.stdout.splitlines():
    if line.startswith('##INFO=<ID='):
      info_fields.add(line[11:].split(',', 1)[0])

  return info_fields

# Check the value of an annotation. Missing values are output as empty. Numeric annotations can have multiple comma
# separated values, and all must be numbers within the bounds that Mosaic can store (1E-37 < |x| < 1E37)
def check_value(value, tag_type, annotation, record):
  if tag_type == 'string':
    return '' if value == '.' else value

  final_values = []
  for value in value.split(','):
    if value == '.':
      final_values.append('')
      continue
    try:
      number = float(value)
    except ValueError:
      fail('Invalid value for annotation ' + str(annotation) + ': ' + record.rstrip())
    if math.isnan(number):
      final_values.append('')
    elif abs(number) <= MIN_VALUE:
      final_values.append('0')
    elif abs(number) >= MAX_VALUE:
      final_values.append('-1e37' if number < 0 else '1e37')
    else:
      final_values.append(value)

  return ','.join(final_values)

# Determine if a checked value contains an annotation. Numeric values of zero are treated as no annotation
def value_present(value, tag_type):
  if tag_type == 'string':
    return value != ''

  return any(x != '' and float(x) != 0. for x in value.split(','))

# If the script fails, provide an error message and exit
def fail(message):
  print(message, file = sys.stderr)
  sys.exit(1)

if __name__ == "__main__":
  main()

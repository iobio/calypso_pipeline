from os.path import exists

import argparse
import glob
import os
import sys

# Merge multiple Calypso annotation tsv files (e.g. the HGVS tsv generated for many projects) into a single tsv for
# upload to Mosaic. All input files must share the same header (i.e. the same annotation uids). Variants present in
# multiple files are written once. If the same variant has different annotation values in different files, the first
# value seen is kept and the conflict is reported
def main():

  # Parse the command line
  args = parse_command_line()

  # Get the list of files to merge
  files = get_input_files(args)
  if not files: fail('No input tsv files were found')
  print('Merging ', len(files), ' files', sep = '')

  header    = None
  variants  = {}
  conflicts = 0
  for i, filename in enumerate(files):
    with open(filename, 'r') as tsv:
      fileHeader = tsv.readline().rstrip('\n')

      # Skip empty files
      if not fileHeader:
        print('WARNING: Empty file skipped: ', filename, sep = '', file = sys.stderr)
        continue

      # All files must have the same columns. The annotation uids in the header define which annotations the values
      # are uploaded to, so files with different uids (e.g. project private annotations) cannot be merged
      if header is None:
        header = fileHeader
        if not header.startswith('CHROM\tSTART\tEND\tREF\tALT'): fail('Unexpected header in ' + filename + ': ' + header)
        noColumns = len(header.split('\t'))
      elif fileHeader != header:
        fail('Header in ' + filename + ' does not match the header in ' + files[0] + ':\n  ' + fileHeader + '\n  ' + header)

      for lineNo, line in enumerate(tsv, start = 2):
        line = line.rstrip('\n')
        if not line: continue
        fields = line.split('\t')

        # Pad lines that are missing trailing empty values
        if len(fields) < noColumns: fields += [''] * (noColumns - len(fields))
        elif len(fields) > noColumns: fail('Too many columns at line ' + str(lineNo) + ' of ' + filename)

        # Skip variants with no annotation values
        values = '\t'.join(fields[5:])
        if not values.strip('\t'): continue

        key = (fields[0], int(fields[1]), int(fields[2]), fields[3], fields[4])
        if key not in variants: variants[key] = values
        elif variants[key] != values:
          conflicts += 1
          if conflicts <= args.max_conflict_warnings:
            print('WARNING: Conflicting values for ', ':'.join(str(x) for x in key), ' in ', filename, ' (kept first value)', sep = '', file = sys.stderr)

    if (i + 1) % 100 == 0: print('  processed ', i + 1, ' files, ', len(variants), ' unique variants', sep = '')

  if header is None: fail('All input files were empty')

  # Write out the merged file, sorted by chromosome and position
  with open(args.output, 'w') as output:
    print(header, file = output)
    for key in sorted(variants, key = lambda x: (chrom_order(x[0]), x[1], x[2], x[3], x[4])):
      print(key[0], key[1], key[2], key[3], key[4], variants[key], sep = '\t', file = output)

  print('Wrote ', len(variants), ' unique variants to ', args.output, sep = '')
  if conflicts > 0: print('WARNING: ', conflicts, ' variants had conflicting values across files', sep = '', file = sys.stderr)

# Input options
def parse_command_line():
  parser = argparse.ArgumentParser(description='Merge annotation tsv files into a single tsv for upload to Mosaic')
  parser.add_argument('--input_files', '-i', nargs = '+', metavar = 'string', help = 'The tsv files to merge (shell globs are expanded)')
  parser.add_argument('--input_list', '-l', metavar = 'string', help = 'A file containing the paths of the tsv files to merge, one per line')
  parser.add_argument('--input_directory', '-d', metavar = 'string', help = 'A directory to search recursively for tsv files to merge')
  parser.add_argument('--file_name', '-f', default = 'hgvs.tsv', metavar = 'string', help = 'The tsv file name to look for when using --input_directory (default: hgvs.tsv)')
  parser.add_argument('--output', '-o', required = True, metavar = 'string', help = 'The merged output tsv file')
  parser.add_argument('--max_conflict_warnings', '-w', type = int, default = 20, metavar = 'integer', help = 'The maximum number of conflict warnings to print (default: 20)')

  return parser.parse_args()

# Collect the input files from all of the provided sources
def get_input_files(args):
  files = []
  if args.input_files:
    for pattern in args.input_files: files += sorted(glob.glob(pattern)) or [pattern]
  if args.input_list:
    if not exists(args.input_list): fail('The file ' + args.input_list + ' was not found')
    with open(args.input_list, 'r') as fileList:
      files += [line.strip() for line in fileList if line.strip()]
  if args.input_directory:
    files += sorted(glob.glob(os.path.join(args.input_directory, '**', args.file_name), recursive = True))

  # Remove duplicates, preserving order, and the output file if it was picked up
  output = os.path.abspath(args.output)
  files  = [f for f in dict.fromkeys(files) if os.path.abspath(f) != output]
  for filename in files:
    if not exists(filename): fail('The file ' + filename + ' was not found')

  return files

# Order chromosomes numerically, then X, Y, M/MT, then anything else
def chrom_order(chrom):
  if chrom.isdigit(): return (0, int(chrom), '')
  order = {'X': 23, 'Y': 24, 'M': 25, 'MT': 25}
  return (0, order[chrom], '') if chrom in order else (1, 0, chrom)

# If the script fails, provide an error message and exit
def fail(message):
  print('ERROR: ', message, sep = '', file = sys.stderr)
  exit(1)

if __name__ == "__main__":
  main()

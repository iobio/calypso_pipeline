import argparse
import os
import sys

# Split a large annotation tsv (e.g. the output of merge_annotation_tsvs.py) into multiple tsv files that are each small
# enough to upload to Mosaic. Mosaic rejects annotation uploads over 1000MB (maxFileSize in
# server/routes/external/v1/annotations/annotations-routes.js), so the default is 900MB to leave headroom. Every output
# file starts with the header from the input file, and lines are never split across files
def main():

  # Parse the command line
  args = parse_command_line()
  if not os.path.exists(args.input): fail('The file ' + args.input + ' was not found')
  maxBytes = int(args.max_size_mb * 1e6)

  # Default the output prefix to the input file name without the .tsv extension
  prefix = args.output_prefix or (args.input[:-4] if args.input.endswith('.tsv') else args.input)

  # Read in binary mode so the byte counts match the size of the files that are written
  with open(args.input, 'rb') as tsv:
    header = tsv.readline()
    if not header: fail('The file ' + args.input + ' is empty')
    if not header.startswith(b'CHROM\t'): fail('Unexpected header in ' + args.input + ': ' + header.decode().rstrip())
    if not header.endswith(b'\n'): header += b'\n'
    if len(header) >= maxBytes: fail('The header alone is larger than the maximum file size')

    outputFiles = []
    output      = None
    for line in tsv:
      if not line.strip(): continue
      if not line.endswith(b'\n'): line += b'\n'

      # Start a new file if this line would take the current one over either limit
      if output is None or size + len(line) > maxBytes or (args.max_lines and lines >= args.max_lines):
        if output: output.close()
        if len(header) + len(line) > maxBytes: fail('A single line is larger than the maximum file size')
        filename = prefix + '_' + str(len(outputFiles) + 1).zfill(3) + '.tsv'
        output   = open(filename, 'wb')
        output.write(header)
        size  = len(header)
        lines = 0
        outputFiles.append([filename, 0])

      output.write(line)
      size  += len(line)
      lines += 1
      outputFiles[-1][1] = lines

    if output: output.close()

  if not outputFiles: fail('The file ' + args.input + ' contains only a header, so no files were written')
  for filename, lines in outputFiles:
    print(filename, '\t', lines, ' variants\t', round(os.path.getsize(filename) / 1e6, 1), 'MB', sep = '')
  print('Wrote ', len(outputFiles), ' files', sep = '')

# Input options
def parse_command_line():
  parser = argparse.ArgumentParser(description='Split an annotation tsv into multiple files small enough to upload to Mosaic')
  parser.add_argument('--input', '-i', required = True, metavar = 'string', help = 'The annotation tsv file to split')
  parser.add_argument('--output_prefix', '-o', metavar = 'string', help = 'The prefix for the output files, which are named <prefix>_001.tsv etc. (default: the input file name without .tsv)')
  parser.add_argument('--max_size_mb', '-s', type = float, default = 900, metavar = 'float', help = 'The maximum size of each output file in MB (default: 900, Mosaic\'s limit is 1000)')
  parser.add_argument('--max_lines', '-n', type = int, metavar = 'integer', help = 'Optionally, also limit the number of variants in each output file')

  return parser.parse_args()

# If the script fails, provide an error message and exit
def fail(message):
  print('ERROR: ', message, sep = '', file = sys.stderr)
  exit(1)

if __name__ == "__main__":
  main()

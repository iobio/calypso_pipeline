#!/usr/bin/python

from __future__ import print_function
from datetime import date
from os.path import exists

import sys
import os
import argparse
import json
import math
import glob

# Add the path of the common functions and import them
from sys import path
path.append(os.path.dirname(__file__) + "/mosaic_commands")
import mosaic_config
import api_samples as api_s

def main():
  global mosaicConfig
  global mosaicInfo

  # Parse the command line
  args = parseCommandLine()

  # Check the supplied parameters are as expected
  #checkArguments(args)

  # Create a directory to store all created files
  setWorkingDir()

  # Parse the Mosaic config file to get the token and url for the api calls
  mosaicRequired = {"token": True, "url": True, "attributeProjectId": False}
  mosaicConfig   = mosaic_config.parseConfig(args, mosaicRequired)

  # Determine the family structure
  getProband(args)

  # Parse the patients HPO terms
  parsePatientHpoTerms(args.hpo_terms)

  # Get HPO / gene associations
  parseHpoGene(args)

  # Parse vcf file
  parseVcf(args)

# Input options
def parseCommandLine():
  global version

  parser = argparse.ArgumentParser(description='Process the command line')

  # Required arguments
  #parser.add_argument('--reference', '-r', required = True, metavar = "string", help = "The reference genome to use. Allowed values: '37', '38'")
  #parser.add_argument('--family_type', '-f', required = True, metavar = "string", help = "The familty structure. Allowed values: 'singleton', 'duo', 'trio'")
  #parser.add_argument('--data_directory', '-d', required = True, metavar = "string", help = "The path to the directory where the resources live")
  #parser.add_argument('--input_vcf', '-i', required = True, metavar = "string", help = "The input vcf file to annotate")
  parser.add_argument('--hpo', '-o', required = True, metavar = "string", help = "The gene:HPO associations file")
  parser.add_argument('--project_id', '-p', required = True, metavar = "string", help = "The project id that variants will be uploaded to")
  parser.add_argument('--hpo_terms', '-r', required = True, metavar = "string", help = "The HPO terms associated with the patient")

  # Optional reanalysis arguments
  #parser.add_argument('--previous_vcf', '-s', required = False, metavar = "string", help = "A previously annotated Calypso vcf file")

  # Optional pipeline arguments
  parser.add_argument('--ped', '-e', required = False, metavar = "string", help = "The pedigree file for the family. Not required for singletons")
  #parser.add_argument('--resource_json', '-j', required = False, metavar = "string", help = "The json file describing the annotation resources")

  # Optional mosaic arguments
  parser.add_argument('--config', '-c', required = False, metavar = "string", help = "The config file for Mosaic")
  parser.add_argument('--token', '-t', required = False, metavar = "string", help = "The Mosaic authorization token")
  parser.add_argument('--url', '-u', required = False, metavar = "string", help = "The base url for Mosaic")
  parser.add_argument('--attributes_project', '-a', required = False, metavar = "integer", help = "The Mosaic project id that contains public attributes")
  parser.add_argument('--mosaic_json', '-m', required = False, metavar = "string", help = "The json file describing the Mosaic parameters")

  # Version
  parser.add_argument('--version', '-v', action="version", version='HPO prioritization pipeline version: ' + str(version))

  return parser.parse_args()

# Check the supplied arguments
def checkArguments(args):
  global allowedReferences
  global allowedFamily
  global isFamily

  # Check the reference is allowed
  if args.reference not in allowedReferences:
    print("The allowed references (--reference, -r) are:")
    for ref in allowedReferences: print("  ", ref, sep = "")
    exit(1)

  # Check that the input files exist and have the correct extensions
  #if not exists(args.input_vcf): fail("The vcf file could not be found")
  #elif not args.input_vcf.endswith(".vcf.gz"): fail("The input vcf file (--vcf, -v) must be a compressed, indexed vcf and have the extension '.vcf.gz'")

# Create a working directory
def setWorkingDir():
  global version
  global workingDir

  workingDir += version + "/"

  # If the directory doesn't exist, create it
  if not os.path.exists(workingDir): os.makedirs(workingDir)
  
#########
#########
######### Handle different delimiters etc.
#########
#########
# Parse the patients HPO terms
def parsePatientHpoTerms(hpoTerms):
  global patientTerms

  patientTerms = hpoTerms.split(",")

# Determine the id of the proband
def getProband(args):
  global mosaicConfig
  global samples
  global proband
  global mother
  global father

  # Open the ped file and get the samples
  try: pedFile = open(args.ped, "r")
  except: fail("Couldn't open the ped file (" + args.ped + ")")

  # Get the samples Mosaic id and store
  mosaicSamples = {}
  for sample in json.loads(os.popen(api_s.getSamples(mosaicConfig, args.project_id)).read()):
    mosaicSamples[sample["name"]] = sample["id"]

  # Get information on all the samples
  noAffected = 0
  for line in pedFile:

    # Ignore the header line
    if not line.startswith("#"):
      fields   = line.rstrip().split("\t")
      sample   = fields[1]

      # Determine if the sample has parents
      father = False if fields[2] == "0" else fields[2]
      mother = False if fields[3] == "0" else fields[3]
      if mother and father: noParents = 2
      elif not mother and not father: noParents = 0
      else: noParents = 1

      # Determine if the sample is affected
      try: affectedStatus = int(fields[5])
      except: fail("Ped file does not contain enough fields")
      if int(fields[5]) == 2:
        isAffected = True
        noAffected += 1
        proband = sample
      else: isAffected = False
      samples[sample] = {"noParents": noParents, "father": father, "mother": mother, "isAffected": isAffected, "relationship": False}
      if isAffected: samples[sample]["relationship"] = "Proband"

      # Attach the Mosaic sample id
      try: samples[sample]["mosaicId"] = mosaicSamples[sample]
      except: fail("Sample " + str(sample) + " is not present in the Mosaic project")

  # Check that the ped file conforms to the supplied family type
#  if args.family_type == "singleton":
#    if len(samples) != 1: fail("Family type was specified as a singleton, but the ped file contains multiple samples")
#  elif args.family_type == "duo":
#    if len(samples) != 2: fail("Family type was specified as a duo, but the ped file doesn't contain two samples")
#  elif args.family_type == "trio":
#    if len(samples) != 3: fail("Family type was specified as a trio, but the ped file doesn't contain three samples")
#  elif args.family_type == "quad":
#    if len(samples) != 4: fail("Family type was specified as a quad, but the ped file doesn't contain three samples")

  # If multiple samples are affected, throw a warning
  if noAffected != 1: fail("Cannot determine proband")

  # Identify the mother and father of the proband
  mother = samples[proband]["mother"]
  father = samples[proband]["father"]
  if samples[proband]["mother"]: samples[mother]["relationship"] = "Mother"
  if samples[proband]["father"]: samples[father]["relationship"] = "Father"

  # Identify siblings
  for sample in samples:
    if samples[sample]["mother"] == mother and samples[sample]["father"] and not samples[sample]["isAffected"]: samples[sample]["relationship"] = "Sibling"
    if not samples[sample]["relationship"]: fail("Sample " + str(sample) + " has an unknown relationship to the proband")

# Get HPO / gene associations
def parseHpoGene(args):
  global hpoGene

  # Open the file containing HPO:gene associations
  hpoFile = open(args.hpo, 'r')
  while True:
    line = hpoFile.readline().rstrip()
    if not line: break

    # Ignore the header
    if line.startswith("#"): continue
    gene = line.split("\t")[1]
    hpo  = line.split("\t")[2]

    # Store the info
    if gene not in hpoGene: hpoGene[gene] = [hpo]
    else: hpoGene[gene].append(hpo)

  # Close the file
  hpoFile.close()

# Parse vcf file
def parseVcf(args):
  global hpoGene
  global patientTerms
  global samples
  global proband
  global mother
  global father
  global weights
  global clinvarSig
  global moiScores
  global consequences
  global spliceAiCutoff

  # Open output files
  hpoOut      = open('hpo.tsv', 'w')
  priorityOut = open('priority.tsv', 'w')
  scoreOut    = open('score.tsv', 'w')
  plotsOut    = open('plots.tsv', 'w')

  # Write the header lines to the output files
  headerBase = "CHROM\tSTART\tEND\tREF\tALT"
  print(headerBase, "hpo_overlaps_fd47afb7", "hpo_terms_2bb3b6d5", sep = "\t", file = hpoOut)
  print(headerBase, "variant_prioritization_b851ec5c", sep = "\t", file = priorityOut)
  print(headerBase, "variant_score_96c63960", sep = "\t", file = scoreOut)
  print("gnomad2AfScore", "gnomad2homaltScore", "gnomad3AfScore", "gnomad3homaltScore", "pLIScore", "hpoScore", "clinvarScore", "moiScores", "impactScore", "variantScore", sep = "\t", file = plotsOut)

  # Read the standard input
  for line in sys.stdin:

    # Skip header lines
    if line.startswith("#"): continue
    line = line.rstrip()

    # Split the line up
    fields = line.split("\t")

    # Get the genotypes
    isSample = True
    currentSample = False
    for info in fields[14:]:

      # If isSample is true, then the bcftools query output is a sample name...
      if isSample:
        if info not in samples: fail(info + " is not a recognised sample")
        isSample = False
        currentSample = info

      # Otherwise it is the genotype associated with the sample name just read
      else:
        isSample = True
        samples[currentSample]["genotype"] = info

    # Get the mode of inheritance
    probandGen = samples[proband]['genotype']
    motherGen  = samples[mother]['genotype']
    fatherGen  = samples[father]['genotype']
    if (probandGen == "0/1" or probandGen == "1/0") and motherGen == "0/0" and fatherGen == "0/0": moi = "denovo"
    elif probandGen == "1/1" and (motherGen == "0/1" or motherGen == "1/0") and (fatherGen == "0/1" or fatherGen == "1/0"): moi == "recessive"
    elif (probandGen == "0/1" or probandGen == "1/0") and (motherGen == "0/1" or motherGen == "1/0" or motherGen == "1/1") and fatherGen == "0/0": moi = "dominant_mother"
    elif (probandGen == "0/1" or probandGen == "1/0") and (fatherGen == "0/1" or fatherGen == "1/0" or fatherGen == "1/1") and motherGen == "0/0": moi = "dominant_father"
    elif probandGen == "1/1" and motherGen == "0/0" and fatherGen == "0/0": moi = "hom_denovo"
    elif (probandGen == "0/1" or probandGen == "1/0") and (motherGen == "0/1" or motherGen == "1/0") and (fatherGen == "0/1" or fatherGen == "1/0"): moi = "skip"
    elif (probandGen == "0/1" or probandGen == "1/0") and (motherGen == "0/1" or motherGen == "1/0") and fatherGen == "1/1": moi = "skip"
    elif (probandGen == "0/1" or probandGen == "1/0") and (fatherGen == "0/1" or fatherGen == "1/0") and motherGen == "1/1": moi = "skip"
    elif (probandGen == "0/1" or probandGen == "1/0") and motherGen == "1/1" and fatherGen == "1/1": moi = "skip"
    elif probandGen == "1/1" and (motherGen == "0/1" or motherGen == "1/0" or motherGen == "1/1") and (fatherGen == "0/1" or fatherGen == "1/0" or fatherGen == "1/1"): moi = "skip"
    elif probandGen == "1/1" and motherGen == "1/1" and fatherGen == "1/1": moi = "skip"
    elif motherGen == "./." or fatherGen == "./.": moi = "unknown"
    elif probandGen == "1/1" and (motherGen == "0/1" or motherGen == "1/0" or motherGen == "1/1") and fatherGen == "0/0": moi = "unlikely"
    elif probandGen == "1/1" and (fatherGen == "0/1" or fatherGen == "1/0" or fatherGen == "1/1") and motherGen == "0/0": moi = "unlikely"
    else: fail("Unknown moi, proband: " + probandGen + ", mother: " + motherGen + ", father: " + fatherGen)

    # Get the gnomAD annotations. If these are '.', the variants are not in gnomAD so should be set to 0, i.e. they are rare
    gnomad2Af     = float(fields[6]) if fields[6] != '.' else 0.
    gnomad2homalt = float(fields[7]) if fields[7] != '.' else 0.
    gnomad3Af     = float(fields[8]) if fields[8] != '.' else 0.
    gnomad3homalt = float(fields[9]) if fields[9] != '.' else 0.

    # If revel score has multiple values, take the max. The REVEL score will be user to modify the score from the consequence
    # (e.g. missense score = missense weight * REVEL). So a REVEL score of 1 will give the impact score the full weight. A
    # very low REVEL score will give a very low impact score to the variant.
    if ',' in fields[10]: revel = float(max(fields[10].split(','))) if fields[10] != '.' else 0.0001
    else: revel = float(fields[10]) if fields[10] != '.' else 0.0001

    # SpliceAI need to be broken up. The scores are fields 2-5, and we want the max
    spliceAi = float( max(fields[11].split("|")[2:5]) ) if fields[11] != '.' else float( -1. )

    # pLI sometimes has multiple values. Takes the max if it does
    if ',' in fields[12]: pLI = float(max(fields[12].split(','))) if fields[12] != '.' else -1
    else: pLI = float(fields[12]) if fields[12] != '.' else -1

    # The ClinVar score is defined by the clinvarSig array
    clinvarScore = clinvarSig[fields[13]] if fields[13] in clinvarSig else fail("Unknown Clinvar significance: " + fields[13])

    # If there is no BCSQ in the record, a '.' will be supplied and these should be ignored
    bcsq  = fields[5]
    genes = []
    impactScore = 0.
    if bcsq != '.':
      finalConsequence = False
      if ',' in bcsq:
        for anno in bcsq.split(","):

          # Get all the genes
          gene = anno.split("|")[1]
          if gene not in genes: genes.append(gene)

          # Get the consequence
          consequence = anno.split("|")[0]
          if "&" in consequence:
            for a in consequence.split("&"):
              impactScore = consequences[a] if consequences[a] > impactScore else impactScore
              finalConsequence = a if consequences[a] > impactScore else finalConsequence
          else:
            impactScore = consequences[consequence]
            finalConsequence = consequence
      else:
        genes       = [bcsq.split("|")[1]]
        consequence = bcsq.split("|")[0]
        if "&" in consequence:
          for a in consequence.split("&"):
            impactScore = consequences[a] if consequences[a] > impactScore else impactScore
            finalConsequence = a if consequences[a] > impactScore else finalConsequence
        else:
          impactScore = consequences[consequence]
          finalConsequence = consequence

      # If this variant is a missense variant, modify the impact score based on the REVEL score
      if finalConsequence == "missense": impactScore = impactScore * revel
      elif finalConsequence == "splice_acceptor" or finalConsequence == "splice_donor": impactScore = impactScore * ( max(spliceAi - spliceAiCutoff, 0.) / (1. - spliceAiCutoff) )

    # Check if any of the genes have an association with an HPO term association with the patient
    hasAssociation    = False
    associatedTerms   = []
    noAssociatedTerms = 0
    for gene in genes:
      if gene in hpoGene:
        for hpoTerm in hpoGene[gene]:
          if hpoTerm in patientTerms:
            hasAssociation = True
            if hpoTerm not in associatedTerms:
              noAssociatedTerms += 1
              associatedTerms.append(hpoTerm)

    # Calculate a score for the variant
    gnomad2AfScore     = weights['gnomad2Af'] * (1. - 10. * gnomad2Af) if gnomad2Af != -1 else 0.
    gnomad2homaltScore = weights['gnomad2homalt'] * max((1 - gnomad2homalt / 10.), 0.) if gnomad2homalt != -1 else 0.
    gnomad3AfScore     = weights['gnomad3Af'] * (1. - 10. * gnomad3Af) if gnomad3Af != -1 else 0.
    gnomad3homaltScore = weights['gnomad3homalt'] * max((1 - gnomad3homalt / 10.), 0.) if gnomad3homalt != -1 else 0.
    pLIScore           = weights['pLI'] * pLI if pLI != -1 else 0.
    hpoScore           = weights['hpo'] * (noAssociatedTerms / len(patientTerms))
    variantScore       = gnomad2AfScore + gnomad2homaltScore + gnomad3AfScore + gnomad3homaltScore
    variantScore      += pLIScore
    variantScore      += weights['clinvar'] * clinvarScore
    variantScore      += weights['moi'] * moiScores[moi]
    variantScore      += weights['consequence'] * impactScore
    if hasAssociation: variantScore += weights['hpoBase'] + hpoScore

    # Get the basic information required for upload to Mosaic as an annotation
    chrom = line.split("\t")[0]
    if chrom.startswith("chr"): chrom = chrom.strip("chr")
    start = line.split("\t")[1]
    end   = int(line.split("\t")[2]) + 1
    ref   = line.split("\t")[3]
    alt   = line.split("\t")[4]

    # If there is an association with a patient HPO term, write out the tsv for upload to Mosaic
    if hasAssociation: print(chrom, start, end, ref, alt, noAssociatedTerms, ",".join(associatedTerms), sep = "\t", file = hpoOut)

    # Print the variant score to file, so it can be uploaded as an annotation
    print(chrom, start, end, ref, alt, variantScore, sep = "\t", file = scoreOut)
    print(gnomad2Af, gnomad2homalt, gnomad3Af, gnomad3homalt, pLI, noAssociatedTerms, clinvarScore, moiScores[moi], impactScore, variantScore, sep = "\t", file = plotsOut)

    # If the variant score is higher than the lowest score in the priority list, remove the lowest score and add this variant
    lowestScore = min(priorityList.keys())
    if float(variantScore) > float(lowestScore):
      del priorityList[lowestScore]
      priorityList[variantScore] = str(chrom) + "\t" + str(start) + "\t" + str(end) + "\t" + str(ref) + "\t" + str(alt) + "\t" + str(float("{:.3f}".format(variantScore)))

  # Write the prioritized variants to file:
  for i in sorted(priorityList): print(priorityList[i], file = priorityOut)

  # Close output files
  hpoOut.close()
  priorityOut.close()
  scoreOut.close()
  plotsOut.close()

# If the script fails, provide an error message and exit
def fail(message):
  print(message, sep = "")
  exit(1)

# Initialise global variables

# Version
version = "0.0.1"
date    = str(date.today())

# Store information related to Mosaic
mosaicConfig = {}
mosaicInfo   = {}

# The working directory where all created files are kept
workingDir = os.getcwd() + "/calypso_v" + version + "r"

# Store info on allowed values
allowedReferences = ['37', '38']

# Store the patients HPO terms
patientTerms = []

# Store the hpo:gene associations
hpoGene = {}

# Mosaic samples
samples = {}
proband = False
mother  = False
father  = False

# Define the cutoff for SpliceAI. Scores below this will eliminate the contribution from the impact score
spliceAiCutoff = 0.7

# Weights for calculating variant scores
weights                  = {}
weights['consequence']   = 3.
weights['clinvar']       = 4.
weights['gnomad2Af']     = 0.2
weights['gnomad2homalt'] = 0.5
weights['gnomad3Af']     = 0.2
weights['gnomad3homalt'] = 0.5
weights['pLI']           = 0.5
weights['hpo']           = 1.
weights['hpoBase']       = 3.
weights['moi']           = 1.

# Store scores for ClinVar significance
clinvarSig = {}
clinvarSig['Pathogenic'] = 1.
clinvarSig['Pathogenic/Likely_pathogenic'] = 0.95
clinvarSig['Likely_pathogenic'] = 0.9
clinvarSig['Conflicting_interpretations_of_pathogenicity'] = 0.4
clinvarSig['Uncertain_significance'] = 0.2
clinvarSig['Likely_benign'] = -0.2
clinvarSig['Benign/Likely_benign'] = -0.4
clinvarSig['Benign'] = -1.
clinvarSig['.'] = 0.

# Scores for moi
moiScores = {}
moiScores['denovo'] = 1
moiScores['hom_denovo'] = 0.8
moiScores['recessive'] = 0.2
moiScores['dominant_mother'] = -1.
moiScores['dominant_father'] = -1.
moiScores['unknown'] = -0.2
moiScores['skip'] = -0.2
moiScores['unlikely'] = 0.

# Impact scores based on consequence
consequences = {}
consequences['frameshift']        = 1.0 # HIGH
consequences['stop_gained']       = 1.0 # HIGH
consequences['stop_lost']         = 1.0 # HIGH
consequences['start_lost']        = 1.0 # HIGH
consequences['splice_donor']      = 1.0 # HIGH
consequences['splice_acceptor']   = 1.0 # HIGH

consequences['missense']          = 0.6 # MODERATE
consequences['inframe_insertion'] = 0.6 # MODERATE
consequences['inframe_deletion']  = 0.6 # MODERATE

consequences['non_coding']        = 0. # MODIFIER
consequences['coding_sequence']   = 0. # MODIFIER
consequences['intron']            = 0. # MODIFIER
consequences['3_prime_utr']       = 0. # MODIFIER
consequences['5_prime_utr']       = 0. # MODIFIER

consequences['synonymous']        = 0. # LOW
consequences['splice_region']     = 0. # LOW
consequences['stop_retained']     = 0. # LOW
consequences['start_retained']    = 0. # LOW

# Store the top 50 scored variants
priorityList = {}
for i in range(1, 51, 1): priorityList[i * 0.0001] = False

if __name__ == "__main__":
  main()

import argparse
import json
import os
import subprocess
import sys

from datetime import date
from os.path import exists
from pprint import pprint
from subprocess import Popen
from sys import path

import read_resource_jsons as read_resources

def main():
  print()
  print('Starting the Calypso pipeline')

  # Parse the command line
  args = parse_command_line()

  # Check the supplied parameters are as expected, then expand the path to enable scripts and api commands
  # to be accessed for Calypso
  root_path = os.path.dirname(__file__)
  checkArguments(args, root_path)

  # Import the api client
  path.append(args.api_client)
  from mosaic import Mosaic, Project, Store
  api_store  = Store(config_file = args.client_config)
  api_mosaic = Mosaic(config_file = args.client_config)

  # Open an api client project object for the defined project
  project = api_mosaic.get_project(args.project_id)

  # Get the project settings and define the reference genome
  project_settings = project.get_project_settings()
  reference = project_settings['reference']
  if reference not in allowed_references:
    fail('The specified reference (' + str(reference) + ') is not recognised. Allowed values are: ' + str(', '.join(allowed_references)))
  print('Using the reference: ', reference, sep = '')

  # Read the resources json file to identify all the resources that will be used in the annotation pipeline
  if args.resource_json:
    resource_info = check_resources(reference, args.data_directory, args.tools_directory, args.resource_json)
  else:
    args.resource_json = str(args.data_directory) + 'resources_' + str(reference) + '.json'
    resource_info = check_resources(reference, args.data_directory, args.tools_directory, args.resource_json)
  print('Using resources file: ', args.resource_json, sep = '')
  resource_info = read_resources.read_resources(reference, root_path, resource_info, args.no_vep)

  # If threads were not set, default to 4
  args.threads = 4 if not args.threads else args.threads

  # Define the tools to be used by Calypso
  resource_info = calypso_tools(resource_info)

  # Define paths to be used by Calypso
  set_working_directory(resource_info['version'])

  # Get the samples in the Mosaic projects and determine the id of the proband, and the family structure
  mosaic_samples = {}
  mosaic_sample_ids = {}
  proband = False
  for sample in project.get_samples():
    relation = False
    for attribute in project.get_attributes_for_sample(sample['id']):
      if attribute['name'] == 'Relation':
        for value in attribute['values']:
          if value['sample_id'] == sample['id']:
            relation = value['value']
            break
        break

    # Fail if this sample has no Relation set, and store the sample name if this is the proband
    if not relation:
      fail('The Relation attribute is not set for sample ' + sample['name'])
    elif relation == 'Proband':
      proband = sample['name']

    # Add this sample to the list of samples
    mosaic_samples[sample['name']] = {'id': sample['id'], 'relation': relation}
    mosaic_sample_ids[sample['id']] = sample['name']

  # Fail if no proband is defined
  if not proband:
    fail('Could not find a proband for this project')

  # Get the vcf file for each sample. Currently, Calypso only works for projects with a single multi-sample vcf.
  # Determine the name of this vcf file, then determine if this file has chromosomes listed as e.g. '1' or 'chr1'.
  # If a vcf was supplied on the command line, use this. This is for occasions where the vcf file has not been
  # linked to samples in Mosaic
  if args.input_vcf:
    vcf = os.path.abspath(args.input_vcf)

    # Check the file exists
    if not exists(vcf):
      fail('Input vcf file does not exist')

    # Check that the index file also exists
    index_file = vcf + '.tbi'
    if not exists(index_file):
      fail('The vcf index file does not exist')

    # Get the vcf header
    header = os.popen(str(resource_info['tools']['bcftools']) + ' view -h ' + str(vcf)).read()
    for line in header.split('\n'):
      if line.startswith('#CHROM'):
        vcf_samples = line.rstrip().split('\t')[9:]
        break

    # Loop over the samples in the vcf file and find the ones in the mosaic_samples list
    for vcf_sample in vcf_samples:

      # This is a hack for UDN where the sample name has the experiment id appended
      vcf_sample_name = vcf_sample
      if args.udn:
        if '-' in vcf_sample:
          vcf_sample = vcf_sample.split('-')[0]
          mosaic_samples[vcf_sample]['vcf_file'] = vcf
          mosaic_samples[vcf_sample]['vcf_sample_name'] = vcf_sample_name

      else:
        if vcf_sample in mosaic_samples:
          mosaic_samples[vcf_sample]['vcf_file'] = vcf
          mosaic_samples[vcf_sample]['vcf_sample_name'] = vcf_sample
      print('  Sample ', vcf_sample, ' appears as "', vcf_sample_name, '" in the header of vcf file: ', vcf, sep = '')

    # Check that all samples have been associated with a vcf file
    samples_with_no_vcf = []
    for sample in mosaic_samples:
      if 'vcf_file' not in mosaic_samples[sample]:
        samples_with_no_vcf.append(sample)
    if samples_with_no_vcf:
      fail('\nERROR: The following samples do not appear in the vcf header:\n  ' + '\n  '.join(samples_with_no_vcf))

  # If the vcf wasn't specified on the command line, get it from the Mosaic project
  else:
    all_vcfs = []
    vcf_files = {}

    # Loop over all the samples
    for sample in mosaic_samples:
      sample_id = mosaic_samples[sample]['id']
      for sample_file in project.get_sample_files(sample_id):
        if sample_file['type'] == 'vcf':
          vcf_files[sample_file['id']] = {'name': sample_file['name'], 'uri': sample_file['uri'], 'vcf_sample_name': sample_file['vcf_sample_name']}
  
      # If there are no vcf files for this sample output a warning, but continue, unless this is the proband
      if len(vcf_files) == 0:
        if mosaic_samples[sample]['relation'] == 'Proband':
          fail('calypso_vcf_files: The proband (' + sample + ') is not associated with any vcf files. The proband must appear in a vcf file (a vcf file can be specified on the command line using the --input_vcf (-i) argument)')
        else:
          print('  WARNING: Sample ' + sample + ', listed as having the relationship ' + mosaic_samples[sample]['relation'] + ' is not associated with a vcf and so will contribute no variants')
          mosaic_samples[sample]['vcf_file'] = False
          mosaic_samples[sample]['vcf_sample_name'] = False
  
      # If there is more than one vcf file for the sample, fail
      elif len(vcf_files) > 1:
        vcf_file_names = []
        for vcf_file_id in vcf_files:
          vcf_file_names.append(vcf_files[vcf_file_id]['name'])
        print()
        fail('ERROR: Multiple vcf files for the samples. Calypso requires a single multi sample vcf. The following vcf files were found:\n  ' + '\n  '.join(vcf_file_names))
      else:
        for vcf_file in vcf_files:
          uri = vcf_files[vcf_file]['uri']
  
          # If the uri begins with 'file://', this needs to be removed
          final_uri = uri[6:] if uri.startswith('file://') else uri
          vcf_name = vcf_files[vcf_file]['vcf_sample_name']
          if final_uri not in all_vcfs:
            all_vcfs.append(final_uri)
        mosaic_samples[sample]['vcf_file'] = final_uri
        mosaic_samples[sample]['vcf_sample_name'] = vcf_name
        print('Sample ', sample, ' appears as "', vcf_name, '" in the header of vcf file: ', final_uri, sep = '')
  
    # If there are multiple VCF files, not all samples are in a single joint-called vcf. This case
    # is not yet handled:
    if len(all_vcfs) > 1:
      fail('ERROR: Multiple vcf files for the samples. Calypso requires a single multi sample vcf. The following vcf files were found:\n  ' + '\n  '.join(all_vcfs))

    # Store the name of the vcf file
    vcf = mosaic_samples[list(mosaic_samples.keys())[0]]['vcf_file']

  # If the ped file has not been defined, generate a bed file from pedigree information in Mosaic
  if not args.ped:

    # Open a ped file in the working directory
    ped_filename = str(working_directory) + str(proband) + '.ped'
    ped_file = open(ped_filename, 'w')
    print('#Family_id', 'individual_id', 'paternal_id', 'maternal_id', 'sex', 'affected_status', sep = '\t', file = ped_file)
 
    # Parse the pedigree associated with the proband
    for pedigree_line in project.get_pedigree(mosaic_samples[proband]['id']):
      kindred_id = pedigree_line['pedigree']['kindred_id']
      sample_id = mosaic_sample_ids[pedigree_line['pedigree']['sample_id']]
      paternal_id = mosaic_sample_ids[pedigree_line['pedigree']['paternal_id']] if pedigree_line['pedigree']['paternal_id'] else 0
      maternal_id = mosaic_sample_ids[pedigree_line['pedigree']['maternal_id']] if pedigree_line['pedigree']['maternal_id'] else 0
      sex = pedigree_line['pedigree']['sex']
      affection_status = pedigree_line['pedigree']['affection_status']

      # The ped file must use the names as they appear in the vcf file, so use the sample_id to get the vcf_sample_name
      sample_vcf_id = mosaic_samples[sample_id]['vcf_sample_name']
      paternal_vcf_id = mosaic_samples[paternal_id]['vcf_sample_name'] if paternal_id != 0 else 0
      maternal_vcf_id = mosaic_samples[maternal_id]['vcf_sample_name'] if maternal_id != 0 else 0
      print(kindred_id, sample_vcf_id, paternal_vcf_id, maternal_vcf_id, sex, affection_status, sep = '\t', file = ped_file)

    # Close the ped file
    ped_file.close()

    # Set args.ped to the created ped file
    args.ped = ped_filename

  # Check that the reference genome associated with the vcf file matches the project reference
  header = os.popen(str(resource_info['tools']['bcftools']) + ' view -h ' + str(vcf)).read()
  for line in header.split('\n'):
    if line.startswith('##contig=<ID=chr1,'):
      chr1_line = line.rstrip()
      break
    elif line.startswith('##contig=<ID=1,'):
      chr1_line = line.rstrip()
      break
    else:
      chr1_line = False

  # Check the length of chr1
  is_correct = False
  if 'length' not in chr1_line or not chr1_line:
    print('Could not verify the reference genome for the vcf file')
  else:
    length = (line.rstrip().split('=')[3]).rstrip('>')
    if reference == 'GRCh38' and str(length).startswith('248956422'):
      is_correct = True
    elif reference == 'GRCh37' and str(length).startswith('249250621'):
      is_correct = True

  if is_correct:
    print('The vcf file matches the reference genome of the project (', str(reference), ')', sep = '')
  else:
    print('The vcf file DOES NOT match the reference genome of the project (', str(reference), ')', sep = '')
    print('Please verify the reference genome in Mosaic and ensure these match')
    exit(1)

  # Check that we have read access to the vcf file
  if not os.access(vcf, os.R_OK):
    fail('\nFAILED: Calypso does not have access to the vcf file which is required to complete the pipeline')
##########
##########
########## DELETE WHEN HANDLED S3 FILES
##########
##########
  if vcf.startswith('s3'):
    vcf = vcf.split('/')[-1]

  # Determine the format of the chromosomes in the vcf file
  header = os.popen(str(resource_info['tools']['bcftools']) + ' view -h ' + str(vcf)).read()
  for line in header.split('\n'):
    if line.startswith('##contig'):
      chr_id = (line.split('=')[2]).split(',')[0]
      chr_format = True if chr_id.startswith('chr') else False
      break

  # Loop over all the samples. For each sample, get the header final header line from its associated vcf file (which
  # contains the list of all samples in the vcf) and determine at which position this sample appears
  for sample in mosaic_samples:
    vcf = mosaic_samples[sample]['vcf_file']
    vcf_sample_name = mosaic_samples[sample]['vcf_sample_name']

    # Only process samples associated with vcf files
    if vcf:

####################
####################
#################### HACK FOR UDN S3 files
#################### DELETE
####################
####################
      if vcf.startswith('s3://udn-joint-analysis/wgs/joint-by-family/'):
        vcf = vcf.replace('s3://udn-joint-analysis/wgs/joint-by-family/', './')

      # Get the final header line
      data = os.popen(str(resource_info['tools']['bcftools']) + ' view -h ' + str(vcf)).readlines()[-1].rstrip().split('\t')

      # Read through the samples and define the sample order
      for index, vcf_sample in enumerate(data[9:]):
        if vcf_sample == vcf_sample_name:
          mosaic_samples[sample]['vcf_position'] = index

  # Check that every sample in samples has vcf_position set
  for sample in mosaic_samples:
    if vcf and 'vcf_position' not in mosaic_samples[sample]:
      fail('Sample ' + str(sample) + ' is listed in the ped file, but does not appear in the vcf header')

  # If the HPO terms are not specified on the command line, get them from the project. If they are specified on the command
  # line, use these
  if not args.hpo:
    hpo_terms = []
    for hpo_term in project.get_sample_hpo_terms(mosaic_samples[proband]['id']):
      if hpo_term['hpo_id'] not in hpo_terms:
        hpo_terms.append(hpo_term['hpo_id'])
    args.hpo = ','.join(hpo_terms)
  if args.hpo:
    print('Using the HPO terms: ', args.hpo, sep = '')
  else:
    print('Using the HPO terms: No terrms available')

  # Get all the project attributes in the target Mosaic project

  # Get all the project attributes from Mosaic for the project
  project_attributes = {}
  for attribute in project.get_project_attributes():
    for attribute_project in attribute['values']:
      value = False
      if int(attribute_project['project_id']) == int(args.project_id):
        value = attribute_project['value']
    project_attributes[attribute['uid']] = {'name': attribute['name'], 'id': attribute['id'], 'value': value}

  # Get all of the public attributes that are available for import
  public_attributes = {}
  for attribute in api_mosaic.get_public_project_attributes():
    public_attributes[attribute['uid']] = {'name': attribute['name'], 'id': attribute['id']}

  # Read the Mosaic json and validate its contents
  if not args.mosaic_json:
    args.mosaic_json = args.data_directory + 'resources_mosaic_' + reference + '.json'
  print('Using Mosaic resources file: ', args.mosaic_json, sep = '')
  mosaic_info = read_resources.read_mosaic_json(args.mosaic_json, reference)

  # Check that the Mosaic resource file does not include resources that are not defined in the resources json
  check_mosaic_resources(mosaic_info, resource_info)

  # Check that a variant filters file has been supplied. The validity of the file will be checked when it is parsed.
  # For now, it is sufficient that a file exists. If no file is supplied, use the one that should exist in the
  # Calypso data directory
  if not args.variant_filters:
    if mosaic_info['variantFilters']:
      args.variant_filters = str(args.data_directory) + str(mosaic_info['variantFilters'])
    else:
      fail('No json file describing the variant filter is specified. This can be specified in the mosaic resources file or on the command line')
  if not os.path.exists(args.variant_filters):
    fail('ERROR: The json file describing the preset variant filters does not exist (' + str(args.variant_filters) + ')')

  # Get all the project annotations already in the Mosaic project
  project_annotations = {}
  for annotation in project.get_variant_annotations():
    project_annotations[annotation['uid']] = {'id': annotation['id'], 'name': annotation['name'], 'type': annotation['value_type']}

  # Get all available public annotations
  public_annotations = {}
  for annotation in project.get_variant_annotations_to_import():
    public_annotations[annotation['uid']] = {'name': annotation['name'], 'id': annotation['id'], 'type': annotation['value_type']}

  # Create all the required private annotations
  private_annotations = {}

  # Get the names (not the uids) of all the annotations in the project
  existing_annotations = []
  for annotation in project_annotations:
    existing_annotations.append(project_annotations[annotation]['name'])

  # Loop over each resource in the Mosaic resources file
  for resource in mosaic_info['resources']:
    annotations = {}

    # Loop over all annotations for this resource
    if mosaic_info['resources'][resource]['annotation_type'] == 'private':

      # The project id for these annotations is the target project as these annotations will be created there. Update
      # mosaic_info to reflect this.
      mosaic_info['resources'][resource]['project_id'] = args.project_id

      # Loop over all annotations for this resource
      for annotation in mosaic_info['resources'][resource]['annotations']:
        value_type = mosaic_info['resources'][resource]['annotations'][annotation]['type']

        # If the annotation has "isSample" set, there needs to be an annotation for every sample in the vcf file for the 
        # project (some projects include samples that were not sequenced and so these will not have annotations). Append the
        # the samples relationship to the proband to the annotation name. 
        if mosaic_info['resources'][resource]['annotations'][annotation]['isSample']:
          for sample in mosaic_samples:

            # Include the sample name in the annotation, otherwise non-unique names could result. For example, if the proband
            # has 2 brothers, the generated annotation name for the 2 brothers will be identical
            if mosaic_samples[sample]['vcf_sample_name']:
              annotations[str(annotation) + ' ' + str(mosaic_samples[sample]['relation']) + ' ' + str(sample)] = value_type
        else:
          annotations[annotation] = value_type

      # Remove this annotation from the resources dictionary. It will be replaced below with the new annotations with the
      # corresponding uids
      mosaic_info['resources'][resource]['annotations'] = {}

    # Loop over the annotations to be created, check that there doesn't already exist an annotation of that name in the
    # project, and if not, create a new, private annotation
    for annotation in annotations:
      value_type = annotations[annotation]

      # If the annotation already exists, get the uid of the annotation, otherwise create the annotation.
      if annotation in existing_annotations:
        for existing_annotation in project_annotations:
          if str(project_annotations[existing_annotation]['name']) == annotation:
            private_annotations[existing_annotation] = {'name': annotation, 'id': project_annotations[existing_annotation]['id']}
            mosaic_info['resources'][resource]['annotations'][annotation] = {'uid': existing_annotation, 'type': value_type, 'id': project_annotations[existing_annotation]['id']}
            break
      else:
        data = project.post_variant_annotation(name = annotation, value_type = value_type, privacy_level = 'private')
        private_annotations[data['uid']] = {'name': annotation, 'id': data['id']}
        mosaic_info['resources'][resource]['annotations'][annotation] = {'uid': data['uid'], 'type': value_type, 'id': data['id']}

        # Add the created private annotation to the project_annotations dictionary
        project_annotations[data['uid']] = {'id': data['id'], 'name': annotation, 'type': value_type}

  # Loop over all of the public annotations defined in the Mosaic resources json file and check that they exist
  # in Mosaic and can be imported
  for resource in mosaic_info['resources']:
    if mosaic_info['resources'][resource]['annotation_type'] == 'public':

      # Loop over all the annotations required for this resource
      for annotation in mosaic_info['resources'][resource]['annotations']:
        uid = mosaic_info['resources'][resource]['annotations'][annotation]['uid']

        # If the annotation uid is not in the publicAnnotations dictionary, then this annotation does not exist for
        # import. Fail and indicate that the resource file may contain errors
        if uid not in public_annotations:
          fail('Annotation with uid ' + str(uid) + ', for resource ' + str(resource) + ', is not available for import. Check that the uid in the Mosaic resources file is correct')

  # Build the toml file for vcfanno. This defines all of the annotations that vcfanno needs to use
  # with the path to the required files. The build the lua file that vcfanno uses to grab lua functions
  # for postannotation
  toml_filename = buildToml(working_directory, resource_info)
  lua_filename = generate_lua_file(working_directory)

  # Generate the bash script to run the annotation pipeline
  bash_filename, bash_file = open_bash_script(working_directory)
  filtered_vcf = bash_resources(resource_info, working_directory, bash_file, vcf, chr_format, args.ped, lua_filename, toml_filename)
  annotate_vcf(resource_info, bash_file, chr_format, args.threads, mosaic_samples)
  filter_vcf(bash_file, mosaic_samples, proband, resource_info, args.threads)

  # Process the filtered vcf file to extract the annotations to be uploaded to Mosaic
  print('# Generate tsv files to upload annotations to Mosaic', file = bash_file)
  print('GENERATE_TSV=', root_path, '/generate_annotation_tsv.py', sep = '', file = bash_file)
  print('GENERATE_HGVS_TSV=', root_path, '/generate_hgvs_annotation_tsv.py', sep = '', file = bash_file)
  print('GENERATE_VARIANT_QUALITY_TSV=', root_path, '/generate_variant_quality_tsv.py', sep = '', file = bash_file)
  print('MOSAIC_JSON=', args.mosaic_json, sep = '', file = bash_file)

  # Loop over all the resources to be uploaded to Mosaic
  tsv_files = []
  for resource in mosaic_info['resources']:

    # The HGVS annotations require additional processing to identify the codes corresponding to the MANE transcript
    # and strip the transcript information from the annotation
    if resource == 'HGVS':
      generate_hgvs_tsv(bash_file, reference)
      tsv_files.append('hgvs.tsv')

    # Extract the Variant Quality scores after getting the uid of this annotation
    elif resource == 'variant_quality':
      for annotation in private_annotations:
        if private_annotations[annotation]['name'] == 'Variant Quality':
          uid = annotation
          break
      generate_variant_quality_tsv(bash_file, uid)
      tsv_files.append('variant_quality.tsv')

    # Annotations to ignore
    elif resource == 'GQ':
      pass

    # Remaining resources
    else:
      tsv_files.append(generate_annotation_tsv(bash_file, resource, reference))

  # Close the bash script and make it executable
  print('echo "Calypso pipeline version ', version, ' completed successfully"', sep = '', file = bash_file)
  bash_file.close()
  make_executable = os.popen('chmod +x ' + bash_filename).read()

  # Prepare the target project. This includes:
  # 1. Remove any unnecessary annotations from the project
  # 2. Import public annotations into the project
  # 3. Set the default annotations for the project
  # 4. Set up required variant filters
  # 5. Add and set project attributes
  if 'remove' in mosaic_info:

    # Loop over the annotations to remove, get the annotation id, and remove them from the project
    for uid in mosaic_info['remove']:
      if uid in project_annotations:
        project.delete_variant_annotation(project_annotations[uid]['id'])

  # Loop over all of the resources in the Mosaic resources file and get all of the annotation uids to import.
  # Only public attributes can be imported
  for resource in mosaic_info['resources']:
    if mosaic_info['resources'][resource]['annotation_type'] == 'public':

      # Loop over all the annotations required for this resource
      for annotation in mosaic_info['resources'][resource]['annotations']:
        uid = mosaic_info['resources'][resource]['annotations'][annotation]['uid']

        # Check if the annotation is already in the target project. If not, import the annotation and update the
        # projectAnnotations dictionary
        if uid not in project_annotations:
          project.post_import_annotation(public_annotations[uid]['id'])
          project_annotations[uid] = {'id': public_annotations[uid]['id'], 'name': public_annotations[uid]['name'], 'type': public_annotations[uid]['type']}

  # Import any annotations specified in the resource file
  if 'import_annotations' in mosaic_info:
    for uid in mosaic_info['import_annotations']:
      if uid not in project_annotations:
        if uid not in public_annotations:
          print('WARNING: resource file requests annotation with uid ', uid, ' be imported, but this is not a public attribut and so is not imported', sep = '')
        else:
          project.post_import_annotation(public_annotations[uid]['id'])
          project_annotations[uid] = {'id': public_annotations[uid]['id'], 'name': public_annotations[uid]['name'], 'type': public_annotations[uid]['type']}

  # Set the projects default annotations. The default annotations are what a new user will see in the table by default
  # and all users will have the option to reset the annotations table to show only the default annotations
  default_annotation_ids = []

  # Get the ids of all the default annotations
  for annotation in mosaic_info['defaultAnnotations']:
    if annotation in project_annotations:
      default_annotation_ids.append(project_annotations[annotation]['id'])
    elif annotation in public_annotations:
      project.post_import_annotation(public_annotations[annotation]['id'])
      project_annotations[uid] = {'id': public_annotations[annotation]['id'], 'name': public_annotations[annotation]['name'], 'type': public_annotations[annotation]['type']}
      default_annotation_ids.append(public_annotations[annotation]['id'])

    # Check if the annotation is a private annotation. In this case, the supplied annotation will not be the uid as
    # this is not known, but will be the annotation name
    else:
      private_annotation_id = False
      for private_annotation in private_annotations:
        if str(private_annotations[private_annotation]['name']) == str(annotation):
          private_annotation_id = private_annotations[private_annotation]['id']
      if private_annotation_id:
        default_annotation_ids.append(private_annotation_id)
      else:
        fail('Default annotation ' + str(annotation) + ' has not been imported or created in the project')

  # Set the default annotations in the Mosaic project. This requires the annotation version ids, not the
  # annotation ids. Use the latest version. If latest isn't set, use default
  annotation_version_ids = []
  for annotation in project.get_variant_annotations():
    if annotation['id'] in default_annotation_ids:
      latest_annotation_version_id = annotation['latest_annotation_version_id']
      if not latest_annotation_version_id:
        for annotation_version in annotation['annotation_versions']:
          if annotation_version['version'] == 'default':
            annotation_version_ids.append(annotation_version['id'])
      else:
        annotation_version_ids.append(latest_annotation_version_id)
  project.put_project_settings(selected_variant_annotation_version_ids = annotation_version_ids)

  # Generate scripts to upload filtered variants to Mosaic
  upload_filename = working_directory + '02_calypso_upload_variants.sh'
  try:
    upload_file = open(upload_filename, 'w')
  except:
    fail('Could not open ' + str(upload_filename) + ' to write to')

  # Write the command to file to upload the filtered variants
  print('# Upload variants to Mosaic', file = upload_file)
  print('API_CLIENT=', args.api_client, sep = '', file = upload_file)
  print('CONFIG=', args.client_config, sep = '', file = upload_file)
  print('VCF=', filtered_vcf, sep = '', file = upload_file)
  print(file = upload_file)
  print('python3 $API_CLIENT/variants/upload_variants.py ', end = '', file = upload_file)
  print('-a $API_CLIENT ', end = '', file = upload_file)
  print('-c $CONFIG ', end = '', file = upload_file)
  print('-p ', str(args.project_id) + ' ', sep = '', end = '', file = upload_file)
  print('-m "no-validation" ', sep = '', end = '', file = upload_file)
  print('-v $VCF ', file = upload_file)
  print(file = upload_file)

  # Close the file
  upload_file.close()

  # Make the annotation script executable
  make_executable = os.popen('chmod +x ' + str(upload_filename)).read()

  # Open a file with the commands to upload all the annotations
  upload_filename = working_directory + '03_calypso_upload_annotations.sh'
  try:
    upload_file = open(upload_filename, 'w')
  except:
    fail('Could not open ' + str(upload_filename) + ' to write to')

  if reference == 'GRCh37':
    annotation_project_id = api_store.get('Project ids', 'annotations_grch37')
  elif reference == 'GRCh38':
    annotation_project_id = api_store.get('Project ids', 'annotations_grch38')
  print('API_CLIENT=', args.api_client, sep = '', file = upload_file)
  print('CONFIG=', args.client_config, sep = '', file = upload_file)
  print('UPLOAD_SCRIPT=$API_CLIENT/variant_annotations/upload_annotations.py', sep = '', file = upload_file)
  for tsv in tsv_files:
    print(file = upload_file)
    print('# Upload ', tsv, ' annotations', sep = '', file = upload_file)
    print('TSV=', working_directory, tsv, sep = '', file = upload_file)
    print('python3 $UPLOAD_SCRIPT ', end = '', file = upload_file)
    print('-a $API_CLIENT ', end = '', file = upload_file)
    print('-c $CONFIG ', end = '', file = upload_file)
    print('-p ', annotation_project_id, ' ', sep = '', end = '', file = upload_file)
    print('-t $TSV', sep = '', file = upload_file)

  # Close the upload file
  upload_file.close()

  # Make the annotation script executable
  make_executable = os.popen('chmod +x ' + str(upload_filename)).read()

  # Determine all of the variant filters (from the calypso_mosaic_filters.json) that are to be added; remove any filters that already
  # exist with the same name; fill out variant filter details not in the json (e.g. the uids of private annotations created by
  # Calypso); create the filters; and finally update the project settings to put the filters in the correct category and sort order.
  # Note that the filters to be applied depend on the family structure. E.g. de novo filters won't be added to projects without parents
  if not args.variant_filters:
    args.variant_filters = mosaic_info['variant_filters']
  if args.api_client.endswith('/'):
    variant_filter_command = 'python3 ' + args.api_client + 'project_setup/set_variant_filters.py'
  else:
    variant_filter_command = 'python3 ' + args.api_client + '/project_setup/set_variant_filters.py'
  variant_filter_command += ' -a ' + args.api_client
  variant_filter_command += ' -c ' + args.client_config
  variant_filter_command += ' -p ' + args.project_id
  variant_filter_command += ' -f ' + args.variant_filters
  os.system(variant_filter_command)

  # Output a summary file listing the actions undertaken by Calypso with all version histories
  #res.calypsoSummary(working_directory, version, resource_info, reference)
  #print('Calypso pipeline version ', version, ' completed successfully', sep = '')

  # Run Exomiser on the original vcf file and process the exomiser outputs as private annotations. Only
  # run if hpo terms are supplied
  if args.hpo:
    print('Generating exomiser scripts...', end = '')
    application_properties(working_directory, args.tools_directory, reference)
    yaml = generate_yml(working_directory, proband, reference, str(working_directory) + str(filtered_vcf), args.ped, args.hpo)
    exomiser_script_name, exomiser_script = generate_exomiser_script(working_directory, args.tools_directory, yaml)
    print('complete')

    # Check the variant filters json is set
    if not args.exomiser_filters_json:
      if 'exomiser_variant_filters' in mosaic_info:
        args.exomiser_filters_json = str(args.data_directory) + str(mosaic_info['exomiser_variant_filters'])
      else:
        fail('No json file describing the exomiser variant filters is specified. This can be specified in the mosaic resources file or on the command line')
    if not os.path.exists(args.exomiser_filters_json):
      fail('ERROR: The json file describing the preset exomiser variant filters does not exist (' + str(args.exomiser_filters_json) + ')')
  
    # The exomiser script currently uses the filtered vcf, so no new variants will need to be uploaded
    #upload_exomiser_variants(working_directory, args.api_client, args.client_config, args.project_id, proband)
    print('Applying exomiser filters...', end = '')
    exomiser_annotations(root_path, working_directory, args.api_client, args.client_config, args.project_id, proband, args.exomiser_filters_json)
    print('complete')

# Input options
def parse_command_line():
  global version
  parser = argparse.ArgumentParser(description='Process the command line')

  # Define the location of the api_client and the ini config file
  parser.add_argument('--api_client', '-a', required = True, metavar = 'string', help = 'The api_client directory')
  parser.add_argument('--client_config', '-c', required = True, metavar = 'string', help = 'The ini config file for Mosaic')

  # Required arguments
  parser.add_argument('--data_directory', '-d', required = False, metavar = 'string', help = 'The path to the directory where the resources live')
  parser.add_argument('--tools_directory', '-s', required = False, metavar = 'string', help = 'The path to the directory where the tools to use live')
  parser.add_argument('--input_vcf', '-i', required = False, metavar = 'string', help = 'The input vcf file to annotate')
  parser.add_argument('--ped', '-e', required = False, metavar = 'string', help = 'The pedigree file for the family. Not required for singletons')
  parser.add_argument('--project_id', '-p', required = True, metavar = 'string', help = 'The project id that variants will be uploaded to')
  parser.add_argument('--variant_filters', '-f', required = False, metavar = 'string', help = 'The json file describing the variant filters to apply to each project')

  # Exomiser arguments
  parser.add_argument('--exomiser_filters_json', '-x', required = False, metavar = 'string', help = 'The json describing the exomiser filters')

  # Optional pipeline arguments
  parser.add_argument('--reference', '-r', required = False, metavar = 'string', help = 'The reference genome to use. Allowed values: ' + ', '.join(allowed_references))
  parser.add_argument('--resource_json', '-j', required = False, metavar = 'string', help = 'The json file describing the annotation resources')
  parser.add_argument('--no_vep', '-n', required = False, action = "store_false", help = "If set, do not run VEP annotation")
  parser.add_argument('--threads', '-t', required = False, metavar = 'integer', help = 'The number of threads to use')

  # Optional argument to handle HPO terms
  parser.add_argument('--hpo', '-o', required = False, metavar = "string", help = "A comma separate list of hpo ids for the proband")

  # Optional mosaic arguments
  parser.add_argument('--mosaic_json', '-m', required = False, metavar = 'string', help = 'The json file describing the Mosaic parameters')

  # Flag for UDN projects as they have a slightly different naming convention
  parser.add_argument('--udn', '-u', required = False, action = 'store_true', help = 'Set for UDN projects to handle the request id in the sample name')

  # Version
  parser.add_argument('--version', '-v', action='version', version='Calypso annotation pipeline version: ' + str(version))

  return parser.parse_args()

# Check the supplied arguments
def checkArguments(args, root_path):

  # If Calypso is being run from the standard directory, the location of the data directory is
  # predictable. If the directory has not been specified on the command line, check for the
  # existence of this known directory
  if not args.data_directory: args.data_directory = '/'.join(root_path.split('/')[0:-1]) + '/data/'
  if not args.tools_directory: args.tools_directory = '/'.join(root_path.split('/')[0:-1]) + '/tools/'

  # Ensure the data path ends with a "/", then add the reference directory and check that the
  # directory exists
  if args.data_directory[-1] != '/': args.data_directory += '/'
  if not os.path.isdir(args.data_directory): fail('ERROR: The data directory does not exist (' + str(args.data_directory) + ')')

  # Check that the tools directory exists and add to the path so scripts from here can be used
  if args.tools_directory[-1] != '/': args.tools_directory += '/'
  if not os.path.exists(args.tools_directory): fail('ERROR: The tools directory does not exist (' + str(args.tools_directory) + ')')

  # Check that the api_client directory exists
  if args.api_client[-1] != '/': args.api_client += '/'
  if not os.path.exists(args.api_client): fail('ERROR: The api client directory does not exist (' + str(args.api_client) + ')')

# Create a directory where all Calypso associated files will be stored
def set_working_directory(version):
  global working_directory

  working_directory += version + "/"

  # If the directory doesn't exist, create it
  if not os.path.exists(working_directory):
    os.makedirs(working_directory)







######
###### Following are routines for reading the resources json
######

# Check the supplied arguments
def check_resources(reference, data_directory, tools_directory, filename):
  resource_info = {}
  resource_info['path'] = data_directory + str(reference) + '/'

  # Store the directory where tools reside if defined
  if tools_directory:
    resource_info['toolsPath'] = str(tools_directory) if tools_directory.endswith('/') else str(tools_directory) + str('/')
  else:
    resource_info['toolsPath'] = False

  if not exists(filename):
    fail('The resource file "' + filename + '" does not exist')
  resource_info['json'] = filename

  # Return the info on resources
  return resource_info

# Define all of the tools to be used in Calypso, check they exist and output their versions
def calypso_tools(resource_info):

  # Define the tools. If no path to the tools is provided, assume the tools are available in the system and don't need a path.
  # The preference is to use the parh to the Calypso tools directory
  resource_info['tools'] = {}
  if resource_info['toolsPath']:
    resource_info['tools']['bcftools'] = str(resource_info['toolsPath']) + 'bcftools/bcftools'
    resource_info['tools']['vcfanno'] = str(resource_info['toolsPath']) + 'vcfanno'
    resource_info['tools']['slivar'] = str(resource_info['toolsPath']) + 'slivar'
    resource_info['tools']['vep'] = str(resource_info['toolsPath']) + 'ensembl-vep/vep'
  else:
    resource_info['tools']['bcftools'] = 'bcftools'
    resource_info['tools']['vcfanno'] = 'vcfanno'
    resource_info['tools']['slivar'] = 'slivar'
    resource_info['tools']['vep'] = 'vep'
  print('Using the following tools:')

  # Bcftools
  bcftoolsVersion = os.popen(resource_info['tools']['bcftools'] +  ' -v').readlines()[0].rstrip()
  print('    ', bcftoolsVersion, sep = '')

  # vcfanno
  proc = Popen(resource_info['tools']['vcfanno'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  print('    ', proc.stderr.readlines()[2].decode().rstrip(), sep = '')

  # slivar
  proc = Popen(resource_info['tools']['slivar'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  print('    ', proc.stderr.readlines()[0].decode().rstrip()[2:], sep = '')

  # vep
  vepInfo = os.popen(resource_info['tools']['vep']).readlines()
  print('    vep versions:')
  for i in range(5, 9, 1): print('    ', vepInfo[i].rstrip(), sep = '')

  # Return the updated resource_info
  return resource_info









######
###### Following are routines for parsing the mosaic json
######

# Check that all Mosaic resources are defined in the resources json
def check_mosaic_resources(mosaic_info, resource_info):
  failed_resources = []
  for resource in mosaic_info['resources']:
    if resource not in resource_info['resources'].keys():
      failed_resources.append(resource)

  # Return the list of failed resources
  return failed_resources

















######
###### Following are routines for building files to be used in the data processing steps
######

# Build the toml file defining all the annotations that vcfanno should use
def buildToml(working_directory, resource_info):

  # Create a toml file
  toml_filename = 'calypso_annotations.toml'
  toml = str(working_directory) + str(toml_filename)
  try:
    toml_file = open(toml, 'w')
  except:
    fail('There was a problem opening a file (calypso_annotations.toml) to write to')

  # Get the path to the data
  datapath = resource_info['path']

  # Add each required resource to the toml file
  for resource in resource_info['resources']:
    if resource_info['resources'][resource]['toml']:
      print('[[annotation]]', file = toml_file)
      print('file="', str(datapath), str(resource_info['resources'][resource]['file']), '"', sep = '', file = toml_file)

      # Some of the annotations include "fields", "columns", and "names". Include these in the toml if they exist
      if 'fields' in resource_info['resources'][resource]:
        text, no_values = toml_info(resource_info, resource, 'fields')
        print(text, file = toml_file)
      if 'columns' in resource_info['resources'][resource]:
        text, no_values = toml_info(resource_info, resource, 'columns')
        print(text, file = toml_file)
      if 'names' in resource_info['resources'][resource]:
        text, no_values = toml_info(resource_info, resource, 'names')
        print(text, file = toml_file)

      # Write out the ops
      ops = 'ops=["'
      no_ops = len(resource_info['resources'][resource]['ops'])
      for i in range(0, no_ops - 1):
        ops += str(resource_info['resources'][resource]['ops'][i]) + '", "'
      ops += str(resource_info['resources'][resource]['ops'][-1]) + '"]'
      print(ops, file = toml_file)
      print(file = toml_file)

      # If there is post annotation information to include in the toml, include it
      if 'post_annotation' in resource_info['resources'][resource]:
        post = resource_info['resources'][resource]['post_annotation']
        print('[[postannotation]]', file = toml_file)
        if 'name' in post:
          print('name="', post['name'], '"', sep = '', file = toml_file)
        if 'fields' in post:
          fieldString = 'fields=["'
          for i in range(0, len(post['fields']) - 1):
            fieldString += str(post['fields'][i]) + '", "'
          fieldString += str(post['fields'][-1]) + '"]'
          print(fieldString, file = toml_file)
        if 'op' in post:
          print('op="', post['op'], '"', sep = '', file = toml_file)
        if 'type' in post:
          print('type="', post['type'], '"', sep = '', file = toml_file)
        print(file = toml_file)

  # Close the toml file
  toml_file.close()

  return toml_filename

# Include information in the toml file
def toml_info(resource_info, resource, infoType):
  no_values = len(resource_info['resources'][resource][infoType])

  # Define the fields to use
  text = infoType + '=['
  for index, value in enumerate(resource_info['resources'][resource][infoType]):
    if infoType == 'columns':
      text += str(value)
    else:
      text += '"' + str(value) + '"'
    if (index + 1) < no_values:
      text += ', '
  text += ']'

  return text, no_values

# Generate a lua file with required functions
def generate_lua_file(working_directory):

  # Create script file for running exomiser
  script_base = 'calypso_vcfanno_lua.lua'
  script_name = str(working_directory) + str(script_base)
  script = open(script_name, 'w')
  print('function hemi(nonpar, xy)', file = script)
  print('  if (xy)', file = script)
  print('  then', file = script)
  print('    if (nonpar == true)', file = script)
  print('    then', file = script)
  print('      return string.format("%d", xy[1])', file = script)
  print('    end', file = script)
  print('  end', file = script)
  print('end', file = script)

  # Close the exomiser script file
  script.close()

  # Return the name of the lua file
  return script_base

# Open a file to write the annotation pipeline script to
def open_bash_script(working_directory):

  # Create a script file
  bash_filename = working_directory + '01_calypso_annotation_pipeline.sh'
  try:
    bash_file = open(bash_filename, "w")
  except:
    fail('There was a problem opening a file (calypso_annotation_pipeline.sh) to write to')

  # Return the file
  return bash_filename, bash_file

# Write all the resource file information to the bash file for the script to use
def bash_resources(resource_info, working_directory, bash_file, vcf, chr_format, ped, lua_filename, toml_filename):

  # Initial information
  print('set -eou pipefail', file = bash_file)
  print(file = bash_file)

  # Define the tools to use
  print('# Define the tools to use while executing the Calypso pipeline', file = bash_file)
  print('DATAPATH=', resource_info['path'], sep = '', file = bash_file)
  if resource_info['toolsPath']:
    print('TOOLPATH=', resource_info['toolsPath'], sep = '', file = bash_file)
    print('BCFTOOLS=$TOOLPATH/bcftools/bcftools', sep = '', file = bash_file)
    print('export BCFTOOLS_PLUGINS=$TOOLPATH/bcftools/plugins', file = bash_file)
    print('VEP=$TOOLPATH/ensembl-vep/vep', sep = '', file = bash_file)
    print('SLIVAR=$TOOLPATH/slivar', sep = '', file = bash_file)
    print('VCFANNO=$TOOLPATH/vcfanno', sep = '', file = bash_file)
  else:
    print('BCFTOOLS=bcftools', file = bash_file)
    print('VEP=vep', file = bash_file)
    print('SLIVAR=slivar', file = bash_file)
    print('VCFANNO=vcfanno', file = bash_file)
  print(file = bash_file)

  # Define the names of the input and output files
  print('# Following are the input VCF and output files created by the pipeline', file = bash_file)
  print('VCF=', vcf, sep = '', file = bash_file)
  print(file = bash_file)

  # Generate the names of the intermediate and final vcf files
  vcf_base = os.path.abspath(vcf).split('/')[-1].rstrip('vcf.gz')
  filtered_vcf = str(vcf_base) + '_calypso_filtered.vcf.gz'
  print('FILEPATH=', working_directory, sep = '', file = bash_file)
  print('ANNOTATEDVCF=$FILEPATH/' + str(vcf_base) + '_annotated.vcf.gz', sep = '', file = bash_file)
  print('FILTEREDVCF=$FILEPATH/' + str(filtered_vcf), sep = '', file = bash_file)
  print('STDOUT=calypso_annotation_pipeline.stdout', file = bash_file)
  print('STDERR=calypso_annotation_pipeline.stderr', file = bash_file)

  # Write the ped file, if necessary
  print('PED=', ped, sep = '', file = bash_file)
  print(file = bash_file)

  # Define the required resources
  print('# Following is a list of required resources', file = bash_file)
  try:
    print('REF=$DATAPATH/', resource_info['resources']['fasta']['file'], sep = '', file = bash_file)
  except:
    fail('The resources json does not define a reference fasta file')
  try:
    print('GFF=$DATAPATH/', resource_info['resources']['gff']['file'], sep = '', file = bash_file)
  except:
    fail('The resources json does not define a gff file')
  print('TOML=$FILEPATH/', toml_filename, sep = '', file = bash_file)
  print('LUA=$FILEPATH/', lua_filename, sep = '', file = bash_file)

  # If the chr map is required, include the path to it. The chromosome format is different depending on whether this is a GRCh37
  # or GRCh38 reference. GRCh37 needs to use the '1' format, but GRCh38 uses 'chr1'. chr_format is False if in the '1' format
  if resource_info['reference'] == 'GRCh37':
    print('CHR_NOCHR_MAP=$DATAPATH/chr_nochr_map.txt', sep = '', file = bash_file)
    print('NOCHR_CHR_MAP=$DATAPATH/nochr_chr_map.txt', sep = '', file = bash_file)
  elif resource_info['reference'] == 'GRCh38':
    print('CHR_NOCHR_MAP=$DATAPATH/chr_nochr_map.txt', sep = '', file = bash_file)
    print('NOCHR_CHR_MAP=$DATAPATH/nochr_chr_map.txt', sep = '', file = bash_file)

  # The slivar js file is based on the family structure. Check that a file exists for the current family
  js_file = 'slivar-functions.js'
  print('JS=$DATAPATH/slivar/', js_file, sep = '', file = bash_file)
  print(file = bash_file)

  # Return the name of the filtered vcf file
  return filtered_vcf

# Annotate the vcf file using bcftools, vcfanno, and VEP
def annotate_vcf(resource_info, bash_file, chr_format, threads, samples):

  # Generate a comma separated list of samples to extract from the vcf file
  sample_list = ''
  for sample in samples:
    if samples[sample]['vcf_sample_name']:
      sample_list += samples[sample]['vcf_sample_name'] + ','
  sample_list = sample_list.rstrip(',')

  # Generate the list of annotations to be output from VEP
  annotation_string = ''
  for annotation in resource_info['resources']['vep']['fields']:
    annotation_string += annotation + ','
  annotation_string = annotation_string.rstrip(',')

  # Print out status messages
  print('# Normalize, subset, and annotate original VCF', file = bash_file)
  print('echo -n "Subsetting, normalizing, and annotating input VCF..."', file = bash_file)

  # Include additional resources for VEP
  if 'export' in resource_info['resources']['vep']:
    for newPath in resource_info['resources']['vep']['export']: print(newPath, file = bash_file)

################
################
################ Limiting to autosome and X was included due to errors, but Y and MT need to be included
################
################
  # Start the script by extracting the required samples from the vcf, and limiting to the autosome and X chromosomes.
  # Ensure the correct format is used for the chromosomes. chr_format is False if in the '1' format
  if resource_info['reference'] == 'GRCh37':
    if not chr_format:
      chroms = '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X'
      print('$BCFTOOLS view -a -c 1 --threads ', threads, ' -s "', sample_list, '" -r "', chroms, '" $VCF 2>> $STDERR \\', sep = '', file = bash_file)
      print('  | $BCFTOOLS annotate -x INFO --threads ', threads, ' 2>> $STDERR \\', sep = '', file = bash_file)
    else:
      chroms = 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX'
      print('$BCFTOOLS view -a -c 1 --threads ', threads, ' -s "', sample_list, '" -r "', chroms, '" $VCF 2>> $STDERR \\', sep = '', file = bash_file)
      print('  | $BCFTOOLS annotate -x INFO --threads ', threads, ' --rename-chrs $CHR_NOCHR_MAP 2>> $STDERR \\', sep = '', file = bash_file)
  elif resource_info['reference'] == 'GRCh38':
    if not chr_format:
      chroms = '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X'
      print('$BCFTOOLS view -a -c 1 --threads ', threads, ' -s "', sample_list, '" -r "', chroms, '" $VCF 2>> $STDERR \\', sep = '', file = bash_file)
      print('  | $BCFTOOLS annotate -x INFO --threads ', threads, ' --rename-chrs $NOCHR_CHR_MAP 2>> $STDERR \\', sep = '', file = bash_file)
    else:
      chroms = 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX'
      print('$BCFTOOLS view -a -c 1 --threads ', threads, ' -s "', sample_list, '" -r "', chroms, '" $VCF 2>> $STDERR \\', sep = '', file = bash_file)
      print('  | $BCFTOOLS annotate -x INFO --threads ', threads, ' 2>> $STDERR \\', sep = '', file = bash_file)

  # Normalize and decompose the vcf, then extract the variants present in the required samples
  print('  | $BCFTOOLS norm -m - -w 10000 -f $REF --threads ', threads, ' 2> $STDERR \\', sep = '', file = bash_file)
  print('  | $BCFTOOLS view -a -c 1 --threads ', threads, ' -s "', sample_list, '" 2>> $STDERR \\', sep = '', file = bash_file)

  # Annotate the vcf file using both bcftools csq and vcfanno
  print('  | $BCFTOOLS csq -f $REF --ncsq 40 -l -g $GFF --threads ', threads, ' 2>> $STDERR \\', sep = '', file = bash_file)
  print('  | $VCFANNO -lua $LUA -p ', threads, ' $TOML /dev/stdin 2>> $STDERR \\', sep = '', file = bash_file)

  # Annotate with VEP unless it is to be ignored
  if not resource_info['resources']['vep']['ignore']:
    print('  | $VEP \\', file = bash_file)
    print('    --assembly ', resource_info['reference'], ' \\', sep = '', file = bash_file)
    print('    --fasta $REF \\', file = bash_file)
    print('    --cache \\', file = bash_file)
    print('    --dir_cache ', resource_info['resources']['vep']['cache'], ' \\', sep = '', file = bash_file)
    print('    --dir_plugins ', resource_info['resources']['vep']['plugins'], ' \\', sep = '', file = bash_file)
    print('    --offline \\', file = bash_file)
    print('    --vcf \\', file = bash_file)
    print('    --force \\', file = bash_file)
    print('    --check_existing \\', file = bash_file)
    print('    --quiet \\', file = bash_file)
    print('    --fork 40 \\', file = bash_file)
    print('    --format vcf \\', file = bash_file)
    print('    --force_overwrite \\', file = bash_file)

    # Add additional options defined in the resource json
    if 'options' in resource_info['resources']['vep']:
      for vep_option in resource_info['resources']['vep']['options']:
        print('    ', vep_option, ' \\', sep = '', file = bash_file)

    # Include all plugins
    if 'active_plugins' in resource_info['resources']['vep']:
      for plugin in resource_info['resources']['vep']['active_plugins']:
        if 'files' in resource_info['resources']['vep']['active_plugins'][plugin]:
          print('    --plugin ', plugin, ',', resource_info['resources']['vep']['active_plugins'][plugin]['files'], ' \\', sep = '', file = bash_file)
        else:
          print('    --plugin ', plugin, ' \\', sep = '', file = bash_file)
    print('    --fields "', annotation_string, '" \\', sep = '', file = bash_file)
    print('    --no_stats \\', file = bash_file)
    print('    --output_file STDOUT \\', file = bash_file)
    print('    2>> $STDERR \\', file = bash_file)

  # Add VAF to the FORMAT and HWE to the info field
  print('  | $BCFTOOLS +fill-tags - -- -t FORMAT/VAF,HWE \\', file = bash_file)

  # Output the final vcf file. If the original file had the '1' format (not 'chr1') for chromosomes,
  # ensure the final vcd file is in this format
  if resource_info['reference'] == 'GRCh37':
    if chr_format:
      print('  | $BCFTOOLS annotate -O z -o $ANNOTATEDVCF --threads ' , threads, ' --rename-chrs $NOCHR_CHR_MAP - \\', sep = '', file = bash_file)
    else:
      print('  | $BCFTOOLS view -O z -o $ANNOTATEDVCF --threads ', threads, ' - \\', sep = '', file = bash_file)
  elif resource_info['reference'] == 'GRCh38':
    if not chr_format:
      print('  | $BCFTOOLS annotate -O z -o $ANNOTATEDVCF --rename-chrs $CHR_NOCHR_MAP --threads ', threads, ' - \\', sep = '', file = bash_file)
    else:
      print('  | $BCFTOOLS view -O z -o $ANNOTATEDVCF --threads ', threads, ' - \\', sep = '', file = bash_file)
  else:
    fail('Unknown reference: ' + str(resource_info['reference']))
  print('  >> $STDOUT 2>> $STDERR', file = bash_file)
  print('echo "complete"', file = bash_file)
  print(file = bash_file)

# Filter the final vcf file to only include variants present in the proband, and based on some
# basic annotations
def filter_vcf(bash_file, samples, proband, resource_info, threads):

  # Create a file containing the name of the proband for use in the filter. Note that there can be multiple probands, e.g.
  # if the family is two affected siblings
  print('# Create file containing the proband(s) only', file = bash_file)
  print('echo -n "Creating proband(s) file..."', file = bash_file)
  print('touch proband.txt', file = bash_file)
  print('echo ', samples[proband]['vcf_sample_name'], ' >> proband.txt', sep = '', file = bash_file)
  print('echo "complete"', file = bash_file)
  print(file = bash_file)

  # Filter the VCF file to generate variants to pass to Mosaic. This is a relatively permissive filter
  # to allow researchers to search variants. Annotate these variants with additional information with VEP
  #   1. Have PASS in the FILTER field
  #   2. Do not have * as the ALT allele
  #   3. Have a popmax AF < 0.01
  #   4. The proband has a genotype containing the alt allele
  print('# Filter the VCF file to generate variants to pass to Mosaic', file = bash_file)
  print('echo -n "Filtering final VCF file..."', file = bash_file)
  print('$BCFTOOLS view --threads ', threads, ' \\', sep = '', file = bash_file)
  print('  -i \'GT[@proband.txt]="alt" & GT[@proband.txt]!="mis"\' \\', file = bash_file)
  print('  $ANNOTATEDVCF \\', file = bash_file)
  print('  2>> $STDERR \\', file = bash_file)
  print('  | $SLIVAR expr --vcf - \\', file = bash_file)
  print('  -p $PED \\', file = bash_file)
  print('  --js $JS \\', file = bash_file)

  # The filter should return all variants that are not in, or have popmax AF < 0.01 in gnomAD. The versions depend on the reference
  if str(resource_info['reference']) == 'GRCh37':
    print('  --info \'( ( (!("ge2_1_1_AF_popmax" in INFO) || ("ge2_1_1_AF_popmax" in INFO && INFO.ge2_1_1_AF_popmax < 0.01)) && ', file = bash_file)
    print('            (  !("gg2_1_1_AF_popmax" in INFO) || ("gg2_1_1_AF_popmax" in INFO && INFO.gg2_1_1_AF_popmax < 0.01)) ) ||', file = bash_file)
  elif str(resource_info['reference']) == 'GRCh38':
    print('  --info \'( ( (!("ge4_0_0_AF_grpmax" in INFO) || ("ge4_0_0_AF_grpmax" in INFO && INFO.ge4_0_0_AF_grpmax < 0.01)) && ', file = bash_file)
    print('            (  !("gg4_0_0_AF_grpmax" in INFO) || ("gg4_0_0_AF_grpmax" in INFO && INFO.gg4_0_0_AF_grpmax < 0.01)) ) ||', file = bash_file)

  # ...or that the variant is present in ClinVar...
  print('  ("CLNSIG" in INFO) ||', file = bash_file)

  # ...or that the CADD_Phred score is over 15...
  print('  ("CADD_Phred" in INFO && INFO.CADD_Phred > 15) ||', file = bash_file)

  # ...or that the REVEL or MutScores are over 0.1...
  print('  ("REVEL" in INFO && INFO.REVEL > 0.1) ||', file = bash_file)
  print('  ("MutScore" in INFO && INFO.MutScore > 0.1) ||', file = bash_file)

  # ...or that the variant has a spliceAI score > 0.1 for any of the acceptor / donor sites.
  print('  ("SpliceAI_DS_AG" in INFO && INFO.SpliceAI_DS_AG > 0.1) ||', file = bash_file)
  print('  ("SpliceAI_DS_AL" in INFO && INFO.SpliceAI_DS_AL > 0.1) ||', file = bash_file)
  print('  ("SpliceAI_DS_DG" in INFO && INFO.SpliceAI_DS_DG > 0.1) ||', file = bash_file)
  print('  ("SpliceAI_DS_DL" in INFO && INFO.SpliceAI_DS_DL > 0.1) ) &&', file = bash_file)

  # ...and that the variant does not have "*" as the alt allele
  print('  variant.FILTER == "PASS" && variant.ALT[0] != "*"\' \\', file = bash_file)

  # Then tag variants based on their quality
  ###### TEMP CHANGE FOR TESTING
  ###### There medium quality are being called high quality for now pending user testing
#  print('  --trio \"het_low_qual:kid.het && ', file = bash_file)
#  print('    ((kid.GQ >= 10 && kid.GQ < 20 && kid.AB >= 0.1 && kid.AB <= 0.9) || ', file = bash_file)
#  print('    (kid.GQ >= 20 && kid.AB >= 0.1 && kid.AB < 0.2) || ', file = bash_file)
#  print('    (kid.GQ >= 20 && kid.AB > 0.8 && kid.AB <= 0.9))\" \\', file = bash_file)
#  print('  --trio "het_med_qual:kid.het && ', file = bash_file)
#  print('    ((kid.GQ >=20 && kid.GQ < 30 && kid.AB >= 0.2 && kid.AB <= 0.8) || ', file = bash_file)
#  print('    (kid.GQ >=30 && kid.AB >= 0.2 && kid.AB < 0.3) || ', file = bash_file)
#  print('    (kid.GQ >=30 && kid.AB > 0.7 && kid.AB <= 0.8))\" \\', file = bash_file)
#  print('  --trio \"het_hi_qual:kid.het && kid.GQ >= 30 && kid.AB >= 0.3 && kid.AB <= 0.7\" \\', file = bash_file)
#  print('  --trio \"hom_low_qual:kid.hom_alt && ', file = bash_file)
#  print('    ((kid.GQ >= 10 && kid.GQ < 20 && kid.AB >= 0.7) || ', file = bash_file)
#  print('    (kid.GQ >= 20 && kid.AB >= 0.7 && kid.AB < 0.8))\" \\', file = bash_file)
#  print('  --trio "hom_med_qual: kid.hom_alt && ', file = bash_file)
#  print('    ((kid.GQ >= 20 && kid.GQ < 30 && kid.AB >= 0.8) || ', file = bash_file)
#  print('    (kid.GQ >= 30 && kid.AB >= 0.8 && kid.AB < 0.9))\" \\', file = bash_file)
#  print('  --trio \"hom_hi_qual: kid.hom_alt && kid.GQ >= 30 && kid.AB >= 0.9\" \\', file = bash_file)
  print('  --trio \"het_low_qual:kid.het && ', file = bash_file)
  print('    ((kid.GQ >= 10 && kid.GQ < 20 && kid.AB >= 0.1 && kid.AB <= 0.9) || ', file = bash_file)
  print('    (kid.GQ >= 20 && kid.AB >= 0.1 && kid.AB < 0.2) || ', file = bash_file)
  print('    (kid.GQ >= 20 && kid.AB > 0.8 && kid.AB <= 0.9))\" \\', file = bash_file)
  print('  --trio \"het_hi_qual:kid.het && kid.GQ >= 20 && kid.AB >= 0.2 && kid.AB <= 0.8\" \\', file = bash_file)
  print('  --trio \"hom_low_qual:kid.hom_alt && ', file = bash_file)
  print('    ((kid.GQ >= 10 && kid.GQ < 20 && kid.AB >= 0.7) || ', file = bash_file)
  print('    (kid.GQ >= 20 && kid.AB >= 0.7 && kid.AB < 0.8))\" \\', file = bash_file)
  print('  --trio \"hom_hi_qual: kid.hom_alt && kid.GQ >= 20 && kid.AB >= 0.8\" \\', file = bash_file)

  # And final compression
  print('  --pass-only \\', file = bash_file)
  print('  2>> $STDERR \\', file = bash_file)
  print('  | $BCFTOOLS view -O z -o $FILTEREDVCF --threads ', threads, ' - \\', sep = '', file = bash_file)
  print('  >> $STDOUT 2>> $STDERR', file = bash_file)
  print('echo "complete"', file = bash_file)
  print('$BCFTOOLS index -t $FILTEREDVCF', file = bash_file)
  print(file = bash_file)








#####
##### Following are scripts to deal with annotation tsvs
#####

# Call the command to extract the HGVS annotations from the filtered vcf into a tsv file
def generate_hgvs_tsv(bash_file, reference):
  print(file = bash_file)
  print('# Resource: HGVS', sep = '', file = bash_file)
  print('echo -n "Creating tsv file for HGVS..."', sep = '', file = bash_file)
  print('python3 $GENERATE_HGVS_TSV ', end = '', file = bash_file)
  print('-i $FILTEREDVCF ', end = '', file = bash_file)
  print('-r "', reference, '" ', sep = '', end = '', file = bash_file)
  print('-m $MOSAIC_JSON ', end = '', file = bash_file)
  print('-s $TOOLPATH ', file = bash_file)
  print('echo "complete"', file = bash_file)

# Call the command to extract the variant quality values
def generate_variant_quality_tsv(bash_file, uid):
  print(file = bash_file)
  print('# Resource: Variant Quality', sep = '', file = bash_file)
  print('echo -n "Creating tsv file for Variant Quality..."', sep = '', file = bash_file)
  print('python3 $GENERATE_VARIANT_QUALITY_TSV ', end = '', file = bash_file)
  print('-i $FILTEREDVCF ', end = '', file = bash_file)
  print('-u "', uid, '" ', sep = '', end = '', file = bash_file)
  print('-t $TOOLPATH ', file = bash_file)
  print('echo "complete"', file = bash_file)

# Generate the command to process general annotations
def generate_annotation_tsv(bash_file, resource, reference):
  output_tsv = str(resource) + '.tsv'

  print(file = bash_file)
  print('# Resource: ', resource, sep = '', file = bash_file)
  print('echo -n "Creating tsv file for ', resource, '..."', sep = '', file = bash_file)
  print('python3 $GENERATE_TSV ', end = '', file = bash_file)
  print('-i $FILTEREDVCF ', sep = '', end = '', file = bash_file)
  print('-o ', output_tsv, ' ', sep = '', end = '', file = bash_file)
  print('-e ', resource, ' ', sep = '', end = '', file = bash_file)
  print('-r ', reference, ' ', sep = '', end = '', file = bash_file)
  print('-m $MOSAIC_JSON ', end = '', file = bash_file)
  print('-s $TOOLPATH', file = bash_file)
  print('echo "complete"', file = bash_file)

  # Return the name of the output file
  return output_tsv

















#####
##### Following are scripts for Exomiser
#####

# Generate the application properties file
def application_properties(working_dir, tools_dir, reference):

  # Create the application.properties file
  properties_file = open(str(working_dir) + 'application.properties', 'w')

  # Write the data versions to the file
  print('exomiser.data-directory=', tools_dir, 'exomiser-cli-14.0.0/data', sep = '', file = properties_file)
  print('remm.version=0.3.1.post1', file = properties_file)
  print('cadd.version=1.4', file = properties_file)
  if reference == 'GRCh37':
    print('exomiser.hg19.data-version=2402', file = properties_file)
  elif reference == 'GRCh38':
    print('exomiser.hg38.data-version=2402', file = properties_file)
  print('exomiser.phenotype.data-version=2402', file = properties_file)

  # Close the application properties file
  properties_file.close()

# Generate the yml file for running exomiser
def generate_yml(working_dir, proband, reference, vcf, ped, hpo):

  # Determine if exomiser should use resources for hg19 or hg38
  if reference == 'GRCh37':
    exomiser_ref = 'hg19'
  elif reference == 'GRCh38':
    exomiser_ref = 'hg38'
  else:
    fail('Unknown reference for use with exomiser: ' + str(reference))

  # Generate a list of hpo terms, with the terms contained in single quotes
  hpo = hpo.split(',') if ',' in hpo else [hpo]
  hpo_string = '\', \''.join(hpo)

  # Open a new yaml file
  yaml_name = 'exomiser_' + str(proband) + '.yml'
  yaml = open(str(working_dir) + str(yaml_name), 'w')

  # Write the information about this project to the yaml
  print('# Exomiser analysis file for proband: ', proband, sep = '', file = yaml)
  print(file = yaml)
  print('analysis:', file = yaml)
  print('  genomeAssembly: ', str(exomiser_ref), sep = '', file = yaml)
  print('  vcf: ', vcf, sep = '', file = yaml)
  print('  ped: ', ped, sep = '', file = yaml)
  print('  proband: ', proband, sep = '', file = yaml)
  print('  hpoIds: [\'', hpo_string, '\']', sep = '', file = yaml)
  print(file = yaml)

  # Include the inheritance modes to consider
  # These are the default settings, with values representing the maximum minor allele frequency in percent (%) permitted for an
  # allele to be considered as a causative candidate under that mode of inheritance.
  # If you just want to analyse a sample under a single inheritance mode, delete/comment-out the others. For AUTOSOMAL_RECESSIVE
  # or X_RECESSIVE ensure *both* relevant HOM_ALT and COMP_HET modes are present.
  # In cases where you do not want any cut-offs applied an empty map should be used e.g. inheritanceModes: {}
  print('  inheritanceModes: {', file = yaml)
  print('    AUTOSOMAL_DOMINANT: 0.1,', file = yaml)
  print('    AUTOSOMAL_RECESSIVE_HOM_ALT: 0.1,', file = yaml)
  print('    AUTOSOMAL_RECESSIVE_COMP_HET: 2.0,', file = yaml)
  print('    X_DOMINANT: 0.1,', file = yaml)
  print('    X_RECESSIVE_HOM_ALT: 0.1,', file = yaml)
  print('    X_RECESSIVE_COMP_HET: 2.0,', file = yaml)
  print('    MITOCHONDRIAL: 0.2', file = yaml)
  print('  }', file = yaml)
  print(file = yaml)

  # Use the PASS_ONLY analysis mode
  print('  analysisMode: PASS_ONLY', file = yaml)
  print(file = yaml)

  # Define the allele frequency sources to be used
  print('  frequencySources: [', file = yaml)
  print('    THOUSAND_GENOMES,', file = yaml)
  print('    TOPMED,', file = yaml)
  print('    UK10K,', file = yaml)
  print('    ESP_AFRICAN_AMERICAN,', file = yaml)
  print('    ESP_EUROPEAN_AMERICAN,', file = yaml)
  print('    ESP_ALL,', file = yaml)
  print('    EXAC_AFRICAN_INC_AFRICAN_AMERICAN,', file = yaml)
  print('    EXAC_AMERICAN,', file = yaml)
  print('    EXAC_SOUTH_ASIAN,', file = yaml)
  print('    EXAC_EAST_ASIAN,', file = yaml)
  print('    EXAC_FINNISH,', file = yaml)
  print('    EXAC_NON_FINNISH_EUROPEAN,', file = yaml)
  print('    EXAC_OTHER,', file = yaml)
  print('    GNOMAD_E_AFR,', file = yaml)
  print('    GNOMAD_E_AMR,', file = yaml)
  print('    GNOMAD_E_EAS,', file = yaml)
  print('    GNOMAD_E_FIN,', file = yaml)
  print('    GNOMAD_E_NFE,', file = yaml)
  print('    GNOMAD_E_SAS,', file = yaml)
  print('    GNOMAD_E_OTH,', file = yaml)
  print('    GNOMAD_G_AFR,', file = yaml)
  print('    GNOMAD_G_AMR,', file = yaml)
  print('    GNOMAD_G_EAS,', file = yaml)
  print('    GNOMAD_G_FIN,', file = yaml)
  print('    GNOMAD_G_NFE,', file = yaml)
  print('    GNOMAD_G_SAS,', file = yaml)
  print('    GNOMAD_G_OTH', file = yaml)
  print('  ]', file = yaml)
  print(file = yaml)

  # Define the resources to be used for determining pathogenicity
  print('  pathogenicitySources: [', file = yaml)
  print('    REVEL,', file = yaml)
  print('    MVP,', file = yaml)
  print('    ALPHA_MISSENSE', file = yaml)
  print('  ]', file = yaml)
  print(file = yaml)

  # Define the order of steps for exomiser to take
  print('  steps: [', file = yaml)
  print('    failedVariantFilter: { },', file = yaml)
  print('    variantEffectFilter: {', file = yaml)
  print('      remove: [', file = yaml)
  print('        FIVE_PRIME_UTR_EXON_VARIANT,', file = yaml)
  print('        FIVE_PRIME_UTR_INTRON_VARIANT,', file = yaml)
  print('        THREE_PRIME_UTR_EXON_VARIANT,', file = yaml)
  print('        THREE_PRIME_UTR_INTRON_VARIANT,', file = yaml)
  print('        NON_CODING_TRANSCRIPT_EXON_VARIANT,', file = yaml)
  print('        NON_CODING_TRANSCRIPT_INTRON_VARIANT,', file = yaml)
  print('        CODING_TRANSCRIPT_INTRON_VARIANT,', file = yaml)
  print('        UPSTREAM_GENE_VARIANT,', file = yaml)
  print('        DOWNSTREAM_GENE_VARIANT,', file = yaml)
  print('        INTERGENIC_VARIANT,', file = yaml)
  print('        REGULATORY_REGION_VARIANT', file = yaml)
  print('      ]', file = yaml)
  print('    },', file = yaml)
  print('    frequencyFilter: {maxFrequency: 2.0},', file = yaml)
  print('    pathogenicityFilter: {keepNonPathogenic: true},', file = yaml)
  print('    inheritanceFilter: {},', file = yaml)
  print('    omimPrioritiser: {},', file = yaml)
  print('    hiPhivePrioritiser: {runParams: \'human\'},', file = yaml)
  print('  ]', file = yaml)
  print(file = yaml)

  # Define the output options
  print('outputOptions:', file = yaml)
  print('  outputContributingVariantsOnly: false', file = yaml)
  print('  numGenes: 0', file = yaml)
  print('  outputDirectory: ', str(working_dir) + 'exomiser_results', file = yaml)
  print('  outputFileName: ', str(proband), file = yaml)
  print('  outputFormats: [HTML, JSON, TSV_GENE, TSV_VARIANT, VCF]', file = yaml)

  # Close the yaml file
  yaml.close()

  # Return the path and name of the yaml file
  return yaml_name

# Generate a script file to run exomiser
def generate_exomiser_script(working_dir, tools_dir, yaml):

  # Create script file for running exomiser
  script_name = str(working_directory) + '04_calypso_exomiser.sh'
  script = open(script_name, 'w')
  print('set -eou pipefail', file = script)
  print(file = script)

  # Create variables for paths
  print('TOOLS_PATH=', tools_dir, sep = '', file = script)
  print('WORKINGPATH=', working_dir, sep = '', file = script)
#########
#########
######### Update with version from resources json
#########
#########
  print('EXOMISER=$TOOLS_PATH/exomiser-cli-14.0.0/exomiser-cli-14.0.0.jar', file = script)
  print('YML=$WORKINGPATH/', yaml, sep = '', file = script)
  print('STDOUT=', working_dir, 'exomiser.stdout', sep = '', file = script)
  print('STDERR=', working_dir, 'exomiser.stderr', sep = '', file = script)
  print(file = script)

#####
##### UTAH SPECIFIC REQUIREMENT
#####
  print('# Module requirement for Utah environment', file = script)
  print('module load openjdk/17.0.1', file = script)
  print(file = script)

  # Include the executable for running exomiser
  print('java -jar $EXOMISER ', end = '', file = script)
  print('--analysis $YML ', end = '', file = script)
  print('> $STDOUT 2> $STDERR', file = script)
  print(file = script)
  print('# Write to screen if an error was encountered', file = script)
  print('if [[ $? != 0 ]];', file = script)
  print('then', file = script)
  print('  echo "Exomiser script failed. Please check the stderr file"', file = script)
  print('fi', file = script)
  print(file = script)

  # Close the script file and make executable
  script.close()
  make_executable = os.popen('chmod +x ' + str(script_name)).read()

  # Return the script
  return script_name, script

# Create an upload variants script to upload the exomiser variants
def upload_exomiser_variants(working_dir, api_client, client_config, project_id, proband):

  # Open a script file
  filename = working_dir + 'calypso_exomiser_upload_variants.sh'
  try:
    upload_file = open(filename, 'w')
  except:
    fail('Could not open ' + str(filename) + ' to write to')

  # Write the command to file to upload the filtered variants
  print('WORKINGPATH=', str(working_dir), sep = '', file = upload_file)
  print('API_CLIENT=', api_client, sep = '', file = upload_file)
  print('CONFIG=', client_config, sep = '', file = upload_file)
  print('VCF=$WORKINGPATH/exomiser_results/', str(proband), '.vcf.gz', sep = '', file = upload_file)
  print(file = upload_file)
  print('python3 $API_CLIENT/variants/upload_variants.py -a $API_CLIENT -c $CONFIG -p ', str(project_id), ' -m "no-validation" -v $VCF', sep = '', file = upload_file)

  # Close the file
  upload_file.close()

  # Make the annotation script executable
  make_executable = os.popen('chmod +x ' + str(filename)).read()

# Create a script to update and upload the exomiser annotations
def exomiser_annotations(calypso_dir, working_dir, api_client, client_config, project_id, proband, filters_json):

  # Open a script file
  script_name = working_dir + '05_calypso_exomiser_annotations.sh'
  try:
    script = open(script_name, 'w')
  except:
    fail('Could not open ' + str(script_name) + ' to write to')

  # Define the name of the output vcf and tsv
  tsv = 'exomiser.tsv'
  vcf = 'exomiser.vcf.gz'

  print(file = script)
  print('CALYPSO_PATH=', calypso_dir, sep = '', file = script)
  print('WORKING_PATH=', working_dir, sep = '', file = script)
  print('API_CLIENT=', api_client, sep = '', file = script)
  print('CONFIG=', client_config, sep = '', file = script)
  print('SCRIPT=$CALYPSO_PATH/parse_exomiser_variants.py', sep = '', file = script)
  print('INPUT_TSV=$WORKING_PATH/exomiser_results/', str(proband), '.variants.tsv', sep = '', file = script)
  print('OUTPUT_TSV=$WORKING_PATH/', tsv, sep = '', file = script)
  print(file = script)
  print('# Parse the resulting output file', file = script)
  print('echo -n "Creating tsv file for exomiser variants..."', file = script)
  print('python3 $SCRIPT ', sep = '', end = '', file = script)
  print('-a $API_CLIENT ', sep = '', end = '', file = script)
  print('-c $CONFIG ', sep = '', end = '', file = script)
  print('-p ', str(project_id), ' ', sep = '', end = '', file = script)
  print('-i $INPUT_TSV ', sep = '', end = '', file = script)
  print('-o  $OUTPUT_TSV', sep = '', file = script)
  print('echo "complete"', file = script)

  # Apply exomiser variant filters
  print(file = script)
  print('# Set exomiser variant filters', file = script)
  print('UPLOAD_SCRIPT=$API_CLIENT/project_setup/set_variant_filters.py', file = script)
  print('FILTERS_JSON=', filters_json, sep = '', file = script)
  print('python3 $UPLOAD_SCRIPT ', end = '', file = script)
  print('-a $API_CLIENT ', end = '', file = script)
  print('-c $CONFIG ', end = '', file = script)
  print('-p ', str(project_id), ' ', sep = '', end = '', file = script)
  print('-f $FILTERS_JSON', file = script)

  # If the parsing completed successfully, upload the annotations to mosaic
  print(file = script)
  print('# Upload exomiser annotations', sep = '', file = script)
  print('UPLOAD_SCRIPT=$API_CLIENT/variant_annotations/upload_annotations.py', sep = '', file = script)
  print('python3 $UPLOAD_SCRIPT ', end = '', file = script)
  print('-a $API_CLIENT ', end = '', file = script)
  print('-c $CONFIG ', end = '', file = script)
  print('-p ', str(project_id), ' ', sep = '', end = '', file = script)
  print('-t $OUTPUT_TSV', sep = '', file = script)

  # Close the exomiser script and make it executable
  script.close()
  make_executable = os.popen('chmod +x ' + script_name).read()








# If the script fails, provide an error message and exit
def fail(message):
  print(message, sep = '')
  exit(1)

# Initialise global variables

# Pipeline version
version = "2.0"
date = str(date.today())

# The working directory where all created files are kept
working_directory = os.getcwd() + "/calypso_v" + version + "r"

# Store the allowed references that can be specified on the command line
allowed_references = ['GRCh37', 'GRCh38']

if __name__ == "__main__":
  main()

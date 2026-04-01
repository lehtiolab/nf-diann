process createNewSpectraLookup {
  container 'quay.io/biocontainers/msstitch:3.16--pyhdfd78af_0'

  input:
  tuple path(mzml), path(dino)
  
  output:
  path('mslookup_db.sqlite')

  script:
  """
  msstitch storespectra --spectra "${mzml}" --setnames 'QC'
  ${dino.baseName != 'NO__FILE' ? "msstitch storequant --dbfile mslookup_db.sqlite --dinosaur \"${dino}\" --spectra \"${mzml}\" --mztol 20.0 --mztoltype ppm --rttol 5.0" : ''}
  """
}


def read_header(info_fn) {
  def header = []
  def info = file(info_fn).eachLine { line, ix ->
      if (ix == 1) {
        header = line.tokenize('\t')
      }
  }
  return header
}


def create_info_map(info_fn, possible_params) {
  /* From possible params this parses a tab separated input file with
  a header (which has some of those params names.
  It returns a map with params and their values, which are defaults for
  those not set in the input file.
  */
  def info = file(info_fn).readLines().collect { it.tokenize('\t') }
  def header = info.pop()

  def params_not_header = possible_params - header

  def info_map = [:]
  info.findAll{ it[0][0] != '#' }.eachWithIndex { it, ix ->
    tmp_fn = [:]
    header.eachWithIndex{ hfield, hix ->
    tmp_fn[hfield] = it[hix]
    }
    fpath = file(tmp_fn.file_path)
    key_fn = ix //tmp_fn.file_path
    info_map[key_fn] = tmp_fn
    info_map[key_fn].id = ix
    info_map[key_fn].file_path = fpath
    info_map[key_fn].filename = "${info_map[key_fn].file_path.baseName}.${info_map[key_fn].file_path.extension}"
    params_not_header.each {
      info_map[key_fn][it] = params[it]
    }
  }
  return info_map
}


def identify_info_map(info_fn) {
  expected_fields = ["file_path", "create_lib", "train_quantums", "searchbatch", "quantfile"]
  def samples = create_info_map(info_fn, expected_fields)
  // Set all files to go to lib create if none specified (and no params.library passed)
  new_samples = samples
  if (!samples.findAll { k,v -> v.create_lib as Integer }) {
    samples.each { k,v ->
      new_samples[k].create_lib = 1
    }
    samples = new_samples
  }
  // Set all files to go to train quantUMS if none specified
  if (!samples.findAll { k,v -> v.train_quantums as Integer }) {
    samples.each { k,v ->
      new_samples[k].train_quantums = 1
    }
    samples = new_samples
  }
  return samples
}


def listify(it) {
  /* This function is useful when needing a list even when having a single item
  - Single items in channels get unpacked from a list
  - Processes expect lists. Even though it would be fine
  without a list, for single-item-lists any special characters are not escaped by NF
  in the script, which leads to errors. See:
  https://github.com/nextflow-io/nextflow/discussions/4240
  */
  return it instanceof java.util.List ? it : [it]
}


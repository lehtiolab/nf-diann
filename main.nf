#!/usr/bin/env nextflow

include { paramsSummaryMap } from 'plugin/nf-schema'

import groovy.io.FileType

include { identify_info_map; listify } from './modules.nf' 
//include { REPORTING } from './workflows/reporting.nf'



process predictFastaLibrary {
cache 'lenient'

  container "ghcr.io/lehtiolab/nfhelaqc:3.2-diann.2.3.1"
  
  input:
  tuple path(fasta), val(diannparams)

  output:
  path('library.predicted.speclib')

  script:
  """
  # Create predicted library from fasta
  diann-linux --threads ${task.cpus} \
    ${fasta.collect { "--fasta $it" }.join(' ')} \
    --cut ${diannparams.cut} \
    --gen-spec-lib \
    --predictor \
    --fasta-search \
    --out-lib library.speclib \
    --missed-cleavages $diannparams.miscleav \
    ${diannparams.varmods.collect { "--var-mod $it" }.join(' ')} \
    ${diannparams.fixmods.collect { "--fixed-mod $it" }.join(' ')} \
    --var-mods $diannparams.maxvarmods \
    ${diannparams.ntermmetex ? '--met-excision' : ''} \
    ${diannparams.ntermac ? '--var-mod UniMod:1,42.010565,*n' : ''} \
    ${diannparams.idstonames ? '--ids-to-names' : ''} \
    --pg-level ${diannparams.pglvl} \
    --min-pr-charge ${diannparams.mincharge} \
    --max-pr-charge ${diannparams.maxcharge} \
    --min-pep-len ${diannparams.minpeplen} \
    --max-pep-len ${diannparams.maxpeplen} \
    --min-pr-mz ${diannparams.minmz} \
    --max-pr-mz ${diannparams.maxmz} \
    --min-fr-mz ${diannparams.minfrmz} \
    --max-fr-mz ${diannparams.maxfrmz} \
    ${diannparams.pglvl} \
    ${diannparams.exclude_contaminants ? "--cont-quant-exclude ${diannparams.exclude_contaminants}" : ''}

    # Is this used in predict from fasta, but maybe test this:
    #${diannparams.window ? "--window $diannparams.window" : ''} \
    #--mass-acc-ms1 $diannparams.ms1acc \
    #--mass-acc $diannparams.ms2acc \
    """
}


process createEmpiricalLibrary {
cache 'lenient'

  container "ghcr.io/lehtiolab/nfhelaqc:3.2-diann.2.3.1"
  
  input:
  tuple path(predlib), path(raws), path(fasta), val(diannparams)

  output:
  path('library.parquet')

  script:
  """
  # Empirical library by running it with raws
  diann-linux --threads ${task.cpus} \
    ${raws.collect { "--f $it" }.join(' ') } \
    ${fasta.collect { "--fasta $it" }.join(' ')} \
    --cut ${diannparams.cut} \
    --lib $predlib \
    --gen-spec-lib \
    --out-lib library \
    --rt-profiling \
    --missed-cleavages $diannparams.miscleav \
    ${diannparams.ms1acc ? "--mass-acc-ms1 ${diannparams.ms1acc}" : ''} \
    ${diannparams.ms2acc ? "--mass-acc ${diannparams.ms2acc}" : ''} \
    ${diannparams.window ? "--window $diannparams.window" : ''} \
    ${diannparams.varmods.collect { "--var-mod $it" }.join(' ')} \
    ${diannparams.fixmods.collect { "--fixed-mod $it" }.join(' ')} \
    --var-mods $diannparams.maxvarmods \
    ${diannparams.ntermmetex ? '--met-excision' : ''} \
    ${diannparams.ntermac ? '--var-mod UniMod:1,42.010565,*n' : ''} \
    --min-pr-charge ${diannparams.mincharge} \
    --max-pr-charge ${diannparams.maxcharge} \
    --min-pep-len ${diannparams.minpeplen} \
    --max-pep-len ${diannparams.maxpeplen} \
    --min-pr-mz ${diannparams.minmz} \
    --max-pr-mz ${diannparams.maxmz} \
    --min-fr-mz ${diannparams.minfrmz} \
    --max-fr-mz ${diannparams.maxfrmz} \
    ${diannparams.indiwin ? "--individual-windows" : ''} \
    ${diannparams.indiacc ? "--individual-mass-acc" : ''} \
    ${diannparams.exclude_contaminants ? "--cont-quant-exclude ${diannparams.exclude_contaminants}" : ''}
  """
// rt-profiling add to the empirical step
}


process RunDiaAnalysis {

  // "Second pass" after creating the library, manual MBR
  // No reports are output, only quant files
  container "ghcr.io/lehtiolab/nfhelaqc:3.2-diann.2.3.1"

  input:
  tuple val(ids), path(raws), path(lib), path(fasta), val(diannparams)
  
  output:
  tuple val(ids), path(raws), path('*.quant')

  script:
  """
  diann-linux --threads ${task.cpus} \
    ${raws.collect { "--f $it" }.join(' ') } \
    --lib $lib \
    ${fasta.collect { "--fasta $it" }.join(' ')} \
    --cut ${diannparams.cut} \
    --temp ./ \
    ${diannparams.ms1acc ? "--mass-acc-ms1 ${diannparams.ms1acc}" : ''} \
    ${diannparams.ms2acc ? "--mass-acc ${diannparams.ms2acc}" : ''} \
    ${diannparams.window ? "--window $diannparams.window" : ''} \
    ${diannparams.varmods.collect { "--var-mod $it" }.join(' ')} \
    ${diannparams.fixmods.collect { "--fixed-mod $it" }.join(' ')} \
    ${diannparams.ntermmetex ? '--met-excision' : ''} \
    ${diannparams.ntermac ? '--var-mod UniMod:1,42.010565,*n' : ''} \
    ${diannparams.nonorm ? '--no-norm' : ''} \
    ${diannparams.idstonames ? '--ids-to-names' : ''} \
    --pg-level ${diannparams.pglvl} \
    --min-pr-charge ${diannparams.mincharge} \
    --max-pr-charge ${diannparams.maxcharge} \
    --min-pep-len ${diannparams.minpeplen} \
    --max-pep-len ${diannparams.maxpeplen} \
    --min-pr-mz ${diannparams.minmz} \
    --max-pr-mz ${diannparams.maxmz} \
    --min-fr-mz ${diannparams.minfrmz} \
    --max-fr-mz ${diannparams.maxfrmz} \
    ${diannparams.indiwin ? "--individual-windows" : ''} \
    ${diannparams.indiacc ? "--individual-mass-acc" : ''} \
    ${diannparams.exclude_contaminants ? "--cont-quant-exclude ${diannparams.exclude_contaminants}" : ''}
  """
}


process TrainQuantUMS {

  container "ghcr.io/lehtiolab/nfhelaqc:3.2-diann.2.3.1"

  input:
  tuple path(raws), path(quants), path(lib), path(fasta), val(diannparams)
  
  output:
  stdout

  script:
  paramline = '.*Quantification parameters:' 
  """
  diann-linux --threads ${task.cpus} \
    ${raws.collect { "--f \$(realpath $it)"}.join(' ') } \
    --lib $lib \
    ${fasta.collect { "--fasta $it" }.join(' ')} \
    --cut ${diannparams.cut} \
    --use-quant \
    --temp ./ \
    --quant-train-runs 0:${listify(raws).size() -1} \
    ${diannparams.ms1acc ? "--mass-acc-ms1 ${diannparams.ms1acc}" : ''} \
    ${diannparams.ms2acc ? "--mass-acc ${diannparams.ms2acc}" : ''} \
    ${diannparams.window ? "--window $diannparams.window" : ''} \
    ${diannparams.varmods.collect { "--var-mod $it" }.join(' ')} \
    ${diannparams.fixmods.collect { "--fixed-mod $it" }.join(' ')} \
    ${diannparams.ntermmetex ? '--met-excision' : ''} \
    ${diannparams.ntermac ? '--var-mod UniMod:1,42.010565,*n' : ''} \
    ${diannparams.nonorm ? '--no-norm' : ''} \
    ${diannparams.idstonames ? '--ids-to-names' : ''} \
    --pg-level ${diannparams.pglvl} \
    --min-pr-charge ${diannparams.mincharge} \
    --max-pr-charge ${diannparams.maxcharge} \
    --min-pep-len ${diannparams.minpeplen} \
    --max-pep-len ${diannparams.maxpeplen} \
    --min-pr-mz ${diannparams.minmz} \
    --max-pr-mz ${diannparams.maxmz} \
    --min-fr-mz ${diannparams.minfrmz} \
    --max-fr-mz ${diannparams.maxfrmz} \
    ${diannparams.indiwin ? "--individual-windows" : ''} \
    ${diannparams.indiacc ? "--individual-mass-acc" : ''} \
    ${diannparams.exclude_contaminants ? "--cont-quant-exclude ${diannparams.exclude_contaminants}" : ''} \
       > tmpstdout.log
    grep '$paramline' tmpstdout.log | sed 's/$paramline/\\-\\-quant-params/'
  """
}

process DiaQuantificationReport {

  container "ghcr.io/lehtiolab/nfhelaqc:3.2-diann.2.3.1"

  input:
  tuple path(raws), path(quants), path(lib), path(fasta), val(diannparams), val(quantparams)
  
  output:
  tuple path('report.parquet'), path('*.tsv'), path('report.log.txt')

  script:
  """
  diann-linux --threads ${task.cpus} \
    ${raws.collect { "--f \$(realpath $it)"}.join(' ') } \
    --lib $lib \
    ${fasta.collect { "--fasta $it" }.join(' ')} \
    --cut ${diannparams.cut} \
    --use-quant \
    --temp ./ \
    ${quantparams.trim()} \
    ${diannparams.ms1acc ? "--mass-acc-ms1 ${diannparams.ms1acc}" : ''} \
    ${diannparams.ms2acc ? "--mass-acc ${diannparams.ms2acc}" : ''} \
    ${diannparams.window ? "--window $diannparams.window" : ''} \
    ${diannparams.varmods.collect { "--var-mod $it" }.join(' ')} \
    ${diannparams.fixmods.collect { "--fixed-mod $it" }.join(' ')} \
    --matrices \
    ${diannparams.ntermmetex ? '--met-excision' : ''} \
    ${diannparams.ntermac ? '--var-mod UniMod:1,42.010565,*n' : ''} \
    ${diannparams.nonorm ? '--no-norm' : ''} \
    ${diannparams.idstonames ? '--ids-to-names' : ''} \
    --pg-level ${diannparams.pglvl} \
    --min-pr-charge ${diannparams.mincharge} \
    --max-pr-charge ${diannparams.maxcharge} \
    --min-pep-len ${diannparams.minpeplen} \
    --max-pep-len ${diannparams.maxpeplen} \
    --min-pr-mz ${diannparams.minmz} \
    --max-pr-mz ${diannparams.maxmz} \
    --min-fr-mz ${diannparams.minfrmz} \
    --max-fr-mz ${diannparams.maxfrmz} \
    ${diannparams.indiwin ? "--individual-windows" : ''} \
    ${diannparams.indiacc ? "--individual-mass-acc" : ''} \
    --qvalue ${diannparams.precfdr} \
    --matrix-qvalue ${diannparams.protfdr} \
    ${diannparams.exclude_contaminants ? "--cont-quant-exclude ${diannparams.exclude_contaminants}" : ''}
  """
}


workflow {
  main:
  // FIXME
//  ms1acc = [timstof: 20, velos: 10, qe: 10, astral: 10][instrument]
//  ms2acc = 20

  if (params.library && params.outputlib) {
    exit 1, 'Cannot output generated library while also being passed a --library'
  }

  cut = ['trypsin': 'K*,R*,!*P', 'trypsinp': 'K*,R*'][params.enzyme]
  diann_params = [
    varmods: params.varmods.tokenize(';'),
    fixmods: params.fixmods.tokenize(';'),
    ntermmetex: params.ntermmetexcision,
    ntermac: params.ntermacetyl,
    precfdr: params.precconflvl,
    protfdr: params.proteinconflvl,
    maxvarmods: params.maxvarmods,
    ms1acc: params.ms1acc,
    ms2acc: params.ms2acc,
    miscleav: params.miscleav,
    window: params.window,
    mincharge: params.mincharge,
    maxcharge: params.maxcharge,
    minmz: params.minmz,
    maxmz: params.maxmz,
    minfrmz: params.minfragmz,
    maxfrmz: params.maxfragmz,
    minpeplen: params.minpeplen,
    maxpeplen: params.maxpeplen,
    nonorm: params.nonorm,
    pglvl: params.proteotypicity,
    idstonames: params.ids_to_names,
    contam: params.contaminants,
    indiacc: params.individual_massacc,
    indiwin: params.individual_windows,
    cut: cut,
  ]
  


  if (params.input) {
    def infiles = identify_info_map(params.input)

    channel.fromList(infiles.collect { k,v -> [k, v]})
    | map { [it[0], it[1].file_path] } // mapkey, file
    | branch { 
      thermo: it[1].extension == 'raw' 
      bruker: it[1].extension == 'd'
      }
    | set { raw_c }
    
    raw_c.thermo
    | concat(raw_c.bruker)
    | set { diann_in }
  
    db_params = channel.fromPath(params.tdb)
      .toList()
      .map { [it, diann_params] }

    diann_in
    | filter { infiles[it[0]].create_lib as Integer == 1 }
    | map { it[1] } // only need file for library, not id
    | toList
    | toList
    | combine(db_params)
    | set { raws_to_emp_lib }

    if (!params.library) {
      passed_lib = channel.empty()

      db_params
      | predictFastaLibrary
      | set { predicted_lib }

      predicted_lib
      | combine(raws_to_emp_lib)
      | createEmpiricalLibrary
      | set { empirical_lib }

    } else {
      predicted_lib = channel.empty()

      passed_lib = channel.fromPath(params.library)
      passed_lib 
      | filter { it.extension == 'speclib' }
      | combine(raws_to_emp_lib)
      | createEmpiricalLibrary
      | concat(passed_lib.filter { it.extension == 'parquet' })
      | set { empirical_lib }

    }

    // If no output params are given, output only the last step, i.e. report
    outputreport = params.outputreport || (!params.outputlib && !params.outputquant && !params.outputreport)
    if (outputreport || params.outputquant) {
  
      // Pre-made quantfiles go into a channel
      if (infiles.findAll { k,v -> v.quantfile }) {
        // First if any quantfiles are specified in inputdef we only take those
       channel.from(infiles.findAll { k,v -> v.quantfile }
           .collect { k,v -> [k, v.file_path, file(v.quantfile)] })
         .concat(infiles.findAll { k,v -> !v.quantfile }
           .collect { k,v -> [k, v.file_path, false] })
         .set { rawquantfiles }
  
      } else if (params.quantdir) {
        // If not and a quantdir is specified, try to match those with the raws by name
        def qfs = []
        file(params.quantdir).traverse(type: FileType.FILES, maxDepth: 0) { qfs.add(it) }
        tmp_q = channel.from(qfs).map { [file(it).baseName, file(it)] }
        channel.from(infiles.collect { k,v ->
              ["${v.file_path.baseName}_${v.file_path.extension}", k, v.file_path] })
          .join(tmp_q, remainder: true)
          .map { [it[1], it[2], it[3]] } // key, rawfile, quantfile/null
          .set { rawquantfiles }
  
      } else {
        // No quantfiles
        channel.from(infiles.collect { k,v ->
              ["${v.file_path.baseName}_${v.file_path.extension}", k, v.file_path] })
          .map { [it[1], it[2], false] } // key, rawfile, false
          .set { rawquantfiles }
      }
  	
      batchsize = params.batchsize ?: infiles.size()
      rawquantfiles
      | filter { !it[2] } // skip rawfiles without quant
      | map { [it[0], it[1]] } // no quant to this proc
      | collate(batchsize)
      | transpose
      | collate(2) // id, raw
      | combine(empirical_lib)
      | combine(db_params)
      | RunDiaAnalysis
      | map { [it[0], listify(it[1]), listify(it[2])] }
      // first sort keys/raw files on raw basename
      | map { it.transpose().sort({a,b -> a[1].baseName <=> b[1].baseName}).transpose() }
      // now sort non-matched quantfiles by basename so they match up
      | map { [it[0], it[1], it[2].sort({a,b -> a.baseName <=> b.baseName})] }
      | transpose
      | set { new_raw_quants }

    } else {
      new_raw_quants = channel.empty()
    }
  
    if (outputreport) {
      // Run training quantUMS and then full experiment
      rawquantfiles
      | filter { it[2] }
      | concat(new_raw_quants)
      | filter { infiles[it[0]].train_quantums as Integer == 1 }
      | map { [it[1], it[2]] }
      | toList
      | transpose
      | toList
      | combine(empirical_lib)
      | combine(db_params)
      | TrainQuantUMS

      rawquantfiles
      | filter { it[2] }
      | concat(new_raw_quants)
      | map { [it[1], it[2]] }
      | toList
      | transpose
      | toList
      | combine(empirical_lib)
      | combine(db_params)
      | combine(TrainQuantUMS.out)
      | DiaQuantificationReport
      | set { reports_out }
    } else {
      reports_out = channel.empty()
    }


  } else if (params.raw) {
    channel.fromPath(params.raw)
    | branch { 
      thermo: it.extension == 'raw' 
      bruker: it.extension == 'd'
      }
    | set { raw_c }
    mzml_c = channel.empty()
  
  } else if (params.mzml) {
    mzml_c = channel.fromPath(mzml)
    raw_c = channel.empty()
  }

// TODO pass emp lib, refine it with another raw file? Is that a usecase?

    predicted_lib.filter { params.output_pred_lib }
    .concat(createEmpiricalLibrary.out.filter { params.output_emp_lib })
    .concat(new_raw_quants.filter { params.outputquant }.map { it[2] })
    .concat(reports_out)
    .set { ch_wfoutputs }

  publish:
  wfoutputs = ch_wfoutputs
}

output {
  wfoutputs { }
}

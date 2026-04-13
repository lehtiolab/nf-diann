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
  path('library.predicted.speclib'), emit: lib
  path('insilico_predict_lib.log'), emit: log

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
    ${diannparams.excl_contam ? "--cont-quant-exclude ${diannparams.excl_contam}" : ''}

    mv library.log.txt insilico_predict_lib.log
    # Is this used in predict from fasta, but maybe test this:
    #${diannparams.window ? "--window $diannparams.window" : ''} \
    #--mass-acc-ms1 $diannparams.ms1acc \
    #--mass-acc $diannparams.ms2acc \
    """
}


process searchWithPredictedLib {
cache 'lenient'

  container "ghcr.io/lehtiolab/nfhelaqc:3.2-diann.2.3.1"
  
  input:
  tuple path(predlib), path(raws), path(fasta), val(diannparams)

  output:
  path('quants/*.quant'), emit: quants
  path('search_predicted_lib.log'), emit: log

  script:
  """
  mkdir quants
  # Empirical library by running it with raws
  diann-linux --threads ${task.cpus} \
    ${raws.collect { "--f \$(realpath $it)" }.join(' ') } \
    ${fasta.collect { "--fasta $it" }.join(' ')} \
    --cut ${diannparams.cut} \
    --lib $predlib \
    --quant-ori-names \
    --temp quants/ \
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
    ${diannparams.excl_contam ? "--cont-quant-exclude ${diannparams.excl_contam}" : ''}

    mv report.log.txt search_predicted_lib.log
  """
}


process combineEmpiricalLibraryRuns {
  container "ghcr.io/lehtiolab/nfhelaqc:3.2-diann.2.3.1"

  input:
  tuple path('quants/*'), path(predlib), path(raws), path(fasta), val(diannparams)
  
  output:
  path('library.parquet'), emit: lib
  path('create_empirical_lib.log'), emit: log
  
  script:
  """
  # Empirical library by running it with raws
  diann-linux --threads ${task.cpus} \
    ${raws.collect { "--f \$(realpath $it)" }.join(' ') } \
    ${fasta.collect { "--fasta $it" }.join(' ')} \
    --use-quant \
    --quant-ori-names \
    --temp quants/ \
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
    ${diannparams.excl_contam ? "--cont-quant-exclude ${diannparams.excl_contam}" : ''}

    mv report.log.txt create_empirical_lib.log
"""
}


process RunDiaAnalysis {

  // "Second pass" after creating the library, manual MBR
  // No reports are output, only quant files
  container "ghcr.io/lehtiolab/nfhelaqc:3.2-diann.2.3.1"

  input:
  tuple val(ids), path(raws), path(lib), path(fasta), val(diannparams)
  
  output:
  tuple val(ids), path(raws), path('quants/*.quant'), emit: rawquants
  path('search_empirical_lib.log'), emit: log

  script:
  """
  mkdir quants
  diann-linux --threads ${task.cpus} \
    ${raws.collect { "--f \$(realpath $it)" }.join(' ') } \
    --lib $lib \
    ${fasta.collect { "--fasta $it" }.join(' ')} \
    --cut ${diannparams.cut} \
    --quant-ori-names \
    --temp quants \
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
    ${diannparams.excl_contam ? "--cont-quant-exclude ${diannparams.excl_contam}" : ''}

    mv report.log.txt search_empirical_lib.log
  """
}


process TrainQuantUMS {
  /* Since we use stdout as output here (to pass parameters to next process),
  we cannot get the entire stdout in the .command.out, because that would put
  the full log into the next process .command.sh instead of only the params
  */

  container "ghcr.io/lehtiolab/nfhelaqc:3.2-diann.2.3.1"

  input:
  tuple path(raws), path('quants/*'), path(lib), path(fasta), val(diannparams)
  
  output:
  stdout emit: params
  path('train_quantums.log'), emit: log

  script:
  paramline = '.*Quantification parameters:' 
  """
  diann-linux --threads ${task.cpus} \
    ${raws.collect { "--f \$(realpath $it)" }.join(' ') } \
    --lib $lib \
    ${fasta.collect { "--fasta $it" }.join(' ')} \
    --cut ${diannparams.cut} \
    --use-quant \
    --quant-ori-names \
    --temp quants \
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
    ${diannparams.excl_contam ? "--cont-quant-exclude ${diannparams.excl_contam}" : ''} \
       | tee stdout.backup | grep '$paramline' | sed 's/$paramline/\\-\\-quant-params/'

    mv report.log.txt train_quantums.log
  """
}

process DiaQuantificationReport {

  container "ghcr.io/lehtiolab/nfhelaqc:3.2-diann.2.3.1"

  input:
  tuple path(raws), path('quants/*'), path(lib), path(fasta), val(diannparams), val(quantparams)
  
  output:
  tuple path('report.parquet'), path('*.tsv'), emit: report
  path('quantify_report.log'), emit: log

  script:
  """
  diann-linux --threads ${task.cpus} \
    ${raws.collect { "--f \$(realpath $it)"}.join(' ') } \
    --lib $lib \
    ${fasta.collect { "--fasta $it" }.join(' ')} \
    --cut ${diannparams.cut} \
    --use-quant \
    --quant-ori-names \
    --temp quants \
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
    ${diannparams.excl_contam ? "--cont-quant-exclude ${diannparams.excl_contam}" : ''}

    mv report.log.txt quantify_report.log
  """
}


process logConcat {

  input:
  path(logs)

  output:
  path('diann-log.txt')
  
  script:
  """
  cat ${logs.join(' ')} > diann-log.txt
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
    excl_contam: params.exclude_contaminants,
    indiacc: params.individual_massacc,
    indiwin: params.individual_windows,
    cut: cut,
  ]
  

  if (params.input) {
    def infiles = identify_info_map(params.input)
    batchsize = params.batchsize && params.batchsize < infiles.size() ? params.batchsize : false

    channel.fromList(infiles.collect { k,v -> [k, v]})
    | map { [it[0], it[1].file_path] } // mapkey, file
    | branch { 
      thermo: it[1].extension == 'raw' 
      bruker: it[1].extension == 'd'
      }
    | set { raw_c }
    
    raw_c.thermo
    | mix(raw_c.bruker)
    | set { diann_in }
  
    db_params = channel.fromPath(params.tdb)
      .toList()
      .map { [it, diann_params] }

    diann_in
    .filter { infiles[it[0]].create_lib as Integer == 1 }
    .map { it[1] } // only need file for library, not id
    .set { all_raws_to_emp_lib }

    if (batchsize) {
      // Collate into batches if applicable
      all_raws_to_emp_lib
      .collate(batchsize)
      .set { list_of_raws }
    } else {
      // Run all in same proc, we need double toList so it will be a single item in the combine
      all_raws_to_emp_lib
      .toList()
      .toList()
      .set { list_of_raws }
    }
    list_of_raws
    .combine(db_params)
    .set { batched_raws_to_emp_lib }

    
    if (!params.library) {
      passed_lib = channel.empty()
      predicted_lib = predictFastaLibrary(db_params)

      predicted_lib.lib
      | combine(batched_raws_to_emp_lib)
      | searchWithPredictedLib

      searchWithPredictedLib.out.quants
      | flatten | toList | toList // combine all quant files in one big list
      | combine(predicted_lib.lib)
      | combine(all_raws_to_emp_lib.toList().toList())
      | combine(db_params)
      | combineEmpiricalLibraryRuns
      combineEmpiricalLibraryRuns.out.lib
      | set { empirical_lib }

    } else {
      predicted_lib = channel.empty().branch {
        lib: it==1
        log: it==2
      }

      passed_lib = channel.fromPath(params.library)
      passed_lib 
      | filter { it.extension == 'speclib' }
      | combine(batched_raws_to_emp_lib)
      | searchWithPredictedLib
      searchWithPredictedLib.out.quants
      | mix(passed_lib.filter { it.extension == 'parquet' })
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
         .mix(infiles.findAll { k,v -> !v.quantfile }
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
  	
      rawquantfiles
      .filter { !it[2] } // skip rawfiles without quant
      .map { [it[0], it[1]] } // no quant to this proc
      .set { raw_without_q }

      if (batchsize) {
      	raw_without_q.collate(batchsize).set { list_of_raw_wo_q }
      } else {
        // toList to make this work in the combine step which will not flatten the raws
      	raw_without_q.toList().set { list_of_raw_wo_q }
      }

      list_of_raw_wo_q
      | transpose
      | collate(2) // id, raw
      | combine(empirical_lib)
      | combine(db_params)
      | RunDiaAnalysis
      RunDiaAnalysis.out.rawquants
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
      | mix(new_raw_quants)
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
      | mix(new_raw_quants)
      | map { [it[1], it[2]] }
      | toList
      | transpose
      | toList
      | combine(empirical_lib)
      | combine(db_params)
      | combine(TrainQuantUMS.out.params)
      | DiaQuantificationReport
      DiaQuantificationReport.out.report
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

  // sort logs
  predicted_lib.log.map { [0, it] }
  | mix(searchWithPredictedLib.out.log.map { [1, it] })
  | mix(combineEmpiricalLibraryRuns.out.log.map { [2, it] })
  | mix(RunDiaAnalysis.out.log.map { [3, it] })
  | mix(TrainQuantUMS.out.log.map { [4, it] })
  | mix(DiaQuantificationReport.out.log.map { [5, it] })
  | toList
  | map { it.sort({a,b -> a[0] <=> b[0]}) }
  | transpose
  | toList
  | map { it[1] } // remove indices after sorting
  | logConcat
  | mix(predicted_lib.lib.filter { params.output_pred_lib })
  | mix(combineEmpiricalLibraryRuns.out.lib.filter { params.output_emp_lib })
  | mix(new_raw_quants.filter { params.outputquant }.map { it[2] })
  | mix(reports_out)
  | set { ch_wfoutputs }

  publish:
  wfoutputs = ch_wfoutputs
}

output {
  wfoutputs { }
}

#!/usr/bin/env nextflow

include { paramsSummaryMap } from 'plugin/nf-schema'

include { identify_info_map; listify; read_header } from './modules.nf' 
include { QC_REPORT } from './workflows/create_report.nf'


process predictFastaLibrary {
cache 'lenient'

  tag 'diann'
  container Containers.containers[task.tag][workflow.containerEngine]
  
  input:
  tuple path(fasta, arity: '1..*'), val(diannparams)

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
    --min-pr-charge ${diannparams.mincharge} \
    --max-pr-charge ${diannparams.maxcharge} \
    --min-pep-len ${diannparams.minpeplen} \
    --max-pep-len ${diannparams.maxpeplen} \
    --min-pr-mz ${diannparams.minmz} \
    --max-pr-mz ${diannparams.maxmz} \
    --min-fr-mz ${diannparams.minfrmz} \
    --max-fr-mz ${diannparams.maxfrmz} \
    --mass-acc-ms1 $diannparams.ms1acc \
    --mass-acc $diannparams.ms2acc \
    ${diannparams.excl_contam ? "--cont-quant-exclude ${diannparams.excl_contam}" : ''} \
      | tee stdout.bak
    grep ERROR stdout.bak && exit 1

    mv library.log.txt insilico_predict_lib.log
    # Is this used in predict from fasta, but maybe test this:
    """
}


process searchWithPredictedLib {
cache 'lenient'

  tag 'diann'
  container Containers.containers[task.tag][workflow.containerEngine]
  
  input:
  tuple path(predlib), path(raws, arity: '1..*'), path(fasta), val(diannparams)

  output:
  path('quants/*.quant'), emit: quants
  path("${raws[0]}_search_predicted_lib.log"), emit: log

  script:
  raws_actual = []
  raws.each { it.isDirectory() ?  it.toRealPath().eachFileMatch('analysis.tdf_bin') { raws_actual << it } : raws_actual << it   }
  
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
    ${diannparams.excl_contam ? "--cont-quant-exclude ${diannparams.excl_contam}" : ''} \
      | tee stdout.bak
    grep ERROR stdout.bak && exit 1

    mv report.log.txt ${raws[0]}_search_predicted_lib.log
  """
}


process combineEmpiricalLibraryRuns {

  tag 'diann'
  container Containers.containers[task.tag][workflow.containerEngine]

  input:
  tuple path('quants/*'), path(predlib), path(raws, arity: '1..*'), path(fasta), val(diannparams)
  
  output:
  path('library.parquet'), emit: lib
  path('create_empirical_lib.log'), emit: log
  
  script:
  raws_actual = []
  raws.each { it.isDirectory() ?  it.toRealPath().eachFileMatch('analysis.tdf_bin') { raws_actual << it } : raws_actual << it }

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
    ${diannparams.excl_contam ? "--cont-quant-exclude ${diannparams.excl_contam}" : ''} \
      | tee stdout.bak
    grep ERROR stdout.bak && exit 1

    mv report.log.txt create_empirical_lib.log
"""
}


process RunDiaAnalysis {

  // "Second pass" after creating the library, manual MBR
  // No reports are output, only quant files
  tag 'diann'
  container Containers.containers[task.tag][workflow.containerEngine]

  input:
  tuple val(ids), path(raws, arity: '1..*'), path(lib), path(fasta), val(diannparams)
  
  output:
  tuple val(ids), path(raws), path('quants/*.quant'), emit: rawquants
  path("${raws[0]}_search_empirical_lib.log"), emit: log

  script:
  raws_actual = []
  raws.each { it.isDirectory() ?  it.toRealPath().eachFileMatch('analysis.tdf_bin') { raws_actual << it } : raws_actual << it }
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
    ${diannparams.excl_contam ? "--cont-quant-exclude ${diannparams.excl_contam}" : ''} \
      | tee stdout.bak
    grep ERROR stdout.bak && exit 1

    mv report.log.txt ${raws[0]}_search_empirical_lib.log
  """
}


process TrainQuantUMS {
  /* Since we use stdout as output here (to pass parameters to next process),
  we cannot get the entire stdout in the .command.out, because that would put
  the full log into the next process .command.sh instead of only the params
  */

  tag 'diann'
  container Containers.containers[task.tag][workflow.containerEngine]

  input:
  tuple path(raws, arity: '1..*'), path('quants/*'), path(lib), path(fasta), val(diannparams)
  
  output:
  stdout emit: params
  path('train_quantums.log'), emit: log

  script:
  raws_actual = []
  raws.each { it.isDirectory() ?  it.toRealPath().eachFileMatch('analysis.tdf_bin') { raws_actual << it } : raws_actual << it }
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

    grep ERROR stdout.backup && exit 1
    mv report.log.txt train_quantums.log
  """
}

process DiaQuantificationReport {

  tag 'diann'
  container Containers.containers[task.tag][workflow.containerEngine]

  input:
  tuple path(raws, arity: '1..*'), path('quants/*'), path(lib), path(fasta), val(diannparams), val(quantparams), val(enzyme)
  
  output:
  tuple path('report.parquet'), path('*.tsv'), emit: report
  path('precursors.txt'), emit: precursors
  path('quantify_report.log'), emit: log

  script:
  raws_actual = []
  raws.each { it.isDirectory() ?  it.toRealPath().eachFileMatch('analysis.tdf_bin') { raws_actual << it } : raws_actual << it }
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
    ${diannparams.excl_contam ? "--cont-quant-exclude ${diannparams.excl_contam}" : ''} \
      | tee stdout.bak
    grep ERROR stdout.bak && exit 1
    parquet_to_tsv.py report.parquet $enzyme

    mv report.log.txt quantify_report.log
  """
}


process logConcat {
  tag 'local'

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

  if (params.library && file(params.library).extension == 'speclib' && params.output_pred_lib) {
    exit 1, 'Cannot output new predicted library while also being passed a --library'
  } else if (params.library && file(params.library).extension == 'parquet' && params.output_emp_lib) {
    exit 1, 'Cannot output new empirical library while also being passed a --library'
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

    raw_c = channel.fromList(infiles.collect { k,v -> [k, v]})
    .map { [it[0], it[1].file_path] } // mapkey, file
    .branch { 
      thermo: it[1].extension == 'raw' 
      bruker: it[1].extension == 'd'
      }
    
    diann_in = raw_c.thermo.mix(raw_c.bruker)
  
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

      emplib_in = predicted_lib.lib.combine(batched_raws_to_emp_lib)
      searchWithPredictedLib(emplib_in)

      all_emp_libs = searchWithPredictedLib.out.quants
      .flatten()
      .toList()
      .toList() // combine all quant files in one big list
      .combine(predicted_lib.lib)
      .combine(all_raws_to_emp_lib.toList().toList())
      .combine(db_params)
      combineEmpiricalLibraryRuns(all_emp_libs)
      empirical_lib = combineEmpiricalLibraryRuns.out.lib

    } else {
      predicted_lib = channel.empty().branch {
        lib: it==1
        log: it==2
      }

      passed_lib = channel.fromPath(params.library)
      predlib_in = passed_lib 
      .filter { it.extension == 'speclib' }
      .combine(batched_raws_to_emp_lib)
      searchWithPredictedLib(predlib_in)
      empirical_lib = searchWithPredictedLib.out.quants
      .mix(passed_lib.filter { it.extension == 'parquet' })
    }

    // If no output params are given, output only the last step, i.e. report
    outputreport = params.outputreport || (!params.output_pred_lib && !params.output_emp_lib && !params.outputquant && !params.outputreport)
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
        file(params.quantdir).traverse(type: groovy.io.FileType.FILES, maxDepth: 0) { qfs.add(it) }
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

      rundia_in = list_of_raw_wo_q
      .transpose()
      .collate(2) // id, raw
      .combine(empirical_lib)
      .combine(db_params)
      RunDiaAnalysis(rundia_in)
      new_raw_quants = RunDiaAnalysis.out.rawquants
      .map { [it[0], listify(it[1]), listify(it[2])] }
      // first sort keys/raw files on raw basename
      .map { it.transpose().sort({a,b -> a[1].baseName <=> b[1].baseName}).transpose() }
      // now sort non-matched quantfiles by basename so they match up
      .map { [it[0], it[1], it[2].sort({a,b -> a.baseName <=> b.baseName})] }
      .transpose()

    } else {
      new_raw_quants = channel.empty()
    }
  
    if (outputreport) {
      // Run training quantUMS and then full experiment
      trainq_in = rawquantfiles
      .filter { it[2] }
      .mix(new_raw_quants)
      .filter { infiles[it[0]].train_quantums as Integer == 1 }
      .map { [it[1], it[2]] }
      .toList()
      .transpose()
      .toList()
      .combine(empirical_lib)
      .combine(db_params)
      TrainQuantUMS(trainq_in)

      qreport_in = rawquantfiles
      .filter { it[2] }
      .mix(new_raw_quants)
      .map { [it[1], it[2]] }
      .toList()
      .transpose()
      .toList()
      .combine(empirical_lib)
      .combine(db_params)
      .combine(TrainQuantUMS.out.params)
      .map { it + [params.enzyme] }
      DiaQuantificationReport(qreport_in)

      raws_ftypes = raw_c.bruker
      .map { [it[1], 'bruker'] }
      .concat(raw_c.thermo
        .map { [it[1], 'thermo'] }
      )
      input_to_qc = channel.fromPath(params.input) 
      QC_REPORT(raws_ftypes, input_to_qc, DiaQuantificationReport.out.precursors, params.proteinconflvl)


      reports_out = DiaQuantificationReport.out.report
    } else {
      reports_out = channel.empty()
    }


  } else if (params.raw) {
    raw_c = channel.fromPath(params.raw)
    .branch { 
      thermo: it.extension == 'raw' 
      bruker: it.extension == 'd'
      }
    mzml_c = channel.empty()
  
  } else if (params.mzml) {
    mzml_c = channel.fromPath(params.mzml)
    raw_c = channel.empty()
  }

  // sort logs
  ch_log_outputs = predicted_lib.log.map { [0, it] }
  .mix(searchWithPredictedLib.out.log.map { [1, it] })
  .mix(combineEmpiricalLibraryRuns.out.log.map { [2, it] })
  .mix(RunDiaAnalysis.out.log.map { [3, it] })
  .mix(TrainQuantUMS.out.log.map { [4, it] })
  .mix(DiaQuantificationReport.out.log.map { [5, it] })
  .toList()
  .map { it.sort({a,b -> a[0] <=> b[0]}) }
  .transpose()
  .toList()
  .map { it[1] } // remove indices after sorting

  ch_wfoutputs = logConcat(ch_log_outputs)
  .mix(QC_REPORT.out)
  .mix(predicted_lib.lib.filter { params.output_pred_lib })
  .mix(combineEmpiricalLibraryRuns.out.lib.filter { params.output_emp_lib })
  .mix(new_raw_quants.filter { params.outputquant }.map { it[2] })
  .mix(reports_out)

  publish:
  wfoutputs = ch_wfoutputs

  onComplete:
  if (workflow.success) {
    def libfile = file("${workflow.outputDir}/libs.js")
    def libs = libfile.readLines()
    def bulma = file("${baseDir}/assets/bulma.js").readLines()
    def psmap = paramsSummaryMap(workflow)

    if (params.input && file(params.input).exists()) {
      inputraws = identify_info_map(params.input)
      def files_header = read_header(params.input)
      files_header[0] = 'file_path'
      infiles = inputraws.collect { k, fn -> files_header.collect { fn[it] }}
      infiles.add(0, files_header)
    } else {
      infiles = [[], []]
    }

    // Get processes tags used from trace to output the software versions used in pipeline
    def x = file("execution_trace.txt").readLines()[1..-1]
      .collect { it.tokenize('\t')[4] } // get process tag
      .grep { it != 'local' }
      .unique()
    def sw_versions = file("execution_trace.txt").readLines()[1..-1]
      .collect { it.tokenize('\t')[4] } // get process tag
      .grep { it != 'local' }
      .unique()
      .collect { [it, Containers.containers[it].version, Containers.containers[it][workflow.containerEngine]] }
    // The above crashes this handler in case the tag (software) is not defined in the lib/containers.groovy
      
    // Set the name of the workflow
    // if not wf.runName (-name or auto) is like "crazy_euler" or other "{adjective}_{scientist}"
    if (!params.name && !(workflow.runName ==~ /[a-z]+_[a-z]+/) ) {
      runname = workflow.runName
    } else if (!params.name) {
      runname = 'untitled'
    } else {
      runname = params.name
    }
    def fields = [runname: runname,
        sw_versions: sw_versions,
        params: psmap['Other parameters'],
        infiles: infiles,
        libs: libs, bulma: bulma]
    def rf = new File("${workflow.outputDir}/report_groovy_template.html")
    def temp_engine = new groovy.text.StreamingTemplateEngine()
    def report_template = temp_engine.createTemplate(rf).make(fields)
    def report_html = report_template.toString()
    def output_rf = new File( "${workflow.outputDir}/report.html" )
    output_rf.withWriter { w -> w << report_html }
    rf.delete()
    libfile.delete()
  }
}

output {
  wfoutputs { }
}

#!/usr/bin/env nextflow

include { paramsSummaryMap } from 'plugin/nf-schema'

import groovy.io.FileType

include { identify_info_map; listify } from './modules.nf' 
//include { REPORTING } from './workflows/reporting.nf'



process createLibrary {
cache 'lenient'

  container "ghcr.io/lehtiolab/nfhelaqc:3.2-diann.2.3.1"
  
  input:
  tuple path(raws), path(fasta), val(diannparams)

  output:
  path('library.parquet')

  script:
  """
  # Create predicted library from fasta
  diann-linux --threads ${task.cpus} \
    --fasta $fasta \
    --gen-spec-lib \
    --predictor \
    --fasta-search \
    --out-lib firstlib \
    --missed-cleavages $diannparams.miscleav \
    --mass-acc-ms1 $diannparams.ms1acc \
    --mass-acc $diannparams.ms2acc \
    --window $diannparams.window \
    ${diannparams.varmods.collect { "--var-mod $it" }.join(' ')} \
    ${diannparams.fixmods.collect { "--fixed-mod $it" }.join(' ')} \
    --var-mods $diannparams.maxvarmods
  
  # Empirical library by running it with raws
  diann-linux --threads ${task.cpus} \
    ${raws.collect { "--f $it" }.join(' ') } \
    --fasta $fasta \
    --lib firstlib.predicted.speclib \
    --gen-spec-lib \
    --out-lib library \
    --missed-cleavages $diannparams.miscleav \
    --mass-acc-ms1 $diannparams.ms1acc \
    --mass-acc $diannparams.ms2acc \
    --window $diannparams.window \
    ${diannparams.varmods.collect { "--var-mod $it" }.join(' ')} \
    ${diannparams.fixmods.collect { "--fixed-mod $it" }.join(' ')} \
    --var-mods $diannparams.maxvarmods
  """
}


process RunDiaAnalysis {

  // "Second pass" after creating the library, manual MBR
  // No reports are output, only quant files
  // q-value cutoff seems 0.01 by default
  container "ghcr.io/lehtiolab/nfhelaqc:3.2-diann.2.3.1"

  input:
  tuple val(ids), path(raws), path(lib), path(tdb), val(diannparams)
  
  output:
  //val('fake'), emit: fake
  tuple val(ids), path(raws), path('*.quant')
  //tuple path('out.txt'), path('out.stats.tsv'), emit: tsv
  //path('out.parquet'), emit: pq

  script:
  """
  diann-linux --threads ${task.cpus} \
    ${raws.collect { "--f $it" }.join(' ') } \
    --lib $lib \
    --fasta $tdb \
    --temp ./ \
    --missed-cleavages $diannparams.miscleav \
    --mass-acc-ms1 $diannparams.ms1acc \
    --mass-acc $diannparams.ms2acc \
    --window $diannparams.window \
    ${diannparams.varmods.collect { "--var-mod $it" }.join(' ')} \
    ${diannparams.fixmods.collect { "--fixed-mod $it" }.join(' ')} \
    --var-mods $diannparams.maxvarmods
  """
}


process TrainQuantUMS {

  // q-value cutoff seems 0.01 by default
  
  container "ghcr.io/lehtiolab/nfhelaqc:3.2-diann.2.3.1"

  input:
  tuple path(raws), path(quants), path(lib), path(tdb), val(diannparams)
  
  output:
  stdout

  script:
  paramline = '.*Quantification parameters:' 
  """
  diann-linux --threads ${task.cpus} \
    ${raws.collect { "--f \$(realpath $it)"}.join(' ') } \
    --lib $lib \
    --fasta $tdb \
    --use-quant \
    --temp ./ \
    --quant-train-runs 0:${listify(raws).size() -1} \
    --missed-cleavages $diannparams.miscleav \
    --mass-acc-ms1 $diannparams.ms1acc \
    --mass-acc $diannparams.ms2acc \
    --window $diannparams.window \
    ${diannparams.varmods.collect { "--var-mod $it" }.join(' ')} \
    ${diannparams.fixmods.collect { "--fixed-mod $it" }.join(' ')} \
    --var-mods $diannparams.maxvarmods \
       > tmpstdout.log
    grep '$paramline' tmpstdout.log | sed 's/$paramline/\\-\\-quant-params/'
  """
}

process DiaQuantificationReport {

  // q-value cutoff seems 0.01 by default
  
  //container 'michelmoser/diann-1.9.2'
  
  container "ghcr.io/lehtiolab/nfhelaqc:3.2-diann.2.3.1"

  input:
  tuple path(raws), path(quants), path(lib), path(tdb), val(diannparams), val(quantparams)
  
  output:
  path('report.parquet')

  script:
  """
  diann-linux --threads ${task.cpus} \
    ${raws.collect { "--f \$(realpath $it)"}.join(' ') } \
    --lib $lib \
    --fasta $tdb \
    --use-quant \
    --temp ./ \
    ${quantparams.trim()} \
    --missed-cleavages $diannparams.miscleav \
    --mass-acc-ms1 $diannparams.ms1acc \
    --mass-acc $diannparams.ms2acc \
    --window $diannparams.window \
    ${diannparams.varmods.collect { "--var-mod $it" }.join(' ')} \
    ${diannparams.fixmods.collect { "--fixed-mod $it" }.join(' ')} \
    --var-mods $diannparams.maxvarmods \
  """
}


    //--var-mod 'UniMod:35,15.994915,M' \
    //--var-mod 'UniMod:4,57.021464,C' \

workflow {
  main:
  // FIXME
//  ms1acc = [timstof: 20, velos: 10, qe: 10, astral: 10][instrument]
//  ms2acc = 20

  if (params.library && params.outputlib) {
    exit 1, 'Cannot output generated library while also being passed a --library'
  }

  diann_params = [
    varmods: params.varmods.tokenize(';'),
    fixmods: params.fixmods.tokenize(';'),
    maxvarmods: params.maxvarmods,
    ms1acc: params.ms1acc,
    ms2acc: params.ms2acc,
    miscleav: params.miscleav,
    window: params.window,
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
  
    db_params = channel.fromPath(params.tdb).map { [it, diann_params] }

    if (!params.library) {
      diann_in
      | filter { infiles[it[0]].create_lib as Integer == 1 }
      | map { it[1] } // only need file for library, not id
      | toList
      | toList
      | combine(db_params)
      | createLibrary
      | set { library }
    } else {
      library = channel.fromPath(params.library)
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
      | combine(library)
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
      | combine(library)
      | combine(db_params)
      | TrainQuantUMS

      rawquantfiles
      | filter { it[2] }
      | concat(new_raw_quants)
      | map { [it[1], it[2]] }
      | toList
      | transpose
      | toList
      | combine(library)
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

  library.filter { !params.library && params.outputlib }
    .concat(new_raw_quants.filter { params.outputquant }.map { it[2] })
    .concat(reports_out)
    .set { ch_wfoutputs }

  publish:
  wfoutputs = ch_wfoutputs
}

output {
  wfoutputs { }
}

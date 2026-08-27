process extractThermoScans {
  /* Used in DIA for Thermo, so we dont need to convert to mzML
  just to get the nr of scans */
  tag 'scanheadsman'
  container Containers.containers[task.tag][workflow.containerEngine]

  input:
  path(raws)

  output:
  path('nrscans')

  script:
  """
  exitcode=0
  ${raws.collect { "timeout --preserve-status 5s dotnet /scanheadsman/ScanHeadsman.dll \$(realpath $it) > tmpfn || exitcode=\$? \
    && if [[ \$exitcode != 0 && \$exitcode != 143 ]] ; then exit \$exitcode ;fi \
    && echo ${it.baseName}\$'\t'\$(grep 'Processing scan [0-9]' tmpfn | sed 's/.* of //') >> nrscans"}.join('\n')}
  """
}


process getBrukerScanNumbers {
  tag 'sqlite'
  container Containers.containers[task.tag][workflow.containerEngine]

  input:
  path(raws)

  output:
  path('nrscans')

  script:
  """
  ${raws.collect { "echo ${it.baseName}\$'\t'\$(sqlite3 ${it}/analysis.tdf 'SELECT COUNT(*) FROM Frames') >> nrscans"}.join(' && ')}
  """
}


process parquetToTsv {
  // Precursor table output, from parquet to tsv, could possibly run directly
  // inside DIA-NN process

  tag 'diann'
  container Containers.containers[task.tag][workflow.containerEngine]

  input:
  tuple path(report), val(enzyme)

  output:
  path('precursors.txt')

  script:
  """
  parquet_to_tsv.py $report $enzyme
  """

}


process precursorPlot {
  tag 'ddamsproteomics'
  container Containers.containers[task.tag][workflow.containerEngine]

  input:
  tuple path('filescans????'), path(precursors), path(inputfn), val(conflvl)
  
  output:
  tuple path('precursorplothtml'), path('*_qc.txt'), path('*__overlap'), path('genesplothtml'), path('proteinsplothtml')
  
  script:
  // FIXME error if not finding these columns!
  """
  cat filescans* > concat_filescans
  mkdir -p precursorplothtml genesplothtml proteinsplothtml
  precursor_qc.R --precursors $precursors --inputfn $inputfn --conflvl $conflvl
  """
}


process summaryReport {
cache false

  tag 'ddamsproteomics'
  container Containers.containers[task.tag][workflow.containerEngine]


  input:
  tuple path('precplothtml'), path(summaries), path(feat_overlaps), path('genesplots'), path('proteinsplots')
  
  //tuple path(platescans), path(plotlibs), path('psmplots'), path(psm_summary), path('psmids'), path('miscleav'), path(featplots), path(feat_summaries), path(feat_overlaps), path('ptmplots'), path(ptmfiles), path('warnings*')
  
  output:
  tuple path('report_groovy_template.html'), path('libs.js')
  
  script:
  """
  # xargs removes trailing whitespace
  report_tables.py --version "${workflow.manifest.version}" --doi "${workflow.manifest.doi}" \
      --templatedir "$baseDir/assets"
  """
}


workflow QC_REPORT {
  take:
  raws_ftypes
  inputfn
  precursors
  proteinfdr
  
  main:
  raws_ftypes
  | filter { it[1] == 'thermo' }
  | map { it[0] }
  | collate(10) // split raw reading in batches of 10
  | filter { it.size() > 0 }
  | extractThermoScans
  
  raws_ftypes
  | filter { it[1] == 'bruker' }
  | map { it[0] }
  | collate(10) // split raw reading in batches of 10
  | filter { it.size() > 0 }
  | getBrukerScanNumbers
  | mix(extractThermoScans.out)
  | combine(precursors)
  | combine(inputfn)
  | map { it + [proteinfdr] }
  | precursorPlot
  | summaryReport
  
  emit:
  summaryReport.out
}

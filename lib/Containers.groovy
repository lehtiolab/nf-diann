class Containers {

  public static containers = [
    diann: [
      version: '2.3.1',
      singularity: 'docker://ghcr.io/lehtiolab/nf-diann:0.2',
      docker: 'ghcr.io/lehtiolab/nf-diann:0.2',
    ],
    ddamsproteomics: [
      version: '3.0',
      docker: 'ghcr.io/lehtiolab/ddamsproteomics:3.0',
      singularity: 'docker://ghcr.io/lehtiolab/ddamsproteomics:3.0',
    ],
    scanheadsman: [
      version: '3.0',
      docker: 'ghcr.io/lehtiolab/kantele_nf_thermoreader:1.4',
      singularity: 'docker://ghcr.io/lehtiolab/kantele_nf_thermoreader:1.4',
    ],
    sqlite: [
      version: '3.33.0',
      docker: 'quay.io/biocontainers/sqlite:3.33.0',
      singularity: 'https://depot.galaxyproject.org/singularity/sqlite:3.33.0',
    ],
  ]
}

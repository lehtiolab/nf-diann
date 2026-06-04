class Containers {

  public static containers = [
    diann: [
      version: '2.3.1',
      singularity: 'oras://ghcr.io/lehtiolab/nfhelaqc:3.2-diann.2.3.1',
      docker: 'ghcr.io/lehtiolab/nfhelaqc:3.2-diann.2.3.1',
    ],
    ddamsproteomics: [
      version: '3.0',
      docker: 'ghcr.io/lehtiolab/ddamsproteomics:3.0',
      singularity: 'docker://ghcr.io/lehtiolab/ddamsproteomics:3.0',
    ],
    scanheadsman: [
      version: '3.0',
      docker: 'ghcr.io/lehtiolab/kantele_thermoreader:latest',
      singularity: 'docker://ghcr.io/lehtiolab/kantele_thermoreader:latest',
    ],
    sqlite: [
      version: '3.33.0',
      docker: 'quay.io/biocontainers/sqlite:3.33.0',
      singularity: 'https://depot.galaxyproject.org/singularity/sqlite:3.33.0',
    ],
  ]
}

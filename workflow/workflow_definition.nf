

params.log_file = "nextflow_process_log.txt"
params.pipe_script = "/call_pipeline_func.py"


process extract_mgf {
    label 'proteoWizard'

    memory '10 GB'
    cpus 4
  
    input:
    path(raw_file)

    output:
    tuple(path(raw_file), path("${raw_file.baseName}.mgf"))

    script:
    """
    /usr/bin/wine64_anyuser msconvert ${raw_file} --mgf --filter "peakPicking true 1-"
    """
}


process extract_mzml {
    label 'proteoWizard'

    memory '10 GB'
    cpus 4

    input:
    path(raw_file)

    output:
    tuple(path(raw_file), path("${raw_file.baseName}.mzML"))

    script:
    """
    /usr/bin/wine64_anyuser msconvert ${raw_file} --mzML --filter "peakPicking true 1-"
    """
}


// It could potentially be better to have sage run on multiple files per VM instance, if it's actually so fast.
process run_sage {
    label 'sageSearch'

    memory '32 GB'
    cpus 8

    input:
    tuple(val(raw_file), path(mzml_file))
    path(fasta_file)
    path(sage_parameters)

    output:
    tuple(val(raw_file.baseName), val(raw_file), path("${raw_file.baseName}.sage.pin"))
    // The duplication of raw_file is needed to make the join of sage and comet results work; the full filepath
    // seen by each process is different, so .baseName is used to get a consistent join key, but the full path
    // needs to be retained somehow so it can be used by the annotation step. It's a bit of a hack, to say the least.

    script:
    """
    /app/sage -f ${fasta_file} -o \$(pwd) ${sage_parameters} --write-pin ${mzml_file}
    mv results.sage.pin ${raw_file.baseName}.sage.pin
    """
}


process run_sage_wide {
    label 'sageSearch'

    memory '32 GB'
    cpus 8

    input:
    tuple(val(raw_file), path(mzml_file))
    path(fasta_file)
    path(sage_parameters)

    output:
    tuple(val(raw_file.baseName), val(raw_file), path("${raw_file.baseName}.widetol_sage.pin"))
    // The duplication of raw_file is needed to make the join of sage and comet results work; the full filepath
    // seen by each process is different, so .baseName is used to get a consistent join key, but the full path
    // needs to be retained somehow so it can be used by the annotation step. It's a bit of a hack, to say the least.

    script:
    """
    /app/sage -f ${fasta_file} -o \$(pwd) ${sage_parameters} --write-pin ${mzml_file}
    mv results.sage.pin ${raw_file.baseName}.widetol_sage.pin
    """
}


process run_comet {
    label 'cometImage'

    memory '32 GB'
    cpus 8

    input:
    tuple(val(raw_file), path(mgf_file))
    path(fasta_file)
    path(comet_parameters)
    path(comet_varmods)

    output:
    tuple(val(raw_file.baseName), path("${mgf_file.baseName}.pin"), path("${mgf_file.baseName}.txt"))

    script:
    """
    python /comet_search.py ${mgf_file} ${fasta_file} --config ${comet_parameters} --varmods ${comet_varmods}
    """
}

process run_dinosaur {
    label 'javaImage'

    memory '10 GB'
    cpus 1  // Dinosaur doesn't seem to be multithreaded, even with the "concurrency" setting.

    time '30m'

    input:
    tuple(path(rawfile), path(mzml_file))

    output:
    tuple(val(rawfile.baseName), path("${mzml_file.baseName}.features.tsv"))
    // There's also .hills.csv, .msInspect.tsv, .mzq, and .features.bin ... whatever those are ... in case they're useful?

    //script:
    //"""
    //filesize=\$(stat --printf="%s" ${mzml_file})
    //filesize_mb=\$((filesize / 1024 / 1024 + 100))  
    //mkdir /ramdisk
    //mount -t tmpfs -o size=\${filesize_mb}M tmpfs /ramdisk 
    //cp ${mzml_file} /ramdisk

    //java -jar /Dinosaur-1.2.1.free.jar --verbose --profiling /ramdisk/${mzml_file}

    //mv /ramdisk/${mzml_file.baseName}.features.tsv ${mzml_file.baseName}.features.tsv
    //"""
    script:
    """
    java -jar /Dinosaur-1.2.1.free.jar --profiling ${mzml_file}
    """
}

process reconcile_search_results {
    label 'reconcileImage'

    memory '6 GB'
    cpus 1
    disk '30 GB'

    input:
    tuple(val(join_key), path(pin_file), path(txt_file), val(raw_file), path(sage_file), val(raw_file_2), path(widesage_file)) 

    output:
    tuple(val(join_key), val(raw_file), path("${pin_file.baseName}.reconciled.pin"))

    script:
    """
    python /search_output_reconciliation.py ${pin_file} ${txt_file} ${sage_file} ${widesage_file}
    """
}

process annotate_psm_file {
    label 'annotateImage'

    memory '32 GB'
    cpus 2  

    input:
    tuple(val(join_key), path(raw_file), path(psm_file), path(dino_file))

    output:
    path("${psm_file.baseName}.annotated.pin")
    path("${psm_file.baseName}.annotated.pin.resource_monitoring.txt")

    script:
    """
    python /annotate_psms_via_models.py ${psm_file} ${raw_file} ${dino_file} ${psm_file.baseName}.annotated.pin
    """
}


process concatenate_and_mokapot {
    label 'annotateImage'

    memory '64 GB'
    cpus 4
    disk '1000 GB'

    input:
    path file_list
    path fasta_file

    output:
    path 'concatenated.pin'
    path 'mokapot.psms.txt'
    path 'mokapot.peptides.txt'
    path 'mokapot.proteins.txt'

    script:
    """
    python /run_mokapot.py --fasta ${fasta_file} --psm ${file_list.join(' ')}
    """
}

process generate_psm_output {
    label 'postprocessImage'

    time '1h'
    memory '16 GB'
    cpus 2

    input:
    path pin_file
    path mokapot_psms

    output:
    path("${pin_file.baseName}.report.csv")

    script:
    """
    echo "foobar"
    python /generate_final_report.py -p ${pin_file} -m ${mokapot_psms}
    """
}

process OLD_dinosaur_annotate {
    label 'postprocessImage'

    memory '6 GB'
    cpus 1

    input:
    path psm_file
    path dino_file_list

    output:
    path("${psm_file}.w_features.csv")

    script:
    """
    python /dinosaur_quantitation.py --psm ${psm_file} --data ${dino_file_list.join(' ')}
    """
}

process extract_isobaric_values {
    label 'postprocessImage'

    memory '6 GB'
    cpus 1

    input:
    tuple(val(raw_file), path(mzml_file))

    output:
    path("${mzml_file}.quants.csv")

    script:
    """
    python /extract_isobaric_quant.py extract --mzml_file ${mzml_file} --output_file ${mzml_file}.quants.csv --channel_mode ${params.isobaric_type} --tol 0.0015 
    """
}

process integrate_isobaric_data {
    label 'postprocessImage'

    memory '6 GB'
    cpus 1

    input:
    path psm_file
    path isobaric_outputs

    output:
    path("${psm_file}.w_isobaric.csv")

    script:
    """
    python /extract_isobaric_quant.py integrate --psm_file ${psm_file} --output_file ${psm_file}.w_isobaric.csv --isobaric_datafiles ${isobaric_outputs.join(' ')}
    """
}

process peptide_protein_aggregation {
    label 'postprocessImage'

    memory '12 GB'
    cpus 2

    publishDir params.run_name

    input:
    path psm_file
    path mokapot_peptides
    path mokapot_proteins

    output:
    path("${psm_file}.peptides.csv")
    path("${psm_file}.peptides.csv.proteins.csv")

    script:
    """
    python /peptide_protein_aggregation.py --psm ${psm_file} --mokapot ${mokapot_peptides} --mokapot_prots ${mokapot_proteins} --run_parse_regexp "${params.run_parse_regexp}"
    """
}


workflow {
    params.comet_varmods
    params.comet_parameters
    params.sage_parameters
    params.sage_widewindow_parameters
    params.fasta_file
    params.data_glob

    rawfiles = Channel.fromPath(params.data_glob)
  
    rawfiles = rawfiles.multiMap{it -> to_comet: to_sage: it}

    mgf_out = extract_mgf(rawfiles.to_comet)
    comet_out = run_comet(mgf_out, params.fasta_file, params.comet_parameters, params.comet_varmods)

    mzml_out = extract_mzml(rawfiles.to_sage)
    mzml_out = mzml_out.multiMap{it -> to_sage: to_isobar: to_dinosaur: to_widesage: it}
    sage_out = run_sage(mzml_out.to_sage, params.fasta_file, params.sage_parameters)
    widesage_out = run_sage_wide(mzml_out.to_widesage, params.fasta_file, params.sage_widewindow_parameters)

    dino_out = run_dinosaur(mzml_out.to_dinosaur)

    full_search_out = comet_out.join(sage_out).join(widesage_out)
    reconciled_out = reconcile_search_results(full_search_out)

    (annotated_out, etc) = annotate_psm_file(reconciled_out.join(dino_out))
    (concat_pin, moka_psms, moka_peps, moka_prots) = concatenate_and_mokapot(annotated_out.collect(), params.fasta_file)
    psm_report = generate_psm_output(concat_pin, moka_psms)

    if (params.isobaric_type == 0) {
      w_isobar_quant = psm_report 
    }
    else {
      isobar_data = extract_isobaric_values(mzml_out.to_isobar)
      w_isobar_quant = integrate_isobaric_data(psm_report, isobar_data.collect())
    }

    (peptide_agg, protein_agg) = peptide_protein_aggregation(w_isobar_quant, moka_peps, moka_prots)
}



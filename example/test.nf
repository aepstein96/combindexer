workflow {

    samples_ch = Channel
        .fromPath("/ru-auth/local/home/aepstein/cursor_projects/combindexer/example/barcodes/Samples.csv")
        .splitText()        .map { it.trim() }
        .filter { it }
        .filter { it != "Sample Name"}

    samples_ch_branched = samples_ch.branch {
        gex: it.startsWith("GEX_")
        targ: it.startsWith("Targeted_")
    }

    samples_ch_branched.gex.view()
}
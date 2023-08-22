library(bisplotti)
library(stringr)

# Creates plot showing reads that overlap both cpg and snp of interest
# fname  - absolute path to epiBED
# chr    - chromosome data is on
# cpg    - name of CpG of interest
# snp    - name of SNP of interest
# oname  - absolute path of output PDF
# region - region to subset plot to (default: NULL)
# width  - width of PDF (default: 24)
# height - height of PDF (default: 24)
make_plot <- function(
    fname,
    chr,
    cpg,
    snp,
    oname,
    region = NULL,
    width = 24,
    height = 24
) {
    # Build table with SNPs and CpGs
    frag <- readEpibed(fname, chr=chr, genome="hg38")
    mats <- tabulateEpibed(frag, region=region)
    mats <- mats$tab_cgvr
    
    # Order reads by methylation in CG of interest
    n_meths <- sum(rowSums(mats == "M", na.rm = TRUE)) + sum(rowSums(mats == "U", na.rm = TRUE))
    srtd <- mats[order(rowSums(mats == "M", na.rm = TRUE) / n_meths, decreasing = TRUE), ]
    srtd_m <- srtd[endsWith(srtd[,cpg], "M"),]
    srtd_m <- srtd_m[!startsWith(srtd_m[,cpg], "N"),]
    srtd_m <- srtd_m[rowSums(is.na(srtd_m)) != ncol(srtd_m),]
    srtd_u <- srtd[endsWith(srtd[,cpg], "U"),]
    srtd_u <- srtd_u[!startsWith(srtd_u[,cpg], "N"),]
    srtd_u <- srtd_u[rowSums(is.na(srtd_u)) != ncol(srtd_u),]
    srtd_comb <- rbind(srtd_u, srtd_m)
    srtd_comb <- srtd_comb[(!is.na(srtd_comb[,snp]) & !is.na(srtd_comb[,cpg])), ]

    # Create plot
    plt <- plotEpibed(srtd_comb, force=TRUE, n_loci=50, show_readnames = FALSE, show_positions = TRUE)
    ggsave(oname, plot=plt$epi, width=width, height=height)
}

dir.loc <- "/path/to/working/directory"

# Chr 15 imprinted CpG
f0 <- paste0(dir.loc, "/epibed/Asample.chr15_24954913_24955009.epibed.gz")
make_plot(
    f0, "chr15", "chr15:24954998-24954999", "chr15:24954924-24954924",
    paste0(dir.loc, "/chr15_asm.pdf"),
    region = "chr15:24954914-24955003", width=11, height=8
)

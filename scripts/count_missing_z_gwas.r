library(data.table)

gwas_dir = "/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/jingling_analysis/IBDverse_data/DeLange_sumstat/"
files = list.files(gwas_dir, pattern = ".*intersected.*\\.txt\\.gz$", full.names = T)

for(f in files){
    df = fread(f)
    df[, z_from_variance := BETA / sqrt(BETA_VAR)]
    df[, z_from_se := BETA / BETA_VAR]
    cat("\nFile:", basename(f), "\n")
    for(col in c("z_from_variance", "z_from_se")){
        x = df[[col]]
        cat(sprintf("  %s -> NA: %d, Inf: %d, -Inf: %d\n",
            col, sum(is.na(x)), sum(x == Inf, na.rm=T), sum(x == -Inf, na.rm=T)))
    }
}
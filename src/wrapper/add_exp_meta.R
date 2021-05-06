meta <- jsonlite::fromJSON('./input/meta.json', simplifyVector = TRUE)
experiment_id <- readLines('./output/experiment_id.txt', warn = FALSE)

n_filtered <- length(list.files('./output', 'pre-emptydrops'))

res <- data.frame(
  date = Sys.time(),
  experiment_id,
  experiment_name = meta$name,
  organism = meta$organism,
  url = sprintf('https://scp.biomage.net/experiments/%s', experiment_id),
  env = Sys.getenv('CLUSTER_ENV', 'staging'),
  n_filtered
)

write.table(res, "exp_meta.csv", sep = ",",
            col.names = !file.exists("exp_meta.csv"),
            row.names = FALSE,
            append = TRUE)

library(R.matlab)

convert_rds_to_mat <- function(directory) {
  # Obtenir tous les fichiers RDS dans les sous-dossiers
  files <- list.files(directory, pattern = "\\.rds$", recursive = TRUE, full.names = TRUE)
  
  if (length(files) == 0) {
    message("Aucun fichier RDS trouvé.")
    return()
  }
  
  for (file in files) {
    tryCatch({
      # Charger le fichier RDS
      data <- readRDS(file)
      
      # Construire le chemin du fichier de sortie .mat
      output_file <- sub("\\.rds$", ".mat", file)
      
      # Sauvegarder au format .mat
      writeMat(output_file, data = data)
      
      message(paste("Converti :", file, "->", output_file))
    }, error = function(e) {
      message(paste("Erreur lors de la conversion de", file, ":", e$message))
    })
  }
}

# Définir le répertoire de travail (le dossier courant)
directory <- "./output/sim_diff"  # Vous pouvez remplacer par un chemin spécifique
convert_rds_to_mat(directory)
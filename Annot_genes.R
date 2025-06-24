#!/bin/Rscript 

## Prend en entrée une liste de miARNs matures, un ou plusieurs fichiers de bases de données, et un répertoire de sortie. ## 
args = commandArgs(trailingOnly=TRUE)

# Vérification des arguments
if (length(args) < 6) {
  stop('Usage: Annot_genes.R <miRNA_list.txt> <databases> <output_directory> <field_mode> <nb_interact> <enrichment_databases>')
}

# Définition des fichiers et dossiers
miRNA_interest <- unlist(strsplit(args[1], ","))  # la liste des miRNA
print(miRNA_interest)
selected_databases <- unlist(strsplit(args[2], ","))  # Liste des bases sélectionnées
print(selected_databases)
output_dir <- args[3]  # Répertoire de sortie
field_mode <- args[4]  # INTER ou UNION
nb_interact <- as.numeric(args[5])  # Seuil de Count
selected_enrichment_dbs <- unlist(strsplit(args[6], ","))
print(selected_enrichment_dbs)
#print(paste("Chaine de ctr en entrée:",selected_databases))

# Création du répertoire de sortie si inexistant
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Lecture du fichier miRNA_list.txt
#if (!file.exists(miRNA_file)) {
#  stop(paste("Fichier miRNA introuvable :", miRNA_file))
#}



# Chargement des bibliothèques
suppressMessages(library(dplyr))
suppressMessages(library(readr))
suppressMessages(library(tidyr))
suppressMessages(library(biomaRt))
suppressMessages(library(clusterProfiler))
suppressMessages(library(ReactomePA))
suppressMessages(library(enrichR))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(DOSE))
suppressMessages(library(Cairo))
suppressMessages(library(purrr))




# Étape 1 : Lecture du fichier et filtrage des miRNA
#filtered_miRNA <- scan(miRNA_file, what = "", sep = "\n")

print(selected_databases)

# Étape 2 : Chargement des bases de données sélectionnées
#Et pré-traitement des bases de données

# Se connecter à la base de données Ensembl
#A voir dans la boucle
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

DBs_dataframe <- list()

for (db_path in selected_databases) {
  print(db_path)
  if (file.exists(db_path)) {
    db <- basename(db_path)
    print(paste("Lecture de la base de données :", db))
    
    #Traitement spécifique à chaque base de données
    if (grepl("ENCORI", db)){
      ENCORI <- read_delim(db_path, delim  =  "\t", comment = "#")
      #Renommer la colonne si nécessaire
      colnames(ENCORI) <- ifelse(colnames(ENCORI) %in% "miRNAname", "miRNA", colnames(ENCORI))
      # Filtre spécifique
      #ENCORI_filtered <- ENCORI %>%
        #filter(clipExpNum >= 2, pancancerNum >= 2, degraExpNum >=1) %>%
        #rowwise() %>%
        #filter(sum(c_across(c(RNA22, miRmap, microT, miRanda, PicTar, TargetScan, PITA)), na.rm = TRUE) >= 2) %>%
        #ungroup()
      
      DBs_dataframe$ENCORI <- ENCORI
      print(head(DBs_dataframe$ENCORI))

    } else if (db == "miRTarBase.txt") {   #Base de données dont l'URL n'est plus disponible
      miRTarBase <- read_delim(db_path,delim=",")
      #Filtre spécifique
      miRTarBase <- miRTarBase[-which(miRTarBase$'Support Type' != 'Functional MTI'),]
      names(miRTarBase)[4]<- 'geneName'
      DBs_dataframe$miRTarBase <- miRTarBase
      print(head(DBs_dataframe$miRTarBase))
      
      
      
    } else if (db == 'miRDB.txt') {
      miRDB <- read.table(db_path, sep =  "\t", header=F)
      colnames(miRDB) = c("miRNA","Gene","Prediction_score")
      miRDB <- miRDB[which(grepl("hsa-miR-",miRDB$miRNA)),]
      #Score
      score_path <- sub("\\.txt$", "_score.txt",db_path)
      ligne_score <- readLines(score_path)
      confidence_score <- as.numeric(sub(".*: ([0-9]+)%.*", "\\1", ligne_score))
      miRDB <- miRDB[miRDB$Prediction_score >= confidence_score,]
      # Convertir les noms des gènes
      converted_genes <- getBM(
        filters = "refseq_mrna",
        attributes = c("refseq_mrna", "hgnc_symbol"),
        values = miRDB$Gene,
        mart = ensembl ) 
      miRDB$Gene <- ifelse(miRDB$Gene %in% converted_genes[,1],converted_genes[,2],NA)
      names(miRDB)[2] <-'geneName'
      DBs_dataframe$miRDB <- miRDB
      print(head(DBs_dataframe$miRDB))
      
      
      
      
    } else if (db == 'TargetScan.txt'){
      TargetScan <- read.table(db_path, sep =  "\t", comment = "#",header = T)
      TargetScan <- TargetScan[which(TargetScan$'Species.ID'== '9606'),]
      #Score
      score_path <- sub("\\.txt$", "_score.txt",db_path)
      ligne_score <- readLines(score_path)
      confidence_score <- (as.numeric(sub(".*: ([0-9]+)%.*", "\\1", ligne_score)))/100
      TargetScan <- TargetScan[as.numeric(TargetScan$PCT) >=confidence_score,]
      # Séparation
      TargetScan <- TargetScan %>%
        separate_rows(miR.Family, sep = "/") %>%  # Séparer les valeurs multiples en lignes distinctes
        mutate(miR.Family = ifelse(grepl("^miR-", miR.Family), miR.Family, paste0("hsa-miR-", miR.Family)))  # Ajouter "miR-" si absent
      names(TargetScan)[1] <-  'miRNA'
      names(TargetScan)[3] <-  'geneName'
      TargetScan <- na.omit(TargetScan)
      DBs_dataframe$TargetScan <- TargetScan
      print(head(DBs_dataframe$TargetScan))
      
      
      
    } else if (db == 'DianaTarbase.csv'){
      DianaTools <- read.table(db_path, sep =  "\t", comment = "#", header=T)
      DianaTools <- DianaTools[which(DianaTools$species == "Homo sapiens"),]
      colnames(DianaTools) <- ifelse(colnames(DianaTools) %in% "mirna", "miRNA", colnames(DianaTools))
      DBs_dataframe$DianaTools <- DianaTools
      print(head(DBs_dataframe$DianaTools))
      
      
      
    } else if (db == 'TarBase.txt') {
      Tarbase <- read.table(db_path, sep =  "\t", header=T)
      #Score
      score_path <- sub("\\.txt$", "_score.txt",db_path)
      ligne_score <- readLines(score_path)
      confidence_score <- (as.numeric(sub(".*: ([0-9]+)%.*", "\\1", ligne_score)))/100
      Tarbase <- Tarbase[which(Tarbase$microt_score >= confidence_score),]
      names(Tarbase)[2] <- 'miRNA'
      #converted_genes <- getBM(
       # filters = "ensembl_gene_id",
        #attributes = c("ensembl_gene_id", "hgnc_symbol"),
        #values = Tarbase$ensembl_gene_id,
        #mart = ensembl ) 
      #Tarbase$ensembl_gene_id <- ifelse(Tarbase$ensembl_gene_id %in% converted_genes[,1],converted_genes[,2],NA)
      names(Tarbase)[4] <- "geneName"
      DBs_dataframe$Tarbase <- Tarbase
      print(head(DBs_dataframe$Tarbase))
      
    }
    
  }else{
    print(paste("Attention : Fichier introuvable pour la base de données", db))
  }
}
print("Toute les bases de données ont été téléchargées ;)")



# Associer les miRNA avec leurs cibles
Interactions <- bind_rows(lapply(DBs_dataframe, function(df) unique(df[, c("miRNA", "geneName")])))

# Appliquer INTER ou UNION
if (field_mode == "INTER") {
  common_mirnas <- Reduce(intersect, lapply(DBs_dataframe, function(df) unique(df$miRNA)))
  Interactions <- Interactions %>% filter(miRNA %in% common_mirnas)
}

# Sélectionner les gènes régulés par plusieurs miRNA uniques
#selected_genes <- Interactions %>%
 # group_by(geneName) %>%
  #summarise(n_unique_miRNA = n_distinct(miRNA)) %>%  # Compter les miRNA uniques
  #filter(n_unique_miRNA >= nb_interact) %>%  # Filtrer selon le seuil
  #pull(geneName)

# Étape 2 : Compter le nombre de miRNAs distincts par gène
gene_counts <- table(Interactions$geneName)

# Étape 3 : Garder uniquement les gènes ciblés par au moins 2 miRNAs
selected_genes <- names(gene_counts[gene_counts >= nb_interact])

# Affiner les interactions en fonction des gènes sélectionnés
#Interactions_genes_filter <- Interactions%>%
#  filter(geneName %in% selected_genes)

# Étape 4 : Filtrer la base de données
Interactions_genes_filter <- Interactions[Interactions$geneName %in% selected_genes, ]

# Filtrer les interactions pour les miRNA d'intérêt
#Interactions_miRNA_filtered <- Interactions_genes_filter %>%
#filter(miRNA %in% unique(miRNA_interest))
Interactions_miRNA_filtered <- Interactions_genes_filter[Interactions_genes_filter$miRNA %in% unique(miRNA_interest), ]

selected_genes_miRNA <- Interactions_miRNA_filtered$geneName

# Générer un tableau miRNA -> liste de gènes
miRNA_gene_table <- Interactions_miRNA_filtered %>%
  group_by(miRNA) %>%
  summarise(genes = paste(unique(geneName), collapse = ", "))


# Conversion en identifiants Entrez pour enrichPathway
selected_genes_entrezIDs <- mapIds(org.Hs.eg.db,
                                   keys = selected_genes,
                                   column = "ENTREZID",
                                   keytype = "SYMBOL",
                                   multiVals = "first")

#print("Gènes cibles régulés par plusieurs miRNA:")
#print(selected_genes)



# Génération des fichiers dans le bon répertoire
output_selected_genes <- file.path(output_dir, "Selected_Genes.txt")
write.table(miRNA_gene_table, file.path(output_dir, "Interactions_miRNA_filtered.txt"), row.names = FALSE, sep = "\t")
write.table(selected_genes_miRNA, output_selected_genes, row.names = FALSE, col.names = FALSE, quote = FALSE)


# Étape 4 : Annotation fonctionnelle avec EnrichR, ReactomePA, ClusterProfiler, Gene Ontology, DAVID, Panther
print("Début des analyses d'enrichissement géniques.\n Cela peut prendre un certain temps...")

# EnrichR
if ("EnrichR" %in% selected_enrichment_dbs) {
  # EnrichR
  print("Lancement de l'enrichissement avec EnrichR...")
  dbs <- c("GO_Biological_Process_2021", "KEGG_2021_Human", "Reactome_2022")
  enrich_results <- enrichr(selected_genes, dbs)
  # Redirection des fichiers de sortie
  output_go <- file.path(output_dir, "GO_Results.txt")
  output_kegg <- file.path(output_dir, "KEGG_Results.txt")
  output_reactome <- file.path(output_dir, "Reactome_Results.txt")
  write.table(enrich_results$GO_Biological_Process_2021, output_go, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = '\t')
  write.table(enrich_results$KEGG_2021_Human, output_kegg, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = '\t')
  write.table(enrich_results$Reactome_2022, output_reactome, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = '\t')
}

#print("Résultats EnrichR:")
#print(enrich_results)

# ClusterProfiler (Gene Ontology)
if ("ClusterProfiler" %in% selected_enrichment_dbs) {
  # ClusterProfiler
  print("Lancement de l'enrichissement avec ClusterProfiler...")
  ego <- enrichGO(gene = selected_genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05)
  print("Enrichissement Gene Ontology...")
  write.csv(ego, file.path(output_dir, "ClusterProfiler_results.csv"))
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    pdf(file.path(output_dir, "ClusterProfiler.pdf"))
    barplot(ego, showCategory=10, title="Gene Ontology (Biological Process)")
    dev.off()
    print("✅ ClusterProfiler PDF généré avec succès")
  } else {
    print("⚠️ Aucune donnée disponible pour ClusterProfiler, pas de PDF généré")
  }
  Cairo::CairoPDF(file.path(output_dir, "ClusterProfiler.pdf"))
  barplot(ego, showCategory=10, title="Gene Ontology (Biological Process)")
  dev.off()
  print("PDF enregistré avec succès")
  
  }#print(ego)

# ReactomePA
if ("ReactomePA" %in% selected_enrichment_dbs) {
  # ReactomePA
  print("Lancement de l'enrichissement avec ReactomePA...")
  reactome_results <- enrichPathway(gene = selected_genes_entrezIDs, organism = "human", pvalueCutoff = 0.05)
  print("Enrichissement Reactome...")
  write.csv(reactome_results, file.path(output_dir, "ReactomePA_results.csv"))
  if (!is.null(reactome_results) && nrow(as.data.frame(reactome_results)) > 0) {
    pdf(file.path(output_dir, "ReactomePA.pdf"))
    barplot(reactome_results, showCategory=10, title="Reactome Pathways")
    dev.off()
    print("✅ ReactomePA PDF généré avec succès")
  } else {
    print("⚠️ Aucune donnée disponible pour ReactomePA, pas de PDF généré")
  }
  Cairo::CairoPDF(file.path(output_dir, "ReactomePA.pdf"))
  barplot(reactome_results, showCategory=10, title="Reactome Pathways")
  dev.off()
  print("PDF enregistré avec succès")
  }#print(reactome_results)








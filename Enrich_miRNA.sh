#!/bin/bash 

## Prend en entrée une liste de miRNAs matures ## 
## Dans R: #
## Télécharge une ou plusieurs DB ou permet à l'utilisateur d'en fournir manuellement ##
## Renvoie une liste de gènes régulés et réalise une analyse d'enrichissement ##

# Vérifier si un fichier contenant les miARN est fourni en argument
if [ -f "$1" ]; then
    miARN_list=($(cat "$1"))  # Lire le fichier et stocker les miARN dans un tableau
else
    echo "Entrez les miARN séparés par des espaces :"
    read -r -a miARN_list
fi

# Demander à l'utilisateur s'il souhaite entrer les bases de données manuellement
echo "Voulez-vous entrer manuellement les chemins vers vos bases de données (oui/non) ?"
read -r manual_input

if [ "$manual_input" == "oui" ]; then
    echo "Entrez les chemins complets vers vos bases de données séparés par des espaces :"
    read -r -a manual_db_paths
    db_files=(${manual_db_paths[@]})
else
    # Exécution normale avec téléchargement des bases de données

    # Bases de données valides
    declare -A database_urls=(
        ["ENCORI"]="https://rnasysu.com/encori/api/miRNATarget/?assembly=hg38&geneType=mRNA&miRNA=all&clipExpNum=2&degraExpNum=1&pancancerNum=2&programNum=2&program=None&target=all&cellType=all"
        ["miRDB"]="https://mirdb.org/download/miRDB_v6.0_prediction_result.txt.gz"
        ["TargetScan"]="https://www.targetscan.org/vert_80/vert_80_data_download/Predicted_Targets_Info.default_predictions.txt.zip"
        ["TarBase"]="https://dianalab.e-ce.uth.gr/tarbasev9/data/Homo_sapiens_TarBase-v9.tsv.gz"
        ["DianaTarbase"]="http://diana.imis.athena-innovation.gr/DianaTools/data/TarBase7data.tar.gz"
    )

    prediction_dbs=("miRDB" "TarBase" "TargetScan")
    declare -A confidence_scores

    while true; do
        echo "Entrez au moins une base de données parmi : ${!database_urls[*]} ou 'All' pour tout télécharger"
        read -r -a user_dbs

        valid_input=true
        for db in "${user_dbs[@]}"; do
            if [[ "$db" != "All" && -z "${database_urls[$db]}" ]]; then
                valid_input=false
                echo "Base de données invalide : $db"
            fi
        done

        if $valid_input; then
            break
        else
            echo "Veuillez entrer des bases valides."
        fi
    done



    if [[ " ${user_dbs[*]} " =~ " All " ]]; then
        user_dbs=("ENCORI" "miRDB" "TargetScan" "TarBase" "DianaTarbase")
    fi
    
    # Si ENCORI est sélectionné, demander les paramètres spécifiques
    if [[ " ${user_dbs[*]} " =~ "ENCORI" ]]; then
    	echo "Vous avez sélectionné ENCORI. Veuillez entrer les paramètres suivants :"

    	while true; do
        	read -p "CLIP-seq_data (2 recommended) : " clip_seq
        	if [[ "$clip_seq" =~ ^[0-9]+$ && "$clip_seq" -ge 0 ]]; then break; fi
        		echo "Valeur invalide. Entrez un nombre >= 0."
    	done

   	while true; do
        	read -p "Degradome_data (>= 1 recommended) : " degradome
        	if [[ "$degradome" =~ ^[0-9]+$ && "$degradome" -ge 0 ]]; then break; fi
        		echo "Valeur invalide. Entrez un nombre >= 0."
    	done

    	while true; do
        	read -p "Pan-Cancer (>= 2 recommended) : " pancancer
        	if [[ "$pancancer" =~ ^[0-9]+$ && "$pancancer" -ge 0 ]]; then break; fi
        		echo "Valeur invalide. Entrez un nombre >= 0."
    	done

    	while true; do
        	read -p "Predicted_targets_programs (>= 2 recommended) : " predicted_targets
        	if [[ "$predicted_targets" =~ ^[0-9]+$ && "$predicted_targets" -ge 0 ]]; then break; fi
        		echo "Valeur invalide. Entrez un nombre >= 0."
    	done

    	# Construire l'URL ENCORI avec les paramètres utilisateur
    	encori_url="https://rnasysu.com/encori/api/miRNATarget/?assembly=hg38&geneType=mRNA&miRNA=all&clipExpNum=${clip_seq}&degraExpNum=${degradome}&pancancerNum=${pancancer}&programNum=${predicted_targets}&program=None&target=all&cellType=all"
    	database_urls["ENCORI"]="$encori_url"
    fi

	

    for db in "${prediction_dbs[@]}"; do
        if [[ " ${user_dbs[*]} " =~ " $db " ]]; then
            while true; do
                read -p "Entrez le score de confiance pour $db (0-100, recommandé : 80) : " confidence
                if [[ "$confidence" =~ ^[0-9]+$ && "$confidence" -ge 0 && "$confidence" -le 100 ]]; then
                    confidence_scores["$db"]="$confidence"
                    break
                fi
                echo "Valeur invalide. Entrez un nombre entre 0 et 100."
            done
        fi
    done

    dest_dir="$(pwd)/DB"
    mkdir -p "$dest_dir"
    chmod u+w "$dest_dir"

    echo "Le téléchargement peut nécessiter plusieurs Mo. Voulez-vous continuer ? (oui/non)"
    read -r confirmation
    if [[ "$confirmation" != "oui" ]]; then
        echo "Téléchargement annulé."
        exit 1
    fi






    for db in "${user_dbs[@]}"; do
        echo "Téléchargement de $db..."
        wget -O "${dest_dir}/${db}" "${database_urls[$db]}"
    done

    for db in "${!confidence_scores[@]}"; do
        echo "Score de confiance : ${confidence_scores[$db]}%" > "${dest_dir}/${db}_score.txt"
    done

    echo "Décompression des fichiers..."
    for compressed_file in "$dest_dir"/*.{gz,zip,tar.gz}; do
        [ -e "$compressed_file" ] || continue
        case "$compressed_file" in
            *.tar.gz) tar -xzf "$compressed_file" -C "$dest_dir" ;;
            *.gz) gunzip "$compressed_file" ;;
            *.zip) unzip -j "$compressed_file" -d "$dest_dir" ;;
        esac
        rm -f "$compressed_file"
    done

    if [ -f "$dest_dir/Predicted_Targets_Info.default_predictions.txt" ]; then
        mv "$dest_dir/Predicted_Targets_Info.default_predictions.txt" "$dest_dir/TargetScan.txt"
    fi

    db_files=()
    for db in "${user_dbs[@]}"; do
        case "$db" in
            "ENCORI") db_file="$dest_dir/ENCORI.txt" ;;
            "miRDB") db_file="$dest_dir/miRDB.txt" ;;
            "TargetScan") db_file="$dest_dir/TargetScan.txt" ;;
            "TarBase") db_file="$dest_dir/TarBase.txt" ;;
            "DianaTarbase") db_file="$dest_dir/DianaTarbase.csv" ;;
        esac
        if [ -f "$db_file" ]; then
            db_files+=("$db_file")
        fi
    done
fi

# Étapes communes dans les deux cas
echo "Choisissez le mode Field (INTER ou UNION) :"
read -r field_mode

echo "Entrez le seuil minimal de nombre d'interactions (nb_interact) :"
read -r nb_interact

echo "Entrez le nom du répertoire de résultats (par défaut Results_enrichment) :"
read -r custom_results_dir
if [ -z "$custom_results_dir" ]; then
    custom_results_dir="Results_enrichment"
fi
results_dir="$(pwd)/$custom_results_dir"
mkdir -p "$results_dir"

echo "Choisissez une ou plusieurs bases d'enrichissement parmi : EnrichR, ClusterProfiler, ReactomePA ou 'All' pour tout inclure."
read -r -a user_enrichment_dbs
if [[ " ${user_enrichment_dbs[*]} " =~ " All " ]]; then
    user_enrichment_dbs=("EnrichR" "ClusterProfiler" "ReactomePA")
fi


##Préparation des arguments finaux pour le script R ##
# Liste de miARNs
miRNA_list_str=$(IFS=,; echo "${miARN_list[*]}")
# Convertir la liste des fichiers en une chaîne séparée par des virgules
db_files_str=$(IFS=,; echo "${db_files[*]}")
# Créer une liste de bases de données d'enrichissement pour le script R
enrichment_dbs_str=$(IFS=,; echo "${user_enrichment_dbs[*]}")

Rscript Annot_genes.R "$miRNA_list_str" "$db_files_str" "$results_dir" "$field_mode" "$nb_interact" "$enrichment_dbs_str"

echo "✅ Annotation terminée. Résultats disponibles dans $results_dir"


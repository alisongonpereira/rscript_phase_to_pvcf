# Load necessary libraries
library(vcfR)
library(dplyr)

# Validation function
validate_inputs <- function(vcf_file, phase_file) {
  
  # Validate VCF file
  cat("Validating VCF file...\n")
  vcf_data <- try(read.vcfR(vcf_file), silent = TRUE)
  if(class(vcf_data) == "try-error" || !("vcfR" %in% class(vcf_data))) {
    cat("Error reading VCF file. Please check the file format.\n")
    return(FALSE)
  }
  
  # Validate PHASE output file
  cat("Validating PHASE output file...\n")
  phase_data <- try(readLines(phase_file), silent = TRUE)
  if(class(phase_data) == "try-error") {
    cat("Error reading PHASE output file. Please check the file format.\n")
    return(FALSE)
  }
  
  # Adjusted validation logic for PHASE output format
  valid_lines <- grepl("^IND: |^[0-9]", phase_data)
  if(!all(valid_lines)) {
    cat("PHASE output file format is incorrect. Please check the formatting.\n")
    return(FALSE)
  }
  
  # Extract sample names and compare between VCF and PHASE output
  vcf_samples <- colnames(vcf_data@gt)[-1]
  cat(length(vcf_samples))
  # Extracting sample names from PHASE output, removing 'IND: ' prefix
  phase_samples <- unique(sapply(strsplit(phase_data[startsWith(phase_data, "IND:")], " "), `[`, 2))
  phase_samples <- gsub("IND: ", "", phase_samples)  # Remove the 'IND: ' prefix
  cat(length(phase_samples))
  # if(length(vcf_samples) != length(phase_samples) || !all(vcf_samples %in% phase_samples)) {
  if(length(vcf_samples) != length(phase_samples)) {
      
    cat("Mismatch in the number of samples or sample names between VCF and PHASE output.\n")
    return(FALSE)
  }
  
  
  # Check for matching variant count
  phase_variant_count <- count_variants_phase(phase_data)
  vcf_variant_count <- nrow(vcf_data@gt)
  if(phase_variant_count != vcf_variant_count) {
    cat("Mismatch in the number of variants between VCF and PHASE output.\n")
    return(FALSE)
  }
  
  cat("Validation completed successfully.\n")
  return(TRUE)
}

# Função modificada para contar variantes no arquivo PHASE
count_variants_phase <- function(phase_data) {
  # Encontrar a primeira linha de dados de haplótipo
  haplotype_line <- phase_data[which(!startsWith(phase_data, "IND:"))[1]]
  
  # Dividir pela vírgula e pegar o primeiro elemento
  haplotype_first_part <- strsplit(haplotype_line, ",")[[1]][1]
  
  # Converter para character e remover espaços
  haplotype_chars <- gsub(" ", "", as.character(haplotype_first_part))
  
  # Contar o número de caracteres
  return(nchar(haplotype_chars))
}


# # Function to parse PHASE output
parse_phase_output <- function(phase_output_file, threshold) {
  cat("Parsing PHASE output file...\n")
  phase_data <- readLines(phase_output_file)
  parsed_data <- data.frame(ID = character(), Haplotype1 = character(), Haplotype2 = character(), Probability = numeric(), stringsAsFactors = FALSE)
  phase_data <- gsub(pattern = " ", replacement = "", x = phase_data)
  current_id <- ""
  for (line in phase_data) {
    if (startsWith(line, "IND:")) {
      current_id <- strsplit(line, ":")[[1]][2]
    } else if (current_id != "" && line != "") {
      # Dividir a linha por vírgula e processar cada haplótipo
      haplotypes_info <- strsplit(line, ",")[[1]]
      for (haplotype_info in haplotypes_info) {
        haplotype1 <- haplotypes_info[1]
        haplotype2 <- haplotypes_info[2]
        probability <- as.numeric(haplotypes_info[3])
        parsed_data <- rbind(parsed_data, data.frame(ID = current_id, Haplotype1 = haplotype1, Haplotype2 = haplotype2, Probability = probability))
      }
    }
  }
  
  
  # Selecionar o haplótipo com maior probabilidade para cada amostra
  parsed_data <- parsed_data %>%
    group_by(ID) %>%
    arrange(desc(Probability)) %>%
    slice(1) %>%
    ungroup()
  
   # print(head(parsed_data))
  # Marcar amostras com probabilidade máxima abaixo do limiar
  low_confidence_samples <- parsed_data %>%
    filter(Probability < threshold) %>%
    select(ID, Probability)
  # print(low_confidence_samples)
  cat("PHASE output file parsed successfully.\n")
  return(list(parsed_data = parsed_data, low_confidence_samples = low_confidence_samples))
}

substitute_genotypes <- function(vcf_data, phase_data) {
  cat("Substituting genotypes in the VCF file...\n")
  
  for (i in 1:nrow(phase_data)) {
    individual <- phase_data$ID[i]
    haplotype1 <- phase_data$Haplotype1[i]
    haplotype2 <- phase_data$Haplotype2[i]
    
    # Encontrar a coluna correspondente ao indivíduo no VCF
    individual_col_index <- which(colnames(vcf_data@gt) == individual)
    
    for (pos in 1:nchar(haplotype1)) {
      # Concatenar os caracteres correspondentes dos haplótipos
      genotype <- paste(substr(haplotype1, pos, pos), substr(haplotype2, pos, pos), sep = "|")
      
      # Substituir o genótipo na posição e coluna apropriadas
      vcf_data@gt[pos, individual_col_index] <- genotype
    }
  }
  
  cat("Genotypes substitution completed.\n")
  return(vcf_data)
}

# Main function to convert PHASE output to VCF format
convert_phase_to_vcf <- function(original_vcf_file, phase_output_file, threshold = 0.8, out_file_name = "modified_vcf.vcf.gz") {
  cat("Validating inputs...\n")
  # Validate inputs
  if(!validate_inputs(original_vcf_file, phase_output_file)) {
    cat("Input validation failed. Terminating process.\n")
    return()
  }
  cat("Inputs validated!\n")
  cat("Starting the conversion process...\n")
  
  # Read original VCF file
  cat("Reading original VCF file...\n")
  vcf_data <- read.vcfR(original_vcf_file)
  cat("Original VCF file read successfully.\n")
  
  # Parse PHASE output
  cat("Reading PHASE file...\n")
  # phase_data <- parse_phase_output(phase_output_file)
  
  
  # Chamada da função parse_phase_output
  phase_data_list <- parse_phase_output(phase_output_file, threshold = threshold)
  
  # Extrair o dataframe 'parsed_data' da lista
  parsed_phase_data <- as.data.frame(phase_data_list$parsed_data)
  # print(head(parsed_phase_data))
  cat("PHASE file read successfully.\n")
  
  
  # Chamar a função substitute_genotypes com o dataframe
  modified_vcf <- substitute_genotypes(vcf_data, parsed_phase_data)

  
  # Generate report for low-confidence haplotypes
  report <- as.data.frame(phase_data_list$low_confidence_samples)
  #print(report)
  
  # Write modified VCF to a file
  cat(paste0("Writing modified VCF to ",out_file_name," file...\n"))
  write.vcf(modified_vcf, out_file_name)
  cat("Modified VCF file written.\n")
  
  # Write the report to a file
  cat("Writing haplotype confidence report to haplotype_low_confidence_report.csv...\n")
  write.csv(report, "haplotype_low_confidence_report.csv")
  cat("Report written. Process completed.\n")
}

# Example usage
# convert_phase_to_vcf("path/to/original.vcf", "path/to/phase_output.txt")

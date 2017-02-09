library(rcellminerPubchem)

pubchem <- read.table("import/tempdir/PubChem.CID_dictionary.csv", sep=",", header=TRUE, row.names=1, stringsAsFactors = FALSE)

smiles <- NULL
inchikeys <- NULL 

for(id in pubchem[,"X"]) {
  smiles <- c(smiles, getIsomericSmilesFromCid(id))
  inchikeys <- c(inchikeys, getInChiKeyFromCid(id))
}

df <- cbind(pubchem, smiles, inchikeys)
write.table(df, "import/tempdir/pubchemStructures.tsv", quote=FALSE, sep="\t", row.names=FALSE)

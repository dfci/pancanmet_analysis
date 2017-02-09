library(treemap)

# Parameters 
hmdbClassificationFilename <- file.path("data", "hmdbClassification.rds")

pdfFilename <- file.path("manuscript", "figures", "Figure7", "metab_treemap.pdf")

# Get saved results or download from HMDB
if(file.exists(hmdbClassificationFilename)) {
  results <- readRDS(hmdbClassificationFilename)
} else {
  # Parse HMDB pages
  library(rvest)
  library(httr)
  library(simpleRCache)
  
  setCacheRootPath()
  GETCached <- addMemoization(GET)
  
  hmdb <- read.table(file.path("import", "tempdir", "Human.Metabolome.Database_dictionary.csv"), sep=",", header=TRUE, stringsAsFactors = FALSE)
  
  results <- NULL
  
  excludeIds <- c("HMDB13123")
  hmdbIds <- setdiff(hmdb$X, excludeIds)
  
  for(hmdbId in hmdbIds) {
    #hmdbId <- "HMDB00122"
    url <- paste0("http://www.hmdb.ca/metabolites/", hmdbId)
    htmlContent <- url %>% GETCached() %>% content("text") %>% read_html()
    
    # Kingdom
    t1 <- html_node(htmlContent, xpath="/html/body/main/table/tbody[1]/tr[23]/td/a")
    kingdom <- html_text(t1, trim=TRUE)
    
    # Super Class
    t1 <- html_node(htmlContent, xpath="/html/body/main/table/tbody[1]/tr[24]/td/a")
    superClass <- html_text(t1, trim=TRUE)
    
    # Class
    class <- NA
  
    tryCatch({
      t1 <- html_node(htmlContent, xpath="/html/body/main/table/tbody[1]/tr[25]/td/a")
      class <- html_text(t1, trim=TRUE)
    }, error = function(e) {
      cat("ERROR: ID: ", hmdbId, "\n")
    })
    
    # Sub Class 
    # NOTE: May be missing
    subClass <- NA 
    
    tryCatch({
      t1 <- html_node(htmlContent, xpath="/html/body/main/table/tbody[1]/tr[26]/td/a")
      subClass <- html_text(t1, trim=TRUE)
    }, error = function(e) {
      cat("ERROR: ID: ", hmdbId, "\n")
    })
    
    tmpResult <- data.frame(hmdbId=hmdbId, kingdom=kingdom, superClass=superClass, class=class, subClass=subClass, stringsAsFactors=FALSE)
    
    results <- rbind(results, tmpResult)  
  }
  
  saveRDS(results, hmdbClassificationFilename)
}

freqCol <- "subClass"

t1 <- as.data.frame(table(results[, freqCol]))
t2 <- data.frame(subClass=t1$Var1, freq=t1$Freq, stringsAsFactors=FALSE)

t3 <- merge(results, t2, by.x=freqCol, by.y=freqCol)
head(t3)

pdfWidth <- 3.4
pdfHeight <- 3.4

pdf(pdfFilename, width=pdfWidth, height=pdfHeight)

treemap(t3,
        index=c("class"),
        vSize = "freq",
        type="index",
        title="",
        align.labels = list(c("center", "center"))
)

dev.off()

treemap(t3,
        index=c("class", "subClass"),
        vSize = "freq",
        type="index",
        title="HMDB Classes and Sub-Classes of Pan-Cancer Metabolites",
        fontsize.title = 10,
        fontsize.labels = c(8, 8),
        overlap.labels = 1,
        align.labels = list(c("left", "top"), c("center", "center")),
        title.legend = "Legend",
        position.legend = "right",
        fontsize.legend = 6,
        force.print.labels = FALSE
)

dev.off()

# treemap(t3, #Your data frame object
#         index=c("subClass"),  #A list of your categorical variables
#         vSize = "freq",  #This is your quantitative variable
#         type="index", #Type sets the organization and color scheme of your treemap
#         title="HMDB Classes and Sub-Classes of Pan-Cancer Metabolites",
#         fontsize.title = 10,
#         fontsize.labels = c(8, 8),
#         overlap.labels = 1,
#         align.labels = list(c("left", "top"), c("center", "center")),
#         title.legend = "Legend",
#         position.legend = "right",
#         fontsize.legend = 6
# )
# 
# unique(results$class)
# sort(table(results$class))






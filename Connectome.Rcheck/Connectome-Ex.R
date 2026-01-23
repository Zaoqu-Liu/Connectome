pkgname <- "Connectome"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('Connectome')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("CreateConnectome")
### * CreateConnectome

flush(stderr()); flush(stdout())

### Name: CreateConnectome
### Title: CreateConnectome
### Aliases: CreateConnectome

### ** Examples

## Not run: 
##D # Basic usage
##D connectome <- CreateConnectome(seurat_obj, species = "human")
##D 
##D # Without p-value calculation (faster)
##D connectome <- CreateConnectome(seurat_obj, species = "human", p.values = FALSE)
##D 
##D # With custom ligand-receptor list
##D my_lr <- data.frame(ligand = c("IL6", "TNF"), receptor = c("IL6R", "TNFRSF1A"),
##D                     mode = c("cytokine", "cytokine"))
##D connectome <- CreateConnectome(seurat_obj, LR.database = "custom", custom.list = my_lr)
## End(Not run)




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')

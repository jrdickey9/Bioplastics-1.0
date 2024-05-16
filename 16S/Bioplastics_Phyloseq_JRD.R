#Bioplastics 2 *preliminary* R script -- preliminary because we are waiting on 4 samples with low read depth to be resequenced. 

#This R script will import the exports from my QIIME2 pipeline script so that they can be combined into a single Phyloseq object for further processing, object extraction, and statistical analysis. 

#installing packages
#install.packages("microbiome")
#BiocManager::install("microbiome") #works better for me

#check R version
R.version.string
#install packages
install.packages("BiocManager")
BiocManager::install("microbiome")
install.packages("C:/Users/smort/Downloads/vegetarian_1.2.tar.gz", repos = NULL, type = "source")
BiocManager::install("phyloseq", force = TRUE)
install.packages("gdata")
install.packages("ecodist")
install.packages("ape")
install.packages("phytools")
install.packages("castor")
install.packages("doParallel")

#load packages; most of these are useful in general
library(microbiome)
library(vegetarian) #might be a tricky install. 
library(phyloseq); packageVersion("phyloseq") #1.26.1; not sure if these are the most updated in comparison to a fresh install, but they will be useful in noting ahead of time for manuscripts. 
library(ggplot2); packageVersion("ggplot2") #3.4.1
library(gdata)
library(ecodist)
library(vegan) #we love her! 
library("car") #very likely will be a tricky install if you haven't experienced this already.  
library(dplyr)
library(biomformat)
library("ape") #phylogenetic tools packages
library(phytools) #phylogenetic tools packages
library(castor)
library(doParallel) #used for UniFrac but meh not super necessary for small data sets.

set.seed(32) #permutational analyses need to be easily reproducable (this includes rarefaction with vegan)

#16S data import
#set your working directory and replace the path
setwd("C:/Users/smort/Documents/R/projects/MP_16S")
#Jonathan code - setwd("~/Desktop/Jonathan_Dickey_LabMac/Bioplastics_Shurin/Pyro2_Prelim/16S/BioP_paired_dadatableR5")

#read in biom table
ASV_reads <- read_biom("feature-table.biom") #uses biomformat package
ASV_table <- as.data.frame(as.matrix(biom_data(ASV_reads))) #shape shifting
otu_tab<-t(ASV_table) #transpose cols and rows
dim(otu_tab) #are the dimensions correct? #150 samps x 6987 ASVS 
colnames(otu_tab) #what about col names? #feature IDs
rownames(otu_tab) #what about row names? #Sample names
df.ASV<-as.data.frame(otu_tab)
df.ASV[1:10,1:10]
rownames(df.ASV) #nothing has changed, gucci
colnames(df.ASV) #nothing has changed, gucci

#reading in meta-data file
setwd("C:/Users/smort/Documents/R/projects/MP_16S")
metadata<-read.csv("16S_BioP_metadata.csv",header=TRUE)
#READ: need to include your own meta data file with treatment groups. Highly suggest you work off this one though as to not create downstream effects when importing into a phyloseq object. In other words, add treatment columns and other necessary columns to this file. Do not change the number of rows as they correspond to the same number of rows in the ASV table, and the same number of rows in the metadata file, and the same number of tips in the phylogenetic tree. Much more complex to independently change those (particularly if you aren't familiar with Newick trees). One can work in an S4 based R object (created with phyloseq) and prune from all four components (ASV table, metadata file, taxonomy table, and tree) all at once. 
dim(metadata) #150 samps x 4 columns 

#fucking around with taxonomy vector to turn into table
#write.csv(taxonomy,"taxa.csv", header=FALSE)
setwd("C:/Users/smort/Documents/R/projects/MP_16S")
taxa2<-read.csv("taxonomy3.csv")
ph.headers<-c("ID", "Kingdom","Phylum","Class","Order","Family","Genus","species")
colnames(taxa2)<-ph.headers
rownames(taxa2) #just numbers, boo & hiss
rownames(taxa2)<-taxa2[,1]
rownames(taxa2) #much better
dim(taxa2) #6,987 x 8
taxa3<-taxa2[,2:8]
colnames(taxa3) #Gucci

#reading in phylogenetic tree
setwd("C:/Users/smort/Documents/R/projects/MP_16S")
dated.16Stree<-read.tree(file="tree.nwk") #reading in tree, uses ape package
is.rooted(dated.16Stree) #TRUE
sample_names(dated.16Stree) #NULL
dated.16Stree$tip.label #for example "bf08d62b32cced86e829cba893bdf318" 

#creating phyloseq object
str(taxa3) #data.frame 291 x 7
taxa3.m<-as.matrix(taxa3) #VERY NECESSARY TO DO, DON'T SKIP. 
str(taxa3.m)
colnames(taxa3.m)
str(df.ASV) #data.frame 150 obs. of 6987 vars
str(metadata) #data.frame 150 obs. of 4 vars

rownames(df.ASV)<-as.character(rownames(df.ASV))
colnames(df.ASV) #accession numbers
rownames(metadata)<-as.character(metadata[,1])
rownames(taxa3.m)<-as.character(rownames(taxa3.m)) #these the accession numbers
samp.names<-as.character(metadata[,1])

#to set up sample names to match (originally marked NULL)
sample_names(df.ASV)<-samp.names
sample_names(metadata)<-samp.names
sample_names(taxa3.m)<-samp.names
sample_names(dated.16Stree)<-samp.names

#Here is the actual phyloseq object
BioP.phylo<-phyloseq(otu_table(df.ASV, taxa_are_rows=FALSE), sample_data(metadata), tax_table(taxa3.m), phy_tree(dated.16Stree))

#replace "taxa names" to short hand so that they are easier to view in R. No more accession numbers which are not helpful here anyhow. 
taxa_names(BioP.phylo) <- paste0("ASV", seq(ntaxa(BioP.phylo)))
BioP.phylo@otu_table[1:10,1:10]
dim(BioP.phylo@otu_table) #150 x 6987
rowSums(BioP.phylo@otu_table) #eye opening
colSums(BioP.phylo@otu_table) #less eye opening

#Examining taxonomic ranks to examine where chloroplasts and mitochondria are nested within
#Need to examine Unassigned / Unknown taxa and do blasts to see if they TRULY are unknown or are chloroplasts/host DNA in disguise. The taxonomic resolution you wish to do this at is dependent on maybe an in person conversation. One needs the .qzv sequences file from that ASV table generation in QIIME2 to acquire the sequence to blast. I think looking at the higher orders is more important as resolution at lower levels isn't as good and is well... full of unknown taxa. This I have not done as it requires extra thought and researchers vary on best practices here. 

table(tax_table(BioP.phylo)[, "Kingdom"], exclude = NULL) #9 archaea, 6923 bacteria, 4 euks, 51 unassigned #i'm removing everything but bacteria
table(tax_table(BioP.phylo)[, "Phylum"], exclude = NULL) #examine 
table(tax_table(BioP.phylo)[, "Class"], exclude = NULL) #examine 
table(tax_table(BioP.phylo)[, "Order"], exclude = NULL) #303 reads assigned as chloroplast at this taxonomic rank
table(tax_table(BioP.phylo)[, "Family"], exclude = NULL) #225 reads assigned as mitochondria at this taxonomic rank
table(tax_table(BioP.phylo)[, "Genus"], exclude = NULL) 
table(tax_table(BioP.phylo)[, "species"], exclude = NULL) #2 unidentified

p1<-subset_taxa(BioP.phylo,  !Kingdom %in% "Unassigned") #requires dplyr (%in% syntax)
p2<-subset_taxa(p1,  !Kingdom %in% "Archaea") #be no more!
p3<-subset_taxa(p2,  !Kingdom %in% "Eukaryota") #sayonara euks!
p4<-subset_taxa(p3,  !Order %in% " Chloroplast") #Chlorplasts be gone! 
p5<-subset_taxa(p4,  !Family %in% " Mitochondria") #Adios amigos!

#verify
table(tax_table(p4)[, "Order"], exclude = NULL) #si
table(tax_table(p5)[, "Family"], exclude = NULL) #si

#removing a true and verified duplicate sample and renaming the R object, but need to examine read depth first for JRD657 and JRD669. Take the higher.
rowSums(p5@otu_table) #JRD657= 3 and JRD669= 155721, note 654 and 666 have zeros because they might have lost reads to high host gene content. Review line 96 in comparison (i.e. before the taxonomic filtering). Looks like the had low reads in general, these are ones to be reseqed. 

# BioPlastics_phylo<-subset_samples(p5, sample_name != "JRD657") #official object
# dim(BioPlastics_phylo@otu_table) #149 x 6,395 one sample did fail to get seqed entirely hence 149 instead of 150. 

BioPlastics_phylo<-subset_samples(p5, sample_name != "JRD654" & sample_name != "JRD657" & sample_name != "JRD666") #official object
dim(BioPlastics_phylo@otu_table) #147 x 6,395 removed three samples on the line above with low read depth to avoid rarefying to 0 or 3 reads. You can see how low these read depths are (hence us wanting to reseq them with the code on L136).

BioPlasticsASVtable<-BioPlastics_phylo@otu_table #  Use this to rarefy on line 162 and then follow to build a Quant. Jaccard Distance for modeling with dbrda (unless Jon wants to model other ways)...
dim(BioPlasticsASVtable) #147 x 6,395 -- dimensions match, need to check if we need to remove suprious zeros found in colSums due to the removal of samples with low read depth. 
min(colSums(BioPlasticsASVtable)) #minimum is 1, we are OK unless you want to be hell bent on removing singletons. No need for the lines below.
#BioPlasticsASVtable1<-BioPlasticsASVtable[,-(which(colSums(BioPlasticsASVtable)==0))] #you can replace the zero with a one here if Jon doesn't like singletons.
#dim(BioPlasticsASVtable1)
#colSums(BioPlasticsASVtable1)
#subset.obj<-colnames(BioPlasticsASVtable1) #shows us which samples to keep for pruning prior to weighted unifrac. 
# BioPlastics_phylo1<-prune_taxa(taxa=subset.obj, x=BioPlastics_phylo@phy_tree) #does to the phylogenetic tree
# str(BioPlastics_phylo1)


#################################################################################
#################################################################################
#################################################################################
#######    data analysis can begin below this line break     ####################
#################################################################################
#################################################################################
#################################################################################

#this is how you pull out an ASV table as a R object
#BioP_ASV_table<-BioPlastics_phylo@otu_table #use the @ sign to navigate between table, metadata, tree, and taxonomy info

#rarefaction
#READ: so many people have so many POV on this, myself included. Doing it once like below is pretty standard, though I could argue since its permutational process that doing it one thousand or more times and taking the mean abundance and comparing it to a non rarefied table is an approach to see if rarefaction is necessary at all. Then some could argue about where to take the mean. Here. or the mean test stat. BUT I digress, for now.  

#find the minimum read depth, notice now its not zero because I removed the samples with low depth that will be reqeuenced. In the final data set you wont remove samples and will use whatever the lowest read depth is on the line below. 30,023.
min(rowSums(BioPlasticsASVtable)) #0 because of shit sample, when that doesn't exist you will use this line to rarefy down to the lowest number of reads. You can eyeball which is the lowest number for now. A quick look would be JRD582.
BioP_ASV_table.df<-as.data.frame(BioP_ASV_table)
rdat<-rrarefy(BioP_ASV_table.df,30023) #the lowest number of reads goes here. If you wish to proceed with stats you will need to also remove the samples (654,657,666) before rarefaction as done on line 123. 

#If you care, we can have a chat about rarefying (doing it once) vs rarefactions (doing it 1000X). I'm team rarefaction, but for a preliminary dataset out of QIIME2. No need to do extra work right now. read: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003531 and the response this year https://journals.asm.org/doi/full/10.1128/msphere.00355-23

#standardize abundances into proportions. 
rowSums(rdat) #rarefaction worked
std.tab<-decostand(rdat,"total")
std.tab[1:10,1:10]

#create quantitative Jaccard distance - takes into account relative abundance of OTUs vs presence/absence
drdat<-vegdist(std.tab,"jaccard")

#to start ... a distance based redundancy analysis on the quantitative jaccard distance matrix predicted by your treats.
mod<-dbrda() #https://uw.pressbooks.pub/appliedmultivariatestatistics/chapter/rda-and-dbrda/ about halfway down the page. 
h<-how(nperm=10000)
anova(mod,permutations = h, by="margin") 
summary(mod)

#UniFrac
#Weighted UniFrac distance
registerDoParallel(cores=4)
# BioPlastics_phylo_obj<-phyloseq(otu_table(rdat, taxa_are_rows=FALSE), phy_tree(BioPlastics_phylo1)) #rebuild for weighted unifrac, use the rarefied otu table rdat here that was made to create the quant jacquard distance… 
# dim(BioPlastics_phylo_obj@otu_table)
# length(BioPlastics_phylo_obj@phy_tree$tip.label)
#use Phylo.weight to make weighted UniFrac distance with UniFrac()
UniFrac_Phylo_Object<-phyloseq(otu_table(rdat, taxa_are_rows=FALSE), phy_tree(BioPlastics_phylo@phy_tree)) #uses rarefied data in OTU table. A MUST!

wdistUni<-UniFrac(physeq=UniFrac_Phylo_Object, weighted=TRUE, parallel=TRUE, fast=TRUE) #weighted #warning about randomly selected root for tree resolved. It's OK.

mod.Uni<-dbrda() #dbrda on weighted unifrac distance predicted by your treats
anova(mod.Uni, permutations =h, by="margin")
summary(mod.Uni)





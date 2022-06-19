library(ggplot2)
library(patchwork)
library(uwot)
library(RColorBrewer)
library(gplots)

# read data

ileum.data <- read.table("ileum/ileum_raw_UMIcounts.txt", header = TRUE, sep = "\t")
row.names(ileum.data) <- ileum.data$GENE
ileum.data <- subset(ileum.data, select = -c(GENE))

ileum.data <- ileum.data[which(rowSums(ileum.data) > 0), which(colSums(ileum.data) > 0)]
ileum.data <- t(ileum.data)
ileum.data <- as.data.frame(ileum.data)

# select genes and randomly sample data

ileum.data.1 <- subset(ileum.data, select = c(ACE2,ANPEP,ENPEP,DPP4))
data.selected <- sample(1:nrow(ileum.data.1), floor(nrow(ileum.data.1)/5), replace = FALSE)
data.selected <- sort(data.selected, decreasing = FALSE)
ileum.data.1 <- ileum.data.1[data.selected,]
#ileum.data.1 <- ileum.data.1[seq(2, nrow(ileum.data.1), 5),]
ileum.data <- subset(ileum.data, select = -c(ACE2,ANPEP,ENPEP,DPP4))
data.selected.2 <- sample(1:ncol(ileum.data), floor(ncol(ileum.data)/5), replace = FALSE)
data.selected.2 <- sort(data.selected.2, decreasing = FALSE)
ileum.data.2 <- ileum.data[data.selected, data.selected.2]
#ileum.data.2 <- ileum.data[seq(2, nrow(ileum.data), 5), seq(4, ncol(ileum.data), 5)]
ileum.data <- cbind.data.frame(ileum.data.1, ileum.data.2)

# scale data

ileum.data[is.na(ileum.data)] <- 0
ileum.data <- scale(ileum.data)

# read cell info

cell.info <- read.csv("ileum/ileum_cell_info.txt", sep = "\t")
rownames(cell.info) <- cell.info$UniqueCell_ID
cell.info <- subset(cell.info, select = -c(UniqueCell_ID))
rownames(cell.info) <- gsub("-", ".", rownames(cell.info))

# UMAP and visualization

ileum.data[is.na(ileum.data)] <- 0
umap.ileum <- umap(ileum.data, scale = "Z")
df.umap.ileum <- data.frame(umap.ileum)
colnames(df.umap.ileum) <- c("UMAP1","UMAP2")
rownames(df.umap.ileum) <- rownames(ileum.data)
cell.info.select <- cell.info[rownames(df.umap.ileum),]
df.umap.ileum <- cbind.data.frame(df.umap.ileum, cell.info.select)
ggplot(df.umap.ileum, aes(x = UMAP1, y = UMAP2, color=CellType, shape=Sample_ID)) + geom_point()
ggplot(df.umap.ileum, aes(x = UMAP1, y = UMAP2, color=Sample_ID)) + geom_point()

df.umap.ileum[["ACE2"]] <- scale(ileum.data.1$ACE2)
df.umap.ileum[["ANPEP"]] <- scale(ileum.data.1$ANPEP)
df.umap.ileum[["ENPEP"]] <- scale(ileum.data.1$ENPEP)
df.umap.ileum[["DPP4"]] <- scale(ileum.data.1$DPP4)
ggplot(df.umap.ileum, aes(x = UMAP1, y = UMAP2, color=ACE2, shape=Sample_ID)) + geom_point() + scale_color_gradient(low = "blue", high = "red")

umap.ileum.1 <- df.umap.ileum[df.umap.ileum$Sample_ID == "Ileum-1" | df.umap.ileum$Sample_ID == "Ileum-2",]
umap.colon <- df.umap.ileum[df.umap.ileum$Sample_ID == "Colon-1" | df.umap.ileum$Sample_ID == "Colon-2",]
umap.rectum <- df.umap.ileum[df.umap.ileum$Sample_ID == "Rectum-1" | df.umap.ileum$Sample_ID == "Rectum-2",]

ggplot(umap.ileum.1, aes(x = UMAP1, y = UMAP2, color=CellType)) + geom_point()
ggplot(umap.ileum.1, aes(x = UMAP1, y = UMAP2, color=ACE2, shape=CellType)) + geom_point() + scale_color_gradient(low = "blue", high = "red")
ggplot(umap.ileum.1, aes(x = UMAP1, y = UMAP2, color=ANPEP, shape=CellType)) + geom_point() + scale_color_gradient(low = "blue", high = "red")
ggplot(umap.ileum.1, aes(x = UMAP1, y = UMAP2, color=ENPEP, shape=CellType)) + geom_point() + scale_color_gradient(low = "blue", high = "red")
ggplot(umap.ileum.1, aes(x = UMAP1, y = UMAP2, color=DPP4, shape=CellType)) + geom_point() + scale_color_gradient(low = "blue", high = "red")

ggplot(umap.colon, aes(x = UMAP1, y = UMAP2, color=CellType)) + geom_point()
ggplot(umap.colon, aes(x = UMAP1, y = UMAP2, color=ACE2, shape=CellType)) + geom_point() + scale_color_gradient(low = "blue", high = "red")
ggplot(umap.colon, aes(x = UMAP1, y = UMAP2, color=ANPEP, shape=CellType)) + geom_point() + scale_color_gradient(low = "blue", high = "red")
ggplot(umap.colon, aes(x = UMAP1, y = UMAP2, color=ENPEP, shape=CellType)) + geom_point() + scale_color_gradient(low = "blue", high = "red")
ggplot(umap.colon, aes(x = UMAP1, y = UMAP2, color=DPP4, shape=CellType)) + geom_point() + scale_color_gradient(low = "blue", high = "red")

ggplot(umap.rectum, aes(x = UMAP1, y = UMAP2, color=CellType)) + geom_point()
ggplot(umap.rectum, aes(x = UMAP1, y = UMAP2, color=ACE2, shape=CellType)) + geom_point() + scale_color_gradient(low = "blue", high = "red")
ggplot(umap.rectum, aes(x = UMAP1, y = UMAP2, color=ANPEP, shape=CellType)) + geom_point() + scale_color_gradient(low = "blue", high = "red")
ggplot(umap.rectum, aes(x = UMAP1, y = UMAP2, color=ENPEP, shape=CellType)) + geom_point() + scale_color_gradient(low = "blue", high = "red")
ggplot(umap.rectum, aes(x = UMAP1, y = UMAP2, color=DPP4, shape=CellType)) + geom_point() + scale_color_gradient(low = "blue", high = "red")
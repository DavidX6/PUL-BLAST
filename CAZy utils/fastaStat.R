install.packages("seqinr")
library("seqinr")

seqs <- read.fasta(file = "CAZyDB.07312019.fa", seqtype = "AA", seqonly = TRUE, whole.header = TRUE)
seqlens <- nchar(seqs)
a <- hist(x = seqlens, breaks = 20)

hist(x = seqlens[seqlens <= 2000], breaks = 20, xlab = "dolžine sekvenc", ylab = "število", main = "Histogram sekvenc krajših od 2000")
hist(x = seqlens[seqlens > 2000 & seqlens <= 5000], breaks = 20, xlab = "dolžine sekvenc", ylab = "število", main = "Histogram sekvenc dolžine med 2000-5000")
hist(x = seqlens[seqlens > 5000 & seqlens <= 10000], breaks = 20, xlab = "dolžine sekvenc", ylab = "število", main = "Histogram sekvenc dolžine med 5000-10000")
hist(x = seqlens[seqlens > 10000], breaks = 20, xlab = "dolžine sekvenc", ylab = "število", main = "Histogram sekvenc daljših od 10000")


min(seqlens)
which(seqlens == min(seqlens))
seqs[[556158]]
# ----------------------------------------
seqs <- read.fasta(file = "../izClankovZAccZimeni.fasta", seqtype = "AA", seqonly = TRUE, whole.header = TRUE)
seqlens <- nchar(seqs)
a <- hist(x = seqlens)
min(seqlens)
which(seqlens == min(seqlens))
seqs[[61]]

# ----------------------------------------
seqs <- read.fasta(file = "../queryTempSequence.txt", seqtype = "AA", seqonly = TRUE, whole.header = TRUE)
seqlens <- nchar(seqs)
hist(x = seqlens[seqlens <= 200], breaks = 20, xlab = "dolžine sekvenc", ylab = "število", main = "Histogram sekvenc krajših od 2000")
a <- hist(x = seqlens)
min(seqlens)
which(seqlens == min(seqlens))
seqs[[61]]

#
seqs <- read.fasta(file = "../PULDB_merged.fasta", seqtype = "AA", seqonly = TRUE, whole.header = TRUE)
seqlens <- nchar(seqs)
a <- hist(x = seqlens)
min(seqlens)
which(seqlens == min(seqlens))
seqs[[138]]

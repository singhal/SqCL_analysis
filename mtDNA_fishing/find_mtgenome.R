library(ape)

t = read.tree("/Users/Sonal/macroevolution/published_data/Pyron2013_squamates.tre")
sp = read.csv("/Volumes/sosi/brazil/samples.csv", stringsAsFactors=F)
mtsp = names(read.dna("/Users/sonal/Desktop/reptile_mitogenomes.fa", format="fasta"))
mtg = unique(gsub("_\\S+", "", mtsp))
tg = gsub("_\\S+", "", t$tip.label)

x = cophenetic(t)
x = x[,f]

d = data.frame(sp$lineage, sp$sample, rep(NA, nrow(sp)), stringsAsFactors=F)
names(d) = c("species", "sample", "genome")

for (i in 1:nrow(d)) {
	sp = d[i, "species"]
	genus = gsub("_\\S+", "", sp)
	
	if (length(grep(genus, mtsp)) > 0) {
		d[i, "genome"] = mtsp[grep(genus, mtsp)][1]
	} else {
	
		tmp = x[grep(genus, dimnames(x)[[1]]),]
		if (length(tmp) > 0) {
			inds = which(tmp == min(tmp), arr.ind=TRUE)
		
			if (length(inds) > 1) {
				keep = colnames(tmp)[inds[,2]]
			} else {
				keep = names(tmp)[inds[1]]
			}
		keep = mtsp[grep(gsub("_\\S+", "", keep), mtsp)]
		d[i, "genome"] = keep[1]
		}
	}	
}
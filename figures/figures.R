library(ape)
library(geosphere)
library(RColorBrewer)
library(scales)
library(zoo)
library(AICcmodavg)
library(ggplot2)
library(gridViewer)
library(RSQLite)
library(grid)
library(gridExtra)
library(vegan)
library(ade4)
library(ggtree)
library(phytools)
library(phangorn)

library(treescape)

######################
# the data dumps
######################

indir = "/Volumes/sosi/brazil/"
outdir =  "/Users/sonal/publications/SqCL/figures/"
moutdir =  "/Users/sonal/publications/SqCL/main_figures/"

c = read.csv(paste(indir, "metadata/collection.csv", sep=""), stringsAsFactors=F)
s = read.csv(paste(indir, "metadata/samples.csv", sep=""), stringsAsFactors=F)
sites = read.csv(paste(indir, "metadata/sites.csv", sep=""), stringsAsFactors=F)
t = read.csv(paste(indir, "metadata/taxonomy.csv", sep=""), stringsAsFactors=F)

cols = brewer.pal(3, "Dark2")

######################
# functions
######################

get_cumulative <- function(cut, vals) {
	return(sum(vals > cut) / length(vals))
}

put.fig.letter <- function(label, location="topleft", x=NULL, y=NULL, 
                           offset=c(0, 0), ...) {
  if(length(label) > 1) {
    warning("length(label) > 1, using label[1]")
  }
  if(is.null(x) | is.null(y)) {
    coords <- switch(location,
                     topleft = c(0.015,0.98),
                     topcenter = c(0.5525,0.98),
                     topright = c(0.985, 0.98),
                     bottomleft = c(0.015, 0.02), 
                     bottomcenter = c(0.5525, 0.02), 
                     bottomright = c(0.985, 0.02),
                     c(0.015, 0.98) )
  } else {
    coords <- c(x,y)
  }
  this.x <- grconvertX(coords[1] + offset[1], from="nfc", to="user")
  this.y <- grconvertY(coords[2] + offset[2], from="nfc", to="user")
  text(labels=label[1], x=this.x, y=this.y, xpd=T, ...)
}

###############################
#  data completeness
###############################

d = read.csv(paste(indir, "metadata/locus_data.csv", sep=""), stringsAsFactors=F)
d$type = gsub("-.*$", "", d$locus)

cuts = seq(0, 1, by=0.01)
all = sapply(cuts, get_cumulative, vals=d$missingness)
types = split(d, d$type)
res = vector("list", length(types))~/publications/SqCL/SqCL_revisions.Md
for (i in 1:length(types)) {
	res[[i]] = sapply(cuts, get_cumulative, vals=types[[i]]$missingness)
}
names(res) = names(types)

pdf(file=paste(outdir, "completeness.pdf", sep=""), width=3, height=3, pointsize=8)
par(mar=c(4,5,1,1), tck=-0.01)
plot(NULL, xlim=c(0,1), ylim=range(0,1), axes=F, xlab="", ylab="")
lines(cuts, all)
lines(cuts, res[[1]], col=cols[1], lwd=2)
lines(cuts, res[[2]], col=cols[2], lwd=2)
lines(cuts, res[[3]], col=cols[3], lwd=2)
axis(1, at=c(0, 0.25, 0.5, 0.75, 1), labels=c(0, 0.25, 0.5, 0.75, 1), cex.axis=1.3)
mtext("proportion complete", side=1, line=2.8, cex = 1.4)
axis(2, at=c(0, 0.25, 0.5, 0.75, 1), labels=c(0, 0.25, 0.5, 0.75, 1), las=1, cex.axis=1.3)
mtext("proportion of loci", side=2, line=3.4, cex=1.4)

labels = c("all loci", names(types))
labcols = c("black", cols)
for (i in 1:length(labels)) {
	text(0.1, 0.3 - 0.075 * (i - 1), labels[i], col=labcols[i], cex=1.5)
}
dev.off()

###############################
# coverage correlation
###############################

d = read.csv(paste(indir, "coverage/coverage_by_locus_sample.csv", sep=""), stringsAsFactors=F)
to_plot = list(c(NA, NA), c(NA, NA), c(NA, NA))
ids = names(d)[2:length(names(d))]
sps = t[match(c[match(ids, c$sample), "brazil_species"], t$brazil_species), "nominal_species"]
ids = data.frame(ids, sps)
# randomly get one per species
ids = do.call( rbind, lapply( split(ids, ids$sps), function(ids) ids[sample(nrow(ids), 1) ,]))
ids = data.frame(lapply(ids, as.character), stringsAsFactors=F)

for (i in 1:length(to_plot)) {
	to_plot[[i]] = sample(ids$ids, 2)
}

get_corr = function(x) {
	return(cor.test(d[, x[1]], d[, x[2]])$estimate)
}

comp = combn(ids$ids, 2, simplify=F)
corr = unlist(lapply(comp, get_corr))

letters = c("A.", "B.", "C.")
pdf(file=paste(outdir, "coverage_correlation.pdf", sep=""), width=7, height=2, pointsize=8)
par(mfrow=c(1,4), mar=c(4,6,1,1), tck=-0.01, las=0, cex.axis=1.4, cex.lab=1.4)
for (i in 1:length(to_plot)) {
	sp1 = to_plot[[i]][1]
	sp2 = to_plot[[i]][2]
		
	xvals = round(exp(pretty(log(d[, sp1]), 3)), -1)
	xvals = c(1, unique(xvals[xvals > 0]))
	atx = log(xvals)
	
	yvals = round(exp(pretty(log(d[, sp2]), 3)), -1)
	yvals = c(1, unique(yvals[yvals > 0]))
	aty = log(yvals)
	
	plot(log(d[, sp1]), log(d[, sp2]), pch=16, cex=0.5, col=alpha("black", 0.4), axes=F, xlab="", ylab="", xlim=c(0, max(atx)), ylim=c(0, max(aty)))
	
	r = cor.test(d[, sp1], d[, sp2])
	
	axis(1, at=atx, labels=xvals)
	mtext(ids[ids$ids == sp1, "sps"], side=1, line=2.8)
	axis(2, at=aty, labels=yvals, las=1)
	mtext(ids[ids$ids == sp2, "sps"], side=2, line=3.4)
	put.fig.letter(letters[i], offset=c(0.025,-0.025), cex=2)
	text(max(atx) * 0.9, 0.6, paste("r=",round(r$estimate, 2)), cex=1.4)
}
hist(corr, xlab="coverage correlation", ylab="frequency", col="gray", border=NA, breaks=20, main="", las=1)
put.fig.letter("D.", offset=c(0.025,-0.025), cex=2)
dev.off()

x = d
d = x[grep("uce", x$locus),]
comp = combn(ids$ids, 2, simplify=F)
corr = unlist(lapply(comp, get_corr))
mean(corr)
d = x[grep("AHE", x$locus),]
comp = combn(ids$ids, 2, simplify=F)
corr = unlist(lapply(comp, get_corr))
mean(corr)
d = x[grep("gene", x$locus),]
comp = combn(ids$ids, 2, simplify=F)
corr = unlist(lapply(comp, get_corr))
mean(corr)

###############################
# coverage uniformity
###############################

d = read.csv(paste(indir, "coverage/coverage_by_locus_sample.csv", sep=""), stringsAsFactors=F)
cv = apply(d[,2:ncol(d)], 2, sd) / apply(d[,2:ncol(d)], 2, mean)

pdf(file=paste(outdir, "coverage_uniformity.pdf", sep=""), width=3, height=3, pointsize=8)
par(mar=c(4,5.5,1,1), tck=-0.01, las=0, cex.axis=1.4, cex.lab=1.4)
xvals = pretty(cv, 3)
yvals = pretty(hist(cv, plot=F, breaks=20)$counts, 3)

hist(cv, xlab="c.v. in coverage", ylab="frequency", col="gray", border=NA, breaks=20, main="", las=1, axes=F, xlim=range(xvals), ylim=range(yvals))
axis(1, at=xvals, labels=xvals)
axis(2, at=yvals, labels=yvals, las=1)
dev.off()

x = d
d = x[grep("uce", x$locus),]
cv = apply(d[,2:ncol(d)], 2, sd) / apply(d[,2:ncol(d)], 2, mean)
mean(cv)
d = x[grep("AHE", x$locus),]
cv = apply(d[,2:ncol(d)], 2, sd) / apply(d[,2:ncol(d)], 2, mean)
mean(cv)
d = x[grep("gene", x$locus),]
cv = apply(d[,2:ncol(d)], 2, sd) / apply(d[,2:ncol(d)], 2, mean)
mean(cv)

###############################
# loci that just failed
###############################

d = read.csv(paste(indir, "coverage/coverage_by_locus_sample.csv", sep=""), stringsAsFactors=F)
count_zero <- function(x) {return(length(x[x < 10]))}
d$zero = apply(d[,2:57], 1, count_zero)
miss = d[d$zero == 56, "locus"]

length(miss[grep("AHE", miss)]) 
length(miss[grep("AHE", miss)]) / nrow(d[grep("AHE", d$locus), ])
length(miss[grep("uce", miss)])
length(miss[grep("uce", miss)]) / nrow(d[grep("uce", d$locus), ])
length(miss[grep("gene", miss)])
length(miss[grep("gene", miss)]) / nrow(d[grep("gene", d$locus), ])

###############################
# loci locations on genome
###############################

width = 10
d = read.csv(paste(indir, "ref_genomes/Anolis_brasiliensis_anoCar2.locations", sep=""), stringsAsFactors=F)
l = read.table(paste(indir, "ref_genomes/chromInfo.txt", sep=""), stringsAsFactors=F)
d$chrs = d$start / max(l$V2) * (width - 0.5) + 0.25
d$col = rep(NA, nrow(d))
d$y = rep(NA, nrow(d))
d[grep("AHE", d$locus), "col"] = cols[1]
d[grep("gene", d$locus), "col"] = cols[2]
d[grep("uce", d$locus), "col"] = cols[3]
d[grep("AHE", d$locus), "y"] = 0.01
d[grep("gene", d$locus), "y"] = 0.01
d[grep("uce", d$locus), "y"] = 0.04

counts = table(d$chr)

pdf(file=paste(outdir, "genome_locations.pdf", sep=""), width=7, height=5, pointsize=8)
chrs = names(counts[counts > 24])
par(mfrow=c(length(chrs), 1), mar=c(1, 0, 1, 0))
for (i in 1:length(chrs)) {
	plot(NULL, xlim=c(0, width + 0.5), ylim=c(0, 0.05), axes=F)
	chrw = (l[l$V1 == chrs[i], "V2"] / max(l$V2)) * (width - 0.5)
	
	text(chrw + 0.4, 0.025, chrs[i], adj=c(0, 0.5), cex=1.2)
	
	x = c(0.25, chrw + 0.25, chrw + 0.25, 0.25)
	y = c(0.01, 0.01, 0.04, 0.04)
	lines(c(0.25, chrw + 0.25), y = c(0.025, 0.025), col=alpha("gray", 0.3), lwd=5)
	# polygon(x, y, border="black", col="white")
	
	tmp = d[d$chr == chrs[i],]
	
	draw_line <- function(x) {
		points(x["chrs"], x["y"], col=alpha(x["col"], 0.5), pch=16)
	}
	
	apply(tmp, 1, draw_line)
}
dev.off()

############################
# capture efficiency etc
############################

get_sample <- function(x) {
	rand = sample(seq(1, nrow(x)), 1)
	return(x[rand, "lineage"])
}

get_vals <- function(d, x, tips, name) {
	df1 = aggregate(d[, x], by=list(d$family), median)
	df2 = aggregate(d[, x], by=list(d$family), quantile, c(0.05, 0.95))
	
	df = cbind(df1, as.data.frame(df2[,2]))
	names(df) = c("family", "mean", "low", "high")
	rownames(df) = df$family
	df = df[tips,]
	
	par(mar=c(5, 2, 1, 2), tck=-0.02)
	xvals = pretty(range(df[,2:4]), 2)
	plot.new()
	plot.window(xlim=range(xvals), ylim=c(1, nrow(df)))
		
	for (i in 1:nrow(df)) {
		pt1 = c(df[i,3], df[i,4])
		pt2 = c(i, i)
		lines(pt1, pt2, col="black")
		points(df[i, 2], i, pch=21, bg="gray", cex=1)
	}
	axis(1, at=xvals, labels=NA, line=1)
	axis(1, at=xvals, labels=xvals, line=0.4, lwd = 0, cex.axis=1)
	mtext(name, side=1, line=3.3, cex=1)
}

d = read.csv(paste(indir, "coverage/summary_statistics.csv", sep=""), stringsAsFactors=F)
x = read.csv(paste(indir, "metadata/assembly_info.csv", sep=""), stringsAsFactors=F)
d = merge(d, x, by="sample")
d$lineage = gsub("_", " ", d$lineage)
d$family = t[match(d$lineage, t$nominal_species), "family"]

tree = read.nexus("~/publications/SqCL/revision_results/div_dating/AHE_L1_sub.tre")
tree$tip.label = gsub("_", " ", tree$tip.label)
x = split(d, d$family)
keep = unlist(lapply(x, get_sample))
drop = setdiff(tree$tip.label, keep)
tree = drop.tip(tree, drop)
tree$tip.label = d[match(tree$tip.label, d$lineage), "family"]
tree = ladderize(tree)
tree = drop.tip(tree, "Gallus_gallus")

tiporder = tree$edge[which(tree$edge[,2] <= length(tree$tip.label)), 2]
tips = tree$tip.label[tiporder]

pdf(paste(outdir, "data_quality.pdf", sep=""), width=7, height=3, pointsize=8)
par(mfrow=c(1,5))
par(mar=c(5,0,1,0))
plot(tree, edge.width=0.6, cex=1.2)
get_vals(d, "map_reads", tips, "capture efficiency")
get_vals(d, "all_num", tips, "# of loci")
get_vals(d, "all_mean_length", tips, "mean locus length")
get_vals(d, "all_mean_cov", tips, "mean coverage")
dev.off()

############################
# capture efficiency etc by locus
############################

get_vals2 <- function(d, x, tips, name) {
	types = c("AHE", "gene", "uce")
	vals = paste(types, x, sep="_")
	counts = c(372, 38, 5032)
	
	df = data.frame(sort(unique(d$family)))
	names(df) = c("family")
	for (i in 1:length(vals)) {
		y = vals[i]
		
		df1 = aggregate(d[, y], by=list(d$family), median)
		if (x == "num") {
			df1[ ,2] = df1[ ,2] / counts[i]
		}
		
		df2 = aggregate(d[, y], by=list(d$family), quantile, c(0.05, 0.95))
		if (x == "num") {
			df2[ ,2] = df2[ ,2]  / counts[i]
		}
		
		df = cbind(df, df1[ , 2], as.data.frame(df2[,2]))
		n1 = c("mean", "low", "high")
		names = paste(types[i], n1, sep="_")
		names(df)[(ncol(df) - 2):ncol(df)] = names	
		}
	
	rownames(df) = df$family
	df = df[tips,]
	
	par(mar=c(6, 2, 2, 2), tck=-0.02)
	xvals = pretty(range(df[,2:ncol(df)]), 3)
	
	for (j in 1:length(types)) {
		plot.new()
		plot.window(xlim=range(xvals), ylim=c(1, nrow(df)))		
		for (i in 1:nrow(df)) {
			pt1 = c(df[i, paste(types[j], "low", sep="_")], df[i, paste(types[j], "high", sep="_")])
			pt2 = c(i, i)
			lines(pt1, pt2, col="black")
			points(df[i,  paste(types[j], "mean", sep="_")], i, pch=21, bg=cols[j], cex=1.3)
			}
		axis(1, at=xvals, labels=NA, line=1)
		axis(1, at=xvals, labels=xvals, line=0.4, lwd = 0, cex.axis=1)
		mtext(paste(name, ", ", types[j], sep=""), side=1, line=3.3)
		}
}

pdf(paste(moutdir, "2data_quality_by_loci.pdf", sep=""), width=7, height=6, pointsize=8)
par(mfrow=c(3,4), mar=c(6,0,2,1))
plot(tree, edge.width=0.6, cex=1.2)
get_vals2(d, "num", tips, "% loci recovered")
plot(tree, edge.width=0.6, cex=1.2)
get_vals2(d, "mean_length", tips, "mean length")
plot(tree, edge.width=0.6, cex=1.2)
get_vals2(d, "mean_cov", tips, "mean cov.")
dev.off()

############################
# corr plot
############################

x = read.csv(paste(indir, "metadata/species_diversity.csv", sep=""), stringsAsFactors=F)
x = x[x$type == "IND",]
d$pi = x[match(d$sample, x$ind), "pi"]
d$X260_280 = s[match(d$sample, s$sample), "X260_280"]

keep = names(d)
keep = keep[grep("AHE", keep, invert=T)]
keep = keep[grep("uce", keep, invert=T)]
keep = keep[grep("gene", keep, invert=T)]
keep = keep[keep != "sample"]
keep = keep[keep != "lineage"]
keep = keep[keep != "family"]
keep = keep[keep != "all_median_length"]
keep = keep[keep != "all_sites_10x"]
keep = keep[keep != "assembled_length"]
keep = keep[keep != "assembled_n50"]
keep = keep[keep != "on_target"]

x = d[,keep]
x$paired = x$paired / x$map_reads
x$duplicates = x$duplicates / x$map_reads
names(x) = c("mean cov.", "mean length", "# of loci", "duplicates", "capture efficiency", "# raw reads", "paired", "# contigs", "heterozygosity", "260/280")

pdf(paste(outdir, "correlations.pdf", sep=""), width=5, height=5, pointsize=8)
M = cor(x)
corrplot.mixed(M, upper="ellipse", tl.col="black", tl.cex=0.5)
dev.off()

############################
# matrix
############################

d = read.csv("/Volumes/sosi/brazil/metadata/heat_map.csv", header=T, row.names=1)
counts = apply(d, 2, sum)
loci = names(counts[order(counts, decreasing=T)])

for (i in 1:ncol(d)) {
	name = gsub("\\.\\S+", "", names(d)[i])
	if (name == "AHE") {
		factor = 1
	}
	if (name == "gene") {
		factor = 2
	}
	if (name == "uce") {
		factor = 3
	}
	d[,i] = d[,i] * factor
}

tree = read.tree("~/publications/SqCL/revision_results/trees/best_trees_tol1e-05_collapse1_AHE_ASTRAL.tre")
tree = root(tree, "Gallus_gallus")
tree2 = read.tree("~/publications/SqCL/data/starting_tree/astral_AHE_partitioned_chrono.tre")
tree = drop.tip(tree, "Gallus_gallus")
tree2 = drop.tip(tree2, "Gallus_gallus")
write.tree(ladderize(tree), "~/Desktop/test.tre")
tree = read.tree("~/Desktop/test.tre")
write.tree(ladderize(tree2), "~/Desktop/test.tre")
tree2 = read.tree("~/Desktop/test.tre")

m = as.matrix(d)
m = m[tree$tip.label,]
m = t(m)
m = m[loci,]

pdf(paste(moutdir, "1heat_map.pdf", sep=""), width=6, height=4, pointsize=8)
par(mfrow=c(1,2))
par(mar = c(3.5, 0, 0, 2))
plot.phylo(tree2, direction = "rightwards", show.tip.label = TRUE, cex=0.6, y.lim = c(1, length(tree2$tip.label)))

tree$node.label = as.numeric(tree$node.label)
for (i in 1:length(tree$node.label)) {
	if (!is.na(tree$node.label[i])) {
		if (tree$node.label[i] >= 0.95) {
			nodelabels(node=i + length(tree$tip.label), bg="gray", pch=21, cex=0.6, lwd=0.5)
		}
	}
}

par(mar = c(4, 0, 0.5, 0))
gl <- 1:(length(tree$tip.label) + 1)
image(1:(dim(m)[1] + 1), 1:(length(tree$tip.label) + 1), m, col = c("lightgray", cols), xlim = c(1, dim(m)[1] - 1), ylim = c(1, length(tree$tip.label) - 1), axes=F, xlab="")
xvals = pretty(c(1,dim(m)))[1:6]
axis(1, at=xvals, labels=xvals)
mtext("number of loci", 1, line=2.5)
dev.off()

############################
# PICs
############################

get_cov <- function(a, b) {
	v = seq(min(a$rel_ungap_pos), max(a$rel_ungap_pos))
	cov = rep(NA, length(v))
	
	get_count <- function(x) {
		tmp = b[b$min <= x,]
		tmp = tmp[tmp$max >= x,]
		return(nrow(tmp))
	}
	
	cov = sapply(v, get_count)
	names(cov) = v
	
	p = table(a$rel_ungap_pos)
	p = p[names(cov)]
	names(p) = names(cov)
	p[is.na(p)] = 0
	
	res = list(cov, p)
	return(res)
}

get_cov_type <- function(a, b, type) {
	a = a[grep(type, a$locus), ]
	b = b[grep(type, b$locus), ]
	v = seq(min(a$rel_ungap_pos), max(a$rel_ungap_pos))
	cov = rep(NA, length(v))
	
	get_count <- function(x) {
		tmp = b[b$min <= x,]
		tmp = tmp[tmp$max >= x,]
		return(nrow(tmp))
	}
	
	cov = sapply(v, get_count)
	names(cov) = v
	
	p = table(a$rel_ungap_pos)
	p = p[names(cov)]
	names(p) = names(cov)
	p[is.na(p)] = 0
	
	res = list(cov, p)
	return(res)
}

get_var <- function(a, b, width) {
	v = seq(min(a$relgap), max(a$relgap))
	cov = rep(NA, length(v))
	
	get_count <- function(x) {
		tmp = b[b$min <= x,]
		tmp = tmp[tmp$max >= x,]
		return(nrow(tmp))
	}
	
	cov = sapply(v, get_count)
	names(cov) = v
	
	p = table(a$relgap)
	p = p[names(cov)]
	names(p) = names(cov)
	p[is.na(p)] = 0
	
	var = p / cov
	var = rollapply(var, width, mean)
}

get_var_type <- function(a, b, type, width) {
	a = a[grep(type, a$locus), ]
	b = b[grep(type, b$locus), ]
	
	v = seq(min(a$relgap), max(a$relgap))
	cov = rep(NA, length(v))
	
	get_count <- function(x) {
		tmp = b[b$min <= x,]
		tmp = tmp[tmp$max >= x,]
		return(nrow(tmp))
	}
	
	cov = sapply(v, get_count)
	names(cov) = v
	
	p = table(a$relgap)
	p = p[names(cov)]
	names(p) = names(cov)
	p[is.na(p)] = 0
	
	var = p / cov
	var = rollapply(var, width, mean)
}

d = read.csv("~/Desktop/number_of_pics.csv", stringsAsFactors=F)
d1 = d[d$type == 'ALL',]
d2 = d[d$type == 'SNAKE',]
x = read.csv(paste(indir, "metadata/locus_lengths_pics.csv", sep=""), stringsAsFactors=F)
x1 = x[x$type == 'ALL',]
x2 = x[x$type == 'SNAKE',]

v1 = read.csv(paste(indir, "metadata/Colobosaura_modesta_variants.csv", sep=""), stringsAsFactors=F)
c1 = read.csv(paste(indir, "metadata/Colobosaura_modesta_locus_lengths_variants.csv", sep=""), stringsAsFactors=F)
var1 = get_var(v1, c1, 10)

v2 = read.csv(paste(indir, "metadata/Bothrops_moojeni_variants.csv", sep=""), stringsAsFactors=F)
c2 = read.csv(paste(indir, "metadata/Bothrops_moojeni_locus_lengths_variants.csv", sep=""), stringsAsFactors=F)
var2 = get_var(v2, c2, 10)

res1 = get_cov(d1, x1)
res2 = get_cov(d2, x2)

pdf("~/Desktop/PIC_variant_density.pdf", width=6.5, height=2.5, pointsize=8)
xlim = c(-500, 500)
ylim=c(0, 0.8)
xvals = pretty(range(xlim), 3)
yvals = pretty(ylim, 4)

par(mar=c(5,5,1,1), tck=-0.02, mfrow=c(1,2))
plot(names(res1[[1]]), res1[[2]] / res1[[1]], pch=16, col="black", cex=0.5, axes=F, xlab="", ylab="", xlim=xlim, ylim=range(yvals))
points(names(res2[[1]]), res2[[2]] / res2[[1]], pch=16, col="gray", cex=0.5)
axis(1, at=xvals, labels=xvals, line=1)
mtext("nucleotide position", side=1, line=3.3)
axis(2, at=yvals, labels=yvals, line=1, las=2)
mtext("# of variable sites / # of loci", side=2, line=4)

xlim = c(-500, 500)
ylim=c(0, 0.003)
xvals = pretty(range(xlim), 3)
yvals = pretty(ylim, 4)

par(mar=c(5,7,1,1))
plot(names(var1), var1, pch=16, cex=0.5, col="black", xlab="", ylab="", axes=F, ylim=range(yvals), xlim=xlim)
points(names(var2), var2, pch=16, cex=0.5, col="gray", xlab="", ylab="")
axis(1, at=xvals, labels=xvals, line=1)
mtext("nucleotide position", side=1, line=3.3)
axis(2, at=yvals, labels=yvals, line=1, las=2)
mtext("# of SNPs / # of loci", side=2, line=5)
dev.off()

################################
# now do the pics by type
################################

types = c("AHE", "gene", "uce")
restype1 = vector("list", length(types))
restype2 = vector("list", length(types))

vartype1 = vector("list", length(types))
vartype2 = vector("list", length(types))

for (i in 1:length(types)) {
	restype1[[i]] = get_cov_type(d1, x1, types[i])
	restype2[[i]] = get_cov_type(d2, x2, types[i])
	vartype1[[i]] = get_var_type(v1, c1, types[i], 10)
	vartype2[[i]] = get_var_type(v2, c2, types[i], 10)
}

pdf(paste(moutdir, "PIC_density_types.pdf", sep=""), width=6.5, height=2.5, pointsize=10)
ylim=c(0, 0.8)
yvals = pretty(ylim, 4)
xlims = list(c(-500, 500), c(-500, 500), c(-500, 500))
par(mfrow=c(1,3), mar=c(5,7,1,1), tck=-0.02)
for (i in 1:length(types)) {
	xlim = xlims[[i]]
	xvals = pretty(range(xlim), 2)

	plot(names(restype1[[i]][[1]]), restype1[[i]][[2]] / restype1[[i]][[1]], pch=16, col="black", cex=0.5, axes=F, xlab="", ylab="", xlim=xlim, ylim=range(yvals))
	points(names(restype2[[i]][[1]]), restype2[[i]][[2]] / restype2[[i]][[1]], pch=16, col="gray", cex=0.5)
	axis(1, at=xvals, labels=xvals, line=1)
	mtext(paste("nucleotide position, ", types[[i]], sep=""), side=1, line=3.3)
	axis(2, at=yvals, labels=yvals, line=1, las=2)
	mtext("# of variable sites / # of loci", side=2, line=4)
	}
dev.off()

pdf("~/Desktop/variant_density_types.pdf", width=6.5, height=2, pointsize=10)
ylims=list(c(0, 0.003), c(0, 0.01), c(0, 0.003))
xlim = c(-500,500)
xvals = pretty(range(xlim), 2)
par(mfrow=c(1,3), mar=c(5,8,1,1), tck=-0.02)

for (i in 1:length(types)) {
	ylim = ylims[[i]]
	yvals = pretty(ylim, 3)

	plot(names(vartype1[[i]]), vartype1[[i]], pch=16, col="black", cex=0.5, axes=F, xlab="", ylab="", xlim=xlim, ylim=range(yvals))
	points(names(vartype2[[i]]), vartype2[[i]], pch=16, col="gray", cex=0.5)
	axis(1, at=xvals, labels=xvals, line=1)
	mtext(paste("nucleotide pos., ", types[[i]], sep=""), side=1, line=3.3)
	axis(2, at=yvals, labels=yvals, line=1, las=2)
	mtext("# of SNPs / # of loci", side=2, line=5)
	}
dev.off()

################################
# coverage
################################

c1 = read.csv(paste(indir, "coverage/Colobosaura_modesta_probes.csv", sep=""), stringsAsFactors=F)
c1$mean = apply(c1[,3:ncol(c1)], 1, mean)

c2 = read.csv(paste(indir, "coverage/Bothrops_moojeni_probes.csv", sep=""), stringsAsFactors=F)
c2$mean = apply(c2[,3:ncol(c2)], 1, mean)

xlim = c(-500,500)
xvals = pretty(range(xlim), 2)
ylim = c(0,150)
yvals = pretty(range(ylim), 2)

pdf("~/Desktop/locus_coverage.pdf", width=6.5, height=3, pointsize=10)
par(mfrow=c(1,2), mar=c(4,5,1,1), tck=-0.02)
plot(c1$loc, c1$mean, col=cols[as.factor(c1$type)], axes=F, xlab="", ylab="", xlim=range(xvals), ylim=range(yvals), pch=16, cex=0.4)
abline(h=10, lty="dashed")
axis(1, at=xvals, labels=xvals, line=0.2)
mtext("nucleotide position", side=1, line=2.3)
axis(2, at=yvals, labels=yvals, line=0.2, las=2)
mtext("coverage", side=2, line=3)
put.fig.letter("A.", cex=1.5, offset=c(0.05, -0.025))
plot(c2$loc, c2$mean, col=cols[as.factor(c2$type)], axes=F, xlab="", ylab="", xlim=range(xvals), ylim=range(yvals), pch=16, cex=0.4)
abline(h=10, lty="dashed")
axis(1, at=xvals, labels=xvals, line=0.2)
mtext("nucleotide position", side=1, line=2.3)
axis(2, at=yvals, labels=yvals, line=0.2, las=2)
mtext("coverage", side=2, line=3)
put.fig.letter("B.", cex=1.5, offset=c(0.05, -0.025))
dev.off()

######################
# mitochondrial data
######################

mt = list.files("/Volumes/sosi/brazil/mitogenomes/", pattern="*mitogenome.fa", full.names=T)
mt = sapply(mt, read.dna, format="fasta", as.character=T)

get_n <- function(x) {
	return(sum(x == 'n'))
}

n_num = lapply(mt, get_n)
names(n_num) = gsub(".*//", "", names(n_num))
names(n_num) = gsub("_.*", "", names(n_num))

d = read.csv(paste(indir, "coverage/summary_statistics.csv", sep=""), stringsAsFactors=F)
rownames(d) = d$sample
d[d$sample == 'CHUNB45362', "lineage"] = 'Corallus_hortulanus'
map = d[names(n_num), "map_reads"]
cor.test(map, unlist(n_num))

tree = read.tree("/Volumes/sosi/brazil/mito_aln/RAxML_bestTree.mito")
tree = drop.tip(tree, "chicken")

tree$tip.label = d[tree$tip.label, "lineage"]
pdf("~/Desktop/mito_gene_tree.pdf", width=5, height=6, pointsize=10)
par(mar=c(1,1,1,1))
plot(ladderize(tree), cex=0.7)
dev.off()

######################
# linear model for individual success
######################

a = read.csv("/Volumes/sosi/brazil/metadata/linear_model_for_inds.csv", stringsAsFactors=F, na.string=c(NA, ""))
a$date = as.Date(a$date, "%m/%d/%y")

iv = c("date", "undiluted_amount", "pool_number", "library_ng", "X260.280", "family", "snake_lizard", "orig_reads", "median_insert_size")
dv = "all_num"

c = a[complete.cases(a[,iv]),]
models = vector("list", length(iv))
for (i in 1:length(iv)) {
	# x = lm(as.formula(paste(dv, " ~ ", paste(iv, collapse= "+"))), data=c)
	models[[i]] = lm(as.formula(paste(dv, " ~ ", iv[i])), data=a)
}
r2 = lapply(models, function(x) {return(summary(x)$adj.r.squared)})
names(r2) = iv

fac2 = combn(iv, 2, simplify=F)
fmla = lapply(fac2, function(x) {return(as.formula(paste(dv, " ~ ", paste(x, collapse= "+"))))})
models2 = vector("list", length(fmla))
for (i in 1:length(models2)) {
	# x = lm(as.formula(paste(dv, " ~ ", paste(iv, collapse= "+"))), data=c)
	models2[[i]] = lm(fmla[[i]], data=a)
}
r2 = lapply(models2, function(x) {return(summary(x)$adj.r.squared)})
r2 = unlist(r2)
names(r2) = fac2

pdf("~/Desktop/model_for_individualsuccess.pdf", width=6.5, height=3, pointsize=10)
ggplot(a, aes(pool_number, all_num)) + geom_boxplot() + theme_bw() + xlab("pool number") + ylab("# loci assembled")
dev.off()

######################
# linear model for loci success
######################

x = read.csv("~/macroevolution/brazil/probe_design/SqCL_v1/develop/locus_divergence.csv", stringsAsFactors=F)
x[x$type == 'genes', "locus"] = paste("gene-", x[x$type == 'genes', "locus"], sep="")
x[x$type == 'AHE', "locus"] = paste("AHE-", x[x$type == 'AHE', "locus"], sep="")
x[x$type == 'UCE', "locus"] = gsub('_\\S+$', "", x[x$type == 'UCE', "locus"])

a$probe_div = x[match(a$locus, x$locus), "lizard_mean"]

x = read.csv("~/Desktop/prg/repeats.csv", stringsAsFactors=F)
a$repeats = rep(0, nrow(a))
a[a$locus %in% x$locus, "repeats"] = 1
a$repeat_len = rep(0, nrow(a))
a[match(x$locus, a$locus), "repeat_len"] = x$repeats

write.csv(a, "~/Desktop/model_for_loci_missingness.csv", row.names=F)

a = read.csv("/Volumes/sosi/brazil/metadata/model_for_loci_missingness.csv", stringsAsFactors=F)
a[a$locus == 'gene-NTF-3', "type"] = "gene" 

iv = c("num_probes", "probe_GC", "probe_density", "probe_length", "seq_GC", "seq_div", "type", "probe_div", "repeats", "repeat_len")
dv = "complete"

c = a[complete.cases(a[,iv]),]
models = vector("list", length(iv))
for (i in 1:length(iv)) {
	# x = lm(as.formula(paste(dv, " ~ ", paste(iv, collapse= "+"))), data=c)
	models[[i]] = lm(as.formula(paste(dv, " ~ ", iv[i])), data=a)
}
r2 = lapply(models, function(x) {return(summary(x)$adj.r.squared)})
names(r2) = iv

fac2 = combn(iv, 2, simplify=F)
fmla = lapply(fac2, function(x) {return(as.formula(paste(dv, " ~ ", paste(x, collapse= "+"))))})
models2 = vector("list", length(fmla))
for (i in 1:length(models2)) {
	# x = lm(as.formula(paste(dv, " ~ ", paste(iv, collapse= "+"))), data=c)
	models2[[i]] = lm(fmla[[i]], data=a)
}
r2 = lapply(models2, function(x) {return(summary(x)$adj.r.squared)})
r2 = unlist(r2)
names(r2) = fac2

pdf("~/Desktop/model_for_locisuccess.pdf", width=4.5, height=3, pointsize=10)
ggplot(a, aes(probe_div, complete)) + geom_point(alpha=0.3, size=0.7) + theme_bw() + xlab("probe divergence") + ylab("# inds assembled") + geom_smooth(method = "lm", se = FALSE, col="red")
dev.off()

pdf("~/Desktop/sequence_divergence_locus.pdf", width=6, height=3, pointsize=10)
x = ggplot(a, aes(type, seq_div)) + geom_boxplot() + theme_bw() + xlab("locus type") + ylab("locus divergence") + ggtitle('B.') + theme(plot.title=element_text(hjust=0))
y = ggplot(a, aes(type, probe_div)) + geom_boxplot() + theme_bw() + xlab("locus type") + ylab("probe divergence") + ggtitle('A.') + theme(plot.title=element_text(hjust=0))
grid.arrange(y, x, ncol=2)
dev.off()

#######################
# complexity
#######################

files = list.files("/Volumes/sosi/brazil/metadata/complexity/", full.names=T)
vals = pretty(range(0,9.5e6),3)

names = gsub(".realigned\\S+", "", files)
names = gsub('^.*\\/', '', names)
p2 = rep(NA, length(files))
names(p2) = names

pdf("~/Desktop/complexity.pdf", width=3, height=3, pointsize=10)
par(mar=c(4,5,1,1), tck=-0.02)
plot(NULL, xlim=range(vals), ylim=range(vals), axes=F, xlab="", ylab="")
for (i in 1:length(files)) {
	d = read.table(files[i], header=T, sep="\t")
	p2[i] = d[d$TOTAL_READS == 2e6, 'EXPECTED_DISTINCT']
	lines(d[,1], d[,2], col=alpha("black", 0.5), lwd=2)
}
abline(a=0,b=1,col="red", lty="dotted")
axis(1, at=vals, labels=vals, line=0)
mtext("sequenced reads", side=1, line=2.3)
axis(2, at=vals, labels=vals, line=0, las=2)
mtext("distinct reads", side=2, line=3.5)
dev.off()

#######################
# subsetting analysis
#######################

d = read.csv('/Volumes/sosi/brazil_sub/coverage/summary_statistics.csv', stringsAsFactors=F)
d$reads = as.numeric(gsub('.*_', '', d$lineage)) * 1e6
d$lineage1 = gsub('_\\d+$', '', d$lineage)

x = ggplot(d, aes(as.factor(reads), all_num)) + geom_boxplot(aes(fill = factor(lineage1))) + theme_bw() + xlab("number of reads") + ylab("# of loci assembled") + theme(legend.position="none") + scale_fill_manual(values=c("gray", "black"))
y = ggplot(d, aes(as.factor(reads), all_sites_10x)) + geom_boxplot(aes(fill = factor(lineage1))) + theme_bw() + xlab("number of reads") + ylab("# of sites with high cov.") + theme(legend.position="none") + scale_fill_manual(values=c("gray", "white"))

pdf(paste(moutdir, "7subset_experiment.pdf", sep=""), width=6.5, height=2.5, pointsize=8)
grid.arrange(x, y, ncol=2)
dev.off()

##################
# variant data
##################

d = read.csv("/Volumes/sosi/brazil/metadata/Colobosaura_modesta_variants.csv", stringsAsFactors=F)
x = read.csv("/Volumes/sosi/brazil/metadata/Colobosaura_modesta_locus_lengths_variants.csv", stringsAsFactors=F)

x$type = gsub('-\\S+', '', x$locus)
x$length = abs(x$min) + abs(x$max)
var = table(d$locus)
rownames(x) = x$locus
x$var = var[rownames(x)]
x$var_n = x$var / x$length
aggregate(x$var_n, by=list(x$type), mean, na.rm=T)

##################
# missingness
##################

d = read.csv("~/Desktop/amount_missing.csv", stringsAsFactors=F)

x = d[grep("uce", d$locus),]
y = d[grep("AHE", d$locus),]
z = d[grep("gene", d$locus),]

x = aggregate(x$per_miss, by=list(x$rel_pos), mean, na.rm=T)
y = aggregate(y$per_miss, by=list(y$rel_pos), mean, na.rm=T)
z = aggregate(z$per_miss, by=list(z$rel_pos), mean, na.rm=T)

x1 = x
y1 = y
z1 = z

x1[,2] = rollmean(x[,2], k=10, na.pad=T)
y1[,2] = rollmean(y[,2], k=10, na.pad=T)
z1[,2] = rollmean(z[,2], k=10, na.pad=T)

xvals = pretty(c(-450,450), 6)[2:6]
yvals = pretty(c(0,0.3), 3)

pdf("~/Desktop/alignment_missingness.pdf", width=3, height=3, pointsize=8)
par(mar=c(4,5,1,1), tck=-0.02)
plot(x1[,1], x1[,2], col=cols[3], axes=F, xlab="", ylab="", xlim=c(-450,450), ylim=range(yvals), lwd=2, type="l")
lines(y1[,1], y1[,2], col=cols[1], lwd=2)
lines(z1[,1], z1[,2], col=cols[2], lwd=2)

axis(1, at=xvals, labels=xvals, line=0.2)
mtext("nucleotide position", side=1, line=2.3)
axis(2, at=yvals, labels=yvals, line=0.2, las=2)
mtext("percent of alignment missing", side=2, line=3)
dev.off()

#################
# Table S1
###############

a = read.csv("/Volumes/sosi/brazil/metadata/linear_model_for_inds.csv", stringsAsFactors=F, na.string=c(NA, ""))
d = read.csv(paste(indir, "coverage/summary_statistics.csv", sep=""), stringsAsFactors=F)

drop = c("lineage", "orig_reads", "map_reads", "median_insert_size", "all_num", "all_mean_cov")
a = a[,setdiff(names(a), drop)]
a = merge(a, d, by="sample")

a$municipality = c[match(a$sample, c$sample), "municipality"]
rownames(sites) = sites$municipality
a = cbind(a, sites[a$municipality, c("latitude", "longitude")])

keep = c("sample", "family", "lineage", "municipality", "longitude", "latitude", "date", "X260.280", "pool_number", "orig_reads", "map_reads", "duplicates", "all_num", "all_mean_length", "all_mean_cov", "all_sites_10x")
a = a[, keep]
a = a[order(a$sample),]
a$duplicates = a$duplicates / a$map_reads
newnames = c("sample", "family", "lineage", "municipality", "longitude", "latitude", "date", "260/280", "capture pool", "raw reads", "per. mapped reads", "per. duplicate reads", "number of assembled loci", "mean locus length", "mean locus coverage", "number sites >10x coverage")
names(a) = newnames
write.csv(a, "/Users/sonal/Desktop/activeWork/SqCL/TableS1.csv", row.names=F)


################################
# phylogenetic informativeness
################################

drv <- dbDriver("SQLite")
con <- dbConnect(drv, dbname = "~/publications/SqCL/revision_results/phylogenetic-informativeness.sqlite")
c = dbGetQuery(con, "select locus, time, pi from loci, discrete where loci.id = discrete.id;")
 
c$type = gsub('-.*', '', c$locus)
loci = unique(c$locus)
lengths = read.csv("~/publications/SqCL/revision_results/loci_lengths.csv", stringsAsFactors=F)

c$len = lengths[match(c$locus, lengths$loci), 'length']
c$pi_nuc = c$pi / c$len

maxpi = data.frame(loci=unique(c$locus))
maxpi$type = gsub('-.*', '', maxpi$loci)
maxpi$time = rep(NA, nrow(maxpi))
for (i in 1:nrow(maxpi)) {
	locus = maxpi[i, "loci"]
	tmp = c[c$locus == locus,]
	maxval = max(tmp$pi)
	maxpi[i, "time"] = tmp[which(tmp$pi == maxval), "time"]
}

# change the time to the right time
# http://www.wienslab.com/Publications_files/Zheng_Wiens_2015b_MPE.pdf
maxpi$time = (maxpi$time / 130) * 200
c$time = (c$time / 130) * 200

pdf(file=paste(outdir, "phyinformativeness.pdf", sep=""), width=6, height=2.5, pointsize=8)
a = ggplot(c, aes(time, pi_nuc, colour=type)) + stat_smooth(span=0.1, alpha=0.3) + theme_bw()  + scale_colour_manual(values=cols) + xlab("time (mya)") + ylab("PI per nucleotide") + theme(legend.position="none")
b = ggplot(maxpi, aes(time, colour = type)) + geom_line(stat="density")+ theme_bw() + xlab("time (mya) of maximum PI") + scale_colour_manual(values=cols)
grid.arrange(a, b, ncol=2)
dev.off()

################################
#  tree certainty
################################

d = read.csv("~/publications/SqCL/revision_results/tree_certainty.csv", stringsAsFactors=F)
d$type = gsub('-.*', '', d$locus)
d$length = lengths[match(d$locus, lengths$loci), 'length']
pdf(file=paste(outdir, "treecertainty.pdf", sep=""), width=2.5, height=2.5, pointsize=8)
ggplot(d, aes(tc, colour = type)) + geom_line(stat="density")+ theme_bw() + xlab("tree certainty") + scale_colour_manual(values=cols)
dev.off()

# linear models suggest this just isn't a function of locus lengths
a = lm(d$tc ~ d$length)
a = lm(d$tc ~ d$length + d$type)
a = lm(d$tc ~ d$length * d$type)

###################################
# correlation of genetic divergence
####################################

get_div <- function(tag) {
	file = paste("~/publications/SqCL/revision_results/", tag, "_divergence.csv", sep="")
	d = read.csv(file, stringsAsFactors=F)

	inds = unique(c(d$ind1, d$ind2))
	dist = matrix(NA, nrow=length(inds), ncol=length(inds), dimnames=list(inds, inds))
	
	for (i in 1:(length(inds) - 1)) {
		for (j in (i+1):length(inds)) {
			ids = sort(c(inds[i], inds[j]))
			div = d[d$ind1 == ids[1] & d$ind2 == ids[2], "divergence"]
			dist[ids[1], ids[2]] = div
			dist[ids[2], ids[1]] = div
		}
	}
	return(dist)
}

ahe_dist = get_div("AHE")
uce_dist = get_div("uce")
gene_dist = get_div("gene")

graph <- function(x, y, xlab, ylab) {
	xtcks = pretty(x, 3)
	ytcks = pretty(y, 3)
	plot(NULL, axes=F, xlab="", ylab="", xlim=range(xtcks), ylim=range(ytcks), xaxs="i", yaxs="i")
	points(x, y, pch=16, col=alpha("gray", 0.2))
	axis(1, at=xtcks, labels=NULL)
	axis(2, at=ytcks, labels=NULL, las=2)
	mtext(xlab, side=1, line=2.75)
	mtext(ylab, side=2, line=3)
	abline(a=0, b=1, col="black", lty=3)
	cor = mantel(x, y)
	text(xtcks[2] * 0.5, max(ytcks) * 0.8, paste("r=", round(cor$statistic, 3), '\np-val=', cor$signif, sep=""), adj=c(0,0))
}

pdf(file=paste(outdir, "genetic_divergence.pdf", sep=""), width=6.5, height=2.5, pointsize=8)
par(mfrow=c(1,3), mar=c(4,4.5,1,1), tck=-0.01)
graph(ahe_dist, uce_dist, "AHE, genetic distance", "UCE, genetic distance")
graph(ahe_dist, gene_dist, "AHE, genetic distance", "gene, genetic distance")
graph(uce_dist, gene_dist, "UCE, genetic distance", "gene, genetic distance")
dev.off()

####################################
# IBD distance
####################################

indir = "/Users/sonal/Dropbox/Cerrado_Australia_NSF_Proposal/UCE_results/"

c = read.csv(paste(indir, "collection.csv", sep=""), stringsAsFactors=F)
s = read.csv(paste(indir, "samples.csv", sep=""), stringsAsFactors=F)
sites = read.csv(paste(indir, "sites.csv", sep=""), stringsAsFactors=F)
t = read.csv(paste(indir, "taxonomy.csv", sep=""), stringsAsFactors=F)

get_distance <- function(lat1, long1, lat2, long2) {
	dist = distm(c(lat1, long1), c(lat2, long2), fun=distHaversine)
	return(dist)
}

get_latlongs <- function(d) {
	d$loc1 = c[match(d$ind1, c$sample), 'municipality']
	d$loc2 = c[match(d$ind2, c$sample), 'municipality']
	d = cbind(d, sites[match(d$loc1, sites$municipality), 5:6])
	d = cbind(d, sites[match(d$loc2, sites$municipality), 5:6])
	names(d)[10:13] = c("lat1", "long1", "lat2", "long2")
	d$dist = mapply(get_distance, d$lat1, d$long1, d$lat2, d$long2)
	
	# how to transform zero distances? drop for now
	d = d[d$dist > 0,]
	d$inv_fst = sapply(d$fst, inv_fst)
	d$ln_m = log(d$dist)

	return(d)
}

inv_fst <- function(fst) {
	if (fst == 0) {
		fst == 0.0001
	}
	fst = fst / (1 - fst)
	return(fst)
}

d1 = read.csv(paste(indir, "Bothrops_moojeni_divergence.csv", sep=""), stringsAsFactors=F)
# CHUNB45362 is a misid? 
d1 = d1[d1$ind1 != 'CHUNB45362',]
d2 = read.csv(paste(indir, "Colobosaura_modesta_divergence.csv", sep=""), stringsAsFactors=F)

d1 = get_latlongs(d1)
d2 = get_latlongs(d2)

xvals = pretty(range(d1$ln_m, d2$ln_m), 3)
yvals = pretty(range(d1$inv_fst, d2$inv_fst), 3)	

bwcols = c(cols[1], cols[2])
pdf(file=paste(moutdir, "6IBD.pdf", sep=""), width=2.5, height=2.5, pointsize=8)

par(mar=c(4,5,1,1), tck=-0.01)

plot(NULL, xlim=range(xvals), ylim=range(yvals), axes=F, xlab="", ylab="", xaxs="i", yaxs="i")
points(d1$ln_m, d1$inv_fst, pch=21, bg=bwcols[1])
abline(lm(d1$inv_fst ~ d1$ln_m), col=alpha(bwcols[1], 0.8))
points(d2$ln_m, d2$inv_fst, pch=23, bg=bwcols[2])
abline(lm(d2$inv_fst ~ d2$ln_m), col=alpha(bwcols[2], 0.8))

axis(1, at=xvals, labels=xvals)
mtext("ln(distance in m)", side=1, line=2.4)
axis(2, at=yvals, labels=yvals, las=1)
mtext(expression(F[ST] / (1 - F[ST])), side=2, line=3)

xloc = (max(xvals) - min(xvals)) * 0.1 + min(xvals)
xloc2 = (max(xvals) - min(xvals)) * 0.05 + min(xvals)
text(xloc, par()$usr[4]*0.75, "Bothrops moojeni", col=bwcols[1], adj=c(0,0.5), font=3)
points(xloc2, par()$usr[4]*0.75,  pch=21, bg=bwcols[1])
text(xloc, par()$usr[4]*0.9, "Colobosaura modesta", col=bwcols[2], adj=c(0,0.5), font=3)
points(xloc2, par()$usr[4]*0.9,  pch=23, bg=bwcols[2])
dev.off()

####################################
# trees across marker types
####################################

getColDiff <- function (tr1, tr2, tipDiff = NULL, vec1 = NULL, vec2 = NULL, 
    baseCol = "grey", col1 = "peachpuff", col2 = "red2", colourMethod = "ramp", 
    palette = lightseasun, ...) 
{
    l <- length(tr1$tip.label)
    if (is.null(tipDiff)) {
        tipDiff <- tipDiff(tr1, tr2, vec1, vec2)
    }
    tipSignificance1 <- sapply(tr1$tip.label, function(x) tipDiff[which(tipDiff[, 
        1] == x), 2])
    tipSignificance2 <- sapply(tr2$tip.label, function(x) tipDiff[which(tipDiff[, 
        1] == x), 2])
    if (colourMethod == "ramp") {
        colfunc <- colorRampPalette(c(col1, col2))
        if (min(tipDiff[, 2]) == 0) {
            tipSignificance1 <- tipSignificance1 + 1
            tipSignificance2 <- tipSignificance2 + 1
            if (max(tipDiff[, 2]) == 0) {
                numCols <- 0
            }
            else {
                numCols <- max(tipDiff[, 2]) - min(tipDiff[, 
                  2][which(tipDiff[, 2] != 0)]) + 1
            }
            pal <- c(baseCol, colfunc(numCols))
        }
        else {
            numCols <- max(tipDiff[, 2]) - min(tipDiff[, 2]) + 
                1
            pal <- colfunc(numCols)
        }
        tipCols1 <- pal[as.factor(tipSignificance1)]
        tipCols2 <- pal[as.factor(tipSignificance2)]
    }
    else {
        tipCols1 <- num2col(tipSignificance1, col.pal = palette)
        tipCols2 <- num2col(tipSignificance2, col.pal = palette)
        if (min(tipDiff[, 2]) == 0) {
            tipCols1[which(tipSignificance1 == 0)] <- baseCol
            tipCols2[which(tipSignificance2 == 0)] <- baseCol
        }
    }
	return(list(tipCols1, tipCols2))
}

ref = paste("~/publications/SqCL/revision_results/trees/best_trees_tol1e-05_collapse1_uce_ASTRAL.tre", sep="")
ref = read.tree(ref)
ref = drop.tip(ref, "Gallus_gallus")
ref = ladderize(ref)
write.tree(ref, file="~/Desktop/test.tre")
ref = read.tree("~/Desktop/test.tre")

types = c("AHE", "gene", "uce")
trees = vector("list", 6)
for (i in 1:length(types)) {
	tree1 = paste("~/publications/SqCL/revision_results/trees/RAxML_bipartitions.concatenated_", types[i], "_partitioned", sep="")
	t1 = read.tree(tree1)
	t1 = root(t1, "Gallus_gallus")
	t1 = drop.tip(t1, "Gallus_gallus")
	t1 = ladderize(t1)
	write.tree(t1, file="~/Desktop/test.tre")
	t1 = read.tree("~/Desktop/test.tre")
	t1 = tipRotate(t1, setNames(1:length(ref$tip.label), rev(ref$tip.label)))
	# plot(t1)
	tree2 = paste("~/publications/SqCL/revision_results/trees/best_trees_tol1e-05_collapse1_", types[i], "_ASTRAL.tre", sep="")
	t2 = read.tree(tree2)
	t2 = root(t2, "Gallus_gallus")
	t2 = drop.tip(t2, "Gallus_gallus")
	write.tree(t2, file="~/Desktop/test.tre")
	t2 = read.tree("~/Desktop/test.tre")
	t2 = tipRotate(t2, setNames(1:length(ref$tip.label),  rev(ref$tip.label)))
	# plot(t2)
	cat(RF.dist(t1, t2), types[i], "\n")
	trees[[i]] = t1
	names(trees)[i] = paste("concat_", types[i], sep="")
	trees[[i + 3]] = t2
	names(trees)[i + 3] = paste("ASTRAL_", types[i], sep="")
}

rfmat = matrix(NA, nrow=length(types), ncol=length(types))
for (i in 1:length(types)) {
	for (j in 1:length(types)) {
		rfmat[i, j] = round(RF.dist(trees[[paste("ASTRAL_", types[i], sep="")]], trees[[paste("ASTRAL_", types[j], sep="")]], normalize=T), 3)
	}
}
dimnames(rfmat) = list(types, types)

unc1 = c("Copeoglossum_nigropunctatum", "Brasiliscincus_heathi", "Notomabuya_frenata")
unc2 = c("Erythrolamprus_reginae", "Erythrolamprus_poecilogyrus", "Erythrolamprus_almadensis", "Xenodon_merremi", "Lygophis_paucidens", "Psomophis_joberti", "Thamnodynastes_hypoconia", "Philodryas_nattereri", "Philodryas_olfersii", "Apostolepis_polylepis", "Apostolepis_cearensis", "Taeniophallus_occipitalis", "Pseudoboa_neuwiedii", "Pseudoboa_nigra", "Oxyrhopus_petolarius", "Phimophis_guerini", "Oxyrhopus_trigeminus")
	
pdf(file=paste(moutdir, "4speciestrees_across_markers.pdf", sep=""), width=6.5, height=3.5, pointsize=8)
par(mfrow=c(1,3), mar=c(0,0,2,0))
for (i in 1:(length(types))) {
	tree= trees[[paste("ASTRAL_", types[i], sep="")]]
	
	plot(tree, cex=0.7)
	
	pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
	
	yvert1 = range(which(tree$tip.label %in% unc1))
	xvert1 = c(30, pp$x.lim[2], pp$x.lim[2], 30)
	yvert1 = c(yvert1[1] - 0.5, yvert1[1] - 0.5, yvert1[2] + 0.5, yvert1[2] + 0.5)
	polygon(xvert1, yvert1, col=alpha("gray", 0.5), border=NA)
	
	yvert1 = range(which(tree$tip.label %in% unc2))
	xvert1 = c(25, pp$x.lim[2], pp$x.lim[2], 25)
	yvert1 = c(yvert1[1] - 0.5, yvert1[1] - 0.5, yvert1[2] + 0.5, yvert1[2] + 0.5)
	polygon(xvert1, yvert1, col=alpha("gray", 0.5), border=NA)
	
	mtext(types[i], side=3, line=0, font=2)
	tree$node.label = as.numeric(tree $node.label)
		for (x in 1:length(tree$node.label)) {
			if (!is.na(tree$node.label[x])) {
				if (tree$node.label[x] < 0.95) {
					nodelabels(node=x + length(tree$tip.label), bg="red", pch=21, cex=1.2, lwd=0.5)
					}
				}
			}
	}
addtable2plot(par()$usr[1] + 1, par()$usr[4] * 0.05, rfmat, bty="o", display.rownames=TRUE, hlines=TRUE, vlines=TRUE)
dev.off()


rfmat = matrix(NA, nrow=length(types), ncol=length(types))
for (i in 1:length(types)) {
	for (j in 1:length(types)) {
		rfmat[i, j] = round(RF.dist(trees[[paste("concat_", types[i], sep="")]], trees[[paste("concat_", types[j], sep="")]], normalize=T), 3)
	}
}
dimnames(rfmat) = list(types, types)

pdf(file=paste(outdir, "concatenatedtrees_across_markers.pdf", sep=""), width=6.5, height=3.5, pointsize=8)
par(mfrow=c(1,3), mar=c(0,0,2,0))
for (i in 1:(length(types))) {
	tree= trees[[paste("concat_", types[i], sep="")]]
	
	plot(tree, cex=0.7)
	
	pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
	
	yvert1 = range(which(tree$tip.label %in% unc1))
	node = getMRCA(tree, unc1)
	xv = nodeheight(tree, node)
	xvert1 = c(xv, pp$x.lim[2], pp$x.lim[2], xv)
	yvert1 = c(yvert1[1] - 0.5, yvert1[1] - 0.5, yvert1[2] + 0.5, yvert1[2] + 0.5)
	polygon(xvert1, yvert1, col=alpha("gray", 0.5), border=NA)
	
	yvert1 = range(which(tree$tip.label %in% unc2))
	node = getMRCA(tree, unc2)
	xv = nodeheight(tree, node)
	xvert1 = c(xv, pp$x.lim[2], pp$x.lim[2], xv)
	yvert1 = c(yvert1[1] - 0.5, yvert1[1] - 0.5, yvert1[2] + 0.5, yvert1[2] + 0.5)
	polygon(xvert1, yvert1, col=alpha("gray", 0.5), border=NA)
	
	mtext(types[i], side=3, line=0, font=2)
	tree$node.label = as.numeric(tree $node.label)
		for (x in 1:length(tree$node.label)) {
			if (!is.na(tree$node.label[x])) {
				if (tree$node.label[x] < 95) {
					nodelabels(node=x + length(tree$tip.label), bg="red", pch=21, cex=1.2, lwd=0.5)
					}
				}
			}
	}
addtable2plot(par()$usr[2] * 0.1, par()$usr[4] * 0.05, rfmat, bty="o", display.rownames=TRUE, hlines=TRUE, vlines=TRUE)
dev.off()


pdf(file=paste(outdir, "speciestrees_across_markers2.pdf", sep=""), width=3.5, height=6.5, pointsize=8)
par(mfrow=c(3,2), mar=c(0,0,0,0))
for (i in 1:(length(types) - 1)) {
	for (j in (i+1):length(types)) {
		tree1 = trees[[paste("ASTRAL_", types[i], sep="")]]
		tree2 = trees[[paste("ASTRAL_", types[j], sep="")]]
		res = getColDiff(tree1, tree2)
		plot(tree1, tip.color=res[[1]], cex=0.7)
		tree1$node.label = as.numeric(tree1$node.label)
		for (x in 1:length(tree1$node.label)) {
			if (!is.na(tree1$node.label[x])) {
				if (tree1$node.label[x] < 0.95) {
					nodelabels(node=x + length(tree1$tip.label), bg="red", pch=21, cex=0.6, lwd=0.5)
					}
				}
			}
		text(par()$usr[2] * 0.1, par()$usr[4] * 0.1, types[i])
		plot(tree2, tip.color=res[[2]], cex=0.7)
		tree2$node.label = as.numeric(tree2$node.label)
		for (x in 1:length(tree2$node.label)) {
			if (!is.na(tree2$node.label[x])) {
				if (tree2$node.label[x] < 0.95) {
					nodelabels(node=x + length(tree2$tip.label), bg="red", pch=21, cex=0.6, lwd=0.5)
					}
				}
			}
		text(par()$usr[2] * 0.1, par()$usr[4] * 0.1, types[j])
	}
}
dev.off()

pdf(file=paste(outdir, "concatenatedtrees_across_markers2.pdf", sep=""), width=3.5, height=6.5, pointsize=8)
par(mfrow=c(3,2), mar=c(0,0,0,0))
for (i in 1:(length(types) - 1)) {
	for (j in (i+1):length(types)) {
		tree1 = trees[[paste("concat_", types[i], sep="")]]
		tree2 = trees[[paste("concat_", types[j], sep="")]]
		res = getColDiff(tree1, tree2)
		plot(tree1, tip.color=res[[1]], cex=0.7)
		tree1$node.label = as.numeric(tree1$node.label)
		for (x in 1:length(tree1$node.label)) {
			if (!is.na(tree1$node.label[x])) {
				if (tree1$node.label[x] < 95) {
					nodelabels(node=x + length(tree1$tip.label), bg="red", pch=21, cex=0.6, lwd=0.5)
					}
				}
			}
		text(par()$usr[2] * 0.1, par()$usr[4] * 0.1, types[i])
		plot(tree2, tip.color=res[[2]], cex=0.7)
		tree2$node.label = as.numeric(tree2$node.label)
		for (x in 1:length(tree2$node.label)) {
			if (!is.na(tree2$node.label[x])) {
				if (tree2$node.label[x] < 95) {
					nodelabels(node=x + length(tree2$tip.label), bg="red", pch=21, cex=0.6, lwd=0.5)
					}
				}
			}
		text(par()$usr[2] * 0.1, par()$usr[4] * 0.1, types[j])
	}
}
dev.off()

#############################################
# genetic differentiation across marker types
#############################################

get_diff <- function(d, tag) {
	tmp = d[d$type == tag,]

	inds = unique(c(tmp$ind1, tmp$ind2))
	dist = matrix(NA, nrow=length(inds), ncol=length(inds), dimnames=list(inds, inds))
	
	for (i in 1:(length(inds) - 1)) {
		for (j in (i+1):length(inds)) {
			ids = c(inds[i], inds[j])
			x = tmp[tmp$ind1 == ids[1] & tmp$ind2 == ids[2],]
			if (nrow(x) < 1) {
				x = tmp[tmp$ind1 == ids[2] & tmp$ind2 == ids[1],]
			}
			dist[ids[1], ids[2]] = x$fst
			dist[ids[2], ids[1]] = x$fst
		}
	}
	return(dist)
}

d = read.csv("~/publications/SqCL/revision_results/Colobosaura_modesta_divergence.csv", stringsAsFactors=F)
ahe_diff1 = get_diff(d, "AHE")
uce_diff1 = get_diff(d, "uce")
gene_diff1 = get_diff(d, "gene")

d = read.csv("~/publications/SqCL/revision_results/Bothrops_moojeni_divergence.csv", stringsAsFactors=F)
ahe_diff2 = get_diff(d, "AHE")
uce_diff2 = get_diff(d, "uce")
gene_diff2 = get_diff(d, "gene")

graph <- function(x, y, xlab, ylab) {
	xtcks = pretty(x, 3)
	ytcks = pretty(y, 3)
	plot(NULL, axes=F, xlab="", ylab="", xlim=range(xtcks), ylim=range(ytcks), xaxs="i", yaxs="i")
	points(x, y, pch=16, col="gray")
	axis(1, at=xtcks, labels=NULL)
	axis(2, at=ytcks, labels=NULL, las=2)
	mtext(xlab, side=1, line=2.75)
	mtext(ylab, side=2, line=3)
	abline(a=0, b=1, col="red", lty=3)
	
	
	x1 = data.frame(x)
	x1 = data.frame( t(combn(names(x1),2)), dist=t(x1)[lower.tri(x1)] )
	y1 = data.frame(y)
	y1 = data.frame( t(combn(names(y1),2)), dist=t(y1)[lower.tri(y1)] )
	z = merge(x1, y1, by=c("X1", "X2"))
	fit = lm(z$dist.y ~ z$dist.x)
	abline(lm(z$dist.y ~ z$dist.x))
	# cat(coef(fit)[2], "\n")
	
	cor = mantel(x, y)
	text((par()$usr[2] - par()$usr[1]) * 0.1 + par()$usr[1], (par()$usr[4] - par()$usr[3]) * 0.8 + par()$usr[3], paste("r=", round(cor$statistic, 3), '\np-val=', round(cor$signif, 3), "\n", 'm=', round(coef(fit)[2], 2), sep=""), adj=c(0,0))
}

pdf(file=paste(outdir, "genetic_diff_Colo.pdf", sep=""), width=6.5, height=2.5, pointsize=8)
par(mfrow=c(1,3), mar=c(4,4.5,1,1), tck=-0.01)
graph(ahe_diff1, uce_diff1, "AHE, Fst", "UCE, Fst")
graph(ahe_diff1, gene_diff1, "AHE, Fst", "gene, Fst")
graph(uce_diff1, gene_diff1, "UCE, Fst", "gene, Fst")
dev.off()

#############################################
# loci filtering for clock-likeness
#############################################

# the Jeremy Brown approach is very sensitive to number of tips
# and locus length

d = read.csv("~/publications/SqCL/revision_results/ultrametric_filtering.csv", stringsAsFactors=F)
d = d[complete.cases(d),]
d$diff = -2 * (d$lnl_ultra  - d$lnl_nonultra)
d$type = gsub("-.*$", "", d$locus)

d = d[d$ntips > 41,]
a = lm(d$diff ~ d$loc_length + d$ntips)
b = lm(d$diff ~ d$loc_length + d$ntips + d$type)
c = lm(d$diff ~ d$loc_length * d$type + d$ntips)


pdf(file=paste(outdir, "clocklikeness.pdf", sep=""), width=3.5, height=2.5, pointsize=8)
ggplot(dsub, aes(loc_length, diff, colour = type)) + geom_point(alpha=0.05) + geom_smooth(se = FALSE, method = "lm") + scale_colour_manual(values=cols) + theme_bw() + xlab("locus length") + ylab("difference in LnL")
dev.off()

# https://github.com/marekborowiec/good_genes

Clocklikeness <- function(file) {
  # read in tree
  tree <- read.tree(file)

  # root tree
  rooted_tr <- midpoint(tree)
  # get matrix diagonal of phylogenetic variance-covariance matrix
  # these are your distances from root
  root_dist <- diag(vcv.phylo(rooted_tr))
  std_dev_root_dist <- sd(root_dist)
  mean_root_dist <- mean(root_dist)
  CV = (std_dev_root_dist/mean_root_dist)*100
 
  return(CV)
}

trees = list.files("~/publications/SqCL/data/trees/", full.names=T)
clocks = sapply(trees, Clocklikeness)
names(clocks) = gsub("/Users/Sonal/publications/SqCL/data/trees//", "", names(clocks))
names(clocks) = gsub(".bestTree.tre", "", names(clocks))
x = data.frame(clock=clocks, locus=names(clocks))
x$type = gsub("-.*$", "", x$locus)

# pdf(file=paste(outdir, "clocklikeness.pdf", sep=""), width=3.5, height=2.5, pointsize=8)
ggplot(x, aes(clock, colour = type)) + geom_line(stat="density")+ theme_bw() + xlab("clocklikeness measure") + scale_colour_manual(values=cols)
# dev.off()

#############################################
# tree clustering
#############################################

# d = read.csv("~/Desktop/dist_matrix.csv", as.is=TRUE, na.string=c("", 'NA'), row.names=1)
# d = as.matrix(d)
# saveRDS(d, "~/Desktop/dist.rds")
td = readRDS("~/publications/SqCL/revision_results/dist.rds")

d = d[d$ntips > 35,]
keep = c(d$locus[grep("AHE", d$locus)], d$locus[grep("gene", d$locus)], sample(d$locus[grep("uce", d$locus)], 300))
keep = sort(keep)
keep = gsub('-', ".", keep)

dimnames(td)[[1]] = gsub("-", ".", dimnames(td)[[1]])
dimnames(td)[[2]] = gsub("-", ".", dimnames(td)[[2]])

tdm1 = td[keep, keep]
tdm2 = as.dist(tdm1, upper=T)
tdm3 = cailliez(tdm2)
pc = dudi.pco(tdm3, nf=4)

pdf(file=paste(outdir, "tree_clustering.pdf", sep=""), width=4.5, height=4.5, pointsize=8)
par(mar=c(4,4,1,1), tck=-0.01)
xtcks = pretty(pc$li[,1], 4)
ytcks = pretty(pc$li[,2], 4)
plot(NULL, axes=F, xlab="", ylab="", xlim=range(xtcks), ylim=range(ytcks), xaxs="i", yaxs="i")
points(pc$li[,1], pc$li[,2], col=alpha(cols[as.factor(gsub("\\..*", "", keep))], 0.7), pch=16)
legend(-1.85, -1, legend=c("AHE", "gene", "uce"), col=alpha(cols, 0.7), pch=16, bty = "n")
axis(1, at=xtcks, labels=NULL)
axis(2, at=ytcks, labels=NULL, las=2)
mtext("PC axis 1", side=1, line=2.5)
mtext("PC axis 2", side=2, line=2.75)
dev.off()

#############################################
# genetic diversity
#############################################

d = read.csv("~/publications/SqCL/revision_results/genetic_diversity_loci.csv", stringsAsFactors=F)
d = d[d$ind != 'CHUNB52120',]

graph_div <- function(x, y, xlab, ylab) {
	xtcks = pretty(x, 3)
	ytcks = pretty(y, 3)
	plot(NULL, axes=F, xlab="", ylab="", xlim=range(xtcks), ylim=range(ytcks), xaxs="i", yaxs="i")	
	points(x, y, pch=21, bg='gray')
	axis(1, at=xtcks, labels=NULL)
	axis(2, at=ytcks, labels=NULL, las=2)
	mtext(xlab, side=1, line=2.5)
	mtext(ylab, side=2, line=3.75)
	# lines(x=c(d[i, "AHE_div"], d[i, "AHE_div"]), y=c(d[i, "uce_div_sub_AHE_95l"], d[i, "uce_div_sub_AHE_95u"]), col="gray")
	fit = lm(y ~ x)
	abline(lm(y ~ x))
	abline(a=0, b=1, col="red", lty=2)
	corval = cor.test(x, y)
	text(par()$usr[2] * 0.2, par()$usr[4] * 0.9, paste("r=", round(corval$estimate, 3), "\n", "p-val=", sprintf("%.1e", corval$p.value), "\n", "m=", round(coefficients(fit)[2], 2), sep=""))
	c
}

pdf(file=paste(outdir, "genetic_diversity.pdf", sep=""), width=6.5, height=2.5, pointsize=8)
par(mar=c(4,7,1,1), tck=-0.01, mfrow=c(1,3))
graph_div(d$uce_div, d$AHE_div, expression(pi ~ ", UCE"), expression(pi ~ ", AHE"))
graph_div(d$uce_div, d$gene_div, expression(pi ~ ", UCE"), expression(pi ~ ", gene"))
graph_div(d$AHE_div, d$gene_div, expression(pi ~ ", AHE"), expression(pi ~ ", gene"))
dev.off()

#############################################
# div dating
#############################################

trees = list.files("~/publications/SqCL/revision_results/div_dating/", pattern="*tre$", full.names=T)
trees = c(trees[1], trees[6], trees[7])
tnames = gsub('/Users/Sonal/publications/SqCL/revision_results/div_dating//', '', trees)
tnames = gsub('_sub.tre', '', tnames)

plot_tree <- function(tree, n) {
	x = read.beast(tree)
	y = read.nexus(tree) 
	
	heights = x@stats$height_0.95_HPD
	heights = lapply(heights, as.numeric)
	xlim = max(nodeHeights(y))

	plot(y)
	lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)

	for (i in 1:length(heights)) {
		node = as.numeric(names(heights)[i])
		range = heights[[i]]
	
		#cat(xlim - range, "\n")
		range = sort(xlim - range)
		if ((range[2] - range[1]) > 0.0001) {
			lines(range, c(lastPP$yy[node], lastPP$yy[node]), col=alpha(cols[n], 0.5), lwd=4)
			}
		}
	}	

# pdf(file=paste(outdir, "divergence_trees.pdf", sep=""), width=6.5, height=4.5, pointsize=8)
par(mfrow=c(1,3), mar=c(1,1,1,1))
plot_tree(trees[1], 1)
plot_tree(trees[2], 2)
plot_tree(trees[3], 3)
# dev.off()

rescaleTree<-function(tree,scale){
tree$edge.length<-
   tree$edge.length/max(nodeHeights(tree)[,2])*scale
return(tree)
}

graph_div <- function(x, y, xlab, ylab) {
	x = read.nexus(x)
	y = read.nexus(y)
	
	rescale = max(nodeHeights(x), nodeHeights(y))
	x = rescaleTree(x, rescale)
	y = rescaleTree(y, rescale)
	
	x = cophenetic.phylo(x) 
	y = cophenetic.phylo(y)
	
	xtcks = pretty(x, 3)
	ytcks = pretty(y, 3)
	plot(NULL, axes=F, xlab="", ylab="", xlim=range(xtcks), ylim=range(ytcks), xaxs="i", yaxs="i")	
	
	abline(a=0, b=1, col="red", lty=2)
	points(x, y, pch=16, col='gray')
	axis(1, at=xtcks, labels=NULL)
	axis(2, at=ytcks, labels=NULL, las=2)
	mtext(xlab, side=1, line=2.5)
	mtext(ylab, side=2, line=3.25)

	corval = mantel(x, y)
	text(par()$usr[2] * 0.2, par()$usr[4] * 0.9, paste("r=", round(corval$statistic, 3), "\n", "p-val=", sprintf("%.1e", corval$signif), sep=""))
}

pdf(file=paste(outdir, "tree_divergences.pdf", sep=""), width=6.5, height=2, pointsize=8)
par(mar=c(4,6,1,1), tck=-0.01, mfrow=c(1,3))
graph_div(trees[1], trees[2], "phenetic distances, AHE", "phenetic distances, gene")
graph_div(trees[1], trees[3], "phenetic distances, AHE", "phenetic distances, UCE")
graph_div(trees[2], trees[3], "phenetic distances, gene", "phenetic distances, UCE")
dev.off()

#############################################
# segregating sites
#############################################

files = list.files("~/publications/SqCL/revision_results/segregating_sites/", pattern="*csv", full.names=T)

ss = list(c(), c(), c())
names(ss) = c("AHE", 'gene', 'uce')

ss_nuc = list(c(), c(), c())
names(ss_nuc) = c("AHE", 'gene', 'uce')

for (file in files) {
	d = read.csv(file, stringsAsFactors=F)
	d$type = gsub("-.*", "", d$locus)
	tmpss = aggregate(ss ~ type, data = d, FUN = function(x) c(mean = mean(x), sd = quantile(x, c(0.1, 0.9)) ) )
	tmpnuc = aggregate(ss_nuc ~ type, data = d, FUN = function(x) c(mean = mean(x), sd = quantile(x, c(0.1, 0.9)) ) )
	
	for (i in 1:nrow(tmpss)) {
		ss[[tmpss[i, "type"]]] = c(ss[[tmpss[i, "type"]]], tmpss[i, "x"])
	}
	
	for (i in 1:nrow(tmpnuc)) {
		ss_nuc[[tmpnuc[i, "Group.1"]]] = c(ss_nuc[[tmpnuc[i, "Group.1"]]], tmpnuc[i, "x"])
	}
}

lapply(ss, mean)
lapply(ss_nuc, mean)

#file = files[grep("Colobosaura_modesta", files)]
file = files[grep("Bothrops_moojeni", files)]
d = read.csv(file, stringsAsFactors=F)
d$type = gsub("-.*", "", d$locus)
tmpss = aggregate(ss ~ type, data = d, FUN = function(x) c(mean = mean(x), sd = quantile(x, c(0.1, 0.9)) ) )
tmpnuc = aggregate(ss_nuc ~ type, data = d, FUN = function(x) c(mean = mean(x), sd = quantile(x, c(0.1, 0.9)) ) )
	
# Colobosaura
#  type  ss.mean ss.sd.10% ss.sd.90%
#1  AHE 3.974696  0.993418  6.953925
#2 gene 1.778341  0.000000  3.774988
#3  uce 3.566598  0.662279  6.953925

#  type ss_nuc.mean ss_nuc.sd.10% ss_nuc.sd.90%
#1  AHE 0.003498118   0.001166800   0.006316400
#2 gene 0.002521407   0.000000000   0.004828000
#3  uce 0.004156315   0.001016000   0.008141200

# Bothrops
#  type   ss.mean ss.sd.10% ss.sd.90%
# 1  AHE 3.3914996 0.6716234 6.7162290
# 2 gene 2.0227242 0.0000000 5.1255435
# 3  uce 1.8402674 0.0000000 4.2418290

#  type ss_nuc.mean ss_nuc.sd.10% ss_nuc.sd.90%
#1  AHE 0.001944875   0.000331000   0.004000400
#2 gene 0.002249472   0.000000000   0.007115000
#3  uce 0.002093903   0.000000000   0.004956600

#############################################
# linkage disequilbrium
#############################################

# we don't have enough data
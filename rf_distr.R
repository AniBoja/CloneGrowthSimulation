library(ape)
library(phangorn)
args <- commandArgs(trailingOnly = TRUE)
gt_tree<-read.tree(args[1])
is.rooted(gt_tree)
gt_tree<-unroot(gt_tree)
recon_tree<-read.tree(args[2])
is.rooted(recon_tree)
recon_tree<-unroot(recon_tree)
rf_dist<-RF.dist(gt_tree, recon_tree)
set.seed(123)
n_permutations<-1000
null_rf<-replicate(n_permutations, {
  permuted_tree<-rtree(n = length(gt_tree$tip.label), tip.label = gt_tree$tip.label)
  RF.dist(gt_tree, permuted_tree)
})
hist(null_rf)
p_value<-sum(rf_dist >= null_rf) / n_permutations
sprintf("Robinson-Foulds distance of %i has a p-value of %i", rf_dist, p_value)

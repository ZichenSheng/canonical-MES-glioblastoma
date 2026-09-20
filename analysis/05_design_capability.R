# Aggregate requalification of frozen recovery counts; does not regenerate worlds.
a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grepl("^--file=", a)])
source(file.path(dirname(dirname(normalizePath(f))),"R/utils.R"));root<-release_root()
x<-read_tsv(file.path(root,"results_reference/design_recovery.tsv"))
stopifnot(nrow(x)==32L,all(x$size==20L),all(x$sum>=0 & x$sum<=x$size))
x$reachable<-x$minimum_attainable_exact_p<=0.05
q<-x[x$q==1,];ok<-q[q$reachable,]
stopifnot(nrow(q)==8L,nrow(ok)==7L,sum(ok$sum==20)==3L,sum(ok$sum==0)==4L,
 q$minimum_attainable_exact_p[q$signature_id=="PF4C-FINAL-053"]==0.0625)
v<-q[q$signature_id %in% c("PF4C-B4-007","PF4C-B4-029"),]
stopifnot(all(v$K_target==7L),all(v$E_s_n==10L),setequal(v$sum,c(0,20)))
cat("PASS: 8 prescribed, 7 reachable, 3 complete recovery, 4 no recovery; FINAL-053 retained but unevaluable.\n")

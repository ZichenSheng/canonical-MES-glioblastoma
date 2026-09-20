# Fail-closed release validator. --integrity-only does not certify publication.
a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grepl("^--file=",a)])
source(file.path(dirname(dirname(normalizePath(f))),"R/utils.R"))
root <- release_root()
if (!requireNamespace("digest",quietly=TRUE)) stop("Install digest 0.6.39 from renv.lock")
checks <- 0L
check <- function(ok,label) {
 if (!isTRUE(ok)) stop(paste("FAIL",label),call.=FALSE)
 checks <<- checks+1L
}
m <- read_tsv(file.path(root,"reproducibility/AUTHORITY_HASHES.tsv"))
for (i in seq_len(nrow(m))) {
 p<-file.path(root,m$file[i]);check(file.exists(p),paste("missing",m$file[i]))
 check(identical(digest::digest(file=p,algo="sha256",serialize=FALSE),m$sha256[i]),paste("SHA256",m$file[i]))
 check(length(readLines(p,warn=FALSE))==m$line_count[i],paste("line count",m$file[i]))
 if (grepl("\\.tsv$",p)) check(identical(readLines(p,n=1,warn=FALSE),m$header[i]),paste("schema",m$file[i]))
}
p <- read_tsv(file.path(root,"results_reference/printed_value_checks.tsv"))
for(i in seq_len(nrow(p))) {
 x<-read_tsv(file.path(root,p$file[i]));v<-as.numeric(x[p$row[i]-1L,p$column[i]])
 check(is.finite(v) && abs(v-p$expected[i])<=p$tolerance[i],paste("printed",p$result_id[i]))
}
c <- read_tsv(file.path(root,"reproducibility/NUMERICAL_COMPARISONS.tsv"))
check(all(is.finite(c$expected)&is.finite(c$observed)&c$tolerance>=0),"numeric receipt schema")
check(all(abs(c$observed-c$expected)<=c$tolerance),"certified numeric outputs")
check(all(c$status=="PASS"),"numeric receipt status")
d<-read_tsv(file.path(root,"results_reference/denominator_lock.tsv"))
e<-read_tsv(file.path(root,"results_reference/estimand_lock.tsv"))
check(nrow(d)==17L&&!anyDuplicated(d$flow_id),"denominator lock")
check(nrow(e)==15L&&!anyDuplicated(e$estimand_id),"estimand lock")
s<-read_tsv(file.path(root,"reproducibility/SCRIPT_CERTIFICATION.tsv"))
check(nrow(s)==8L&&all(s$status=="PASS"),"eight scoped principal certifications")
w<-read_tsv(file.path(root,"docs/RESULT_CODE_CROSSWALK.tsv"))
check(!anyDuplicated(w$result_id)&&all(paste0("R",1:6)%in%w$results_section),"six Results sections")
check(all(w$reproduction_level%in%c("FULL_FROM_SOURCE_DATA","FROM_PREPARED_PUBLIC_INPUT","AGGREGATE_VALIDATION","SOFTWARE_ONLY","NOT_PUBLICLY_REPRODUCIBLE")),"reproduction vocabulary")
check(all(w$status%in%c("CERTIFIED","CERTIFIED_WITH_LIMITATION","NOT_REQUIRED_FOR_MAIN_RELEASE","BLOCKED")),"status vocabulary")
for(script in c("01_inspect_signatures.R","05_design_capability.R","06_crc04_exact_inference.R")) {
 out<-system2(file.path(R.home("bin"),"Rscript"),shQuote(file.path(root,"analysis",script)),stdout=TRUE,stderr=TRUE)
 check(is.null(attr(out,"status"))||attr(out,"status")==0,paste("runnable",script))
 cat(paste(out,collapse="\n"),"\n")
}
python <- Sys.getenv("MANUSCRIPT_PYTHON", "python3")
pyout <- system2(python,shQuote(file.path(root,"analysis/08_aggregate_evidence.py")),stdout=TRUE,stderr=TRUE)
check(is.null(attr(pyout,"status"))||attr(pyout,"status")==0,"runnable aggregate Python module")
cat(paste(pyout,collapse="\n"),"\n")
x<-read_tsv(file.path(root,"results_reference/main_text_numbers.tsv"))
check(nrow(x)==78L && identical(x$result_id,w$result_id),"current78-claim number index")
check(identical(as.character(x$reported_value),as.character(w$reported_value)),"current manuscript reported-value lock")
check(all(file.exists(file.path(root,w$public_reproduction_module))),"all public module mappings exist")
cat(sprintf("INTEGRITY PASS: %d checks; %d printed values; %d frozen-output comparisons.\n",checks,nrow(p),nrow(c)))
b<-read_tsv(file.path(root,"docs/RELEASE_BLOCKERS.tsv"))
open<-b[b$status=="OPEN",]
cat(sprintf("PUBLICATION GATE: %d open blockers; %d BLOCKED claim rows.\n",nrow(open),sum(w$status=="BLOCKED")))
if (!"--integrity-only" %in% commandArgs(TRUE) && (nrow(open)>0L || any(w$status=="BLOCKED"))) {
 cat(paste(open$blocker_id,open$requirement,sep=": "),sep="\n")
 quit(save="no",status=2L)
}

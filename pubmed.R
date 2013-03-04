library("XML")

#ディレクトリ移動
x <- commandArgs()
name <- as.character(x[5])
setwd(paste(getwd(),"/Pubmed/",name,sep=""))
paper <- list.files()
for(i in 1:length(paper)){
print(paste(i," / ", length(paper)," at ",name,sep=""))
d <- xmlToList(paper[i])
#雑誌名
A <- name

#年代
B <- d$front$'article-meta'$`pub-date`$year

#タイトル
C <- ""
l <- length(d$front$'article-meta'$`title-group`)
for(i in 1:l){
	C <- paste(C,gsub(pattern = "\n",replacement=" ",d$front$'article-meta'$`title-group`[i],perl=TRUE),sep="")
}

#アブスト
D <- ""
l <- length(d$front$'article-meta'$abstract)
for(i in 1:l){
	D <- paste(D,d$front$'article-meta'$abstract[i],sep=" ")
}

#PMCかどうか
PMC <- d$front$'article-meta'$`article-id`$.attrs=="pmc"
print(PMC)
#PMCIDとURL
E <- "NO_PMCID"
F <- "NO_URL"
if(PMC){
E <- paste("PMC",d$front$`article-meta`$'article-id'$text,sep="")
F <- paste("http://www.ncbi.nlm.nih.gov/pmc/articles/",E,"/pdf/",sep="")
}
A <- gsub("\n"," ",A)
B <- gsub("\n"," ",B)
C <- gsub("\n"," ",C)
D <- gsub("\n"," ",D)
E <- gsub("\n"," ",E)
F <- gsub("\n"," ",F)
A <- gsub("\t"," ",A)
B <- gsub("\t"," ",B)
C <- gsub("\t"," ",C)
D <- gsub("\t"," ",D)
E <- gsub("\t"," ",E)
F <- gsub("\t"," ",F)
A <- paste(A,"\t",sep="")
B <- paste(B,"\t",sep="")
C <- paste(C,"\t",sep="")
D <- paste(D,"\t",sep="")
E <- paste(E,"\t",sep="")
F <- paste(F,"\n",sep="")
result <- c(A,B,C,D,E,F)
sink(file="../../pubmed.txt",append=T)
cat(result)
sink()
}
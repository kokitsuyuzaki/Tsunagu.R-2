# 「第2回Rでつなぐ次世代オミックス情報統合解析研究会」
での発表資料とソースコード
=======

# イントロ
```r
# Rのパッケージを根こそぎダウンロード (Rのコンソール画面で)
source("http://bioconductor.org/biocLite.R")
options("BioC_mirror" = "http://bioconductor.jp/")
biocLite(all_group(),type="source")
options(repos="http://cran.md.tsukuba.ac.jp")
packs<- available.packages(contriburl=contrib.url("http://cran.md.tsukuba.ac.jp/"))
install.packages(packs[,1],type="source")
```
SQLiteが使われているパッケージを閲覧（Macの場合、ターミナルで）

    ls /Library/Frameworks/R.framework/Resources/library/*/extdata/*.sqlite
    ls /Library/Frameworks/R.framework/Resources/library/*/extdata/*.sql
    ls /Library/Frameworks/R.framework/Resources/library/*/extdata/*.db

# 1. SQLiteの基礎
```r


# 2. RSQLiteの利用
```r
# ロード
library("RSQLite")
library("DBI")

# コネクション
driver <- dbDriver("SQLite")
db <- "pubmed.sqlite"
con <- dbConnect(driver, db)

# どんなテーブルがあるか
dbListTables(con)

# 行数を見てみる
dbGetQuery(con, "SELECT COUNT(*) FROM pubmed;")

# 一番古い年の論文は?
dbGetQuery(con, "SELECT MIN(year) FROM pubmed;")

# 年代別
y <- c()
for(i in 1928:2012){
command <- paste("SELECT COUNT(*) FROM pubmed WHERE year =", i, ";", sep="")
prey <- dbGetQuery(con, command)
y <- cbind(y, as.numeric(prey))
}

# 年代別プロット
jpeg(file="year.jpeg")
plot(1928:2012, y, "l", ylab="Frequency", xlab="Year")
dev.off()

# アブストにRNA-Seqを含むものを表示
dbGetQuery(con, "SELECT title FROM pubmed WHERE abst like '%RNA-Seq%';")

# 年度毎にRNA-Seqを含む論文が何件あるか集計
paper <- rep(0:0, length=85)
for(i in 1928:2012){
command <- paste("SELECT COUNT(*) FROM pubmed WHERE abst like '%RNA-Seq%' AND year = ", i, ";", sep="")
prepaper <- dbGetQuery(con, command)
	if(as.numeric(prepaper) != 0){
 	paper[i-1928] <- as.numeric(prepaper)
	}
}
 
# ヒストグラムで俯瞰
jpeg(file = "paper.jpeg")
plot(1928:2012, paper, "l", ylab="Frequency", xlab="Year")
dev.off()

# RNA-Seq, Microarray, RT-PCR, Western, two-hybrid, GFP 
# ワードの頻度 → ノードの大きさ
# 共起回数 → エッジの長さ

# 実験自体の頻度
freq <- c()
exp <- c("RNA-Seq", "microarray", "RT-PCR", "Western blot", "two-hybrid", "GFP")
for(i in 1:6){
	command <- paste("SELECT COUNT(*) FROM pubmed WHERE abst like '%", exp[i], "%';", sep="")
	freq[i] <- as.numeric(dbGetQuery(con, command))
}

# 実験同士の共起
result <- c()
for(i in 1:5){
n <- i+1
	for(j in n:6){
		command1 <- paste("SELECT COUNT(*) FROM pubmed WHERE abst like '%", exp[i], "%", exp[j], "%';", sep="")
		cofreq1 <- as.numeric(dbGetQuery(con, command1))
		command2 <- paste("SELECT COUNT(*) FROM pubmed WHERE abst like '%", exp[j], "%", exp[i], "%';", sep="")
		cofreq2 <- as.numeric(dbGetQuery(con, command2)) 
		cofreq <- cofreq1 + cofreq2
		result <- rbind(result, paste(exp[i], "	", cofreq, "	", exp[j], sep=""))
	}
}

# CytoScape用結果を出力
# arrtibute.naとcollocation.sifをCytoScapeに入力
write.table(cbind(paste(exp, "="), freq), "attribute.na", row.names=F, col.names=F, quote=F)
write.table(result, "collocation.sif", row.names=F, col.names=F, quote=F)

# PDFファイルをダウンロードする関数を定義
pubmed <- function(input){
	filename <- paste0(input[3], "/", input[2], ".pdf")
	download.file(url = input[1], destfile = filename, quiet=T)
}

pubmed.download <- function(keyword){
preurl <- c()
url <- c(url, as.matrix(dbGetQuery(con, paste("SELECT url FROM pubmed WHERE abst like '%", keyword[1], "%';", sep=""))))

	for(i in 1:length(keyword)){
		preurl <- c(url, as.matrix(dbGetQuery(con, paste("SELECT url FROM pubmed WHERE abst like '%", keyword[i], "%';", sep=""))))
		url <- intersect(url, preurl)
	}

keywords <- c()
	if(length(keyword) == 1){
		keywords <- keyword
	}
	if(length(keyword) > 1){
		keywords <- keyword[1]
		for(j in 2:length(keyword)){
			keywords <- paste0(keywords, "_", keyword[j])
		}
	}

url <- sub(" ", "", url)
pmd <- sub("http://www.ncbi.nlm.nih.gov/pmc/articles/", "", url)
pmd <- sub("/pdf/", "", pmd)
file <- keywords
input <- cbind(url, pmd, file)

d <- getwd()
dir.create(paste0(d, "/", keywords))
apply(input, 1, try(pubmed))
}

# RNA-Seqというキーワードが含まれるPDFファイルを全てダウンロード
keyword <- "RNA-Seq"
pubmed.download(keyword)

# 複数キーワードでもOK
keyword <- c("RNA-Seq", "microarray")
pubmed.download(keyword)

```

# 3. MeSHパッケージの紹介
```r
## 以下は現状Developer版Rでのみ利用可能(2012.3.2)
## 安定して動くのはBioConductor2.13以降

## バイナリ版（クリックしていくだけ）
## Mac用 R-devel : http://R.research.att.com/
## Windows用 R-devel : http://cran.r-project.org/bin/windows/base/rdevel.html

## ソース版（自分でmakeしないといけない）
## http://cran.r-project.org/ > R Sources > R-devel.tar.gzをダウンロード 
tar xvf R-devel.tar.gz
cd R-devel
sudo ./configure --enable-R-shlib
sudo make
sudo make install


## パッケージダウンロード
source("http://bioconductor.org/biocLite.R")
biocLite("GO.db",type="source")
biocLite("hgu95av2.db",type="source")
biocLite("fdrtools",type="source")
biocLite("Category",type="source")
biocLite("cummeRbund",type="source")
biocLite("AnnotationForge",type="source")
biocLite("AnnotationDbi",type="source")
biocLite("DBI",type="source")
biocLite("RSQLite",type="source")

## Githubで以下の3パッケージをダウンロード
## MeSH.db : https://github.com/kokitsuyuzaki/MeSH.db からMeSH.db_1.0.tar.gz
## gendoo.Hs.db : https://github.com/dritoshi/gendoo.Hs.db からgendoo.Hs.db_0.99.0.tar.gz
## meshr : https://github.com/morota/meshr から meshr_0.99.0.tar.gz

# ビルド (ターミナルから)
R CMD INSTALL MeSH.db_1.0.tar.gz
R CMD INSTALL gendoo.Hs.db_0.99.0.tar.gz
R CMD INSTALL meshr_0.99.0.tar.gz

# ロード (R起動後)
library("GO.db")
library("hgu95av2.db")
library("fdrtools")
library("cummeRbund")
library("MeSH.db")
library("gendoo.Hs.db")
library("meshr")

# CummeRbund内のテストデータ用意
cuff <- readCufflinks(dir=system.file("extdata", package="cummeRbund"))
gene.symbols <- annotation(genes(cuff))[,4]
mySigGeneIds <- getSig(cuff,x='hESC',y='iPS',alpha=0.05,level='genes')
mySigGenes <- getGenes(cuff,mySigGeneIds)
sig.gene.symbols <- annotation(mySigGenes)[,4]
gene.symbols <- gene.symbols[!is.na(gene.symbols)]
sig.gene.symbols <- sig.gene.symbols[!is.na(sig.gene.symbols)]
geneid <- select(org.Hs.eg.db, keys=gene.symbols, keytype="SYMBOL", cols="ENTREZID")
sig.geneid <- select(org.Hs.eg.db, keys=sig.gene.symbols, keytype="SYMBOL", cols="ENTREZID")
na.index1 <- which(is.na(geneid[,2]))
for (i in na.index1){
	s <- unlist(strsplit(as.character(geneid[i,][1]), ","))[1]
	sym <- get(s, org.Hs.egALIAS2EG)[1]
	geneid[i,2] <- as.integer(sym)
}
na.index2 <- which(is.na(sig.geneid[,2]))
for (i in na.index2){
	s <- unlist(strsplit(as.character(sig.geneid[i,][1]), ","))[1]
	sym <- get(s, org.Hs.egALIAS2EG)[1]
	sig.geneid[i,2] <- as.integer(sym)
}
geneid <- geneid[!duplicated(geneid[,2]), ]
sig.geneid <- sig.geneid[!duplicated(sig.geneid[,2]), ]


## GO.dbによる遺伝子アノテーション

### BP(Biological Process)によるエンリッチメント解析
paraBP <- new("GOHyperGParams", geneIds=sig.geneid[,2], universeGeneIds=geneid[,2],annotation="hgu95av2.db",ontology="BP",pvalueCutoff=0.05,conditional=F,testDirection="over")
BP <- hyperGTest(paraBP)

### MF(Molecular Function)によるエンリッチメント解析
paraMF <- new("GOHyperGParams", geneIds=sig.geneid[,2], universeGeneIds=geneid[,2],annotation="hgu95av2.db",ontology="MF",pvalueCutoff=0.05,conditional=F,testDirection="over")
MF <- hyperGTest(paraMF)

### CC(Cellular Component)によるエンリッチメント解析
paraCC <- new("GOHyperGParams", geneIds=sig.geneid[,2], universeGeneIds=geneid[,2],annotation="hgu95av2.db",ontology="CC",pvalueCutoff=0.05,conditional=F,testDirection="over")
CC <- hyperGTest(paraCC)

### 結果を集計
summary(BP)
summary(MF)
summary(CC)

### 結果を保存
write.table(summary(BP),"GO_BP.txt")
write.table(summary(MF),"GO_MF.txt")
write.table(summary(CC),"GO_CC.txt")

## MeSH.dbによる遺伝子アノテーション

### A(Anatomy)によるエンリッチメント解析
paraA <- new("MeSHHyperGParams", geneIds=sig.geneid[,2], universeGeneIds=geneid[,2],annotation="GendooMeSHA", pvalueCutoff=0.05, pAdjust="none")
A <- meshHyperGTest(paraA)

### B(Organisms)によるエンリッチメント解析
paraB <- new("MeSHHyperGParams", geneIds=sig.geneid[,2], universeGeneIds=geneid[,2],annotation="GendooMeSHB", pvalueCutoff=0.05, pAdjust="none")
B <- meshHyperGTest(paraB)

### C(Diseases)によるエンリッチメント解析
paraC <- new("MeSHHyperGParams", geneIds=sig.geneid[,2], universeGeneIds=geneid[,2],annotation="GendooMeSHC", pvalueCutoff=0.05, pAdjust="none")
C <- meshHyperGTest(paraC)

### D(Chemicals and Drugs)によるエンリッチメント解析
paraD <- new("MeSHHyperGParams", geneIds=sig.geneid[,2], universeGeneIds=geneid[,2],annotation="GendooMeSHD", pvalueCutoff=0.05, pAdjust="none")
D <- meshHyperGTest(paraD)

### G(Phenomena and Processes)によるエンリッチメント解析
paraG <- new("MeSHHyperGParams", geneIds=sig.geneid[,2], universeGeneIds=geneid[,2],annotation="GendooMeSHG", pvalueCutoff=0.05, pAdjust="none")
G <- meshHyperGTest(paraG)

### 結果を集計
summary(A)
summary(B)
summary(C)
summary(D)
summary(G)

### 結果を保存
write.table(summary(A),"MeSH_A.txt")
write.table(summary(B),"MeSH_B.txt")
write.table(summary(C),"MeSH_C.txt")
write.table(summary(D),"MeSH_D.txt")
write.table(summary(G),"MeSH_G.txt")
```
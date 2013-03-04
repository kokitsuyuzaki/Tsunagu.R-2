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

# 2. RSQLiteの利用

# 3. MeSHパッケージの紹介



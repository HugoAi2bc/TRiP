Sys.getenv("RSTUDIO_PANDOC")
pandoc=Sys.getenv("RSTUDIO_PANDOC")
print(paste0(pandoc))
print(pandoc)
print(Sys.getenv("RSTUDIO_PANDOC"))
print(paste0(Sys.getenv("RSTUDIO_PANDOC")))

Sys.getenv("RSTUDIO_PANDOC") > /data/tmp1.txt
Sys.getenv("RSTUDIO_PANDOC") > ./tmp2.txt
pandoc=Sys.getenv("RSTUDIO_PANDOC")
print(paste0(pandoc)) > /data/tmp3.txt
print(pandoc) > /data/tmp4.txt
print(Sys.getenv("RSTUDIO_PANDOC")) > /data/tmp5.txt
print(paste0(Sys.getenv("RSTUDIO_PANDOC"))) > ./tmp6.txt

ips:
	Rscript -e "knitr::spin('bench-ips.R')"

glasso:
	Rscript -e "knitr::spin('bench-glasso.R')"

grid:
	Rscript -e "knitr::spin('bench-grid.R')"

table:
	Rscript -e "rmarkdown::render('bench-report.Rmd')"

move:
	Rscript -e "source('bench-move.R')"


all:
	Rscript -e "knitr::spin('bench-ips.R')"
	Rscript -e "knitr::spin('bench-glasso.R')"
	Rscript -e "knitr::spin('bench-grid.R')"
	Rscript -e "rmarkdown::render('bench-report.Rmd')"
	Rscript -e "source('bench-move.R')"










# all:
# 	Rscript -e "DATA='PROS'; knitr::spin('bench-ips.R')"
# 	Rscript -e "DATA='N01';  knitr::spin('bench-ips.R')"
# 	Rscript -e "DATA='PROS'; knitr::spin('bench-glasso.R')"
# 	Rscript -e "DATA='N01';  knitr::spin('bench-glasso.R')"
# 	Rscript -e "knitr::spin('bench-grid.R')"
# 	Rscript -e "rmarkdown::render('bench-report.Rmd')"

# table:
# 	Rscript -e "rmarkdown::render('bench-report.Rmd')"




# ipsp:
# 	Rscript -e "DATA='PROS'; knitr::spin('bench-ips.R')"

# ipsn:
# 	Rscript -e "DATA='N01';  knitr::spin('bench-ips.R')"

# ips:
# 	Rscript -e "DATA='PROS'; knitr::spin('bench-ips.R')"
# 	Rscript -e "DATA='N01';  knitr::spin('bench-ips.R')"


# glassop:
# 	Rscript -e "DATA='PROS'; knitr::spin('bench-glasso.R')"

# glasson:
# 	Rscript -e "DATA='N01';  knitr::spin('bench-glasso.R')"

# glasso:
# 	Rscript -e "DATA='PROS'; knitr::spin('bench-glasso.R')"
# 	Rscript -e "DATA='N01';  knitr::spin('bench-glasso.R')"

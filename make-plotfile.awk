# ls -d *.smcdir | awk -f make-plotfile.awk > plotfile

BEGIN{OFS="\t"}

/Gbb/{print "psmcdir", $1, "-", "#B94100"}
/Gbg/{print "psmcdir", $1, "-", "green"}
/Ggg/{print "psmcdir", $1, "-", "blue"}
/Ggd/{print "psmcdir", $1, "-", "magenta"}

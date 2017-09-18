# computational.biology

This app will take a text file with GWAS summary statistics as input and output Manhattan plot and QQ-plot.
Run it as an executable from command line to see all the available options. 

Usage: ./manqq.R [options]

Options:
        -f CHARACTER, --file=CHARACTER
                Text file with the columns P, CHROM and POS, it can also be compressed (.gz)

        -r CHARACTER, --chromosome=CHARACTER
                Column name for the chromosome, [default = CHR]

        -b CHARACTER, --basepair=CHARACTER
                Column name for the genomic position, [default = BP]

        -p CHARACTER, --pvalue=CHARACTER
                Column name for the association p-value, [default = P]

        -c LOGICAL, --cut=LOGICAL
                To speed up plotting, cut out snps with P > 0.05, [default = FALSE]

        -m CHARACTER, --manout=CHARACTER
                output file name [default = out(.manhattan.png)]

        -t NUMERIC, --threshold=NUMERIC
                threshold line for suggestive significance [default = 5e-06]

        -q CHARACTER, --qqout=CHARACTER
                output file name [default = out(qq.png)]

        -h, --help
                Show this help message and exit

for chromNumber in range(1, 23):
    dosageFile = open("hg38.chr" + str(chromNumber) + ".maf0.01.R20.8.dosage.txt", "r")
    annoFile = open("pred_db_anno_hg38.chr" + str(chromNumber) + ".maf0.01.R20.8.dosage.txt", "w")
    genoFile = open("pred_db_geno_hg38.chr" + str(chromNumber) + ".maf0.01.R20.8.dosage.txt", "w")

    for line in dosageFile:
        if line[:30] == "chr snp_ID pos ref alt AA_freq":
            header_genoFile = line[31:len(line)]
            genoFile.write("snp_ID " + header_genoFile)
            header_annoFile = "chr pos varID refAllele effectAllele rsid"
            annoFile.write(header_annoFile)
        else:
            column = line.split(" ")
            chrom = column[0]
            snpID = column[1]
            pos = column[2]
            refAllele = column[3]
            effectAllele = column[4].strip("\n")

            varID = chrom + "_" + pos + "_" + refAllele + "_" + effectAllele + "_b38"
            rsid = snpID

            anno = chrom + " " + pos + " " + varID + " " + refAllele + " " + effectAllele + " " + rsid
            annoFile.write("\n" + anno)

            # writing the geno file

            geno = snpID + " "
            for i in range(6, len(column)):
                if i == len(column) - 1:
                    geno += column[i]
                else:
                    geno += column[i] + " "
                    ++i
            genoFile.write(geno)

    dosageFile.close()
    annoFile.close()
    genoFile.close()
++chromNumber
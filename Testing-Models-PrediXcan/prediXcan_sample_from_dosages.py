predixcanSampleFile = open("prediXcan_samples.txt", "w")

with open("hg38.chr1.maf0.01.R20.8.dosage.txt", "r") as dosageFileChr1:
    header_dosage = dosageFileChr1.readline().strip("\n").split(" ")

    AA_freq = header_dosage.index("AA_freq") + 1

    for i in range(AA_freq, len(header_dosage)):
        if i == len(header_dosage):
            predixcanSampleFile.write(header_dosage[i] + "\t" + header_dosage[i])
        else:
            predixcanSampleFile.write(header_dosage[i] + "\t" + header_dosage[i] + "\n")

dosageFileChr1.close()
predixcanSampleFile.close()

import logging


def normalize_vcf(vcf_in: str, vcf_out: str):
    outText = ""
    print("Reading file " + vcf_in + "...")
    vcf = open(vcf_in, "r")
    j = 0
    print("Modifying contents...")
    for line in vcf:
        
        j = j+1
        line = line.strip()

        # preserve headers
        if line.startswith("#"):
            outText = outText + line + "\n"
            continue

        colValues = line.split("\t")

        # Remove any records that don't have a FILTER of "PASS" or "LowGQX"
        if(colValues[6] != "PASS" and colValues[6] != "LowGQX"):
            continue

        # The final file should only have variant rows with an alternate allele.
        if(colValues[4] == "."):
            continue

        # The INFO column should only contain this variant's specific allele frequency, i.e. "AF={actual value}"

        # This dataset does not contain AF values. Using AF1000G instead as per email clarification.
        infoParts = colValues[7].split(";")
        AFs = list(filter(lambda x: "AF1000G=" in x, infoParts))

        
        assert len(AFs) <= 1, "Unexpected values in INFO column on line " + str(j) + ": " + line

        if len(AFs) == 0:
            colValues[7] = "."
        else:
            AF = AFs[0].split("=")
            colValues[7] = "AF=" + AF[1] # Technically "AF" is a reserved subfield but following given instructions


        # The FORMAT column should only include the following values "GT:AD:DP"
        formatCodes = colValues[8].split(":")
        formatValues = colValues[9].split(":") # NA12878 column
        codes = list()
        values = list()
        for i in range(len(formatCodes)):
            if(formatCodes[i] == "GT" or formatCodes[i] == "AD" or formatCodes[i] == "DP"):
                codes.append(formatCodes[i])
                #get the corresponding values from NA12878 column
                values.append(formatValues[i])

        colValues[8] = ":".join(codes)
        colValues[9] = ":".join(values)



        newline = "\t".join(colValues)
        outText = outText + newline + "\n"

    #last item will have a trailing \n. Remove it
    outText = outText.strip()



    vcf.close()
    return outText




# Remove any records that don't have a FILTER of "PASS" or "LowGQX" X
# The INFO column should only contain this variant's specific allele frequency, i.e. "AF={actual value}" (GMAF)
# The FORMAT column should only include the following values "GT:AD:DP"
# The final file should only have variant rows with an alternate allele. X
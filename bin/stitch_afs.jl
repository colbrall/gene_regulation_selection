# stitch_afs.jl
#
# @author laura colbran
# super-ugly code to stitch vcftools frequency files together. not generalizable.

function thousGensMAF()
    chrs = 1:1:22
    for CHR in chrs
        # files for each population
        ACB = eachline(open("chr$(CHR).ACB.frq"))
        ASW = eachline(open("chr$(CHR).ASW.frq"))
        BEB = eachline(open("chr$(CHR).BEB.frq"))
        CDX = eachline(open("chr$(CHR).CDX.frq"))
        CEU = eachline(open("chr$(CHR).CEU.frq"))
        CHB = eachline(open("chr$(CHR).CHB.frq"))
        CHS = eachline(open("chr$(CHR).CHS.frq"))
        CLM = eachline(open("chr$(CHR).CLM.frq"))
        ESN = eachline(open("chr$(CHR).ESN.frq"))
        FIN = eachline(open("chr$(CHR).FIN.frq"))
        GBR = eachline(open("chr$(CHR).GBR.frq"))
        GIH = eachline(open("chr$(CHR).GIH.frq"))
        GWD = eachline(open("chr$(CHR).GWD.frq"))
        IBS = eachline(open("chr$(CHR).IBS.frq"))
        ITU = eachline(open("chr$(CHR).ITU.frq"))
        JPT = eachline(open("chr$(CHR).JPT.frq"))
        KHV = eachline(open("chr$(CHR).KHV.frq"))
        LWK = eachline(open("chr$(CHR).LWK.frq"))
        MSL = eachline(open("chr$(CHR).MSL.frq"))
        MXL = eachline(open("chr$(CHR).MXL.frq"))
        PEL = eachline(open("chr$(CHR).PEL.frq"))
        PJL = eachline(open("chr$(CHR).PJL.frq"))
        PUR = eachline(open("chr$(CHR).PUR.frq"))
        STU = eachline(open("chr$(CHR).STU.frq"))
        TSI = eachline(open("chr$(CHR).TSI.frq"))
        YRI = eachline(open("chr$(CHR).YRI.frq"))

        files_to_stitch = zip(ACB,ASW,BEB,CDX,CEU,CHB,CHS,CLM,ESN,FIN,GBR,GIH,GWD,IBS,ITU,
                            JPT,KHV,LWK,MSL,MXL,PEL,PJL,PUR,STU,TSI,YRI)
        POPS = join(["ACB","ASW","BEB","CDX","CEU","CHB","CHS","CLM","ESN","FIN","GBR","GIH","GWD","IBS","ITU",
                            "JPT","KHV","LWK","MSL","MXL","PEL","PJL","PUR","STU","TSI","YRI"],"\t")
        open("chr$(CHR).ALL.frq","w") do outf
            write(outf,"CHROM\tPOS\tN_ALLELES\tN_CHR\t$(POPS)\n")
            for line in files_to_stitch
                if startswith(line[1],"CHROM") continue end
                base_line = split(line[1],"\t")
                base_line[5] = join(base_line[5:end],",")
                deleteat!(base_line,6:length(base_line))
                for pop in line[2:end]
                    base_line = vcat(base_line,join(split(pop,"\t")[5:end],","))
                end
                full_line = join(base_line,"\t")
                write(outf,"$(full_line)\n")
            end
        end

    end
end

function hgdpMAF()
    chrs = 1:1:22
    for CHR in chrs
        AFR = eachline(open("chr$(CHR).AFRICA.frq"))
        OCE = eachline(open("chr$(CHR).OCEAN.frq"))
        AMER = eachline(open("chr$(CHR).AMER.frq"))
        EURO = eachline(open("chr$(CHR).EURO.frq"))
        CS_ASIA = eachline(open("chr$(CHR).CS_ASIA.frq"))
        M_EAST = eachline(open("chr$(CHR).M_EAST.frq"))
        E_ASIA = eachline(open("chr$(CHR).E_ASIA.frq"))

        files_to_stitch = zip(AFR,AMER,EURO,M_EAST,CS_ASIA,E_ASIA,OCE)
        POPS = join(["AFRICA","AMER","EURO","M_EAST","CS_ASIA","E_ASIA","OCEAN"],"\t")
        
        open("chr$(CHR).ALL.frq","w") do outf
            write(outf,"CHROM\tPOS\tN_ALLELES\tN_CHR\t$(POPS)\n")
            for line in files_to_stitch
                if startswith(line[1],"CHROM") continue end
                base_line = split(line[1],"\t")
                base_line[5] = join(base_line[5:end],",")
                deleteat!(base_line,6:length(base_line))
                for pop in line[2:end]
                    base_line = vcat(base_line,join(split(pop,"\t")[5:end],","))
                end
                full_line = join(base_line,"\t")
                write(outf,"$(full_line)\n")
            end
        end
    end
end

function main()
    # thousGensMAF()
    hgdpMAF()
end

main()

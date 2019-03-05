import argparse
import re
import cairo
import random
def filelist():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta', help = "Fasta file with genes, capitalized exons and lower case introns",
        required = True, type = str)
    parser.add_argument('-m', '--motif', help = "Motif list, separated with newline characters",
        required = True, type = str)
    parser.add_argument('-F', '--folder', help = "Output folder for svg files",
        required = False, type = str, default = "")
    parser.add_argument('-W', '--width', help = "Width of svg output image, defaults to 256",
        required = False, type = int, default = 256)
    parser.add_argument('-H', '--height', help = "Height of svg output image, defaults to 256",
        required = False, type = int, default = 256)
    return parser.parse_args()
args = filelist()
fasta = args.fasta
motifs = args.motif
width = args.width
height = args.height
folder = args.folder
def motifgen(file):
    '''Generates list of motifs from a .txt file'''
    motiflist = []
    with open(motifs, "r") as m:
        for line in m:
            line = line.strip()
            motiflist.append(line)
    nummotif = len(motiflist)
    return(motiflist, nummotif)
def motifconvert(motif):
    '''Creates regex search phrase list from motif ordered as [main strand intron, main strand exon, reverse complement intron, reverse complement exon]'''
    regexmotif = []
    mainin = []
    mainex = []
    revin = []
    revex = []
    count = 0
    motif = motif.lower()
    for x in motif:
        if x == "y": #ambiguous base
            mainin.append("[ct]")
            revin = ["[ag]"] + revin
            if count == 0:
                mainex.append("[CT]")
                revex = ["[AG]"] + revex
            else:
                mainex.append("[ctCT]")
                revex = ["[agAG]"] + revex
        elif x == "a":
            mainin.append("a")
            revin = ["t"] + revin
            if count == 0:
                mainex.append("A")
                revex = ["T"] + revex
            else:
                mainex.append("[aA]")
                revex = ["[tT]"] + revex
        elif x == "t":
            mainin.append("t")
            revin = ["a"] + revin
            if count == 0:
                mainex.append("T")
                revex = ["A"] + revex
            else:
                mainex.append("[tT]")
                revex = ["[aA]"] + revex
        elif x == "c":
            mainin.append("c")
            revin = ["g"] + revin
            if count == 0:
                mainex.append("C")
                revex = ["G"] + revex
            else:
                mainex.append("[cC]")
                revex = ["[gG]"] + revex
        elif x == "g":
            mainin.append("g")
            revin = ["c"] + revin
            if count == 0:
                mainex.append("G")
                revex = ["C"] + revex
            else:
                mainex.append("[gG]")
                revex = ["[cC]"] + revex
        count += 1
    regexmotif = ["".join(mainin), "".join(mainex), "".join(revin), "".join(revex)]
    return(regexmotif)
def findmotif(regexmotif, sequence, direction):
    '''Creates lists of all introninc and exonic positions of motifs in a given gene'''
    intronlist = []
    exonlist = []
    if direction == "F":
        intron = regexmotif[0]
        exon = regexmotif[1]
    else:
        intron = regexmotif[2]
        exon = regexmotif[3]
    for x in re.finditer(intron, sequence):
        intronlist.append(x.start())
    for x in re.finditer(exon, sequence):
        exonlist.append(x.start())
    return(intronlist, exonlist)
def cairoconvert(motiflist, genes_dict, width):
    '''Scales exon location/lengths and motif locations to svg file dimensions'''
    W = width
    motifdict = dict()
    for gene in genes_dict.keys(): #Structure is "Gene name" : [[exonlocs and len], [[[intron motif][exon motif]]]
        motifdict[gene] = [] #Key for new gene
        mlist = genes_dict[gene]
        scalar = W/mlist[1] #Scaling everything based on length of gene and svg scale
        exonloc = mlist[2]
        for x in range(0,len(exonloc)): #scaling exon locations and lengths
            for y in range(0,2):
                exonloc[x][y] = exonloc[x][y]*scalar
        motifdict[gene].append(exonloc)
        for y in range(0, len(motiflist)):
            motifloc = mlist[y+3]
            for a in range(0,2): #First loop covers introns, second covers exons
                for b in range(0,len(motifloc[a])):
                    motifloc[a][b] = motifloc[a][b]*scalar
            motifdict[gene].append(motifloc)
        motifdict[gene].append(scalar)
    return(motifdict)
def cairoout(gene,motifdictlist, motiflist, w, h, folder):
    '''Generates SVG files for each gene'''
    motifs = motiflist
    if folder[-1] == "/":
        name = folder + gene + ".svg"
    else:
        name = folder + "/" + gene + ".svg"
    exons = motifdictlist[0]
    scalar = motifdictlist[-1]
    width, height = w, h
    surface = cairo.SVGSurface(name, width, height)
    context = cairo.Context(surface)
    #Gene name
    context.move_to(10,10)
    context.show_text(gene)
    context.stroke()


    #Full line for gene
    context.set_source_rgb(0.0, 0.0, 0.0)
    context.set_line_width(2)
    context.move_to(5, height/2)
    context.line_to(width-5, height/2)
    context.stroke()

    #Exons
    context.set_line_width(12)
    for x in exons:
        context.move_to(x[0]+5, height/2)
        context.line_to(x[0]+5+x[1], height/2)
        context.stroke()
    #Motifs label
    for n in range(0,len(motifs)):
        motiflen = len(motifs[n])
        revn = len(motifs)-n
        num = n+1
        inmotifs = motifdictlist[1+n][0]
        exmotifs = motifdictlist[1+n][1]
        context.set_source_rgb(random.randint(0,1)*(num+.4),random.randint(0,1)*(num-.3),random.randint(0,1)*(num + .15))
        context.set_line_width(10)
        for i in inmotifs:
            context.move_to(i+5, height/2)
            context.line_to(i+5 + motiflen * scalar, height/2)
            context.stroke()
        context.set_line_width(15)
        for e in exmotifs:
            context.move_to(e+5, height/2)
            context.line_to(e+5 + motiflen * scalar, height/2)
            context.stroke()
        context.set_line_width(4)
        context.move_to(10, height - (8*revn))
        context.line_to(15, height - (8*revn))
        context.stroke()
        context.line_to(17, height - (8*revn)+2)
        context.set_font_size(8)
        context.set_source_rgb(0,0,0)
        context.show_text(motifs[n])
        context.stroke()
    surface.finish()
    return()
motiflist, nummotif = motifgen(motifs)
count = 0
regexmotiflist = []
for x in motiflist:
    regexmotiflist.append(motifconvert(motiflist[count]))
    count +=1
genes = dict()
with open(fasta, "r") as f:
    sequence = ""
    for line in f:
        if line[0] == ">": #Header line
            if len(genes.keys()) > 0: #Not first time through loop
                genes[gene].append(len(sequence)) # Dictionary list is now [direction, length]
                genes[gene].append([])
                for x in  re.finditer(r"[A-Z]+", sequence): #Dictionary list is now [direction, length, exon_locations/len]
                    genes[gene][2].append([x.start(), x.end() - x.start()])
                for n in range(0, nummotif): #For each motif, search the sequence for intronic and exonic instances
                    intron, exon = findmotif(regexmotiflist[n], sequence, genes[gene][0])
                    motifloc = [intron, exon]
                    genes[gene].append(motifloc) #Dictionary list is now [direction, length, exon_locations/len, [motif_n_intron_locs, motif_n_exon_locs],...]
                sequence = ""
            else: #First time through the loop
                pass
            if "reverse complement" in line: # New gene header
                rev = "T"
                gene = line[1:-22]
                genes[gene] = [rev]
            else:
                rev = "F"
                gene = line[1:-1]
                genes[gene] = [rev]
        else:
            line = line.strip()
            sequence = sequence + line
genes[gene].append(len(sequence))
genes[gene].append([])
for x in  re.finditer(r"[A-Z]+", sequence):
    genes[gene][2].append([x.start(), x.end() - x.start()])
for n in range(0, nummotif): #Final sequence search
    intron, exon = findmotif(regexmotiflist[n], sequence, genes[gene][0])
    motifloc = [intron, exon]
    genes[gene].append(motifloc)
cairoready = cairoconvert(motiflist, genes, width)
for x in cairoready:
    cairoout(x, cairoready[x], motiflist, width, height, folder)

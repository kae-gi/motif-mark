#!/usr/bin/env python
"""
Motif marker using object oriented programming.
Author: Kaetlyn Gibson
BI 625
"""
import cairo
import argparse
import re
import random

# == CLASSES ==
class Gene:
    def __init__(self, name, sequence):
        """Initializer"""
        ## Data ##
        self.name = name
        self.sequence = sequence
    # == METHODS ==
    def getSeqLoc(self) -> int:
        """ Returns end location of sequence (start is always 0). """
        return len(self.sequence)
    
    def getELoc(self) -> list:
        """ 
        Gets the start, end location of exon parts of the sequence
        and returns them in tuple form within a list.
        """
        exon = list(re.finditer(r'[A-Z]*', self.sequence))
        exon = [x.span() for x in exon if x.group(0) != '']
        return exon

    def getMotifLocs(self, motif_regex:dict) -> dict:
        """ 
        Gets the start, end location of exon or intron motifs and 
        returns in the form of a dictionary.
        """
        loc_dict = {}     
        for motif in motif_regex:
            motifs = list(filter(None, re.finditer(r'%s' %motif_regex[motif], gene.sequence)))
            motifs = [x.span() for x in motifs]
            loc_dict[motif] = motifs
        return loc_dict


# == HELPER FUNCTIONS ==
def getArgs():
    """ Global variable flags. """
    parser = argparse.ArgumentParser(description='A program for...')
    parser.add_argument('-f', '--fasta', help='Specify the input fasta file name.', type=str, required=True)
    parser.add_argument('-m', '--motif', help='Specify the input motif file name.', type=str, required=True)
    return parser.parse_args()	
args = getArgs()

def oneLineFasta(fastaFile: str) -> dict:
    """
    Takes in a string for a file name as input. Returns a dictionary
    with the fasta record headers as keys and fasta record sequences
    as values.
    """
    fasta_dict = {}
    with open(args.fasta, 'r') as fa_file:
        headers = []
        seqs = []
        seq = ''
        for i, line in enumerate(fa_file):
            line = line.strip()
            if '>' in line[0]:
                header = line
                if i > 1:
                    seqs.append(seq)
                    seq = ''
                headers.append(header)
            else:
                seq += line
        seqs.append(seq)
        for i in range(len(headers)):
            fasta_dict[headers[i]] = seqs[i]
    return fasta_dict

def getRegexStr(IorE:str, motif_list:list) -> dict:
    """
    Given choice of "intron" or "exon", takes the list of motifs 
    and returns a dictionary consisting of motif as key
    and the regex string for that motif as value.
    """
    motif_regex = {}
    regex_str = ''
    for motif in motif_list:
        motif = motif.upper()
        base_count = 0
        for base in motif:
            if base == 'Y':
                base_count += 1
                if base_count > 1:
                    try:
                        motif[base_count]
                    except:
                        regex_str += '{' + str(base_count) + '}'
                else:
                    regex_str += '[C,T]'
            else:
                base_count = 0
                regex_str += base
        if IorE == "intron":
            motif = motif.lower()
            regex_str = regex_str.lower()
        motif_regex[motif] = regex_str
        regex_str = ''
    return motif_regex

def extract_motifs(motif_file:str) -> tuple:
    """
    Extracts the exon and intron motifs from given motif file and
    formulates the regex string for each motif. Then
    returns a tuple that consists of 3 dictionaries -
    the exon motif regex, intron motif regex, and color dictionary.
    """
    intron_motifs = []
    exon_motifs = []
    color_dict = {}
    with open(motif_file, 'r') as mf:
        for motif in mf:    
            exon_motifs += list(filter(None, re.findall(r'[A-Z]*', motif)))
            intron_motifs += list(filter(None, re.findall(r'[a-z]*', motif)))
            color_dict[motif.strip()] = ''
    exon_motif_regex = getRegexStr('exon', exon_motifs)
    intron_motif_regex = getRegexStr('intron', intron_motifs)
    color_dict = generate_colors(color_dict)
    return exon_motif_regex, intron_motif_regex, color_dict

def extract_genes(fasta_file:str) -> list:
    """ 
    Extracts the gene name and sequence from given fasta file and creates
    an object of class Gene based on the gene. Then returns the list
    of Gene objects. 
    """
    genes = []
    fasta_dict = oneLineFasta(fasta_file)
    for record in fasta_dict:
        name = record.strip('>')
        sequence = fasta_dict[record]
        genes.append(Gene(name, sequence))
    return genes

def draw_rectangle(cr, x:float, y:float, width:float, height:float, color:tuple) -> None:
    """ Draws a rectangle. """
    # rectangle
    cr.rectangle(x, y, width, height) # x, y coord of top-left corner, width, height
    cr.set_source_rgba(color[0], color[1], color[2]) # color
    cr.fill() # fill rectangle with color
    return None

def draw_line(cr, base_x:float, x:float, y:float) -> None:
    """ Draws a line. """
    # line
    cr.set_line_width(5) # line thickness
    cr.move_to(base_x, y) # x, y (start)
    cr.line_to(x, y) # x, y (end)
    cr.set_source_rgb(0, 0, 0) # black
    cr.stroke() # draw the line
    return None

def write(cr, x:float, y:float, text:str, size:int, color:tuple) -> None:
    """ Writes text. """
    cr.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
    cr.set_font_size(size)
    cr.set_source_rgb(color[0], color[1], color[2])
    cr.move_to(x,y)
    cr.show_text(text)
    return None

def generate_colors(color_dict:dict)-> dict:
    """ 
    Chooses a color from the color palette. Limited colors
    since it is assumed that only 5 motifs max are provided.
    """
    color_palette = [(226/255,124/255,181/255),(230/255,159/255,0),(86/255,180/255,233/255),(0,158/255,115/255),(240/255,228/255,66/255),(0,103/255,160/255),(213/255,72/255,0)]
    for motif in color_dict:
        value = random.randint(0,len(color_palette)-1)
        color_dict[motif] = color_palette[value]
        color_palette.pop(value)  
    return color_dict 


# == MAIN ==
# get data
intron_motif_regex, exon_motif_regex, color_dict = extract_motifs(args.motif)
genes = extract_genes(args.fasta)

# canvas width: width of sequences <= 1000, add 200px padding
# canvas height: height adjusts based on number of genes and motifs
pad = 175
c_width, c_height = 1000+pad, ((len(genes))*60)+len(color_dict)*60+60
surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, c_width, c_height)
ctx = cairo.Context(surface)
draw_rectangle(ctx, 0, 0, c_width, c_height, (1,1,1)) # canvas

# drawing
end_y_loc = 0
for i, gene in enumerate(genes):
    seq_loc = gene.getSeqLoc()
    e_loc = gene.getELoc()
    e_motif_locs = gene.getMotifLocs(exon_motif_regex)
    i_motif_locs = gene.getMotifLocs(intron_motif_regex)
    y_loc = ((c_height-(len(color_dict)*30)-30)/len(genes))*(i+1)+20
    name = gene.name.split()[0]
    info = ' '.join(gene.name.split()[1:])
    write(ctx, 50, y_loc+8, name, 20, (0,0,0)) # gene name
    write(ctx, pad, y_loc-22, info, 20, (0,0,0)) # gene info
    draw_line(ctx, pad, pad+seq_loc, y_loc) # sequence line
    draw_rectangle(ctx, pad+e_loc[0][0], y_loc-14, (e_loc[0][1]-e_loc[0][0]), 28, (0,0,0)) # exon
    for motif in e_motif_locs:
        if e_motif_locs[motif] != []:
            for loc in e_motif_locs[motif]:
                draw_rectangle(ctx, pad+loc[0], y_loc-10, (loc[1]-loc[0]), 20, color_dict[motif]) # exon motif
    for motif in i_motif_locs:
        if i_motif_locs[motif] != []:
            for loc in i_motif_locs[motif]:
                draw_rectangle(ctx, pad+loc[0], y_loc-10, (loc[1]-loc[0]), 20, color_dict[motif]) # intron motif
    end_y_loc = y_loc+8

# color key
end_y_loc += 60
write(ctx, pad, end_y_loc, 'Key:', 20, (0,0,0))
write(ctx, pad+75, end_y_loc, 'EXON', 20, (0,0,0))
write(ctx, pad+250, end_y_loc, 'intron', 20, (0,0,0))
exon_i, intron_i = 1, 1
for motif in color_dict:
    exon_y_loc, intron_y_loc  = end_y_loc, end_y_loc
    if motif == motif.upper():
        exon_y_loc += 25*(exon_i)
        write(ctx, pad+75, exon_y_loc, motif, 20, color_dict[motif])
        exon_i+=1
    else:
        intron_y_loc += 25*(intron_i)
        write(ctx, pad+250, intron_y_loc, motif, 20, color_dict[motif])
        intron_i+=1
        
# write image to png 
img_name = re.match(r'^(.*?(?=\.fasta))',args.fasta).group(0)
write(ctx, 50, 50, img_name, 40, (0,0,0)) # gene name
surface.write_to_png(img_name+".png")

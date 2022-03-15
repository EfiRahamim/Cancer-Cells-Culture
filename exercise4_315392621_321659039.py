"""
Author: Rivka Bouhnik & Efraim Rahamim
ID: 321659039, 315392621
"""
import json
import math
import random
import sys
import re
import csv
import os
import Bio
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO


class Polymerase:
    # constructor for the Polymerase object
    def __init__(self, type, error_rate=0):
        # validation check
        assert type == 'DNA' or type == 'RNA', type
        self.type = type
        # validaition check
        assert 0 <= error_rate <= 1, error_rate
        self.error_rate = error_rate

    def check_dna_input(self, seq):
        """ this function checks if a given nucleotides sequence is valid"""
        nucleotides = ['A', 'a', 'T', 't', 'C', 'c', 'G', 'g']
        for nuc in seq:
            if nuc not in nucleotides:
                return False
        return True

    def insert_mutation(self, translated_seq):
        "this function puts mutations randomly in K places in the given seq"
        # define K = N * error_rate
        k = math.ceil(len(translated_seq) * self.error_rate)
        # create list of random indexes
        # random.seed(1)
        indexes = random.sample(range(0, len(translated_seq) - 1), k)
        # create lists of DNA and RNA nucleotides
        dna_nuc = ['A', 'T', 'C', 'G']
        rna_nuc = ['A', 'U', 'C', 'G']
        # run on each random index
        for index in indexes:
            original_nuc = translated_seq[index]  # save the original nucleotide
            # choose random nucleotide of DNA or RNA beside the original nuc
            if self.type == 'DNA':
                new_dna_nuc = dna_nuc.copy()
                new_dna_nuc.remove(original_nuc)
                mutation_nuc = random.choice(new_dna_nuc)
            elif self.type == 'RNA':
                new_rna_nuc = rna_nuc.copy()
                new_rna_nuc.remove(original_nuc)
                mutation_nuc = random.choice(new_rna_nuc)
            # replace nucleotides
            translated_seq = translated_seq[:index] + mutation_nuc + translated_seq[index + 1:]
        # return the new sequence with the mutations
        return translated_seq

    def transcribe(self, dna_seq):
        """this function gets a sequence and returns the RNA/DNA seq that would be transcripted from it"""
        # Input integrity check
        assert self.check_dna_input(dna_seq), dna_seq
        # declare the RNA/DNA seq that would be printed
        trans_seq = ''
        # initialize matching dictionary
        if self.type == 'RNA':
            match_dict = {'A': 'U', 'a': 'U',
                          'T': 'A', 't': 'A',
                          'G': 'C', 'g': 'C',
                          'C': 'G', 'c': 'G'}
        else:
            match_dict = {'A': 'T', 'a': 'T',
                          'T': 'A', 't': 'A',
                          'G': 'C', 'g': 'C',
                          'C': 'G', 'c': 'G'}
        # running over the given seq
        for i in range(1, len(dna_seq) + 1):
            # the current nucleotide that we need to transcript
            # It goes from the end of the seq to the start because of the directionality from 5' to 3'
            nuc = dna_seq[-i]
            # replace the nucleotide to it's compatible nucleotide
            trans_seq = trans_seq + match_dict[nuc]
        # add mutation to the seq according to the error rate
        trans_seq = self.insert_mutation(trans_seq)
        # return the translated seq
        return trans_seq


class Ribosome:
    # constructor for the Ribosome object
    def __init__(self, genetic_code, start_codons):
        self.genetic_code = genetic_code
        self.start_codons = start_codons

    def check_rna_input(self, seq):
        """ this function checks if a given nucleotides sequence is valid"""
        nucleotides = ['A', 'a', 'U', 'u', 'C', 'c', 'G', 'g']
        for nuc in seq:
            if nuc not in nucleotides:
                return False
        return True

    def find_longest_frame(self, rna_seq):
        """this function returns the first index of the longest reading frame in the given RNA seq"""
        # initialize the maximum frame size
        max_frame = 0
        # initialize the index of the longest reading frame
        best_index = None
        # running on the RNA seq
        for i in range(0, len(rna_seq)):
            # check for start codon
            if rna_seq[i:i + 3] in self.start_codons:
                # initialize the index of the Start codon
                first_index = i
                # initialize the index of the stop codon, if exist
                last_index = first_index
                # running on the rest of the RNA seq, in threes, to find stop codon
                for j in range(i + 3, len(rna_seq), 3):
                    # check for stop codon
                    if j + 2 < len(rna_seq) and self.genetic_code[rna_seq[j:j + 3]] is None:
                        # update the index of the stop codon
                        last_index = j
                        # stop the loop and continue to the next reading frame, if excist
                        break
                # check if we found stop codon or the sequence is over
                if last_index == first_index:
                    # there is no stop codon, update the index of the 'stop codon' to be the last index of the seq
                    last_index = len(rna_seq)
                # check if this reading frame is the longest
                if last_index - first_index > max_frame:
                    # update the max frame  length
                    max_frame = last_index - first_index
                    # update the index of the best reading frame
                    best_index = first_index
        # return the first index of the longest reading frame
        return best_index

    def translate(self, rna_seq):
        """ this function translate the RNA seq into it's possible longest reading frame"""
        # check validation of input
        assert self.check_rna_input(rna_seq), rna_seq
        # make sequence uppercase
        rna_seq = rna_seq.upper()
        # initialize the string that the function returns
        reading_frame = []
        # find the index of the first AUG in the longest reading frame
        index_of_longest_frame = self.find_longest_frame(rna_seq)
        # check if there is a reading frame at all
        if index_of_longest_frame is None:
            return None
        else:
            # print all the three codon in the reading frame
            for i in range(index_of_longest_frame, len(rna_seq), 3):
                # # make sure not to exceed the indexes of the RNA seq and not to include stop codons
                if (i + 2) < len(rna_seq) and self.genetic_code[rna_seq[i:i + 3]] is not None:
                    reading_frame.append(rna_seq[i:i + 3])
                else:
                    break
        # returns the longest reading frame that was found (or none otherwise)
        return reading_frame

    def synthesize(self, rna_seq):
        """this function synthesize a protein from a given rna sequence"""
        # get the longest ORF of the rna seq
        orf = self.translate(rna_seq)
        if orf is None:
            return None
        # translate the codons to amino acids
        protein = ""
        for codon in orf:
            protein = protein + self.genetic_code[codon]
        # return the amino acids sequence -> the protein
        return protein


class Cell:
    def __init__(self, name, genome, num_copies, genetic_code, start_codons, division_rate):
        self.name = name
        # validation check for the genome
        flag = True
        for seq in genome:
            for nuc in seq:
                if nuc not in ['A', 'a', 'T', 't', 'G', 'g', 'C', 'c']:
                    flag = False
        assert flag, genome
        self.genome = genome
        # validation check for the num_copies
        assert isinstance(num_copies, int) and num_copies > 0, num_copies
        self.num_copies = num_copies
        self.genetic_code = genetic_code
        self.start_codons = start_codons
        assert isinstance(division_rate, int) and division_rate > 1, division_rate
        self.division_rate = division_rate
        # initialize polymerases and ribosom of the cell
        self.dna_poly = Polymerase('DNA')
        self.rna_poly = Polymerase('RNA')
        self.ribosome = Ribosome(self.genetic_code, self.start_codons)

    # print object as requested
    def __repr__(self):
        return "<" + self.name + "," + str(self.num_copies) + "," + str(self.division_rate) + ">"

    # implement mitosis by copying the cell in times of divion rate
    def mitosis(self):
        return self * self.division_rate

    # implement mitosis by the operator * (operator overloading)
    def __mul__(self, num):
        return [self] * num

    # implement meiosis of the cell
    def meiosis(self):
        # odd dna copies cannot undergo meiosis
        if self.num_copies % 2 != 0:
            return None
        # create first cell with exact genome and half of num copies
        cell1 = Cell(self.name, self.genome, int(self.num_copies / 2), self.genetic_code, self.start_codons,
                     self.division_rate)
        # second cell should have the complementary strands of the genome
        # create the complementary strands
        comp = []
        for sence in self.genome:
            comp.append(self.dna_poly.transcribe(sence))
        # create the second cell
        cell2 = Cell(self.name, comp, int(self.num_copies / 2), self.genetic_code, self.start_codons,
                     self.division_rate)
        # return the two cells in a list
        return [cell1, cell2]

    def check_repeats(self, repeat_seq, repeat_seq_len, long_seq):
        """ the function checks if a short given sequence repeats in a long given sequence and counts the times"""
        counter = 0
        # running on the long sequence
        for i in range(0, len(long_seq), repeat_seq_len):
            # checks if the sequence appears
            if repeat_seq == long_seq[i:i + repeat_seq_len]:
                counter = counter + 1
            else:
                # no need to check the rest of the long sequence
                break
        return counter

    def find_srr(self, dna_seq):
        """ the function finds the maximum repeats of each SSR in a DNA sequence """
        # validation check for the sequence
        flag = True
        for seq in dna_seq:
            if seq not in ['A', 'a', 'T', 't', 'G', 'g', 'C', 'c']:
                flag = False
        assert flag, dna_seq
        # change seq to uppercase
        dna_seq = dna_seq.upper()
        # declare dictionary of SRR repeats that the function returns
        # the keys are the repeats and the values are the numbers
        srr_dict = {}
        # dictionary only for the maximum number of repeats
        max_srr_dict = {}
        # sorte dict by keys
        sort_srr_dict = {}
        # flag is an indicator to know if we found repeats
        flag = 0
        # finding sequences in sizes 1 to 6
        for seq_size in range(1, 7):
            #
            for i in range(len(dna_seq)):
                # define the current repeat seq that we are looking for
                repeat_seq = dna_seq[i:i + seq_size]
                # define the rest of the sequence, where we'll look for the repeated seq
                rest_of_seq = dna_seq[i: len(dna_seq)]
                # checks if the current sequence repeats in the rest of the sequence
                counter = self.check_repeats(repeat_seq, seq_size, rest_of_seq)
                # checks if the repeats if more then 2
                if counter > 2:
                    # update the flag to know that we found repeats
                    flag = 1
                    # check if the sequence was already found
                    if repeat_seq in srr_dict.keys():
                        # add the counter as a value
                        srr_dict[repeat_seq] = srr_dict[repeat_seq] + [counter]
                    else:
                        # add the sequence as a new key and the counter as a value
                        srr_dict[repeat_seq] = [counter]
        # check if we found any repeats
        if flag == 0:
            return None
        else:
            # get only the maximum number of repeats for each repeat
            for key, value in srr_dict.items():
                max_srr_dict[key] = max(value)
        # sort the max dict by the repeats
        for key in sorted(max_srr_dict):
            sort_srr_dict[key] = max_srr_dict[key]
        # return srr's that were found
        return sort_srr_dict

    # change the dictionary of srr to string as requested in the output format
    def srr_to_string(self, srr):
        srr_str = ""
        isFirst = True
        # run over each repeat and number in the srr dictionary
        for key, value in srr.items():
            # don't add ";" to the first srr
            if isFirst:
                srr_str = key + "," + str(value)
                isFirst = False
            else:
                srr_str += ";" + key + "," + str(value)
        return srr_str

    # return the repertoire genome of the cell
    def repertoire(self):
        # initialize list of tupels
        rep_list = []
        # running over every dna sequence in the cell genome
        for seq in self.genome:
            # get srr
            if self.find_srr(seq) is None:
                srr = "No simple repeats in DNA sequence"
            else:
                srr = self.srr_to_string(self.find_srr(seq))
            # get transcribe of RNA
            rna = self.rna_poly.transcribe(seq)
            # get protein from the rna seq
            protein = self.ribosome.synthesize(rna)
            if protein is None:
                protein = "Non-coding RNA"
            # make tupel and add it to the list
            tupel = (srr, rna, protein)
            rep_list.append(tupel)
        return rep_list


# class EukaryoticCell that inherits Cell class
class EukaryoticCell(Cell):
    # standard genetic code as a dictionary
    standard_genetic_code = {
        'AUA': 'I', 'AUC': 'I', 'AUU': 'I', 'AUG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACU': 'T',
        'AAC': 'N', 'AAU': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGU': 'S', 'AGA': 'R', 'AGG': 'R',
        'CUA': 'L', 'CUC': 'L', 'CUG': 'L', 'CUU': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P',
        'CAC': 'H', 'CAU': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGU': 'R',
        'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A',
        'GAC': 'D', 'GAU': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G',
        'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S',
        'UUC': 'F', 'UUU': 'F', 'UUA': 'L', 'UUG': 'L',
        'UAC': 'Y', 'UAU': 'Y', 'UAA': None, 'UAG': None,
        'UGC': 'C', 'UGU': 'C', 'UGA': None, 'UGG': 'W'}

    # set start codons of Eukaryotic cell
    start_codons = {'AUG'}

    def __init__(self, name, genome, division_rate):
        super().__init__(name, genome, 2, self.standard_genetic_code, self.start_codons, division_rate)


# class ProkaryoticCell that inherits Cell class
class ProkaryoticCell(Cell):
    # prokaryotic genetic code as a dictionary
    prokaryotic_genetic_code = {
        'AUA': 'I', 'AUC': 'I', 'AUU': 'I', 'AUG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACU': 'T',
        'AAC': 'N', 'AAU': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGU': 'S', 'AGA': 'R', 'AGG': 'R',
        'CUA': 'L', 'CUC': 'L', 'CUG': 'L', 'CUU': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P',
        'CAC': 'H', 'CAU': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGU': 'R',
        'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A',
        'GAC': 'D', 'GAU': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G',
        'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S',
        'UUC': 'F', 'UUU': 'F', 'UUA': 'L', 'UUG': 'L',
        'UAC': 'Y', 'UAU': 'Y', 'UAA': None, 'UAG': None,
        'UGC': 'C', 'UGU': 'C', 'UGA': 'U', 'UGG': 'W'}

    # start codons of Prokaryotic cell
    start_codons = {"AUG", "GUG", "UUG"}

    def __init__(self, genome):
        super().__init__("ProkaryoticCell", genome, 1, self.prokaryotic_genetic_code, self.start_codons, 4)


# define class StemCell
class StemCell(EukaryoticCell):
    def __init__(self, genome):
        EukaryoticCell.__init__(self, "StemCell", genome, 3)


# define class NeuronCell
class NeuronCell(EukaryoticCell):
    def __init__(self, genome):
        EukaryoticCell.__init__(self, "NeuronCell", genome, 2)


# define class MutantCell
class MutantCell(StemCell):
    def __init__(self, genome, num_mutations=0, error_rate=0.05):
        StemCell.__init__(self, genome)
        self.dna_poly.error_rate = error_rate
        self.mutations_num = num_mutations

    # this function calculates the number of mutations in a given mutant cell
    def update_mutation_num(self, mutant_call):
        k = 0
        # print(self.dna_poly.error_rate)
        for seq in mutant_call.genome:  # run over each sequence in the genome and counts it's mutations
            k += (math.ceil(len(seq) * mutant_call.dna_poly.error_rate))
        return k

    def mitosis(self):
        list_of_cells = []
        mutant_cells = self * (self.division_rate - 1)  # create mutant cells with mutations
        for cell in mutant_cells:  # run over each mutant cell and update it's mutation num field
            cell.mutations_num = self.mutations_num + self.update_mutation_num(cell)
        # check for over 10 mutation in each mutant cell and replace it with cancer cell
        for i in range(len(mutant_cells)):
            if mutant_cells[i].mutations_num > 10:
                mutant_cells[i] = CancerCell(mutant_cells[i].genome, mutant_cells[i].mutations_num)
                mutant_cells[i].dna_poly.error_rate = self.dna_poly.error_rate  # keep error rate as the father's
        list_of_cells.append(self)  # create 1 regular cell with no mutations
        list_of_cells += mutant_cells  # add the mutant cells to the list of all cells
        # list_of_cells.append(self)  # create 1 regular cell with no mutations
        return list_of_cells

    # override the operator * by implement mitosis only on cells with mutations
    def __mul__(self, num):
        # create list of the mutant cells
        list_of_mutant_cells = []
        # for i in range(num):
        # create list of the mutant genome for the mutant cell
        mutant_genome = []
        # variable to save the sequence after mutation
        mutant_seq = ""
        # keep the original error rate before changing it to zero
        original_error_rate = self.dna_poly.error_rate
        for seq in self.genome:
            # transcribe the seq with mutation
            mutant_seq = self.dna_poly.transcribe(seq)
            # transcribe again but without mutation
            self.dna_poly.error_rate = 0
            mutant_genome.append(self.dna_poly.transcribe(mutant_seq))
            # resat the error_rate
            self.dna_poly.error_rate = original_error_rate
        # create mutant cell with the mutant genome
        mutant_cell = MutantCell(mutant_genome, error_rate=original_error_rate)
        # duplicate the mutant cell
        for i in range(num):
            list_of_mutant_cells.append(mutant_cell)

        # return list of the mutant cells
        return list_of_mutant_cells

    # print object as requested
    def __repr__(self):
        return "<" + "MutantCell" + ", " + str(self.num_copies) + ", " + str(self.division_rate) + ">"


# define Cancer cell
class CancerCell(MutantCell):
    def __init__(self, genome, num_mutations):
        MutantCell.__init__(self, genome, num_mutations)
        self.division_rate = 10
        self.mutations_num = num_mutations

    # print object as requested
    def __repr__(self):
        return "<" + "CancerCell" + "," + str(self.num_copies) + "," + str(self.division_rate) + ">"


class SequenceClassifier:
    def __init__(self, pattern_file):
        assert os.path.isfile(pattern_file), pattern_file  # check if file exist
        self.pattern_file = pattern_file
        self.dict = self.__patterns_to_domains(self.pattern_file)

    def prosite_to_regex(self, pattern):
        """this function gets a pattern in prosite and convert it to a regex in python.
        if it is not possiable - return None"""

        # define bad regexes for the PROSITE pattern
        bad_reg1 = '[^-x{}/\[/\]/()<>ARNDBCQEZGHILKMFPSTWYV0-9,]'  # for all not valid characyers
        bad_reg2 = '[(][A-Za-z]+[)]'  # for round brackets with letters inside
        bad_reg3 = '[^-A-Z/\[/(/\{/<][A-Za-z]'  # for missing '-' between chars
        bad_reg4 = '[A-Z][^-A-Z/\]/\}/>]'  # for missing '-' between chars
        bad_reg5 = '[/\]/\{/\}/\[/()][/\]/\{/\}/\[/]'  # for two pairs of brackets without '-' between them
        bad_reg6 = '[>].|.[<]'  # for <> not in the right place
        bad_reg7 = '[(][0-9]+[/\]/|/\}/|/\[/|/\{/|(]|[(][0-9]+[,][0-9]+[/\]/|/\}/|/\[/|/\{/|(]'  # for wrong closing () brackets
        bad_reg8 = '[/\[/][A-Z]+[/\}/|/\[/|/\{/|()]'  # for wrong closing [] brackets
        bad_reg9 = '[/\{/][A-Z]+[/\]/|/\[/|/\{/|()]'  # for wrong closing {} brackets
        bad_reg10 = '[(][)]|[/\[/][/\]/]|[/\{/][/\}/]'  # for empty brackets

        # check if the given PROTISTE pattern is valid
        if re.search(bad_reg1, pattern) or re.search(bad_reg2, pattern) or re.search(bad_reg3, pattern) or re.search(
                bad_reg4, pattern) or re.search(bad_reg5, pattern) or re.search(bad_reg6, pattern) or re.search(bad_reg7
            , pattern) or re.search(bad_reg8, pattern) or re.search(bad_reg9, pattern) or re.search(bad_reg10, pattern):
            # raise ValueError('Invalid PROSITE pattern')
            return None  # the pattern is not Valid

        # change prosite pattern to python regex
        regex = pattern
        regex = regex.replace('-', '')
        regex = regex.replace('x', '.')
        regex = regex.replace('{', '[^')
        regex = regex.replace('}', ']')
        regex = regex.replace('(', '{')
        regex = regex.replace(')', '}')
        regex = regex.replace('<', '^')
        regex = regex.replace('>', '$')
        re.compile(regex)
        return regex

    def __prosite_to_python(self, pattern_dict):
        # create new dictionary of python regexes
        regex_dict = {}
        # run over each prosite key and replace it to python regex
        for pattern in pattern_dict.keys():
            regex = self.prosite_to_regex(pattern)
            if regex is not None:  # check if pattern was valid
                regex_dict[regex] = pattern_dict[pattern]
        return regex_dict

    def __patterns_to_domains(self, pattern_file):
        # check if file exist
        assert os.path.isfile(pattern_file), pattern_file
        csv_file = open(pattern_file, 'r')  # open file
        csv_reader = csv.reader(csv_file)
        next(csv_reader)  # skip title
        patterns_dict = {}
        for line in csv_reader:  # convert each line to key and value
            patterns_dict[line[0]] = line[1]

        return self.__prosite_to_python(patterns_dict)  # convert the prosit dictionary to python regexes dictionary

    def classify(self, seq_list, csv_file):
        assert os.path.isfile(csv_file), csv_file  # check if path of output exist
        output_file = open(csv_file, 'w')  # open output file
        output_file.write("Sequence, Domains\n")  # first line - titles
        for seq in seq_list:  # run over each seq and look for domains
            domains = []
            for pattern in self.dict.keys():  # check for math patterns
                if re.search(pattern, seq) and self.dict[pattern] not in domains:  # avoid duplicate domains
                    domains.append(self.dict[pattern])
            if not domains:  # check if any domain was found
                domains.append("NA")
            output_file.write(seq + "," + ";".join(domains) + "\n")


# helper function - function to get genome sequences from txt file
def file_to_genomes(genomes_file):
    assert os.path.isfile(genomes_file), genomes_file  # check if file exist
    genomes = []
    file = open(genomes_file, 'r')
    for row in file:
        row = row[:-1]  # remove \n from the string
        genomes.append(row)  # add the gene to the list
    return genomes


# helper function - this function gets all the unique proteins in a culture and return it as a list
def get_proteins(tarbit):
    list_of_proteins = []  # list of all uniqe proteins in eac cell in the tarbit
    for cell in tarbit:  # run over each cell in the culture
        cell_repertoire = cell.repertoire()
        for repertoire in cell_repertoire:  # check if protein exist and not in the list
            if repertoire[2] != "Non-coding RNA" and repertoire[2] not in list_of_proteins:
                list_of_proteins.append(repertoire[2])  # add protein to the list
    return list_of_proteins


# helper function - find the cell with the max mutation number and return the number
def find_max_mutation(tarbit):
    max_mutations = 0
    # run over each cell (mutant or cancer) and find the maximum mutations
    for cell in tarbit:
        if cell.mutations_num > max_mutations:
            max_mutations = cell.mutations_num
    return max_mutations


def simulate(repeats, division_num, error_rate, gene_id, gene_seq):
    """this function simulates culture for given parameters and returns the results as a list"""
    # define counters for mutant and cancer cells
    mutant_cells_sum = 0
    cancer_cells_sum = 0
    proteins_sum = 0
    original_division_num = division_num  # keep the division num before substract it
    # repeat the simulate as much as requested
    for i in range(repeats):
        # mitosis
        mutant_cell = MutantCell(gene_seq, 0, error_rate)
        tarbit = [mutant_cell]
        new_tarbit = []
        while division_num > 0:
            # count number of mitosis to calculate the number of cells in the cultur
            split = 0
            # mitosis for each cell in the cultur
            for cells in tarbit:
                # check for enough space
                # if len(tarbit) - split + len(new_tarbit) + cells.division_rate - 1 <= max_cells:
                new_tarbit += cells.mitosis()
                # update mitosis counter
                split = split + 1
            # update tarbit list to be the rest of cells in the old tarbit (in case not every one undergo mitosis)
            # and add the new tarbit
            tarbit = tarbit[split:] + new_tarbit.copy()
            # initialize the new tarbit for the next run
            new_tarbit = []
            division_num = division_num - 1

        # run over the culture and count the mutant and cancer cells
        for cell in tarbit:
            if type(cell) is MutantCell:
                mutant_cells_sum += 1
            elif type(cell) is CancerCell:
                cancer_cells_sum += 1
        # count number of proteins in the culture
        proteins_sum += len(get_proteins(tarbit))
        division_num = original_division_num  # resat division number before iterating again

    # calculate averages
    mutant_cells_avg = mutant_cells_sum / repeats
    cancer_cells_avg = cancer_cells_sum / repeats
    proteins_avg = proteins_sum / repeats

    # return list of the results
    results = [error_rate, original_division_num, gene_id, cancer_cells_avg, mutant_cells_avg, proteins_avg]
    return results


def main(fasta_file):
    # check the input is a fasta file
    assert fasta_file.endswith('.fa'), fasta_file
    # create genes dictionary from given fasta file
    # key: ID of gene, Value: sequence of gene
    genes_dict = {}
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        genes_dict[seq_record.id] = str(seq_record.seq)  # create key with id and assign sequence as value

    keys = list(genes_dict.keys())
    values = list(genes_dict.values())

    # define
    repeats = 3
    division_num_list = [1, 2, 3, 4, 5]
    error_rates_list = [0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5]

    # create csv file
    file = open("exercise4_315392621_321659039.csv", 'w', newline='')
    csv_file = csv.writer(file)
    # add titles to the file
    titles = ["Error Rate", "Division Cycles", "Gene ID", "Average Cancer cells", "Average Mutant cells",
              "Average Proteins Num"]
    csv_file.writerow(titles)

    # print("File was being opened")

    max_division_list = []  # list of results of all the experiments with the maximum division cycles (5)
    # run simulates and add results to the csv file
    # counter = 1
    for error_rate in error_rates_list:  # run for every error rate
        for division_num in division_num_list:  # run for every division cycle
            for i in range(len(genes_dict)):  # run for every gene
                # write to the csv file the results of the simulation
                new_row = simulate(repeats, division_num, error_rate, keys[i], [values[i]])
                csv_file.writerow(new_row)

                # print(f'row {counter} was inserted to the file')
                # counter +=1
                # save rows of division cycle 5 for the plots
                if division_num == 5:
                    max_division_list.append(new_row)
    file.close()

    # print("Finish making the csv file")

    # create list for each gene that contains 3 lists of the results of each average
    gene_1_avgs_list = [[], [], []]
    gene_2_avgs_list = [[], [], []]
    gene_3_avgs_list = [[], [], []]
    for row in max_division_list:
        if row[2] == keys[0]:  # first gene
            # add the averages to the lists of the first gene
            gene_1_avgs_list[0].append(row[3])
            gene_1_avgs_list[1].append(row[4])
            gene_1_avgs_list[2].append(row[5])
        elif row[2] == keys[1]:  # second gene
            # add the averages to the lists of the second gene
            gene_2_avgs_list[0].append(row[3])
            gene_2_avgs_list[1].append(row[4])
            gene_2_avgs_list[2].append(row[5])
        elif row[2] == keys[2]:  # third gene
            # add the averages to the lists of the third gene
            gene_3_avgs_list[0].append(row[3])
            gene_3_avgs_list[1].append(row[4])
            gene_3_avgs_list[2].append(row[5])

    # print("Finish making genes list")

    # collect the three lists to one list
    all_genes_avgs = [gene_1_avgs_list, gene_2_avgs_list, gene_3_avgs_list]
    # create plots of seif A
    prefix_file = "exercise4_315392621_321659039"
    avg_list = ["Average Cancer cells", "Average Mutant cells", "Average Proteins Num"]
    # run on each average list of every gene
    for i, gene_id in enumerate(keys):
        for j, avg in enumerate(avg_list):
            # create plot of current average for current gene
            plt.subplot(1, 3, j + 1)
            plt.plot(error_rates_list, all_genes_avgs[i][j])
            plt.xlabel("Error Rates")
            plt.ylabel(avg)
            # save plots as png and pdf files
            plt.savefig(prefix_file + "_" + gene_id + "_seifA" + ".png")
            plt.savefig(prefix_file + "_" + gene_id + "_seifA" + ".pdf")
            plt.suptitle(gene_id)
        plt.close()

    # print("Finish making plots for seif A")

    # lists of all proteins averages in order of division cycles
    # 3 lists for 3 genes, each list containes 5 list for averages in each division cycle
    all_prot_avg = [[[], [], [], [], []],
                    [[], [], [], [], []],
                    [[], [], [], [], []]]
    # get data from the csv file that was generated before
    file2 = open('exercise4_315392621_321659039.csv', 'r')
    csv_reader = csv.reader(file2)
    next(csv_reader)
    # run on every row in file and collect data about protein average
    for row in csv_reader:
        current_gene = row[2]
        division_cycle = int(row[1])  # for use as index
        avg_prot = float(row[5])
        # assign each gene's data to it's appropriate list
        if current_gene == keys[0]:
            all_prot_avg[0][division_cycle - 1].append(avg_prot)
        elif current_gene == keys[1]:
            all_prot_avg[1][division_cycle - 1].append(avg_prot)
        elif current_gene == keys[2]:
            all_prot_avg[2][division_cycle - 1].append(avg_prot)
    file2.close()
    # print("Finish making proteins list")

    # run on the lists of all gene's proteins averages and make multi-bars plots
    for gene_prot, gene_name in zip(all_prot_avg, keys):
        data = gene_prot  # assign current data
        x = np.arange(len(error_rates_list))
        width = 0.25  # width of each bar
        ax = plt.subplot()
        # define bars - 5 bars for 5 division cycles
        ax.bar(x - width, data[0], width, label='Division cycle #1')
        ax.bar(x - width / 2, data[1], width, label='Division cycle #2')
        ax.bar(x, data[2], width, label='Division cycle #3')
        ax.bar(x + width / 2, data[3], width, label='Division cycle #4')
        ax.bar(x + width, data[4], width, label='Division cycle #5')
        # titles
        ax.set_ylabel('Average Proteins Num')
        ax.set_xlabel('Error Rates')
        ax.set_title(gene_name)
        # x axies
        ax.set_xticks(x)
        ax.set_xticklabels(error_rates_list)
        ax.legend()
        # save plots as png and pdf files
        plt.savefig(prefix_file + '_' + gene_name + '_seifB' + '.pdf')
        plt.savefig(prefix_file + '_' + gene_name + '_seifB' + '.png')
        plt.close()
        # print("Finish making plot - Seif B")

    # print("Finish all plots for seif B")


if __name__ == '__main__':
    main(sys.argv[1])

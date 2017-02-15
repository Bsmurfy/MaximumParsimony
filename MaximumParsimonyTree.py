'''
Brendan Murphy

A program to find a consensus tree for a set of four sequences,
using bootstrapping and maximum parsimony.
'''

import random

def main():
    # read information from input file
    sequenceinfo = readsequences('FinalProblemCS4.txt')

    sequencelist = []
    labellist = []
    # make a list of the sequences in the file
    # and a separate list of their labels
    for i in range(0, len(sequenceinfo)): 
        sequence = sequenceinfo[i][1]
        sequencelist.append(sequence)
        label = sequenceinfo[i][0]
        labellist.append(label)

    # trim the sequences down to only informative sites
    sequencelist = showinform(sequencelist)
    n = len(sequencelist[0]) # the number of informative sites

    # for a set of four sequences A, B, C, and D...
    cladepair1 = 0 # counts the number of A-B clades (and thus also C-D clades)
    cladepair2 = 0 # counts the number of A-C clades (and thus also B-D clades)
    cladepair3 = 0 # counts the number of A-D clades (and thus also B-C clades)

    # run bootstrapping 1000 times
    numtrials = 1000
    for trial in range(numtrials):
        sampleseqlist = []
        for seq in sequencelist: # for each sequence...
            sample = ''
            for x in range(n): # and for each informative site...
                site = random.randrange(n) # sample with replacement
                sample += seq[site]
            sampleseqlist.append(sample) # make a list of sequences after replication

        initialtree = build_random_tree(sampleseqlist) # choose a random start tree
        length = computeL(initialtree) # find the total branchlength
        besttree = initialtree

        # count the number of times trees are compared with no change to the best tree
        nochange = 0
        # 10 is an arbitrarily chosen threshhold, but should reflect the best tree fine
        while (nochange < 10):
            treetocompare = build_random_tree(sampleseqlist) # build a new random tree...
            lengthtocompare = computeL(treetocompare) # find its branchlength...
            if (lengthtocompare < length): # if this is the best branchlength so far,
                length = lengthtocompare # make it the new value to be compared with
                besttree = treetocompare # save the tree as the best candidate so far
            else:
                nochange += 1

        # print topologies for the first three trials 
        if (trial == 0 or trial == 1 or trial == 2):
            # treeprint is used to rid the tree of tags before printing
            treeprint = ((besttree[0][0][0], besttree[0][1][0]),
                         (besttree[1][0][0], besttree[1][1][0]))
            # change the prefix for printing according to the trial
            if (trial == 0):
                prefix = 'first' 
            elif (trial == 1):
                prefix = 'second'
            elif (trial == 2):
                prefix = 'third'
            print('The ' + prefix + ' bootstrap replication yielded a topology of ')
            print(treeprint)
            print('') # for reader-friendliness

        # determine the tags of the sequences in the first clade
        # (i.e. 0 and 1, representing a clade of the first sequence and the second sequence)
        tag1 = besttree[0][0][1]
        tag2 = besttree[0][1][1]
        # increment clade counters when the clades are observed in a sample topology
        if ((tag1 == 0 and tag2 == 1) or
            (tag1 == 1 and tag2 == 0)):
            cladepair1 += 1
        elif ((tag1 == 0 and tag2 == 2) or
              (tag1 == 2 and tag2 == 0)):
            cladepair2 += 1
        elif ((tag1 == 0 and tag2 == 3) or
              (tag1 == 3 and tag2 == 0)):
            cladepair3 += 1 
        elif ((tag1 == 1 and tag2 == 2) or
              (tag1 == 2 and tag2 == 1)):
            cladepair3 += 1
        elif ((tag1 == 1 and tag2 == 3) or
              (tag1 == 3 and tag2 == 1)):
            cladepair2 += 1 
        elif ((tag1 == 2 and tag2 == 3) or
              (tag1 == 3 and tag2 == 2)):
            cladepair1 += 1

    # determine which clade pair appears the most
    consensus = max(cladepair1, cladepair2, cladepair3)
    # consensustree has tags removed for printing
    # clade counters are divided by total number of trials to obtain bootstrap value
    if (consensus == cladepair1):
        consensustree = [[labellist[0], labellist[1]],
                         [labellist[2], labellist[3]]]
        bootstrap = cladepair1/numtrials
    elif (consensus == cladepair2):
        consensustree = [[labellist[0], labellist[2]],
                         [labellist[1], labellist[3]]]
        bootstrap = cladepair2/numtrials
    elif (consensus == cladepair3):
        consensustree = [[labellist[0], labellist[3]],
                         [labellist[1], labellist[2]]]
        bootstrap = cladepair3/numtrials
    # print the consensus tree and bootstrap value in a reader-friendly format
    print('The consensus tree for the given sequences was found to be: ')
    print('') 
    print(str(consensustree[0][0]) + '         ' + str(consensustree[1][0]))
    print('                 \________/') # nice output for reader
    print('                 /        \\') # extra \ needed to print \
    print(str(consensustree[0][1]) + '         ' + str(consensustree[1][1]))
    print('')
    print('with a bootstrap value of ' + str(bootstrap) + '.')
    
'''
showinform - A function that checks for informative sites in a list of sequences
and trims the sequences to only show informative sites.
Parameter: seqlist - a list of sequences to be evaluated
Return: 
'''
def showinform(seqlist):

    seqlength = len(seqlist[0])   
    informtracker = [] # tracks the locations of informative sites in the sequences
    for x in range(seqlength): # for each index of the sequences...
        countA = 0
        countC = 0
        countG = 0
        countT = 0
        # iterate over the sequences and count the number of appearances for each nucleotide
        for seq in seqlist: 
            if (seq[x] == 'A'):
                countA += 1
            elif (seq[x] == 'C'):
                countC += 1
            elif (seq[x] == 'G'):
                countG += 1
            elif (seq[x] == 'T'):
                countT += 1
        # an informative site has two or more appearances of at least two different nucleotides
        if ((countA >= 2 and countC >= 2) or
            (countA >= 2 and countG >= 2) or
            (countA >= 2 and countT >= 2) or
            (countC >= 2 and countG >= 2) or
            (countC >= 2 and countT >= 2) or
            (countG >= 2 and countT >= 2)):
            informtracker.append(1) # represent informative sites with a 1 in the tracker
        else:
            informtracker.append(0) # represent non-informative sites with a 0

    informseqlist = []
    # make a list of the new, trimmed down sequences
    for seq in seqlist:
        informseq = ''
        for x in range(seqlength):
            if (informtracker[x] == 1): 
                informseq += seq[x]
        informseqlist.append(informseq)

    return informseqlist    

'''
build_random_tree - A function the builds a random tree from a list of four sequences
Parameter: seqlist - a list of four sequences (only informative sites)
Return: tree - a four-sequence tree, randomly generated
'''
def build_random_tree(seqlist):
    tree = []
    seqnum = 0
    toadd = []
    # assign a tag to each of the sequences to help determine clades later
    for seq in seqlist:
        tag = seqnum
        seqnum += 1
        toadd.append([seq, tag])
    # choose a random order for the sequences 
    a = random.choice(toadd)
    toadd.remove(a)
    b = random.choice(toadd)
    toadd.remove(b)
    c = random.choice(toadd)
    toadd.remove(c)
    d = random.choice(toadd)
    # create a tree reflecting the random order
    tree.append([a,b])
    tree.append([c,d])
    return tree

'''
computeL - A function that computes the total branchlength of a given tree.
Parameter: tree - a tree of four sequences (and their tags)
Return: l - the total branchlength of the tree
'''
def computeL(tree):

    # only look at the sequences, not the tags
    a = tree[0][0][0]
    b = tree[0][1][0]
    c = tree[1][0][0]
    d = tree[1][1][0]

    # length is determined by number of substitutions between two given seqeunces
    ab_length = countsubs(a, b)
    cd_length = countsubs(c, d)
    ac_length = countsubs(a, c)
    bd_length = countsubs(b, d)
    # the total branchlength needs to be calculated from its overlapping parts
    l = (ac_length + bd_length + ab_length + cd_length)/2
    
    return l

'''
countsubs - A function that counts the number of substitutions between two strings
(and thus, the branchlength for the two sequences)
Parameter: seq1 - the first sequence
Parameter: seq2 - the second sequence
Return: subcount - the branchlength, represented as the number of substitutions between the two strings
'''
def countsubs(seq1, seq2):
    subcount = 0
    for x in range(len(seq1)):
        # a substitution occurs whenever the nucleotides at an index do not match
        if (seq1[x] != seq2[x]):
            subcount += 1
    return subcount

'''
readsequences - A function that reads files with a format similar to FASTA.
Parameter: filename - the file to be read as input 
Return: resultList - a list with elements of the form [label, sequence] from the input file
'''
def readsequences(filename):
    resultList = []
    infile = open(filename, 'r')

    line = infile.readline()
    header = line.rstrip() # remove excess white space from end
    label = header[1:] # remove the first character (>)

    sequence = ''

    for line in infile:
        line = line.rstrip()

        # ignore blank lines
        if line != '':

            if line[0] == '>':
                resultList.append([label, sequence])
                header = line.rstrip()
                label = header[1:]
                sequence = ''

            else:
                sequence += line

    infile.close()
    resultList.append([label, sequence])
    return resultList

main()

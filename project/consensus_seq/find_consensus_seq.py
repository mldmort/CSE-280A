from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np

tran_map = {'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C',
            'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C',
            'TTA': 'L', 'TCA': 'S', 'TAA': '|', 'TGA': '|',
            'TTG': 'L', 'TCG': 'S', 'TAG': '|', 'TGG': 'W',

            'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R',
            'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',
            'CTA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',
            'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',

            'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S',
            'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',
            'ATA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',
            'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',

            'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G',
            'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',
            'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',
            'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'}

def translate(seq,start_vntr):
  start_tran = start_vntr%3
  protein = ''

  for prot_ind in range(start_tran,len(seq)-3+int(len(seq)%3==0),3):
    protein += tran_map[seq[prot_ind:prot_ind+3]]
  #print('start_vntr: ', start_vntr)
  #print('start_tran: ', start_tran)
  #print('rnage: ', range(start_tran,len(seq)-3,3))
  #print('len(protein): ', len(protein))
  #print('int((len(seq)-start_tran)/3): ', int((len(seq)-start_tran)/3))
  assert(len(protein)==int((len(seq)-start_tran)/3))
  return protein

def translate_s(seq,start):
  protein = ''
  for prot_ind in range(start,len(seq)-3+int(len(seq)%3==0),3):
    protein += tran_map[seq[prot_ind:prot_ind+3]]
  #print('start_vntr: ', start_vntr)
  #print('start_tran: ', start_tran)
  #print('rnage: ', range(start_tran,len(seq)-3,3))
  #print('len(protein): ', len(protein))
  #print('int((len(seq)-start_tran)/3): ', int((len(seq)-start_tran)/3))
  #assert(len(protein)==int((len(seq)-start_tran)/3))
  return protein

def translate_m_vntr(seq,start,vntr_start,vntr_len):
  protein = ''
  new_seq = seq[:vntr_start] + seq[vntr_start+vntr_len:]
  for prot_ind in range(start,len(new_seq)-3+int(len(new_seq)%3==0),3):
    protein += tran_map[new_seq[prot_ind:prot_ind+3]]
  #print('start_vntr: ', start_vntr)
  #print('start_tran: ', start_tran)
  #print('rnage: ', range(start_tran,len(seq)-3,3))
  #print('len(protein): ', len(protein))
  #print('int((len(seq)-start_tran)/3): ', int((len(seq)-start_tran)/3))
  #assert(len(protein)==int((len(seq)-start_tran)/3))
  return protein

cDNA_map = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}

def make_cDNA(seq):
  cDNA_seq = ''
  for let in seq[(len(seq)-1):0:-1]:
    cDNA_seq += cDNA_map[let]
  cDNA_seq += cDNA_map[seq[0]]
  assert(len(cDNA_seq)==len(seq))
  return cDNA_seq
 
def getStartVNTRList(file_name,verbose=False):
  print("reading state file: "+file_name)
  start_ind_list = []
  with open(file_name,"r") as fh:
    lines = fh.readlines()
    counter = 0
    while (counter<len(lines)):
      if (verbose):
        print(lines[counter])
      line = lines[counter].strip().split()
      read_num = line[2]
      read_len = int(line[4])
      counter += 1
      read_states = lines[counter:counter+read_len]
      counter += read_len
      count_nonMatch = 0
      for i, item in enumerate(read_states):
        word = item.strip().split("_")
        if (word[0][0] != "M"):
          count_nonMatch += 1
        if ( (word[0]=="unit") and (word[1]=="start") ):
          unit_num = int(word[2])
          # count_nonMatch caracters are before unit_start_* which are not matches and should be subtracted from i,
          # also the first VNTR match is 1 index after unit_start_*
          start_ind_list.append(i-unit_num*12-count_nonMatch+1)
          break
  if (verbose):
    print("read "+str(len(start_ind_list))+" starts!")
  return start_ind_list

def getMajorityVote(dic):
  '''
  gets a dictionary and return a string
  dic should have A, C, G, T, -  as keys and a list of qualities as values
  '''
  err = {'A':.5, 'C':.5, 'G':.5, 'T':.5}
  for let in ['A','C','G','T']:
    if (len(dic[let])>0):
      err[let] = 10.**(-0.1*np.mean(np.array(dic[let])))
  Prob = {'A': (1.-err['A'])**(len(dic['A'])) * err['C']**(len(dic['C'])) * err['G']**(len(dic['G'])) * err['T']**(len(dic['T'])), 
          'C': (1.-err['C'])**(len(dic['C'])) * err['A']**(len(dic['A'])) * err['G']**(len(dic['G'])) * err['T']**(len(dic['T'])), 
          'G': (1.-err['G'])**(len(dic['G'])) * err['A']**(len(dic['A'])) * err['C']**(len(dic['C'])) * err['T']**(len(dic['T'])), 
          'T': (1.-err['T'])**(len(dic['T'])) * err['A']**(len(dic['A'])) * err['C']**(len(dic['C'])) * err['G']**(len(dic['G']))}
  major = max(Prob,key=lambda k: Prob[k])
  seq_error = '-'
  for let in (set(['A','C','G','T']).difference(set([major]))):
    if (len(dic[let]) > 0):
      seq_error = let
  return major, seq_error

def getCoverage(dic):
  '''
  gets a dictionary and returns a float
  dic should have A, C, G, T, -  as keys
  '''
  nbp = len(dic['A']) + len(dic['C']) + len(dic['G']) + len(dic['T'])
  return float(nbp)
  

if __name__ == "__main__":
  working_dir = "./"

  fq_files=["class_i_seqs_30x_MSv3.fq","class_i_seqs_20x_MSv3.fq","class_i_seqs_10x_MSv3.fq","class_i_seqs_5x_MSv3.fq",
            "class_j_seqs_30x_MSv3.fq","class_j_seqs_20x_MSv3.fq","class_j_seqs_10x_MSv3.fq","class_j_seqs_5x_MSv3.fq"]

  state_files=["class_i_states_30x_MSv3.dat","class_i_states_20x_MSv3.dat","class_i_states_10x_MSv3.dat","class_i_states_5x_MSv3.dat",
               "class_j_states_30x_MSv3.dat","class_j_states_20x_MSv3.dat","class_j_states_10x_MSv3.dat","class_j_states_5x_MSv3.dat"]

  RU_count=[4,4,4,4,6,6,6,6]

  labels = ['30x class i','20x class i','10x class i','5x class i',
            '30x class j','20x class j','10x class j','5x class j']

  #fq_files = ["class_i_seqs_20x_MSv3_indel.fq","class_j_seqs_20x_MSv3_indel.fq"]

  #state_files = ["class_i_states_20x_MSv3_indel.dat","class_j_states_20x_MSv3_indel.dat"]

  #fq_files = ["class_i_seqs_20x_MSv3_swidel.fq","class_j_seqs_20x_MSv3_swidel.fq"]

  #state_files = ["class_i_states_20x_MSv3_swidel.dat","class_j_states_20x_MSv3_swidel.dat"]

  #working_dir = "/Users/miladmortazavi/Documents/Dev/BIOINF/project/working_dir3/"

  #fq_files = ["class_i_seqs_30x_GP1BA.fq", "class_j_seqs_30x_GP1BA.fq"]

  #state_files = ["class_i_states_30x_GP1BA.dat", "class_j_states_30x_GP1BA.dat"]

  #RU_count=[4,6]

  #labels = ['30x class i','30x class j']

  read_length = 150
  #read_length = 250
  RU_length = 12
  #RU_length = 39

  for ind, item in enumerate(fq_files):
    fq_files[ind] = working_dir+item
  for ind, item in enumerate(state_files):
    state_files[ind] = working_dir+item

  ######## process class files ########
  #plt.figure()
  for fl_ind, fl in enumerate(state_files):
    print()
    print("+"*250)
    start_list = getStartVNTRList(fl)
    print("start loc: ", start_list)

    records = list(SeqIO.parse(fq_files[fl_ind],"fastq"))

    #plt.plot()
    #for i, record in enumerate(records[:5]):
    #  plt.plot(record.letter_annotations['phred_quality'])
    #  print("record: ", i)
    #  print("annot: ",record.letter_annotations)
    #plt.show()
    #break
    ##print("start list len: ", len(start_list))
    #print("records len: ", len(records))

    ind_min = 1e20
    ind_max = -1e20
    for start in start_list:
      ind_min = min(ind_min,start)
      ind_max = max(ind_max,start)

    record_collection_length = ind_max + (read_length - ind_min)
    print('sequence length: ', record_collection_length)
    record_collection = ['-'*record_collection_length]*len(records)
    record_collection_qual = [[0]*record_collection_length]
    for r_ind in range(len(records)-1):
      record_collection_qual.append([0]*record_collection_length)

    for r_ind, record in enumerate(records):
      record_collection[r_ind]  = record_collection[r_ind][:ind_max-start_list[r_ind]] + record.seq + record_collection[r_ind][ind_max-start_list[r_ind]+read_length:]
      record_collection_qual[r_ind][(ind_max-start_list[r_ind]):(ind_max-start_list[r_ind]+read_length)]  = record.letter_annotations['phred_quality']
      #print(str(r_ind)+": \t", record.seq[start_list[r_ind]:start_list[r_ind]+RU_length*i_repeat])

    #for r_ind, record_qual in enumerate(record_collection_qual):
    #  print(str(r_ind)+": \t",record_qual[9])

    for r_ind, record in enumerate(record_collection):
      print(str(r_ind)+": \t", record)

    # find union sequence
    sequence = ''
    sequence_err = ''
    sequence_cov = []
    sequence_vote = []
    for seq_ind in range(record_collection_length):
      dic = {'A':[], 'C':[], 'G':[], 'T':[], '-':[]}
      for r_ind, record in enumerate(record_collection):
        dic[record[seq_ind]].append(record_collection_qual[r_ind][seq_ind])
      sequence_vote.append(dic)
      
    for seq in sequence_vote:
      major, error = getMajorityVote(seq)
      sequence += major
      sequence_err += error
      sequence_cov.append(getCoverage(seq))

    mean_coverage = 0.0
    for cov in sequence_cov:
      mean_coverage += float(cov)
    mean_coverage /= len(sequence_cov)

    # translate mRNA seq to protein
    #protein = translate(sequence,ind_max)

    #new_sequence = sequence[:ind_max] + sequence[ind_max+RU_count[fl_ind]*RU_length+1:]
    #protein_0 = translate_s(new_sequence,0)
    #protein_1 = translate_s(new_sequence,1)
    #protein_2 = translate_s(new_sequence,2)

    #new_sequence_cDNA = make_cDNA(new_sequence)
    #protein_cDNA_0 = translate_s(new_sequence_cDNA,0)
    #protein_cDNA_1 = translate_s(new_sequence_cDNA,1)
    #protein_cDNA_2 = translate_s(new_sequence_cDNA,2)

    protein_0 = translate_s(sequence,0)
    protein_1 = translate_s(sequence,1)
    protein_2 = translate_s(sequence,2)

    cDNA_seq = make_cDNA(sequence)
    #protein_cDNA = translate(cDNA_seq,len(cDNA_seq)-(ind_max+RU_count[fl_ind]*RU_length-1))
    protein_cDNA_0 = translate_s(cDNA_seq,0)
    protein_cDNA_1 = translate_s(cDNA_seq,1)
    protein_cDNA_2 = translate_s(cDNA_seq,2)

    print('err: \t', sequence_err)
    print('seq: \t', sequence)
    print('mean coverage: ', mean_coverage)
    #print('coverage: ', sequence_cov)
    print('sequence length: ', len(sequence))
    print('sequence error: \t', sequence_err)
    print('sequence: \t\t', sequence)
    #print('new sequence: \t\t', new_sequence)
    #print('protein: \t\t', protein)
    print()
    print('protein 0: \t\t', protein_0)
    print('protein 1: \t\t', protein_1)
    print('protein 2: \t\t', protein_2)
    print()
    #print('cDNA seq: \t\t', cDNA_seq)
    #print('protein cDNA: \t\t', protein_cDNA)
    print('protein cDNA 0: \t', protein_cDNA_0)
    print('protein cDNA 1: \t', protein_cDNA_1)
    print('protein cDNA 2: \t', protein_cDNA_2)

    #if (fl_ind <4):
    #  plt.plot(sequence_cov,'o', markersize=3,label=labels[fl_ind])

    print("+"*250)
    print()
  #for loc in spades_loc_diff:
  #  plt.plot([loc, loc], [4,6], linewidth=1, color='black')
  #plt.xlabel('Locus')
  #plt.ylabel('Coverage')
  #plt.title('Coverage vs loci')
  #plt.legend()
  #plt.savefig('Coverage_vs_loci.png')
  ######################################

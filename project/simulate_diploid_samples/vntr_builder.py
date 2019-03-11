import argparse
import os.path
import numpy as np

def read_core(fn):
  '''
    reads a file with file name fn and returns a list of the sequece in the first file.
    only the first line is read.
  '''
  if (not os.path.isfile(fn)):
    print('file '+fn+' does not exist!')
    return []
  with open(fn) as fh:
    seq = fh.readline()
    seq_list = [seq[i] for i in range(len(seq))]
  
  if (seq_list[-1] == '\n'):
    seq_list.pop(-1)

  for l in seq_list:
    if (l not in ['A','T','C','G']):
      print('unrecognized letter in input file: '+fn+'. Check input file!')

  return seq_list

def writeFile(fn,ar):
  with open(fn,'w') as fh:
    for i in ar:
      fh.write(i)

def let2num(seq):
  dic = {'A':0, 'T':1, 'C':2, 'G':3, 'N':4}
  pass

def list2str(lst):
  return ''.join(lst)

parser = argparse.ArgumentParser(description='Generate a short read consisting a vntr with a RU-count')
parser.add_argument('-i','--input',type=str,help='input file containing vntr seq in the first line',default='NO_FILE_SPECIFIED')
parser.add_argument('-o','--output',type=str,help='output file containing the desired read',default='output.dat')
parser.add_argument('-c','--count',type=int,help='RU-count desired',default=1)
parser.add_argument('-l','--length',type=int,help='length of read desired',default=100)
args = parser.parse_args()

print('args:', args)
print('input: ',args.input)
print('count: ',args.count)
print('length: ',args.length)

if __name__ == "__main__":

  # get arguments
  fileName = args.input
  count = args.count
  read_len = args.length
  outFileName = args.output

  # read a list of the sequence from the file
  seq_list_orig = read_core(fileName)

  # make the VNTR as a list
  seq_list = seq_list_orig*count
  seq_len = len(seq_list)

  if (seq_len > read_len):
    print('sequence length is larger than the read length. Aborting...')
    print('original seq length: ', len(seq_list_orig))
    print('VNTR length: ', seq_len)
    print('read length: ', read_len)
    raise(0)
  
  random_len = read_len - seq_len
  random_seq = np.random.choice(['A','T','C','G'],size=random_len,replace=True)
  
  # the location where the VNTR starts from
  random_pos = np.random.choice(random_len+1,size=1)[0]

  read = [random_seq[i] for i in range(random_pos)]
  read += seq_list
  for i in range(random_pos,random_len):
    read += random_seq[i]

  writeFile(outFileName, read)

  print('read length: ', len(read))
  print('random position: ', random_pos)

  print('short seq: ', list2str(seq_list_orig))
  print('VNTR seq: ', list2str(seq_list))
  print('VNTR from read: ', list2str(read[random_pos:(random_pos+seq_len)]))
  print('total read: ', list2str(read))









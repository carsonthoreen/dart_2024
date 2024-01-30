from Bio import SeqIO
import numpy as np
import sys

def process_seq(sequence):
	AUG_location = sequence.find('ATG')
	if AUG_location >= 10:
		Seq_flank_AUG = sequence[AUG_location-10:AUG_location+13]
		Kozak_score = similarity_score(Seq_flank_AUG)
	else:
		Kozak_score = -1	
	return(Kozak_score)

#Kozak Consensus Scoring System

#0=A, 1=T, 2=G, 3=C, 4=N (Missing)
weights = np.array([
       [0.04210526, 0.        , 0.03157895, 0.05263158, 0.        ],
       [0.04210526, 0.05263158, 0.10526316, 0.0625    , 0.        ],
       [0.03157895, 0.04210526, 0.05263158, 0.07368421, 0.        ],
       [0.03157895, 0.01052632, 0.04210526, 0.05263158, 0.        ],
       [0.08421053, 0.07368421, 0.18947368, 0.10526316, 0.        ],
       [0.04210526, 0.05263158, 0.05263158, 0.08421053, 0.        ],
       [0.12631579, 0.0625    , 0.12631579, 0.21052632, 0.        ],
       [0.83157895, 0.12631579, 0.65263158, 0.16842105, 0.        ],
       [0.15789474, 0.06315789, 0.11578947, 0.2       , 0.        ],
       [0.21052632, 0.09473684, 0.31578947, 0.51578947, 0.        ],
       [0.        , 0.        , 0.        , 0.        , 0.        ],
       [0.        , 0.        , 0.        , 0.        , 0.        ],
       [0.        , 0.        , 0.        , 0.        , 0.        ],
       [0.24210526, 0.16666667, 0.53684211, 0.13684211, 0.        ],
       [0.15789474, 0.09473684, 0.09473684, 0.24210526, 0.        ],
       [0.05263158, 0.08421053, 0.14736842, 0.09473684, 0.        ],
       [0.07216495, 0.05263158, 0.10526316, 0.06315789, 0.        ],
       [0.        , 0.        , 0.        , 0.05263158, 0.        ],
       [0.05263158, 0.05263158, 0.10526316, 0.09473684, 0.        ],
       [0.04210526, 0.03157895, 0.05263158, 0.04210526, 0.        ],
       [0.        , 0.        , 0.        , 0.        , 0.        ],
       [0.04210526, 0.04210526, 0.08421053, 0.07368421, 0.        ],
       [0.0625    , 0.04210526, 0.09473684, 0.05263158, 0.        ]
])
#Below function scores using consensus kozak motif scores
def similarity_score(sequence):
    
    assert len(sequence)==23,'Sequence must be 23 bases long. Codon of interest must be centered, with 10 bases flanking both sides.'

            
    #We need consistency and flexibility:
    sequence = sequence.upper()
    for i in np.arange(len(sequence)):
        if sequence[i] =='U':
            sequence = sequence[0:i]+'T'+sequence[i+1:len(sequence)]
    
    numbers=[0]*len(sequence)
    
    for k in np.arange(len(sequence)):
        if sequence[k]=='A':
            numbers[k] = 0
        elif sequence[k]=='T':
            numbers[k] = 1
        elif sequence[k]=='G':
            numbers[k] = 2
        elif sequence[k]=='C':
            numbers[k] = 3
        else:
            numbers[k]=4
                
    
    score = 0
    for k in np.arange(len(numbers)):
        score += weights[k][numbers[k]]
            
    max_score = np.sum(weights.max(axis=1))
    
    score = score/max_score
    
    #Final scoring value: we take the maximum possible score 
    #calculated, and return our score divided by the maximum (to normalize from range 0 to 1) 
    
    return(score)

if __name__ == '__main__':
	FastaFile = sys.argv[1]
    Output = sys.argv[2]
	ID = []
	Kozak_score = []
	with open(FastaFile) as handle:
		for i in SeqIO.parse(handle, "fasta"):
			ID.append(i.id)
			Kozak_score.append(process_seq(str(i.seq)))

	Outfile = open(Output, "w")
	for index in range(len(ID)):
    	Outfile.write(str(ID[index]) + "\t" + str(Kozak_score[index]) + "\n")
	Outfile.close()


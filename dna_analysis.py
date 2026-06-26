# importing the moodule for connectiong a secure network to connect to outdide server
import ssl
ssl._create_default_https_context = ssl._create_unverified_context

# importing all the necessary modules
from Bio import Entrez , SeqIO
import matplotlib.pyplot as my_plot

# time delay effect
import time
def time_print(text , delay = 0.02):
    for char in text:
        print(char, end= '', flush = True)
        time.sleep(delay)
    print()    

# setting up email to connect to NCBI
Entrez.email = "ss4541801@gmail.com"

# taking the acession id as input
accn_id_user_input = input("Enter the acession number: ")

# Extracting data from NCBI
handle = Entrez.efetch(db = "nucleotide",
              id = accn_id_user_input,
              rettype  = "fasta",
              retmode = "text"
              )

# Saving the extracted data
record = SeqIO.read(handle, "fasta")
SeqIO.write(record, "DNA.fasta", "fasta")
handle.close() 


time_print(f"\t Connecting to NCBI database........................")
print(f"\t Connection SUCESSFULL :)")
time_print(f"\t Logging into the NCBI database through provided email....................")
print(f"\t Enttered into the database.")
time_print(f"\t Fetching the data for requested gene({accn_id_user_input})....................")
print(f"\t Data retrived SUCESSFULL :) ")
#------------------------------------------------------------------------------------------------------------------------
# calculating the gc contennt
gc_count = record.seq.count("G")  + record.seq.count("C")
gc_content = (gc_count / len(record.seq)) * 100
time_print("\t \t \t \t \t DNA extraction SUCESSFULL :)")
print("\t FASTA file generated in your current folder ")
time_print(f"\t \t \t \t \t GC count generated SUCESSFULLY :) ")
print(f"\t GC content--> {gc_content:.2f}%")

#--------------------------------------------------------------------------------------------------------------------------
# creating a bar chart to show  amount of nucleotides in the sequence 
a_count = record.seq.count("A")
t_count = record.seq.count("T")
g_count = record.seq.count("G")
c_count = record.seq.count("C")

nucleotides = ['A', 'T', 'G', 'C']
count = [a_count, t_count, g_count, c_count]

my_plot.bar(nucleotides, count, color = ['red', 'blue', 'green', 'pink'])
my_plot.title("Nucleotide composition of  BRACA1 mRNA")
my_plot.xlabel("Nucleotide")
my_plot.ylabel("Count")
graph_name = my_plot.savefig("nucleotide_composition.png") #  saving the plot 
time_print(f"\t \t \t \t \t BAR chart for nucleotides created SUCESSFULLY :)")
print(f"\t BAR chart file generated in your current folder")

#----------------------------------------------------------------------------------------------------------------------------------
# now finding the ORF(open reading frame(for all the three frames))
def ORF_finder(sequence, start_frame):

    all_orf = []
    current_orf = ''
    recording = False
    for start in range(start_frame, len(sequence), 3):
        codon = sequence[start : start + 3]
        if codon == "ATG" and recording == False:
            recording = True
            current_orf = ''
        if recording == True:
            current_orf = current_orf + codon
        if (codon == 'TAA' or codon == 'TAG' or codon == 'TGA') and recording == True:
            all_orf.append(current_orf)
            recording = False
            current_orf = ''
    return all_orf        

sequence = record.seq  # extracted dna sequence is renamed as "sequence"

ORF_list_at_posn_1 = ORF_finder(sequence, 0)
ORF_list_at_posn_2 = ORF_finder(sequence, 1)
ORF_list_at_posn_3 = ORF_finder(sequence, 2)

final_ORF_list = ORF_list_at_posn_1 + ORF_list_at_posn_2 + ORF_list_at_posn_3
    
    
longest_orf = max(final_ORF_list, key = len)
time_print(f"\t \t \t \t \t Codon generation over all three coading frames completed SUCESSSFULLY :)")
with open("Forward_ORF.txt", 'w')as f:
    f.write(f"ORF present in forward direction \n {final_ORF_list}")
print(f"\t Forward ORF file gennerated in current folder ")
print(f"\t Largest ORF in forward direction is --> {len(longest_orf)}")                                           #    NM_007294

#----------------------------------------------------------------------------------------------------------------
# converting all this ORF into a protien sequence (translation)
protien_seq = longest_orf.translate(to_stop = True)
time_print(f"\t \t \t  \t \t Translation completed SUCESSFULLY :) ")
# print(f"\t This is protien sequence --> {protien_seq}")
print(f"\t Length of protien sequence --> {len(protien_seq)} ")
with open("Translation_seq.txt", "w") as f:
    f.write(str(protien_seq))
print(f"\t Protien sequence file generated in your current folder")    


#----------------------------------------------------------------------------------------------------------------------
# Reverse complament
reverse_complment_seq = sequence.reverse_complement()
time_print(f"\t \t \t \t \t Reverse complament completed SUCESSSFULLY :) ")
with open("Reverse_compl_seq.txt", "w") as f:
    f.write(str(reverse_complment_seq))
print("\t Reverse complemennt sequence file generated in your current folder.")    

ORF_of_reversed_seq_1 = ORF_finder(reverse_complment_seq, 0)                         # making ORF of reversed_complamet_seq
ORF_of_reversed_seq_2 = ORF_finder(reverse_complment_seq, 1)
ORF_of_reversed_seq_3 = ORF_finder(reverse_complment_seq, 2)
reversed_complament_ORF_list = ORF_of_reversed_seq_1 + ORF_of_reversed_seq_2 + ORF_of_reversed_seq_3
with open("Backward_ORF.txt", "w") as f:
    f.write(str(reversed_complament_ORF_list))
print(f"\t Largest ORF in reverse complement seq --> {len(max(reversed_complament_ORF_list, key= len))} ")


complete_both_side_ORF_list = final_ORF_list + reversed_complament_ORF_list            # final ORF list for both sidr, strainght and reversed

final_largest_ORF = max(complete_both_side_ORF_list, key = len)
print(f"\t Largest ORF across whole genome is --> {len(final_largest_ORF)}")









    











    


    
    


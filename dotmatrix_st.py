import streamlit as st
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from io import StringIO

def dotmatrix(file1,file2,window=10,width=500):
    seq1={}
    for rec in SeqIO.parse(file1,"fasta"):
        seq1len=len(rec.seq)        
        for i in range(seq1len-window):
            subseq=rec.seq[i:i+window].upper()
            if subseq not in seq1:
                seq1[subseq]=[]
            seq1[subseq].append(i)

        break

    for rec in SeqIO.parse(file2,"fasta"):
        height=int(len(rec.seq)*width/seq1len)
        image=np.zeros((height,width),dtype=np.uint8)
        seq2len=len(rec.seq)
        for j in range(seq2len-window):
            subseq=rec.seq[j:j+window].upper()
            if subseq in seq1:
                y=int(j/seq2len*height)
                for i in seq1[subseq]:
                    x=int(i/seq1len*width)
                    image[y,x]=255
        break

    return image,seq1len,seq2len



st.write("# Dot matrix")

file1=st.sidebar.file_uploader("Sequence file 1:")
file2=st.sidebar.file_uploader("Sequence file 2:")
window=st.sidebar.slider("Window size:",4,100,10)

if file1 and file2:
    with StringIO(file1.getvalue().decode("utf-8")) as f1, StringIO(file2.getvalue().decode("utf-8")) as f2:
        image,seq1len,seq2len=dotmatrix(f1,f2,window,1000)
        plt.imshow(image,extent=(1,seq1len,seq2len,1),cmap="Grays")
        plt.xlabel(file1.name)
        plt.ylabel(file2.name)

st.pyplot(plt)

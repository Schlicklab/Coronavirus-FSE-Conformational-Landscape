library('seqinr')
library('Biostrings')

COVIDseq = 'UUUAAACGGGUUUGCGGUGUAAGUGCAGCCCGUCUUACACCGUGCGGCACAGGCACUAGUACUGAUGUCGUAUACAGGGCUUUU'

PEDV = read.fasta(file='PEDV.fasta')
PEDVseq = PEDV[[1]]
PEDVseq = toupper(PEDVseq)
PEDVseq = chartr('T', 'U', PEDVseq)
PEDVseq = c2s(PEDVseq)
PEDVlocal = pairwiseAlignment(PEDVseq, 'UUUAAACGAG', gapOpening = -2, gapExtension = -1, scoreOnly = FALSE, type="local")
PEDVseq = s2c(PEDVseq)
pos = 12610
start = pos-108
end = pos+113
PEDV222 = PEDVseq[start:end]
PEDV222 = c2s(PEDV222)


PRCV = read.fasta(file='PRCV.fasta')
PRCVseq = PRCV[[1]]
PRCVseq = toupper(PRCVseq)
PRCVseq = chartr('T', 'U', PRCVseq)
PRCVseq = c2s(PRCVseq)
PRCVlocal = pairwiseAlignment(PRCVseq, 'UUUAAACGAG', gapOpening = -2, gapExtension = -1, scoreOnly = FALSE, type="local")
PRCVseq = s2c(PRCVseq)
pos = 12320
start = pos-108
end = pos+113
PRCV222 = PRCVseq[start:end]
PRCV222 = c2s(PRCV222)


TGEV = read.fasta(file='TGEV.fasta')
TGEVseq = TGEV[[1]]
TGEVseq = toupper(TGEVseq)
TGEVseq = chartr('T', 'U', TGEVseq)
TGEVseq = c2s(TGEVseq)
TGEVlocal = pairwiseAlignment(TGEVseq, 'UUUAAACGAG', gapOpening = -2, gapExtension = -1, scoreOnly = FALSE, type="local")
TGEVseq = s2c(TGEVseq)
pos = 12332
start = pos-108
end = pos+113
TGEV222 = TGEVseq[start:end]
TGEV222 = c2s(TGEV222)


Feline = read.fasta(file='Feline.fasta')
Felineseq = Feline[[1]]
Felineseq = toupper(Felineseq)
Felineseq = chartr('T', 'U', Felineseq)
Felineseq = c2s(Felineseq)
Felinelocal = pairwiseAlignment(Felineseq, 'UUUAAACGAG', gapOpening = -2, gapExtension = -1, scoreOnly = FALSE, type="local")
Felineseq = s2c(Felineseq)
pos = 12371
start = pos-108
end = pos+113
Feline222 = Felineseq[start:end]
Feline222 = c2s(Feline222)


Canine = read.fasta(file='Canine.fasta')
Canineseq = Canine[[1]]
Canineseq = toupper(Canineseq)
Canineseq = chartr('T', 'U', Canineseq)
Canineseq = c2s(Canineseq)
Caninelocal = pairwiseAlignment(Canineseq, 'UUUAAACGAG', gapOpening = -2, gapExtension = -1, scoreOnly = FALSE, type="local")
Canineseq = s2c(Canineseq)
pos = 13335
start = pos-108
end = pos+113
Canine222 = Canineseq[start:end]
Canine222 = c2s(Canine222)


Sc_Bat = read.fasta(file='Sc_Bat.fasta')
Sc_Batseq = Sc_Bat[[1]]
Sc_Batseq = toupper(Sc_Batseq)
Sc_Batseq = chartr('T', 'U', Sc_Batseq)
Sc_Batseq = c2s(Sc_Batseq)
Sc_Batlocal = pairwiseAlignment(Sc_Batseq, 'UUUAAACGAG', gapOpening = -2, gapExtension = -1, scoreOnly = FALSE, type="local")
Sc_Batseq = s2c(Sc_Batseq)
pos = 12644
start = pos-108
end = pos+113
Sc_Bat222 = Sc_Batseq[start:end]
Sc_Bat222 = c2s(Sc_Bat222)


Mi_Bat = read.fasta(file='Mi_Bat.fasta')
Mi_Batseq = Mi_Bat[[1]]
Mi_Batseq = toupper(Mi_Batseq)
Mi_Batseq = chartr('T', 'U', Mi_Batseq)
Mi_Batseq = c2s(Mi_Batseq)
Mi_Batlocal = pairwiseAlignment(Mi_Batseq, 'UUUAAACGAG', gapOpening = -2, gapExtension = -1, scoreOnly = FALSE, type="local")
Mi_Batseq = s2c(Mi_Batseq)
pos = 12929
start = pos-108
end = pos+113
Mi_Bat222 = Mi_Batseq[start:end]
Mi_Bat222 = c2s(Mi_Bat222)


Ferret = read.fasta(file='Ferret.fasta')
Ferretseq = Ferret[[1]]
Ferretseq = toupper(Ferretseq)
Ferretseq = chartr('T', 'U', Ferretseq)
Ferretseq = c2s(Ferretseq)
Ferretlocal = pairwiseAlignment(Ferretseq, 'UUUAAACGAG', gapOpening = -2, gapExtension = -1, scoreOnly = FALSE, type="local")
Ferretseq = s2c(Ferretseq)
pos = 12195
start = pos-108
end = pos+113
Ferret222 = Ferretseq[start:end]
Ferret222 = c2s(Ferret222)


Mink = read.fasta(file='Mink.fasta')
Minkseq = Mink[[1]]
Minkseq = toupper(Minkseq)
Minkseq = chartr('T', 'U', Minkseq)
Minkseq = c2s(Minkseq)
Minklocal = pairwiseAlignment(Minkseq, 'UUUAAACGAG', gapOpening = -2, gapExtension = -1, scoreOnly = FALSE, type="local")
Minkseq = s2c(Minkseq)
pos = 12319
start = pos-108
end = pos+113
Mink222 = Minkseq[start:end]
Mink222 = c2s(Mink222)


HCov_NL63 = read.fasta(file='HCov_NL63.fasta')
HCov_NL63seq = HCov_NL63[[1]]
HCov_NL63seq = toupper(HCov_NL63seq)
HCov_NL63seq = chartr('T', 'U', HCov_NL63seq)
HCov_NL63seq = c2s(HCov_NL63seq)
HCov_NL63local = pairwiseAlignment(HCov_NL63seq, 'UUUAAACGAG', gapOpening = -2, gapExtension = -1, scoreOnly = FALSE, type="local")
HCov_NL63seq = s2c(HCov_NL63seq)
pos = 12387
start = pos-108
end = pos+113
HCov_NL63222 = HCov_NL63seq[start:end]
HCov_NL63222 = c2s(HCov_NL63222)


HCov_229E = read.fasta(file='HCov_229E.fasta')
HCov_229Eseq = HCov_229E[[1]]
HCov_229Eseq = toupper(HCov_229Eseq)
HCov_229Eseq = chartr('T', 'U', HCov_229Eseq)
HCov_229Eseq = c2s(HCov_229Eseq)
HCov_229Elocal = pairwiseAlignment(HCov_229Eseq, 'UUUAAACGAG', gapOpening = -2, gapExtension = -1, scoreOnly = FALSE, type="local")
HCov_229Eseq = s2c(HCov_229Eseq)
pos = 12428
start = pos-108
end = pos+113
HCov_229E222 = HCov_229Eseq[start:end]
HCov_229E222 = c2s(HCov_229E222)


write.fasta(PRCV222, 'PRCV', 'Alpha_222nt.fasta', open='a', nbchar=60)
write.fasta(TGEV222, 'TGEV', 'Alpha_222nt.fasta', open='a', nbchar=60)
write.fasta(Feline222, 'Feline', 'Alpha_222nt.fasta', open='a', nbchar=60)
write.fasta(Canine222, 'Canine', 'Alpha_222nt.fasta', open='a', nbchar=60)
write.fasta(HCov_NL63222, 'HCov_NL63', 'Alpha_222nt.fasta', open='a', nbchar=60)
write.fasta(Sc_Bat222, 'Sc_Bat', 'Alpha_222nt.fasta', open='a', nbchar=60)
write.fasta(PEDV222, 'PEDV', 'Alpha_222nt.fasta', open='a', nbchar=60)
write.fasta(Mi_Bat222, 'Mi_Bat', 'Alpha_222nt.fasta', open='a', nbchar=60)
write.fasta(Mink222, 'Mink', 'Alpha_222nt.fasta', open='a', nbchar=60)
write.fasta(Ferret222, 'Ferret', 'Alpha_222nt.fasta', open='a', nbchar=60)
write.fasta(HCov_229E222, 'HCov_229E', 'Alpha_222nt.fasta', open='a', nbchar=60)



# 84-nt
PEDV = read.fasta(file='PEDV.fasta')
PEDVseq = PEDV[[1]]
PEDVseq = toupper(PEDVseq)
PEDVseq = chartr('T', 'U', PEDVseq)
PEDVseq = c2s(PEDVseq)
PEDVlocal = pairwiseAlignment(PEDVseq, 'UUUAAACGAG', gapOpening = -2, gapExtension = -1, scoreOnly = FALSE, type="local")
PEDVseq = s2c(PEDVseq)
pos = 12610
start = pos
end = pos+83
PEDV84 = PEDVseq[start:end]
PEDV84 = c2s(PEDV84)


PRCV = read.fasta(file='PRCV.fasta')
PRCVseq = PRCV[[1]]
PRCVseq = toupper(PRCVseq)
PRCVseq = chartr('T', 'U', PRCVseq)
PRCVseq = c2s(PRCVseq)
PRCVlocal = pairwiseAlignment(PRCVseq, 'UUUAAACGAG', gapOpening = -2, gapExtension = -1, scoreOnly = FALSE, type="local")
PRCVseq = s2c(PRCVseq)
pos = 12320
start = pos
end = pos+83
PRCV84 = PRCVseq[start:end]
PRCV84 = c2s(PRCV84)


TGEV = read.fasta(file='TGEV.fasta')
TGEVseq = TGEV[[1]]
TGEVseq = toupper(TGEVseq)
TGEVseq = chartr('T', 'U', TGEVseq)
TGEVseq = c2s(TGEVseq)
TGEVlocal = pairwiseAlignment(TGEVseq, 'UUUAAACGAG', gapOpening = -2, gapExtension = -1, scoreOnly = FALSE, type="local")
TGEVseq = s2c(TGEVseq)
pos = 12332
start = pos
end = pos+83
TGEV84 = TGEVseq[start:end]
TGEV84 = c2s(TGEV84)


Feline = read.fasta(file='Feline.fasta')
Felineseq = Feline[[1]]
Felineseq = toupper(Felineseq)
Felineseq = chartr('T', 'U', Felineseq)
Felineseq = c2s(Felineseq)
Felinelocal = pairwiseAlignment(Felineseq, 'UUUAAACGAG', gapOpening = -2, gapExtension = -1, scoreOnly = FALSE, type="local")
Felineseq = s2c(Felineseq)
pos = 12371
start = pos
end = pos+83
Feline84 = Felineseq[start:end]
Feline84 = c2s(Feline84)


Canine = read.fasta(file='Canine.fasta')
Canineseq = Canine[[1]]
Canineseq = toupper(Canineseq)
Canineseq = chartr('T', 'U', Canineseq)
Canineseq = c2s(Canineseq)
Caninelocal = pairwiseAlignment(Canineseq, 'UUUAAACGAG', gapOpening = -2, gapExtension = -1, scoreOnly = FALSE, type="local")
Canineseq = s2c(Canineseq)
pos = 13335
start = pos
end = pos+83
Canine84 = Canineseq[start:end]
Canine84 = c2s(Canine84)


Sc_Bat = read.fasta(file='Sc_Bat.fasta')
Sc_Batseq = Sc_Bat[[1]]
Sc_Batseq = toupper(Sc_Batseq)
Sc_Batseq = chartr('T', 'U', Sc_Batseq)
Sc_Batseq = c2s(Sc_Batseq)
Sc_Batlocal = pairwiseAlignment(Sc_Batseq, 'UUUAAACGAG', gapOpening = -2, gapExtension = -1, scoreOnly = FALSE, type="local")
Sc_Batseq = s2c(Sc_Batseq)
pos = 12644
start = pos
end = pos+83
Sc_Bat84 = Sc_Batseq[start:end]
Sc_Bat84 = c2s(Sc_Bat84)


Mi_Bat = read.fasta(file='Mi_Bat.fasta')
Mi_Batseq = Mi_Bat[[1]]
Mi_Batseq = toupper(Mi_Batseq)
Mi_Batseq = chartr('T', 'U', Mi_Batseq)
Mi_Batseq = c2s(Mi_Batseq)
Mi_Batlocal = pairwiseAlignment(Mi_Batseq, 'UUUAAACGAG', gapOpening = -2, gapExtension = -1, scoreOnly = FALSE, type="local")
Mi_Batseq = s2c(Mi_Batseq)
pos = 12929
start = pos
end = pos+83
Mi_Bat84 = Mi_Batseq[start:end]
Mi_Bat84 = c2s(Mi_Bat84)


Ferret = read.fasta(file='Ferret.fasta')
Ferretseq = Ferret[[1]]
Ferretseq = toupper(Ferretseq)
Ferretseq = chartr('T', 'U', Ferretseq)
Ferretseq = c2s(Ferretseq)
Ferretlocal = pairwiseAlignment(Ferretseq, 'UUUAAACGAG', gapOpening = -2, gapExtension = -1, scoreOnly = FALSE, type="local")
Ferretseq = s2c(Ferretseq)
pos = 12195
start = pos
end = pos+83
Ferret84 = Ferretseq[start:end]
Ferret84 = c2s(Ferret84)


Mink = read.fasta(file='Mink.fasta')
Minkseq = Mink[[1]]
Minkseq = toupper(Minkseq)
Minkseq = chartr('T', 'U', Minkseq)
Minkseq = c2s(Minkseq)
Minklocal = pairwiseAlignment(Minkseq, 'UUUAAACGAG', gapOpening = -2, gapExtension = -1, scoreOnly = FALSE, type="local")
Minkseq = s2c(Minkseq)
pos = 12319
start = pos
end = pos+83
Mink84 = Minkseq[start:end]
Mink84 = c2s(Mink84)


HCov_NL63 = read.fasta(file='HCov_NL63.fasta')
HCov_NL63seq = HCov_NL63[[1]]
HCov_NL63seq = toupper(HCov_NL63seq)
HCov_NL63seq = chartr('T', 'U', HCov_NL63seq)
HCov_NL63seq = c2s(HCov_NL63seq)
HCov_NL63local = pairwiseAlignment(HCov_NL63seq, 'UUUAAACGAG', gapOpening = -2, gapExtension = -1, scoreOnly = FALSE, type="local")
HCov_NL63seq = s2c(HCov_NL63seq)
pos = 12387
start = pos
end = pos+83
HCov_NL6384 = HCov_NL63seq[start:end]
HCov_NL6384 = c2s(HCov_NL6384)


HCov_229E = read.fasta(file='HCov_229E.fasta')
HCov_229Eseq = HCov_229E[[1]]
HCov_229Eseq = toupper(HCov_229Eseq)
HCov_229Eseq = chartr('T', 'U', HCov_229Eseq)
HCov_229Eseq = c2s(HCov_229Eseq)
HCov_229Elocal = pairwiseAlignment(HCov_229Eseq, 'UUUAAACGAG', gapOpening = -2, gapExtension = -1, scoreOnly = FALSE, type="local")
HCov_229Eseq = s2c(HCov_229Eseq)
pos = 12428
start = pos
end = pos+83
HCov_229E84 = HCov_229Eseq[start:end]
HCov_229E84 = c2s(HCov_229E84)


write.fasta(PRCV84, 'PRCV', 'Alpha_84nt.fasta', open='a', nbchar=60)
write.fasta(TGEV84, 'TGEV', 'Alpha_84nt.fasta', open='a', nbchar=60)
write.fasta(Feline84, 'Feline', 'Alpha_84nt.fasta', open='a', nbchar=60)
write.fasta(Canine84, 'Canine', 'Alpha_84nt.fasta', open='a', nbchar=60)
write.fasta(HCov_NL6384, 'HCov_NL63', 'Alpha_84nt.fasta', open='a', nbchar=60)
write.fasta(Sc_Bat84, 'Sc_Bat', 'Alpha_84nt.fasta', open='a', nbchar=60)
write.fasta(PEDV84, 'PEDV', 'Alpha_84nt.fasta', open='a', nbchar=60)
write.fasta(Mi_Bat84, 'Mi_Bat', 'Alpha_84nt.fasta', open='a', nbchar=60)
write.fasta(Mink84, 'Mink', 'Alpha_84nt.fasta', open='a', nbchar=60)
write.fasta(Ferret84, 'Ferret', 'Alpha_84nt.fasta', open='a', nbchar=60)
write.fasta(HCov_229E84, 'HCov_229E', 'Alpha_84nt.fasta', open='a', nbchar=60)


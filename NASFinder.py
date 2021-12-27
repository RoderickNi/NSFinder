
#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
import os
print(
'''
 __          __         _                                       _             
 \ \        / /        | |                                     | |            
  \ \  /\  / /    ___  | |   ___    ___    _ __ ___     ___    | |_    ___    
   \ \/  \/ /    / _ \ | |  / __|  / _ \  | '_ ` _ \   / _ \   | __|  / _ \   
    \  /\  /    |  __/ | | | (__  | (_) | | | | | | | |  __/   | |_  | (_) |  
     \/  \/      \___| |_|  \___|  \___/  |_| |_| |_|  \___|    \__|  \___/   
                                                                              
                                                                              
          _   _    _____   ______   _               _                         
         | \ | |  / ____| |  ____| (_)             | |                        
         |  \| | | (___   | |__     _   _ __     __| |   ___   _ __           
         | . ` |  \___ \  |  __|   | | | '_ \   / _` |  / _ \ | '__|          
         | |\  |  ____) | | |      | | | | | | | (_| | |  __/ | |             
         |_| \_| |_____/  |_|      |_| |_| |_|  \__,_|  \___| |_|             
                                                                              
                                                                            
                                          -----Author: Ruoyao Ni
''')
print('''
#  NASFinder 1.0  (Nonsynonymous Substitution Finder)
#  
#  Copyright 2021 Ruoyao Ni,Zoology,CAS <niruoyao@ioz.ac.cn>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
''')


def read_ref(RefFile):
    '''Read ref seq from a FASTA file'''
    REF=open(RefFile,'r').read().split('\n')
    while '' in REF:
        REF.remove('')
    return REF

#=========================================================    
def output(OutputFile):
    '''Define output file path and name'''
    out=open(OutputFile,'a',encoding='utf-8')
    return out

#=========================================================    
def Base_annotation(RefFile):
    '''Exon or not'''
    REF=read_ref(RefFile)
    Seq=REF[1]
    REF_position=0
    cds_position=0
    num=0
    Ex_or_not=''
    Base_index = []
    for i in range(len(Seq)):
        REF_position+=1
        if Seq[i].islower():
            Ex_or_not="."
            Base_index.append(str(REF_position)+'\t'+r'.'+'\t'+Seq[i]+'\t'+Ex_or_not)
        elif Seq[i].isupper():
            if Ex_or_not == ".":
                num+=1
            Ex_or_not='Exon%d'%(num)
            cds_position+=1
            Base_index.append(str(REF_position)+'\t'+str(cds_position)+'\t'+Seq[i]+'\t'+Ex_or_not)
    return Base_index                    # list of base info

#=========================================================    
def read_vcf(VcfFile):
    '''
    Read vcf file as list
    '''
    vcf=open(VcfFile,'r').read().split('\n')
    while '' in vcf:
        vcf.remove('')
    SNP=[]
    for line in vcf:
        if r'#' not in line:
            SNP.append(line)
    return SNP
   
def get_exon_base_pos(RefFile):
    '''
    Getting exon base positon in DNAseq
    '''
    Ann=Base_annotation(RefFile)
    exon_base_pos = []
    for line in Ann:
        if r'Exon' in line:
            exon_base_pos.append(int(line.split('\t')[0]))
    return exon_base_pos
    

def vcf_exon(RefFile,VcfFile):
    '''
    Getting exon-snp info lst
    '''
    SNP_info=read_vcf(VcfFile)
    Exon_pos_lst=get_exon_base_pos(RefFile)
    vcf_exon_lst=[]
    for line in SNP_info:
        pos = int(line.split('\t')[1])
        ref = line.split('\t')[3]
        alt = line.split('\t')[4]
        if len(ref) == 1:
            if pos in Exon_pos_lst:
                vcf_exon_lst.append(line)    
        else:
            pos_list = list(range(pos,int(pos+len(ref))))
            for element in pos_list:
                if element in Exon_pos_lst:
                    vcf_exon_lst.append(line)
                    break
    return vcf_exon_lst

# =============================================================================
Codon_table =   {'UUU': 'F', 'CUU': 'L', 'AUU': 'I', 'GUU': 'V', 
                 'UUC': 'F', 'CUC': 'L', 'AUC': 'I', 'GUC': 'V', 
                 'UUA': 'L', 'CUA': 'L', 'AUA': 'I', 'GUA': 'V', 
                 'UUG': 'L', 'CUG': 'L', 'AUG': 'M', 'GUG': 'V', 
                 'UCU': 'S', 'CCU': 'P', 'ACU': 'T', 'GCU': 'A', 
                 'UCC': 'S', 'CCC': 'P', 'ACC': 'T', 'GCC': 'A', 
                 'UCA': 'S', 'CCA': 'P', 'ACA': 'T', 'GCA': 'A', 
                 'UCG': 'S', 'CCG': 'P', 'ACG': 'T', 'GCG': 'A', 
                 'UAU': 'Y', 'CAU': 'H', 'AAU': 'N', 'GAU': 'D', 
                 'UAC': 'Y', 'CAC': 'H', 'AAC': 'N', 'GAC': 'D', 
                 'UAA': '*', 'CAA': 'Q', 'AAA': 'K', 'GAA': 'E', 
                 'UAG': '*', 'CAG': 'Q', 'AAG': 'K', 'GAG': 'E', 
                 'UGU': 'C', 'CGU': 'R', 'AGU': 'S', 'GGU': 'G', 
                 'UGC': 'C', 'CGC': 'R', 'AGC': 'S', 'GGC': 'G', 
                 'UGA': '*', 'CGA': 'R', 'AGA': 'R', 'GGA': 'G', 
                 'UGG': 'W', 'CGG': 'R', 'AGG': 'R', 'GGG': 'G', 'AUG': 'M'
                 }

def DB_cDNA(RefFile):
    '''
    Build a index for REF
   [[ref base positon],[codon base position],[codon(nt)],[exon_num],{"AA_position":"Amino Acide"}]   
    '''
    RefAnn=Base_annotation(RefFile)
    cnt = 0
    AA_pos = 0
    Ref_base_pos=[]
    codon_pos=[]
    codon=[]
    AA=()
    Data_base=[]
    for line in RefAnn:
        if r'Exon' in line:
            cnt+=1
            if cnt < 3:
                Ref_base_pos.append(line.split('\t')[0])
                codon_pos.append(line.split('\t')[1])
                codon.append(line.split('\t')[2])
                
            elif cnt == 3:
                AA_pos+=1
                Ref_base_pos.append(line.split('\t')[0]) # list
                codon_pos.append(line.split('\t')[1]) # list
                codon.append(line.split('\t')[2]) # list
                a = ''.join(codon)
                b = a.replace('T','U')
                AA=(str(AA_pos),Codon_table[b]) # dic
                EXON=line.split('\t')[3] # str
                
                Data_base.append([Ref_base_pos,codon_pos,codon,EXON,AA])
                
                Ref_base_pos=[]
                codon_pos=[]
                codon=[]
                AA=()                
                cnt = 0
    return Data_base # list

def DB_vcf(RefFile,VcfFile):
    '''
    Build a index for SNP
    [[ref_base_pos],[ref_nt],[alt_nt],(AO,DP,Freq)]
    '''
    SNP_info=vcf_exon(RefFile,VcfFile)
    Data_base=[]
    ref_base_pos=[]
    ref_nt=[]
    alt_nt=[]
    depth=()
    for line in SNP_info:
        if r',' not in line.split('\t')[4]:      # Single varient
            if len(line.split('\t')[3]) == 1 and len(line.split('\t')[4]) == 1: ## Single base
                ref_base_pos.append(line.split('\t')[1])
                ref_nt.append(line.split('\t')[3])
                alt_nt.append(line.split('\t')[4])
                DP=int(line.split('\t')[9].split(r':')[1])
                AO=int(line.split('\t')[9].split(r':')[4])
                depth=(DP,AO,round(AO/DP,5))
                Data_base.append([ref_base_pos,ref_nt,alt_nt,depth])
                ref_base_pos=[]
                ref_nt=[]
                alt_nt=[]
                depth=()
            elif len(line.split('\t')[3]) != 1 and len(line.split('\t')[3]) == len(line.split('\t')[4]):  ## Multiple base
                for pos in range(int(line.split('\t')[1]),int(line.split('\t')[1])+len(line.split('\t')[3])):
                    ref_base_pos.append(str(pos))
                for nt in list(line.split('\t')[3]):
                    ref_nt.append(nt)
                for nt in list(line.split('\t')[4]):
                    alt_nt.append(nt)
                DP=int(line.split('\t')[9].split(r':')[1])
                AO=int(line.split('\t')[9].split(r':')[4])
                depth=(DP,AO,round(AO/DP,5))
                Data_base.append([ref_base_pos,ref_nt,alt_nt,depth])
                ref_base_pos=[]
                ref_nt=[]
                alt_nt=[]
                depth=()                 
        elif r',' in line.split('\t')[4]:       # Multiple varient
            if len(line.split('\t')[3]) == 1:   ## Single base
                ref_base_pos.append(line.split('\t')[1])
                ref_nt.append(line.split('\t')[3])
                var_num = 0
                for var in line.split('\t')[4].split(','):
                    if len(var) == 1:
                        nt = var
                        alt_nt.append(nt)
                        DP=int(line.split('\t')[9].split(r':')[1])
                        AO=int(line.split('\t')[9].split(r':')[4].split(',')[var_num])                        
                        depth=(DP,AO,round(AO/DP,5))
                        Data_base.append([ref_base_pos,ref_nt,alt_nt,depth])
                        alt_nt=[]
                        depth=()
                        var_num += 1
                    else:
                        var_num += 1
                        continue
                ref_base_pos=[]
                ref_nt=[]
                alt_nt=[]
                depth=()                 
            elif len(line.split('\t')[3]) != 1: ## Multiple base
                for pos in range(int(line.split('\t')[1]),int(line.split('\t')[1])+len(line.split('\t')[3])):
                    ref_base_pos.append(str(pos))
                for nt in list(line.split('\t')[3]):
                    ref_nt.append(nt)
                var_num = 0
                for var in line.split('\t')[4].split(','):
                    if len(var) == len(line.split('\t')[3]):
                        for nt in list(var):
                            alt_nt.append(nt)
                        DP=int(line.split('\t')[9].split(r':')[1])
                        AO=int(line.split('\t')[9].split(r':')[4].split(',')[var_num])
                        depth=(DP,AO,round(AO/DP,5))
                        Data_base.append([ref_base_pos,ref_nt,alt_nt,depth])
                        alt_nt=[]
                        depth=()
                        var_num += 1
                    else:
                        var_num += 1
                        continue
                ref_base_pos=[]
                ref_nt=[]
                alt_nt=[]
                depth=()
    return Data_base

def NASF(RefFile,VcfFile,OTPT,T=0.005):
    '''
    Nonsynonymous amino acid substitution finding
    '''
    cDNA_DB = DB_cDNA(RefFile)
    SNP_DB = DB_vcf(RefFile,VcfFile)     
    output = open(OTPT,'a',encoding='utf-8')
    for SNP_line in SNP_DB:
        for REF_line in cDNA_DB:
            SPIR_lst=[]              # SPIR(SNP position in Ref)
            for i in range(len(SNP_line[0])):        # 'i' recording the index
                if SNP_line[0][i] in REF_line[0]:
                    if SNP_line[1][i] != SNP_line[2][i]:
                        SPIR_lst.append(SNP_line[0][i])
            if SPIR_lst != []:
                ref_codon=''.join(REF_line[2])
                int_codon=[]
                for base in REF_line[2]:
                    int_codon.append(base)
                for alt_pos in SPIR_lst:
                    alt_base_pos = REF_line[0].index(alt_pos)
                    alt_base = SNP_line[2][SNP_line[0].index(alt_pos)]
                    int_codon[int(alt_base_pos)]=alt_base
                alt_codon=''.join(int_codon)
                ref_AA=Codon_table[ref_codon.replace('T','U')]
                alt_AA=Codon_table[alt_codon.replace('T','U')]
                if ref_AA != alt_AA:
                    AA_pos = REF_line[4][0]
                    EXON = REF_line[3]
                    DP = str(SNP_line[3][0])
                    AO = str(SNP_line[3][1])
                    Frequency = str(SNP_line[3][2])
                    if float(Frequency) >= T:
                        print('\t'.join([ref_AA+AA_pos+alt_AA,ref_codon+'->'+alt_codon,'AO/DP'+'('+AO+'/'+DP+')',Frequency,EXON]),
                              file = output)
    output.close()


def main():
    while True:
        print("<Enter q to Quit!>")
        RefFile=input('#  Please input the path of reference sequence file (.fasta): ')
        if RefFile == 'q':
            break
        VcfFile=input('#  Please input the path of SNP file (.vcf): ')
        if VcfFile == 'q':
            break
        T=input('#  Frequency greater than or equal to : ')
        if T == 'q':
            break
        OTPT=input('#  Result Output: ')
        NASF(RefFile,VcfFile,OTPT,float(T))
        print("#  Completed!\n")
    print("THANKS FOR USING NASFinder!")

if __name__ == "__main__":
    main()


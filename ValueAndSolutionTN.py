#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 19:11:13 2020

@author: margarida
"""


import numpy as np
import math
import sys
import random
#from sympy import init_printing, latex
#import sympy as sp
#from sympy.printing.latex import print_latex
# importing networkx  
import networkx as nx
# TN is the transportation network

def ValueSolutionMinPairingProblem(i,b=None):
    # i should be an integer >= 2
    if type(i) != type(1):
        print("Error: argument i must be an integer")
        sys.exit()
    if i < 2:
        print("Error: argument i must be greater then or equal to 2")

    
    # Solution vectors
    jcd = np.zeros(4*i)
    rfrak = np.zeros(i-1) # corresponds to arc (u_j,t_j) for j<i
    dfrak = np.zeros(2*i-1)
    esse = np.zeros(2*i)
    
    h = np.zeros(2*i)
    
    '''
    arctu maps arc (t_l,u_j) to appropriate flow variable in s
    '''
    def arctu(l,j):
        alfa = np.sign(j-l)
        if alfa > 0 and l < i:
            if l%2 == 1: return (l-1) # should be j, but python starts at zero
            else: return(2*i+1-l-1)
        if alfa == 0:
            return(i-1)
        if alfa < 0 and l%2 == 1: return (l-2)
        if alfa < 0 and l%2 == 0: return(2*i+1-l)
    
    '''
    arcft maps arc (f_l,t_j) to appropriate flow variable in jcd
    '''
    def arcft(l,j):
        if j> l+1 or j < l:
            print("Error: no arc dfrak from f_",l," to t_",j)
            sys.exit()
        if j == l:
            if l%2 == 1:
                return(2*j-1) # index of flow j^d injcd array
            else: return(2*(2*i+1-j)-2) # index of flow (n-j)^c
        elif l%2 == 0:
            return(2*(l+1)-2) #flow (j+1)^c
        else: return(2*(2*i+1-l-1)-1) # flow (n-j-1)^d
    
    '''
    arcug maps arc (u_l,g_j) to appropriate flow variable in jcd
    '''
    def arcug(l,j):
        if j> l or j < -l:
            print("Error: no arc from u_",l," to g_",j)
            sys.exit()
        if j == l:
            if l%2 == 0:
                return(2*l-1) # index of flow j^d injcd array
            else: return(2*(2*i+1-l)-2) # index of flow (n-j)^c
        elif l%2 == 0:
            return(2*(l-1)) #flow j^c
        else: return(2*(2*i+1-l)-1) # flow (n-j)^d        
    
    '''
    arcut maps arc (u_l,t_j) to appropriate flow variable in dfrak
    '''
    def arcut(l,j):
        alfa = np.sign(j-l)
        if alfa == 0:
            if l != i:
                print("Error: arc (u_j,t_j), for j!=i doesn't exist ")
                sys.exit()
            else: return (i-1)
        if l%2 == 0:
            if alfa > 0:
                return (l-1) # seria j, mas python comeca em zero
            else: return (l-2) # idem
        elif alfa > 0:
            return(2*i+1-l-2) # idem
        else: return(2*i+1-l-1) # idem
    
    '''
    pathtu(l,j) list the nodes in minimum path from t_l to u_j
    '''
    def pathtu(l,j):
        if min(l,j) < 1:
            print("Error: t and u indices must be >= 1")
            sys.exit()
        if max(l,j)>i:
            print("Error: t and u indices must be <= i")
            sys.exit()
        if abs(l-j) == 0:
            if l<i:
                return [l,l+1,l+1,l]
            else: return [l,l-1,l-1,l]
        alfa = np.sign(j-l)
        if abs(l-j)%2 == 1:
            p=[item for item in range(l,j+alfa,alfa)]
            return p
        else:
            p=[l,l+alfa]
            p = p + [item for item in range(l+alfa,j+alfa,alfa)]
            return p
    
    def addflowpathtu(t,u,x):
        print("t =",t," u = ",u)
        path = pathtu(t,u)
        aux = int(len(path)/2)
        for k in range(aux):
            # Arco t -> u
            esse[arctu(path[2*k],path[2*k+1])]+=x
            if k < aux-1:
                # arco u -> t
                if path[2*k+2]==path[2*k+1]:
                    rfrak[path[2*k+1]-1]+=x
                else: dfrak[arcut(path[2*k+1],path[2*k+2])]+=x
    
    '''
    addflow(l,j) increases flow in original variables corresponding to flow (f_l,g_j)
    '''
    def addflow(l,j,x):
        # f0 -> g0
        if l == 0 and j == 0:
            jcd[arcft(l,l+1)]+=x
            jcd[arcug(l+1,j)]+=x
            esse[arctu(1,2)]+=x
            if i == 2: dfrak[arcut(2,2)]+=x
            else: rfrak[1]+=x
            esse[arctu(2,1)]+=x
        # f_0 -> g_1
        if l == 0 and j == 1:
            jcd[arcft(l,l+1)]+=x
            jcd[arcug(l+2,j)]+=x
            esse[arctu(1,2)]+=x
        #f_0 -> j, j>=2
        if l == 0 and j >= 2:
            jcd[arcft(l,l+1)]+=x
            jcd[arcug(j,j)]+=x
            addflowpathtu(l+1,j,x)
        # f_l -> g_l,  1<= l < i
        if l >= 1 and l < i and j == l:
            jcd[arcft(l, l)]+=x
            esse[arctu(l,l+1)]+=x
            jcd[arcug(l+1,j)]+=x
        # f_1 -> g_0
        if l==1 and j == 0:
            jcd[arcft(1,2)]+=x
            esse[arctu(2,1)]+=x
            jcd[arcug(1,0)]+=x
        # f_l -> g_j  j=l+1, 1<= l <i
        if l >= 1 and l < i and j == l+1:
            jcd[arcft(l, l)]+=x
            esse[arctu(l,j)]+=x
            jcd[arcug(j,j)]+=x
        # f_l -> g_{l+2}, 1<= l <= i-2
        if l >= 1 and l < i-1 and j == l+2:
            jcd[arcft(l, l+1)]+=x
            jcd[arcug(j,j)]+=x
            esse[arctu(l+1,j)]+=x 
        # f_l -> g_j, 1<= l,j <i, |l - j|>= 3
        if abs(j-l)>=3 and 1 <= l and l < i and 1<= j and j <=i:
            if j>l:
                jcd[arcft(l,l+1)]+=x
                jcd[arcug(j,j)]+=x
                addflowpathtu(l+1, j, x)
            else:
                jcd[arcft(l,l)]+=x
                jcd[arcug(j+1,j)]+=x
                addflowpathtu(l, j+1, x)            
        if l >=2:
            # f_l -> g_0, l>=2
            if j == 0:
                jcd[arcft(l,l)]+=x
                jcd[arcug(1,0)]+=x
                addflowpathtu(l,1,x)
            # f_l -> g_{l-2} l>=2, j>=0
            if j == l-2 and j > 0:
                jcd[arcft(l,l)]+=x
                esse[arctu(l,l-1)]+=x
                jcd[arcug(l-1,l-2)]+=x
            # f_l -> g_j, j=l-1 2<= l <i
            if l < i and j == l-1:
                jcd[arcft(l, l)]+=x
                esse[arctu(l,j)]+=x
                jcd[arcug(j,j)]+=x
            if l == j and j == i:
                jcd[arcft(l,l)]+=x
                esse[arctu(l,l)]+=x
                jcd[arcug(l,l)]+=x
            if l == i and j <= i-1 and (i-j)%2 == 1 and j > 0:
                jcd[arcft(l,l)]+=x
                jcd[arcug(j,j)]+=x
                addflowpathtu(l,j,x)
            if l == i and j <= i-4 and j > 0 and (i-j)%2 == 0:
                jcd[arcft(l,l)]+=x
                jcd[arcug(j+1,j)]+=x
                addflowpathtu(l,j+1,x)    
    
    def gama0(k):
        if (k<0):
            print("ERRO: argumento de gama < 0")
            sys.exit()
        value = int(np.floor(k/2))
        return value
    def tilgama(k):
        if (k<0):
            print("ERRO: argumento de gama < 0")
            sys.exit()
        if (k==0):
            value = 0
        else: value = gama0(k-1)
        return value
    
    # Random generation of mathcal{B} 
    # b = B^-_0, B^+_0, B_1, ..., B_i
    # satisfies -B^-_0 + B^+_0 - B_1 + B_2 - ... =0 
    # Instances already tested:
    # dim 9 b = [3,0,0,0,0,-3]
    # dim 7 b=[7,2,-1,7,3]
    # dim 7 b=[5, 5, 2, -7, -9]
    # dim 9 b = [5, 5, 5, 3, 0, 2]
    # dim 11 b=[0,7,2,-1,-3,-2,5]
    # dim 17 b=[0,7,3,2,5,7,3,5,0,-10]
    # dim 17 b=[1,0,-2,6,-7,4,1,-4,-4,-17] 
    # dim 17 b=[0,3,-1,6,-4,1,0,-3,5,-7]
    # dim 19 b=[4,1,-4,2,-6,1,-1,-4,3,-6,-2]
    # dim 19 b=[1,3,-5,6,1,-2,6,3,-4,6,17]
    # dim 19 b=[0,6,-5,-6,4,-6,1,-1,0,4,-3]
    
    
    if b == None:
        b = np.array([random.randint(0,7) for j in range(2)])
        b = np.hstack([b,np.array([random.randint(-10,10) for j in range(i-1)]) ])
        soma = sum([(-1)**j * aux for j,aux in enumerate(b)])
        b = np.hstack([b,(-1)**i * soma])
    elif len(b) != i+2 or sum([(-1)**j * aux for j,aux in enumerate(b)]) != 0:
        print("Error: b vector entered is not valid")
        sys.exit()
    
    
    if sum([(-1)**j * aux for j,aux in enumerate(b)]) != 0:
        print("ERROR: alternate sum of mathcal{B} isn't zero")
        sys.exit()
    
    bmenos= [b[0]] + [max((-1)*b[j],0) for j in range(2,i+2)]
    bmais= [b[1]] + [max(b[j],0) for j in range(2,i+2)]
    
    '''
    TN digraph nodes
    0 to i are supply nodes with supplies (-B^-_0,-B^+_1,-B^-_2, etc)
    i+1 to 2*i+1 are demand nodes
    (B^+_0,B^-_1,B^+_2, etc)
    cost of arc (i,j) = tilgama(abs(i-(j-i-1)))
    except for arc ()
    '''
    TN = nx.DiGraph()  
    
    
    TN.add_node(0,demand=-b[0])
    TN.add_node(i+1,demand=b[1])
    for j in range(1,i+1):
        TN.add_node(j,demand=-max((-1)**(j+1) * b[1+j],0))
        TN.add_node(i+j+1,demand=max((-1)**j * b[1+j],0))
    
    for k in range(i+1):
        for j in range(i+1):
            TN.add_edge(k,i+j+1,weight=tilgama(abs(k-j)))
    TN[0][i+1]['weight']=1
    
    flowCostTN, flowDictTN = nx.network_simplex(TN)
    #print("Valor otimo do problema de transporte TN = ", flowCostTN)
    
    
    print("--------","\n","\\fbox{$i = ",i,"\\quad   n = ",2*i+1,"$}\\\\","\n")
    string="$(\\mathcal{B}^-_0,\\mathcal{B}^+_0,\\mathcal{B}) = ("
    for j in range(i+1):
        string+=str(b[j])+","
    string+=str(b[i+1])+")$\\\\[10pt]"
    print(string,"\n")
    
    # First upper bound
    first_upper_bound = 0
    for j in range(1,i+1):
        soma = 0
        for k in range(1,j+1):
            soma += (-1)**(k+j)*(bmenos[k-1]+bmais[k])
        soma = max((-1)*soma,0)
        first_upper_bound += soma
    first_upper_bound += int(np.floor((i-1)/2)*bmais[0])
    # Second upper bound
    second_upper_bound = 0
    for j in range(1,i+1):
        soma = 0
        for k in range(1,j+1):
            soma += (-1)**(k+j)*(bmais[k-1]+bmenos[k])
        soma = max((-1)*soma,0)
        second_upper_bound += soma
    second_upper_bound += int(np.floor(i/2)*bmenos[0])
    
    
    totaldemand = 0
    string = "$ \\begin{array}{|c|*{"+str(i+1)+"}{l@{\\hspace{2pt}}c|}c|}"
    #string += "\cline{1-"+str(2*i+3)+"}\n"
    #string+= "\multicolumn{"+str(2*i+3)
    #string+="}{|c|}{\mbox{Min cost flow for $n="+str(2*i+1)+"$}}\\\\"
    string += "\\hline\n"
    #menosmais=["-","+"]
    for j in range(i+1):
        string+="& \multicolumn{2}{c|}{\,g_{"+str(j)+"\,}}"
    string += "& \\textcolor{blue}{~\\theta~} \\\\ \hline"
    # loop em no de supply k
    for k in range(i+1):
        string+="\n f_{"+str(k)+"}"
        # loop em no de demanda node
        for j in range(i+1,2*i+2):
            value = TN[k][j]['weight']
            if value !=0:
                string+="& {}^"+str(value)
            else: string+="&"
            value = flowDictTN[k][j]
            if value !=0:
                addflow(k,j-i-1,value)
                string+="&\\textcolor{red}{"+str(value)+"} "
            else: string += "&"
        string+="& "
        if TN.nodes[k]['demand'] < 0:
            string+="\\textcolor{blue}{"+str(-TN.nodes[k]['demand'])+"}"
        string += "\\\\ \hline"
    string+="\n \\,\\textcolor{blue}{\\delta}\\,"
    for j in range(i+1,2*i+2):
        totaldemand += TN.nodes[j]['demand'] 
        string += "&\multicolumn{2}{c|}{ "
        if TN.nodes[j]['demand'] > 0:
            string+="\\textcolor{blue}{"+str(TN.nodes[j]['demand'])+"}}"
        else: string += "}"
    string += "&\\\\ \\hline \\end{array}$"
    print(string)
    
    
    #Calculo de lower bound pelas colunas
    lower_bound_col = 0
    # Contribuicao coluna 0
    excesso = 0
    for k in range(int(np.ceil(i/2))-1,0,-1):
        soma = 0
        for linha in range(max(0,1-2*k),min(2*k,i)+1):
            soma += TN.nodes[linha]['demand'] #supply of node linha
        lower_bound_col += k * max(0,TN.nodes[i+1]['demand'] - excesso + soma)
        excesso += max(0,TN.nodes[i+1]['demand']  + soma - excesso)
        #print("k = ",k,"lower bound = ",lower_bound_col," excesso = ",excesso,"\n")
    # Contribuicao da coluna j para j>=1
    for j in range (i+2,2*i+2):
        excesso=0
    #    print("\n ******* coluna ",j-3*i-1)
        for k in range(int(np.ceil(i/2))-1,0,-1):
            soma = 0
            for linha in range(max(0,j-i-1-2*k),min(j-i-1+2*k,i)+1):
                soma += TN.nodes[linha]['demand']
    #            print("linha = ",linha,"  oferta[linha] = ",
    #                 -TN.nodes[linha]['demand'])
            lower_bound_col += k*max(0,TN.nodes[j]['demand'] - excesso + soma)
            excesso += max(0,TN.nodes[j]['demand']  + soma - excesso)
    #        print("j = ",j,"   k = ",k,"   excesso = ",excesso," soma =",soma,
    #              "Demanda = ", TN.nodes[j]['demand']," col_lb = ",lower_bound_col)
    
    
    #Calculo de lower bound pelas linhas
    lower_bound_row = 0
    # Contribuicao linha 0
    excesso = 0
    for k in range(int(np.ceil(i/2))-1,0,-1):
        soma = 0
        for coluna in range(i+1+max(0,1-2*k),i+1+min(2*k,i)+1):
            soma += TN.nodes[coluna]['demand'] #demand of node coluna
        lower_bound_row += k * max(0,-TN.nodes[0]['demand'] - excesso - soma)
    #    print("k = ",k," excesso = ",excesso)
        excesso = max(0,-TN.nodes[0]['demand']  - soma)
    # Contribuicao da linha j para j>=1
    for j in range (1,i+1):
    #    print("\n ******* linha ",j)
        excesso=0
        for k in range(int(np.ceil(i/2))-1,0,-1):
            soma = 0
            for coluna in range(max(0,j-2*k),min(j+2*k,i)+1):
                soma += TN.nodes[i+1+coluna]['demand']
    #            print("coluna = ",coluna,"  demanda[coluna] = ",
    #                 TN.nodes[3*i+1+coluna]['demand'])
            lower_bound_row += k*max(0,-TN.nodes[j]['demand'] - excesso - soma)
            excesso = max(0,-TN.nodes[j]['demand']  - soma)
    #        print("j = ",j,"   k = ",k,"   excesso = ",excesso," soma =",soma,
    #              "Demanda = ", TN.nodes[j]['demand'])
    
    print('''\\hspace{1.5cm}\\parbox{5cm}{%
    \\textbf{Optimal value} \\dotfill \\textbf{''',flowCostTN,'''}\\\\[10pt]
    Column lower bound \dotfill''',lower_bound_col,'''\\\\ 
    Row lower bound  \dotfill ''',lower_bound_row,'''\\\\ 
    First upper bound  \dotfill ''',first_upper_bound,'''\\\\ 
    Second upper bound  \dotfill ''',second_upper_bound,"}\\\\[10pt]\n")
    
    hcd = [item for item in jcd]
    # hjc + = dfrak_{j-1} j=2,...,i
    for j in range(i):
        hcd[2*j]+=dfrak[j-1]
    hcd[2*i]+=dfrak[i-1]
    #print("1) hcd ",hcd)
    for j in range(i+1,2*i):
        hcd[2*j]+=rfrak[2*i-j-1]
        hcd[2*j]+=dfrak[j-1]
    #print("2) hcd ",hcd)
    for j in range(i-1):
        hcd[2*j+1]+=rfrak[j]+dfrak[j]
    hcd[2*i-1]+=dfrak[i-1]
    #print("3) hcd ",hcd)
    for j in range(i,2*i-1):
        hcd[2*j+1]+=dfrak[j]
    #print("4) hcd ",hcd)
        
    '''
    Computing h
    '''
    
    for j in range(2*i):
        h[j] = hcd[2*j]+hcd[2*j+1]
    
    header = ["c","d"]
    
    print("$\\begin{array}{|c||*{",4*i,"}{l@{\hspace{2pt}}c|}}")
    n = np.count_nonzero(jcd)
    print("\\cline{1-",2*n+1,"}\n j^{cd}")
    for j in range(4*i):
        if jcd[j] != 0:
            jj = int(j/2)
            hh = j%2
            print("&{}^{",jj+1,"^",header[hh],"}&\\textcolor{red}{",int(jcd[j]),"}") 
    print("\\\\ \\cline{1-",2*n+1,"}")
    n = np.count_nonzero(rfrak)
    print(" \\cline{1-",2*n+1,"}\n \mathfrak{r}")
    for j in range(i-1):
        if rfrak[j] != 0:
            print("&{}^{",j+1,"}&\\textcolor{red}{",int(rfrak[j]),"}") 
    print("\\\\ \\cline{1-",2*n+1,"}")
    n = np.count_nonzero(dfrak)
    print("\\cline{1-",2*n+1,"}\n \mathfrak{d}")
    for j in range(2*i-1):
        if dfrak[j] != 0:
            print("&{}^{",j+1,"}&\\textcolor{red}{",int(dfrak[j]),"}") 
    print("\\\\ \\cline{1-",2*n+1,"}")
    esse[i] = esse[i-1]
    n = np.count_nonzero(esse)
    print(" \\cline{1-",2*n+1,"}\n s")
    for j in range(2*i):
        if esse[j] != 0:
            print("&{}^{",j+1,"}&\\textcolor{red}{",int(esse[j]),"}") 
    print("\\\\ \\cline{1-",2*n+1,"}")
    
    
    n = np.count_nonzero(hcd)
    print("\\cline{1-",2*n+1,"}\n h^{cd}")
    for j in range(4*i):
        if hcd[j] != 0:
            jj = int(j/2)
            hh = j%2
            print("&{}^{h_{",jj+1,"}^",header[hh],"}&\\textcolor{red}{",int(hcd[j]),"}") 
    print("\\\\ \\cline{1-",2*n+1,"}")
    
    n = np.count_nonzero(h)
    print("\\cline{1-",2*n+1,"}\n h")
    for j in range(2*i):
        if h[j] != 0:
            print("&{}^{h_{",j+1,"}}&\\textcolor{red}{",int(h[j]),"}") 
    print("\\\\ \\cline{1-",2*n+1,"}\\end{array}$\n\n\\vspace{1cm}\n")
    
    '''
    Checking values
    '''
    #print("totaldemand = ",totaldemand)
    
    if sum(jcd) != 2*totaldemand:
        print("Error: jcd imbalance")
        print("sum(jcd) = ",sum(jcd))
        print("2*totaldemand = ",2*totaldemand)
    if sum(esse)-sum(dfrak)-sum(rfrak) != totaldemand:
        print("Error: wrong s or r or d")
        print("sum(s) = ",sum(esse)," sum(r) = ",sum(rfrak),
              " sum(d) = ",sum(dfrak)," totaldemand =",totaldemand)
    
    alternatesum = np.zeros(2*i)
    alternatesum[0] = h[0]
    alternatesum[2*i-1] = h[2*i-1]
    if alternatesum[0] < 0:
        print("Error: alternate sum ",0," is negative, equal to ",alternatesum[0])
    if alternatesum[0] != esse[0]:
        print("Error: s_0 = ",esse[0],"but the alternate sum is ",alternatesum[0])
    if alternatesum[2*i-1] < 0:
        print("Error: alternate sum ",2*i-1," is negative, equal to ",
              alternatesum[2*i-1])
    if alternatesum[2*i-1] != esse[2*i-1]:
        print("Error: s_",2*i-1," = ",esse[2*i-1],"but the alternate sum is ",
              alternatesum[2*i-1])
    
    for j in range(1,i):
        alternatesum[j] = h[j] - alternatesum[j-1]
        if alternatesum[j] < 0:
            print("Error: alternate sum ",j+1," is negative, equal to ",
                  alternatesum[j])
        if alternatesum[j] != esse[j]:
            print("Error: s_",j," = ",esse[j],"but the alternate sum is ",
                  alternatesum[j])
    for j in range(1,i):
        alternatesum[2*i-j-1] = h[2*i-j-1] - alternatesum[2*i-j]
        if alternatesum[2*i-j-1] < 0:
            print("Error: alternate sum ",2*i-j-1+1," is negative, equal to ",
                  alternatesum[2*i-j-1])
        if alternatesum[2*i-j-1] != esse[2*i-j-1]:
            print("Error: s_",2*i-j-1+1," = ",esse[2*i-j-1],"but the alternate sum is ",
                  alternatesum[2*i-j-1])

# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 15:02:12 2022

@author: user
"""
#P='1011'
#K=['011', '101']
#f=lambda x, y: str(int(x[0])^int(y[1]))+str(((int(x[1])&int(y[0]))|int(y[2])))

def E(P,K,f):#코드 짤 때, P를 문자열, K를 문자열들의 리스트들로 받아야 작동되도록 만들었음!
    C=P
    Lp=''
    for i in range(len(K)):
        if i < len(K)-1:
            L=C[0:(len(C)//2)]#len(C)/2의 type은 len(C)가 나눠떨어져도 무조건 float이므로 len(C)//2로 하기
            R=C[(len(C)//2):len(C)]#암호화되고 있는 문장을 절반씩으로 쪼개는 과정(16,17번째 문장)
            for j in range(len(R)):#C=R+(f(R,K[i])^L)->f와 L을 이진수로 변환한 다음 bitwise XOR연산해야 됨!!
                Lp=Lp+str(int(f(R,K[i])[j])^int(L[j]))
            C=R+Lp#이전 round의 R과 이번 라운드에서 암호화된 L(=Lp)를 서로 왼쪽, 오른쪽으로 순서를 바꿔 집어넣어야 됨.
            Lp=''#Lp는 초기화해야 다음 loop에서도 다시 새로운 L을 Lp에 입력받을 수 있음. 초기화하지 않으면 이전의 암호화된 Lp에 추가되는 것이기에 L의 길이가 계속 길어져서 C가 이상하게 나온다!!!
        else:#마지막 round(i=len(K)-1일 때)
            L=C[0:len(C)//2]
            R=C[len(C)//2:len(C)]
            #C=(f(R,K[i])^L)+R
            for j in range(len(R)):#C=R+(f(R,K[i])^L)->f와 L을 이진수로 변환한 다음 bitwise XOR연산해야 됨!!
                Lp=Lp+str(int((f(R,K[i]))[j])^int(L[j]))
            C=Lp+R#마지막 round일 때는 이번 round에 암호화된 L(i.e. Lp)과 바로 전 round에서 암호화된 R의 순서를 바꾸지 않고 그대로 왼쪽, 오른쪽에 두기
    return C

def D(C,K,f):
    K.reverse()
    return E(C,K,f)

#내가 만든 새로운 f와 새로운 P, 새로운 K로도 잘 작동하는지 확인하기(지금까진 잘 됐음!)
#def f(R,Ki):
#    R=list(R)
#    for i in range(len(Ki)):
#        R[i]=str(int(R[i])^int(Ki[i]))
#    R="".join(R)
#    return R

#P='10011101';K=['010','110','101']
#C='10110001'
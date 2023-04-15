# -*- coding: utf-8 -*-
"""
Created on Sat Nov 19 16:10:32 2022

@author: user
"""
#4bit = 16진수(4bit로 총 0부터 15까지의 수(16가지)를 나타낼 수 있으므로)
#AES-128 : P='00041214120412000C00131108231919'; K='2475A2B33475568831E2120013AA5487'; C=='BC028BD3E0E3B195550D6DF8E6F18241'
#AES-192 : P='00041214120412000C00131108231919'; K='2475A2B33475568831E2120013AA54872475A2B334755688'; C=='DAA4DEC0D4FF4F070438BF784B099520'
#AES-256 : P='00041214120412000C00131108231919'; K='2475A2B33475568831E2120013AA54872475A2B33475568831E2120013AA5487'; C=='79CC54F82FE8DFCA1A9F8A4723E5C02C'
#지금 입력값 P,K,C를 16진수로 주어져서 그렇지,
#이진코드로 바꾸면 AES-128은 128bit(16진수로는 32자리)짜리 길이의 key를 쓰고,
#                AES-192은 192bit(16진수로는 48자리)짜리 길이의 key를 쓰고,
#                AES-256은 256bit(16진수로는 64자리)짜리 길이의 key를 쓰고,

############################알고리즘 코딩 시작###########################
# AES에서 쓰이는 각 함수와 알고리즘은 bit단위로 연산이 이뤄지기 때문에,

#1. 주 함수 6가지 안에서 쓰이는 재료들 먼저 정의(RotWord, subbyte, invsubbyte, Subword)
#2. RoundKeyGenerator정의 -> BlockToState정의 -> AddRoundKey정의 -> Subbytes정의 -> ShiftRows정의 -> MixColumns정의
#3. AES코딩
#(알고리즘 처음에 P,K를 이진코드로 바꾸고 시작하기(BlockToState함수는 P를 처음 입력값으로 받을 때, 128bit단위로 받기 때문))
#(RoundKeyGenerator도 K를 이진코드로 입력받아 출력값을 이진코드의 4개의 word단위(총128bit)로 끊어 Key를 만들기 때문)

################################
#1.
#Rotword(32bit의 가장 왼쪽 한개 byte를 왼쪽 순환이동)(word는 4byte짜리 단위임)
def Rotword(word):#문자열 1byte씩 슬라이싱 후, 순환이동
    word1 = word[0:8];
    word2 = word[8:16];
    word3 = word[16:24];
    word4 = word[24:];
    word = word2+word3+word4+word1
    return word

#subbyte, invsubbyte
#subbyte, invsubbyte 계산이 모두 되어있는 table(16x16)이 존재하므로,
# 두 table을 각각 matrix로 만들어 입력값마다 대응되게끔 하기
import sympy as sp
def subbyte(byte):#1byte(8bit)짜리 이진코드 문자열을 입력값으로 받음
#hex(int(subbyte(str(bin((int('6D',16))))),2))이러면, '6D'에 대응되는 수가 나오지 않는다.(bin은 '0b'가 앞에 붙어있고, 128보다 작은 이진수는 자리수가 적기 때문 -> hex(int(subbyte(str('0'+bin((int('6D',16)))[2:])),2))을 입력하면, 맨 앞의 '0b'빼고, '0'을 추가하여 여덟자리 만들 수 있음
    subbyte_table = sp.Matrix([[0x63,0x7C,0x77,0x7B,0xF2,0x6B,0x6F,0xC5,0x30,0x01,0x67,0x2B,0xFE,0xD7,0xAB,0x76],
                               [0xCA,0x82,0xC9,0x7D,0xFA,0x59,0x47,0xF0,0xAD,0xD4,0xA2,0xAF,0x9C,0xA4,0x72,0xC0],
                               [0xB7,0xFD,0x93,0x26,0x36,0x3F,0xF7,0xCC,0x34,0xA5,0xE5,0xF1,0x71,0xD8,0x31,0x15],
                               [0x04,0xC7,0x23,0xC3,0x18,0x96,0x05,0x9A,0x07,0x12,0x80,0xE2,0xEB,0x27,0xB2,0x75],
                               [0x09,0x83,0x2C,0x1A,0x1B,0x6E,0x5A,0xA0,0x52,0x3B,0xD6,0xB3,0x29,0xE3,0x2F,0x84],
                               [0x53,0xD1,0x00,0xED,0x20,0xFC,0xB1,0x5B,0x6A,0xCB,0xBE,0x39,0x4A,0x4C,0x58,0xCF],
                               [0xD0,0xEF,0xAA,0xFB,0x43,0x4D,0x33,0x85,0x45,0xF9,0x02,0x7F,0x50,0x3C,0x9F,0xA8],
                               [0x51,0xA3,0x40,0x8F,0x92,0x9D,0x38,0xF5,0xBC,0xB6,0xDA,0x21,0x10,0xFF,0xF3,0xD2],
                               [0xCD,0x0C,0x13,0xEC,0x5F,0x97,0x44,0x17,0xC4,0xA7,0x7E,0x3D,0x64,0x5D,0x19,0x73],
                               [0x60,0x81,0x4F,0xDC,0x22,0x2A,0x90,0x88,0x46,0xEE,0xB8,0x14,0xDE,0x5E,0x0B,0xDB],
                               [0xE0,0x32,0x3A,0x0A,0x49,0x06,0x24,0x5C,0xC2,0xD3,0xAC,0x62,0x91,0x95,0xE4,0x79],
                               [0xE7,0xC8,0x37,0x6D,0x8D,0xD5,0x4E,0xA9,0x6C,0x56,0xF4,0xEA,0x65,0x7A,0xAE,0x08],
                               [0xBA,0x78,0x25,0x2E,0x1C,0xA6,0xB4,0xC6,0xE8,0xDD,0x74,0x1F,0x4B,0xBD,0x8B,0x8A],
                               [0x70,0x3E,0xB5,0x66,0x48,0x03,0xF6,0x0E,0x61,0x35,0x57,0xB9,0x86,0xC1,0x1D,0x9E],
                               [0xE1,0xF8,0x98,0x11,0x69,0xD9,0x8E,0x94,0x9B,0x1E,0x87,0xE9,0xCE,0x55,0x28,0xDF],
                               [0x8C,0xA1,0x89,0x0D,0xBF,0xE6,0x42,0x68,0x41,0x99,0x2D,0x0F,0xB0,0x54,0xBB,0x16]])
    t=bin(subbyte_table[int(byte[0:4],2),int(byte[4:],2)])#subbyte정의 자체는 1byte 이진 코드 문자열 출력이므로, 이진 코드 출력되게 한 것-subbyte_table은 모두 entry들이 int형태이다.
    t=t[2:]#맨 앞의 0b는 제거해주기
    if len(t)==8:#subbyte_table의 수가 이진수로 변환된 값이 t에 저장되어 있기에, 자릿수가 8보다 적을 수가 있다.
        t=t;#그러므로,8보다 자릿수가 작은 이진수의 경우, 8bit 맞춰줘서 내보내기
    else:
        t='0'*(8-len(t))+t
    return t#subbyte(str('0'+bin((int('6D',16)))[2:]))를 입력하면, '6D'에 대응되는 이진코드 출력, hex(int(subbyte(str('0'+bin((int('6D',16)))[2:])),2))는 6D에 대응되는 table의 수 3C를 그대로 출력해줌!
            #subbyte(str(bin((int('B1',16)))[2:]))를 입력하면, 'B1'에 대응되는 이진코드 출력, hex(int(subbyte(str(bin((int('B1',16)))[2:])),2))는 B1에 대응되는 table의 수 c8을 그대로 출력해준다.
 
def invsubbyte(byte):
    invsubbyte_table = sp.Matrix([[0x52,0x09,0x6A,0xD5,0x30,0x36,0xA5,0x38,0xBF,0x40,0xA3,0x9E,0x81,0xF3,0xD7,0xFB],
                                  [0x7C,0xE3,0x39,0x82,0x9B,0x2F,0xFF,0x87,0x34,0x8E,0x43,0x44,0xC4,0xDE,0xE9,0xCB],
                                  [0x54,0x7B,0x94,0x32,0xA6,0xC2,0x23,0x3D,0xEE,0x4C,0x95,0x0B,0x42,0xFA,0xC3,0x4E],
                                  [0x08,0x2E,0xA1,0x66,0x28,0xD9,0x24,0xB2,0x76,0x5B,0xA2,0x49,0x6D,0x8B,0xD1,0x25],
                                  [0x72,0xF8,0xF6,0x64,0x86,0x68,0x98,0x16,0xD4,0xA4,0x5C,0xCC,0x5D,0x65,0xB6,0x92],
                                  [0x6C,0x70,0x48,0x50,0xFD,0xED,0xB9,0xDA,0x5E,0x15,0x46,0x57,0xA7,0x8D,0x9D,0x84],
                                  [0x90,0xD8,0xAB,0x00,0x8C,0xBC,0xD3,0x0A,0xF7,0xE4,0x58,0x05,0xB8,0xB3,0x45,0x06],
                                  [0xD0,0x2C,0x1E,0x8F,0xCA,0x3F,0x0F,0x02,0xC1,0xAF,0xBD,0x03,0x01,0x13,0x8A,0x6B],
                                  [0x3A,0x91,0x11,0x41,0x4F,0x67,0xDC,0xEA,0x97,0xF2,0xCF,0xCE,0xF0,0xB4,0xE6,0x73],
                                  [0x96,0xAC,0x74,0x22,0xE7,0xAD,0x35,0x85,0xE2,0xF9,0x37,0xE8,0x1C,0x75,0xDF,0x6E],
                                  [0x47,0xF1,0x1A,0x71,0x1D,0x29,0xC5,0x89,0x6F,0xB7,0x62,0x0E,0xAA,0x18,0xBE,0x1B],
                                  [0xFC,0x56,0x3E,0x4B,0xC6,0xD2,0x79,0x20,0x9A,0xDB,0xC0,0xFE,0x78,0xCD,0x5A,0xF4],
                                  [0x1F,0xDD,0xA8,0x33,0x88,0x07,0xC7,0x31,0xB1,0x12,0x10,0x59,0x27,0x80,0xEC,0x5F],
                                  [0x60,0x51,0x7F,0xA9,0x19,0xB5,0x4A,0x0D,0x2D,0xE5,0x7A,0x9F,0x93,0xC9,0x9C,0xEF],
                                  [0xA0,0xE0,0x3B,0x4D,0xAE,0x2A,0xF5,0xB0,0xC8,0xEB,0xBB,0x3C,0x83,0x53,0x99,0x61],
                                  [0x17,0x2B,0x04,0x7E,0xBA,0x77,0xD6,0x26,0xE1,0x69,0x14,0x63,0x55,0x21,0x0C,0x7D]])
    s=bin(invsubbyte_table[int(byte[0:4],2),int(byte[4:],2)])
    s=s[2:]
    if len(s)==8:
        s=s;
    else:
        s='0'*(8-len(s))+s
    return s#invsubbyte(str(bin((int('A9',16)))[2:]))를 입력하면, 'A9'에 대응되는 invsubbyte_table의 정수가 이진수로 출력
            #hex(int(invsubbyte(str(bin((int('A9',16)))[2:])),2))를 입력하면, 'A9'에 대응되는 16진수 B7이 잘 출력된다!

#Subword
def Subword(word):#32bit짜리 문자열 이진코드를 입력값으로 받고 똑같은 길이 출력
    word1 = word[0:8];
    word2 = word[8:16];
    word3 = word[16:24];
    word4 = word[24:];
    word=subbyte(word1)+subbyte(word2)+subbyte(word3)+subbyte(word4)
    return word

#2. 
#RoundKeyGenerator 만들기 전에, Roundconstant먼저 정의하기(AES에서는 기약다항식을 x^8+x^4+x^3+x+1을 사용하기로 약속했음)
from copy import deepcopy
#RoundKeyGenerator
def RoundKeyGenerator(K,Nr):#이진 코드로 받게 하기,Nr=round수(AES에서 이 함수 실행 전에 Nr 결정됨!!)
    #Roundconstant먼저 정의!
    Roundconstant=[bin(int('01000000',16)),bin(int('02000000',16)),bin(int('04000000',16)),bin(int('08000000',16)),bin(int('10000000',16)),bin(int('20000000',16)),bin(int('40000000',16)),bin(int('80000000',16)),bin(int('1B000000',16)),bin(int('36000000',16))]
    for i in range(10):
        Roundconstant[i]=Roundconstant[i][2:]#각각 0b는 빼주기
        #for i in range(10):print(len(Roundconstant[i]))해보면, 모든 원소의 길이가 32bit보다 짧다는 것을 알 수 있다. RoundKeyGenerator에서 각각이 32bit 이진수 문자열이어야 연산이 가능하므로, 32bit맞춰주기
    #Roundconstant(이진코드 형태) 길이 32bit로 맞추기
    for i in range(10):
        if len(Roundconstant[i])==32:
            Roundconstant[i]=Roundconstant[i]
        else:#32보다 길이가 짧은 경우밖에 없으니, else만으로 충분!!
            Roundconstant[i]='0'*(32-len(Roundconstant[i]))+Roundconstant[i]
#    for i in range(10):Roundconstant 잘 만들어졌는지 확인차 쓴 것!
#        print(len(Roundconstant[i]))
#    for i in range(10):
#        print(Roundconstant[i])
    #RoundKey(총 Nr+1개) 만들기
    w=[];  #각 RoundKey는 무조건 w[0]~w[4*Nr+1]까지의 word를 4개씩 묶은 128bit가 되므로, AES버전에 따라 w크기가 결정된다. 
    if len(K) == 128:#Nr=10
        Nk=4
        w=list(range(44))#w초기값 설정(w길이는 총 4*Nr+4이다)
    elif len(K) == 192:#Nr=12
        Nk=6
        w=list(range(52))
    elif len(K) == 256:#key는 AES 안에서 전처리 되고 난 후에 이 Generator에 들어오므로 무조건 3개중 하나!
        Nk=8#Nr=14
        w=list(range(60))
    
    Kp=deepcopy(K)
    for i in range(Nk):#print(w)해보고,
    #AES-128일 땐, K==w[0]+w[1]+w[2]+w[3], #AES-192일 땐, K==w[0]+w[1]+w[2]+w[3]+w[4]+w[5], # AES-256일 땐, K==w[0]+w[1]+w[2]+w[3]+w[4]+w[5]+w[6]+w[7]모두 true라는 것을 위의 예시로부터 확인 가능!
    #이 때, 새로운 K를 16진수 예시로 받을 때마다 당연히, AES E함수 첫부분에서 K를 2진수로 바꾼 다음에 하기!
        w[i]=Kp[:32]
        Kp=Kp[32:]
    
    for i in range(Nk,4*Nr+4):
        if i%Nk==0:
            w[i]=bin(((int(Subword(Rotword(w[i-1])),2))^int(Roundconstant[i//Nk-1],2))^int(w[i-Nk],2))
            w[i]=w[i][2:]#정수끼리 Xor연산 시 이 의미는 각 정수의 이진수에서 bitwise XOR연산 하고 그 결과를 다시 정수의 형태로 돌려주는 것이기에, 각 항이 이진코드 문자열이니, 각각 int형으로 바꿔줘서 XOR 연산 후, 마지막 bin만 씌워주면 됨!(그럼, 이진코드를 문자열로 받음!)
            if len(w[i])<32:
                w[i]='0'*(32-len(w[i]))+w[i]#int에서 bitwise연산 후 결과가 int인 상태에서 bin(무조건 가장 큰 자릿수는 1이 되니,)으로 바꿨기에, 맨앞의 '0b'를 빼더라고, 자릿수가 32보다 적을 수 있다. 적은 만큼 0으로 채워 32bit맞추기!
            
        elif Nk==8 and (i-4)%Nk==0:
            w[i]=bin(int(Subword(w[i-1]),2)^int(w[i-Nk],2))
            w[i]=w[i][2:]
            if len(w[i])<32:
                w[i]='0'*(32-len(w[i]))+w[i]
            
        else:
            w[i]=bin(int(w[i-1],2)^int(w[i-Nk],2))
            w[i]=w[i][2:]
            if len(w[i])<32:
                w[i]='0'*(32-len(w[i]))+w[i]
                
    Roundkey=list(range(Nr+1));#Roundkey들 초기값 설정 후, 다음 줄부터 바로 바꿀 것
    for i in range(Nr+1):
        Roundkey[i]=w[4*i]+w[4*i+1]+w[4*i+2]+w[4*i+3]
    return Roundkey
            
#BlockToState, StateToBlock(BlockToState의 역함수임)
def BlockToState(B):#입력값으로 128bit짜리 이진코드를 받는데 이는 같은 말로, 8byte짜리 block임.
    B_bytes=sp.Matrix([list(range(16))]);
    for i in range(16):#1byte단위로 끊어야되니, 8자리씩 끊기
        B_bytes[i]=int(B[:8],2)#이진코드 문자형은 sympy의 Matrix에서는 지원하지 않고, 정수형만 받으니, 이 이진코드를 정수형으로 바꾸기!
        B=B[8:]
    State=B_bytes.reshape(4,4)
    State=State.transpose()
    return State

def StateToBlock(C):#BlockToState의 inverse이므로, 입력값으로 4x4행렬(sympy용)을 받고, 출력값으로 128bit 문자열 이진코드를 받음
    Block=''
    C_Block=C.transpose()
    C_Block=C_Block.reshape(1,16)
    C_Block=list(C_Block)#C_Block이 list가 되었으므로, 문자열 자료형 받을 수 있음!
    for i in range(16):
        C_Block[i]=bin(C_Block[i])
        C_Block[i]=C_Block[i][2:]
        if len(C_Block[i])<8:
            C_Block[i]='0'*(8-len(C_Block[i]))+C_Block[i]
        Block=Block+C_Block[i]
    return Block

#AddRoundKey
def AddRoundKey(S,Sp):#state들을 입력값으로 갖는 함수들의 state는 모두 entry가 int임을 염두에 두고 프로그래밍!(by sympy)
    #S, Sp모두 4x4 matrix
    T=sp.ones(4,4)#출력값의 초깃값 설정한 것!
    for i in range(4):
        for j in range(4):
            T[i,j]=S[i,j]^Sp[i,j]
    return T
# AddRoundKey 오류 없는지 확인해보기 : S=sp.Matrix([[23,45,22,12],[11,25,34,23],[11,9,6,8],[12,3,33,23]]);Sp=sp.Matrix([[78,92,30,5],[55,34,25,21],[4,98,95,65],[36,23,4,5]])

#Subbytes, InvSubbytes
def Subbytes(S):#state형태의 출력값도 모두 entry가 int이기에, 만약 여기서 entry를 일정 bit수의 이진 문자열 코드로 뽑아내야 한다면, bin을 entry에 씌운후, 원하는 bit수자리보다 부족하면, 그만큼의 '0'을 맨 앞에다가 붙여넣기!
    for i in range(4):
        for j in range(4):
            t=bin(S[i,j])
            t=t[2:]#이진코드의 맨앞 '0b'는 무조건 빼주기!!!
            if len(t)<8:
                t='0'*(8-len(t))+t
            S[i,j]=int(subbyte(t),2)#subbyte의 출력값은 8bit(=1byte)짜리 문자열 이진코드이니, int(,2)로 정수형으로 바꿔야 sympy의 Matrix에 입력 가능
    return S

def invSubbytes(S):
    for i in range(4):
        for j in range(4):
            t=bin(S[i,j])
            t=t[2:]#이진코드의 맨앞 '0b'는 무조건 빼주기!!!
            if len(t)<8:
                t='0'*(8-len(t))+t
            S[i,j]=int(invsubbyte(t),2)#invsubbyte의 출력값은 8bit(=1byte)짜리 문자열 이진코드이니, int(,2)로 정수형으로 바꿔야 sympy의 Matrix에 입력 가능
    return S

#ShiftRows, InvShiftRows
def ShiftRows(S):
    T=sp.ones(4,4)#출력값의 초기값
    for i in range(4):#행
        for j in range(4):#열(각 행, 열의 index는 list의 index처럼 음수 index를 가질 수 있음을 알기(순환성!!!))
            T[i,j-i]=S[i,j]#i번째 행은 각 entry들이 좌측으로 i칸씩 순환이동
    return T

def InvShiftRows(S):
    T=sp.ones(4,4)#출력값의 초기값
    for i in range(4):#행
        for j in range(4):#열(각 행, 열의 index는 list의 index처럼 음수 index를 가질 수 있음을 알기(순환성!!!))but 행, 열의 길이보다 큰 index는 오류남! -> 순환성에 맞추기 위해선, %쓰기! 
            T[i,(j+i)%4]=S[i,j]#i번째 행은 각 entry들이 좌측으로 i칸씩 순환이동
    return T

#MixColumns, InvMixColumns
def MixColumns(S):#입력값 S는 state로, 각 entry들은 8bit이진코드가 정수형의 형태로 저장되어 있음
    import sympy as sp
    from sympy.abc import x
    f=sp.Poly(x**8+x**4+x**3+x+1,x,modulus=2)#Z2[x]상에서의 기약 다항식으로, AES에서는 이 기약다항식을 쓰기로 약속했음. 
    C=sp.Matrix([[0x02,0x03,0x01,0x01],
                 [0x01,0x02,0x03,0x01],
                 [0x01,0x01,0x02,0x03],
                 [0x03,0x01,0x01,0x02]])
    T=sp.ones(4,4)
    for i in range(4):#T의 i번째 행
        for j in range(4):#T의 j번째 열
            c=list(range(4))#C[i,k]의 8bit 이진코드 형태의 문자열을 c[k]에 갖고있는 list -> 이 이진코드들이 대응되는 Z2[x]상의 다항식으로 변화할 것
            s=list(range(4))#S[k,j]의 8bit 이진코드 형태의 문자열을 s[k]에 갖고 있는 list -> 이 이진코드들이 대응되는 Z2[x]상의 다항식으로 변화할 것
            
            for k in range(4):#C*S를 할 때, T[i,j]=sigma(C[i,k]*S[k,j],k=0~3)하는 부분
                c[k]=bin(C[i,k])[2:]#위처럼 C를 hex로 입력해도, sympymatrix에서는 무조건 정수형으로 저장되기에, 그냥 이진수로 변환해주기;bin으로 바꾸면, 무조건 처음에 '0b'가 붙는데, 그 부분 제거
                if len(c[k])<8:#모두 길이가 숫자 크기는 변하지 않으면서, 8bit가 되도록 맞추기
                    c[k]='0'*(8-len(c[k]))+c[k]
                
                #S에 대해서도 똑같이 해주기
                s[k]=bin(S[k,j])[2:]
                if len(s[k])<8:
                    s[k]='0'*(8-len(s[k]))+s[k]
                
                h_ck=list(range(8))#c[k] 이진코드의 각 자릿수들에 대응되는 단항식을 순서대로 갖는 list
                h_sk=list(range(8))#s[k] 이진코드의 각 자릿수들에 대응되는 단항식을 순서대로 갖는 list
                for l in range(8):#c[k], s[k]의 8bit 이진코드들을 다항식으로 변화(H_ck,h_sk는 각각 x**7부터 x**0까지 하나씩만 갖고 있음)
                    h_ck[l]=sp.Poly(int(c[k][l])*x**(7-l),x,modulus=2)
                    h_sk[l]=sp.Poly(int(s[k][l])*x**(7-l),x,modulus=2)
                
                c[k]=sp.Poly(0,x,modulus=2)#c[k],s[k] 초기화->이후, 이전 이진코드에 대응되는 다항식으로 변형할 것
                s[k]=sp.Poly(0,x,modulus=2)
                
                for l in range(8):#각 c[k],s[k]를 Z2[x]의 다항식으로 변환 완료
                    c[k]+=h_ck[l]
                    s[k]+=h_sk[l]
            
            #c, s의 원소가 모두 다항식으로 변환 완료 -> 행렬 곱 연산 시작(entry별 곱셈은 다항식 mod 연산, 덧셈은 bitwise연산이지만, 이것도 다항식 덧셈과 같은 의미이니 상관x) 
            t=(c[0]*s[0]%f)+(c[1]*s[1]%f)+(c[2]*s[2]%f)+(c[3]*s[3]%f)#sigma(C[i,k]*S[k,j],k=0~3) 다항식 version으로 하는 부분
            t=t.all_coeffs()#다항식의 각 단항식의 계수를 최고차항부터 순서대로 1개씩을 원소로 갖는 list가 t가 됨.(그냥 int형태이니 각 원소를 str로 변환)
            
            for n in range(len(t)):
                t[n]=str(t[n])
            t="".join(t)#합치기까지 이진코드 문자열이 되지만, x**7의 계수가 0이면, 8bit가 안 되므로, 8bit 맞춰주기-어차피 T[i,j]로 들어갈 때는 int가 되버리니 필요없을 수도?
            if len(t)<8:
                t='0'*(8-len(t))+t
            T[i,j]=int(t,2)
    return T    

def invMixColumns(S):
    import sympy as sp
    from sympy.abc import x
    f=sp.Poly(x**8+x**4+x**3+x+1,x,modulus=2)#Z2[x]상에서의 기약 다항식으로, AES에서는 이 기약다항식을 쓰기로 약속했음. 
    C=sp.Matrix([[0x0E,0x0B,0x0D,0x09],#Mixcolumns의 C의 역행렬임.
                 [0x09,0x0E,0x0B,0x0D],
                 [0x0D,0x09,0x0E,0x0B],
                 [0x0B,0x0D,0x09,0x0E]])
    T=sp.ones(4,4)
    for i in range(4):#T의 i번째 행
        for j in range(4):#T의 j번째 열
            c=list(range(4))#C[i,k]의 8bit 이진코드 형태의 문자열을 c[k]에 갖고있는 list -> 이 이진코드들이 대응되는 Z2[x]상의 다항식으로 변화할 것
            s=list(range(4))#S[k,j]의 8bit 이진코드 형태의 문자열을 s[k]에 갖고 있는 list -> 이 이진코드들이 대응되는 Z2[x]상의 다항식으로 변화할 것
            
            for k in range(4):#C*S를 할 때, T[i,j]=sigma(C[i,k]*S[k,j],k=0~3)하는 부분
                c[k]=bin(C[i,k])[2:]#위처럼 C를 hex로 입력해도, sympymatrix에서는 무조건 정수형으로 저장되기에, 그냥 이진수로 변환해주기;bin으로 바꾸면, 무조건 처음에 '0b'가 붙는데, 그 부분 제거
                if len(c[k])<8:#모두 길이가 숫자 크기는 변하지 않으면서, 8bit가 되도록 맞추기
                    c[k]='0'*(8-len(c[k]))+c[k]
                
                #S에 대해서도 똑같이 해주기
                s[k]=bin(S[k,j])[2:]
                if len(s[k])<8:
                    s[k]='0'*(8-len(s[k]))+s[k]
                
                h_ck=list(range(8))#c[k] 이진코드의 각 자릿수들에 대응되는 단항식을 순서대로 갖는 list
                h_sk=list(range(8))#s[k] 이진코드의 각 자릿수들에 대응되는 단항식을 순서대로 갖는 list
                for l in range(8):#c[k], s[k]의 8bit 이진코드들을 다항식으로 변화(H_ck,h_sk는 각각 x**7부터 x**0까지 하나씩만 갖고 있음)
                    h_ck[l]=sp.Poly(int(c[k][l])*x**(7-l),x,modulus=2)
                    h_sk[l]=sp.Poly(int(s[k][l])*x**(7-l),x,modulus=2)
                
                c[k]=sp.Poly(0,x,modulus=2)#c[k],s[k] 초기화->이후, 이전 이진코드에 대응되는 다항식으로 변형할 것
                s[k]=sp.Poly(0,x,modulus=2)
                
                for l in range(8):#각 c[k],s[k]를 Z2[x]의 다항식으로 변환 완료
                    c[k]+=h_ck[l]
                    s[k]+=h_sk[l]
            
            #c, s의 원소가 모두 다항식으로 변환 완료 -> 행렬 곱 연산 시작(entry별 곱셈은 다항식 mod 연산, 덧셈은 bitwise연산이지만, 이것도 다항식 덧셈과 같은 의미이니 상관x) 
            t=(c[0]*s[0]%f)+(c[1]*s[1]%f)+(c[2]*s[2]%f)+(c[3]*s[3]%f)#sigma(C[i,k]*S[k,j],k=0~3) 다항식 version으로 하는 부분
            t=t.all_coeffs()#다항식의 각 단항식의 계수를 최고차항부터 순서대로 1개씩을 원소로 갖는 list가 t가 됨.(그냥 int형태이니 각 원소를 str로 변환)
            
            for n in range(len(t)):
                t[n]=str(t[n])
            t="".join(t)#합치기까지 이진코드 문자열이 되지만, x**7의 계수가 0이면, 8bit가 안 되므로, 8bit 맞춰주기-어차피 T[i,j]로 들어갈 때는 int가 되버리니 필요없을 수도?
            if len(t)<8:
                t='0'*(8-len(t))+t
            T[i,j]=int(t,2)
    return T    
       
#####################################################실제 AES에 복붙x


#P,K이진코드로 변환(맨 위 3가지 예시중 하나 고르기)-AES초반에 적용
#K=int(K,16)#16진수로 표현된 문자열을 10진수 정수형의 int로 변환
#K=bin(K)#10진수로 표현된 K를 이진수로 변환(bin은 정수형만 입력값으로 받기에 이처럼 16진수 문자열을 int로 바꾼 후, bin쓰기!)
#K=str(K)
#P도 같은 방식으로!!!
#P=str(bin(int(P,16)))


######################################################
#4. AES
#이진코드문자열로 변환된 K, P가 128bit or 192 or 256bit보다 짧을 수도 있기 때문에
#AES코딩시 처음에 K, P길이부터 각각 맞춰줘야 된다!#'0'필요한 만큼 K,P문자열 앞쪽에 붙여주면, K,P자릿수만 변하고, 실제 값은 차이x
def E(P,K):#P,K는 16진수의 문자열로 입력하기
    #P,K이진코드로 변환(맨 위 3가지 예시중 하나 고르기)
    K=int(K,16)#16진수로 표현된 문자열을 10진수 정수형의 int로 변환
    K=bin(K)#10진수로 표현된 K를 이진수로 변환(bin은 정수형만 입력값으로 받기에 이처럼 16진수 문자열을 int로 바꾼 후, bin쓰기!)
    K=str(K)
    #P도 같은 방식으로!!!
    P=str(bin(int(P,16)))
    #K,P를 이진코드로 변환하면, 앞에 '0b가 붙어있기에, 떼버리기!'
    K=K[2:]
    P=P[2:]
    if len(P)==128:
        P=P;
    else:#무조건 P의 길이는 128bit보다 작거나 같으므로, else만으로 충분
        P='0'*(128-len(P))+P
    if len(K)<=128:#AES-128에 해당, 이때, AES버전도 결정되므로, round수 정의 하기!
        K='0'*(128-len(K))+K
        Nr=10#round수
    elif len(K)<=192:#AES-192에 해당(바로 위의 if문을 거치고 오므로 128<len(K)=<192에 해당)
        K='0'*(192-len(K))+K
        Nr=12
    else:#AES-256에 해당
        K='0'*(256-len(K))+K
        Nr=14
##################################이진 코드 변환 끝->RoundKeyGenerator부터 시작!
    K=RoundKeyGenerator(K, Nr)#Nr+1개의 roundkey 생성
    state=BlockToState(P)
    state=AddRoundKey(state, BlockToState(K[0]))
    for i in range(1,Nr+1):
        state=Subbytes(state)
        state=ShiftRows(state)
        if i!=Nr:
            state=MixColumns(state)
        state=AddRoundKey(state, BlockToState(K[i]))
    C=StateToBlock(state)
    C=hex(int(C,2))
    C=C[2:].upper()
    return C

def D(C,K):
    K=int(K,16)#16진수로 표현된 문자열을 10진수 정수형의 int로 변환
    K=bin(K)#10진수로 표현된 K를 이진수로 변환(bin은 정수형만 입력값으로 받기에 이처럼 16진수 문자열을 int로 바꾼 후, bin쓰기!)
    K=str(K)
    #C도 같은 방식으로!!!
    C=str(bin(int(C,16)))
    #K,C를 이진코드로 변환하면, 앞에 '0b가 붙어있기에, 떼버리기!'
    K=K[2:]
    C=C[2:]
    if len(C)==128:
        C=C;
    else:#무조건 P의 길이는 128bit보다 작거나 같으므로, else만으로 충분
        C='0'*(128-len(C))+C
    if len(K)<=128:#AES-128에 해당, 이때, AES버전도 결정되므로, round수 정의 하기!
        K='0'*(128-len(K))+K
        Nr=10#round수
    elif len(K)<=192:#AES-192에 해당(바로 위의 if문을 거치고 오므로 128<len(K)=<192에 해당)
        K='0'*(192-len(K))+K
        Nr=12
    else:#AES-256에 해당
        K='0'*(256-len(K))+K
        Nr=14
##################################이진 코드 변환 끝->RoundKeyGenerator부터 시작!
    K=RoundKeyGenerator(K, Nr)#Nr+1개의 roundkey 생성
    state=BlockToState(C)
    state=AddRoundKey(state, BlockToState(K[Nr]))
    for i in list(range(Nr-1,-1,-1)):
        state=InvShiftRows(state)
        state=invSubbytes(state)
        state=AddRoundKey(state,BlockToState(K[i]))
        if i!=0:
            state=invMixColumns(state)
    P=StateToBlock(state)#128bit 정확하게 생성됨!
    P=hex(int(P,2))
    P=P[2:].upper()#but, 이진수 크기에 따라, 16진수로 변환시, 자릿수가 32bit보다 작을 수 있다.
    if len(P)<32:
        P='0'*(32-len(P))+P    
    return P
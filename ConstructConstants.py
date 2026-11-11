from matlab_helpers import Struct


def ConstructConstants(ConstructName):
    """
    190116 Steve Yeo
    Output ConstructConst constants like:
    hairpin : Is the contsruct hairpin or not ( 1: T , 0: F )
    dna1label = name of dna1 label 'Hairpin Closed Model';
    ssDNA1 : length of ssDNA1 in closed configuration
    dsDNA1 : length of dsDNA1 in closed configuration
    dna2label = name of dna2 label 'Hairpin Closed Model';
    ssDNA2 : length of ssDNA in open configuration
    dsDNA2 : length of dsDNA in open configuration
    fitoffset : ???
    fitclosed : ???
    fitrange : range of force to fit model and FEC
    fakeoffset : ???
    stdF: ???
    stdF2: ???
    """
    C = Struct()

    if ConstructName == '35dT-SSB-gapped':
        C.hairpin = 0
        C.dna1label = 'Hairpin Closed Model'
        C.ssDNA1 = 35  # 4dT spacer and 10dT tail;
        C.dsDNA1 = 1710+1555  # LH % Other half of the dsDNA is in hairpin
        # ConstructConst.dna2label = 'Hairpin Open Model';
        # ConstructConst.ssDNA2 = 0; % Hairpin is open.
        # ConstructConst.dsDNA2 = 0; % LH
        C.fitoffset = 1
        # fitrange = 100:150;
        C.fitclosed = 1
        C.fitrange = [2, 10]  # Force in pN
        C.fakeoffset = 0
        C.stdF = [2, 10]
        C.stdF2 = [15, 16]
        C.fkorhp = 'hairpin'
    elif ConstructName == '5Fork-Hairpin-10dT':
        C.hairpin = 1
        C.dna1label = 'Hairpin Closed Model'
        C.ssDNA1 = 4+10  # 4dT spacer and 10dT tail;
        C.dsDNA1 = 1710  # LH % Other half of the dsDNA is in hairpin
        C.dna2label = 'Hairpin Open Model'
        C.ssDNA2 = 4+10+1719*2  # Hairpin is open.
        C.dsDNA2 = 1710  # LH
        C.fitoffset = 1
        # fitrange = 100:150;
        C.fitclosed = 1
        C.fitrange = [2, 8]  # Force in pN
        C.fakeoffset = 0
        C.stdF = [2, 10]
        C.stdF2 = [15, 16]
        C.fkorhp = 'hairpin'
    elif ConstructName == '5Fork-Hairpin-30dT':
        C.hairpin = 1
        C.dna1label = 'Hairpin Closed Model'
        C.ssDNA1 = 4+30  # 4dT spacer and 20dT tail;
        C.dsDNA1 = 1710  # LH % Other half of the dsDNA is in hairpin
        C.dna2label = 'Hairpin Open Model'
        C.ssDNA2 = 4+30+1719*2  # Hairpin is open.
        C.dsDNA2 = 1710  # LH
        C.fitoffset = 1
        # fitrange = 100:150;
        C.fitclosed = 1
        C.fitrange = [2, 8]  # Force in pN
        C.fakeoffset = 0
        C.stdF = [2, 10]
        C.stdF2 = [15, 16]
        C.fkorhp = 'hairpin'
    elif ConstructName == '5Fork-Hairpin-20dT':
        C.hairpin = 1
        C.dna1label = 'Hairpin Closed Model'
        C.ssDNA1 = 4+20  # 4dT spacer and 20dT tail;
        C.dsDNA1 = 1710  # LH % Other half of the dsDNA is in hairpin
        C.dna2label = 'Hairpin Open Model'
        C.ssDNA2 = 4+20+1719*2  # Hairpin is open.
        C.dsDNA2 = 1710  # LH
        C.fitoffset = 1
        # fitrange = 100:150;
        C.fitclosed = 1
        C.fitrange = [2, 8]  # Force in pN
        C.fakeoffset = 0
        C.stdF = [2, 10]
        C.stdF2 = [15, 16]
        C.fkorhp = 'hairpin'
    elif ConstructName == '5Fork-Hairpin-70dT':
        C.hairpin = 1
        C.dna1label = 'Hairpin Closed Model'
        C.ssDNA1 = 4+70  # 4dT spacer and 70dT tail;
        C.dsDNA1 = 1710  # LH % Other half of the dsDNA is in hairpin
        C.dna2label = 'Hairpin Open Model'
        C.ssDNA2 = 4+70+1719*2  # Hairpin is open.
        C.dsDNA2 = 1710  # LH
        C.fitoffset = 1
        # fitrange = 100:150;
        C.fitclosed = 1
        C.fitrange = [2, 8]  # Force in pN
        C.fakeoffset = 0
        C.stdF = [2, 10]
        C.stdF2 = [15, 16]
        C.fkorhp = 'hairpin'
    elif ConstructName == '5Fork-Hairpin-40dT':
        C.hairpin = 1
        C.dna1label = 'Hairpin Closed Model'
        C.ssDNA1 = 4+40  # 4dT spacer and 20dT tail;
        C.dsDNA1 = 1710  # LH % Other half of the dsDNA is in hairpin
        C.dna2label = 'Hairpin Open Model'
        C.ssDNA2 = 4+40+1719*2  # Hairpin is open.
        C.dsDNA2 = 1710  # LH
        C.fitoffset = 1
        # fitrange = 100:150;
        C.fitclosed = 1
        C.fitrange = [2, 8]  # Force in pN
        C.fakeoffset = 0
        C.stdF = [2, 10]
        C.stdF2 = [15, 16]
        C.fkorhp = 'hairpin'
    elif ConstructName == '5Fork-Hairpin-50dT':
        C.hairpin = 1
        C.dna1label = 'Hairpin Closed Model'
        C.ssDNA1 = 4+50  # 4dT spacer and 20dT tail;
        C.dsDNA1 = 1710  # LH % Other half of the dsDNA is in hairpin
        C.dna2label = 'Hairpin Open Model'
        C.ssDNA2 = 4+50+1719*2  # Hairpin is open.
        C.dsDNA2 = 1710  # LH
        C.fitoffset = 1
        # fitrange = 100:150;
        C.fitclosed = 1
        C.fitrange = [2, 8]  # Force in pN
        C.fakeoffset = 0
        C.stdF = [2, 10]
        C.stdF2 = [15, 16]
        C.fkorhp = 'hairpin'
    elif ConstructName == '5Fork-Hairpin-60dT':
        C.hairpin = 1
        C.dna1label = 'Hairpin Closed Model'
        C.ssDNA1 = 4+60  # 4dT spacer and 20dT tail;
        C.dsDNA1 = 1710  # LH % Other half of the dsDNA is in hairpin
        C.dna2label = 'Hairpin Open Model'
        C.ssDNA2 = 4+60+1719*2  # Hairpin is open.
        C.dsDNA2 = 1710  # LH
        C.fitoffset = 1
        # fitrange = 100:150;
        C.fitclosed = 1
        C.fitrange = [2, 8]  # Force in pN
        C.fakeoffset = 0
        C.stdF = [2, 10]
        C.stdF2 = [15, 16]
        C.fkorhp = 'hairpin'
    elif ConstructName == 'Hairpin-LambdaRH-dT-100-86bp':
        C.hairpin = 1
        C.dna1label = 'Hairpin Closed Model'
        C.ssDNA1 = 100  # 0;
        C.dsDNA1 = 3.265e3
        C.dna2label = 'Hairpin Open Model'
        C.dsDNA2 = 3.265e3
        C.ssDNA2 = 2*86+4+C.ssDNA1  # 2*50+4+10; %
        C.fitoffset = 1
        # fitrange = 100:150;
        C.fitclosed = 1
        # ConstructConst.fitrange = [3 9]; % Force in pN
        C.fitrange = [5, 9]  # Force in pN
        C.fakeoffset = 0
        C.stdF = [2, 10]
        C.stdF2 = [15, 16]
        C.fkorhp = 'hairpin'
    elif ConstructName == 'Hairpin-LambdaRH-dT-10-86bp':
        C.hairpin = 1
        C.dna1label = 'Hairpin Closed Model'
        C.ssDNA1 = 10  # 0;
        C.dsDNA1 = 3.265e3
        C.dna2label = 'Hairpin Open Model'
        C.dsDNA2 = 3.265e3
        C.ssDNA2 = 2*86+4+C.ssDNA1  # 2*50+4+10; %
        C.fitoffset = 1
        # fitrange = 100:150;
        C.fitclosed = 1
        C.fitrange = [4, 11]  # Force in pN
        C.fakeoffset = 0
        C.stdF = [2, 10]
        C.stdF2 = [15, 16]
    elif ConstructName == '5fork-10dT':
        C.hairpin = 0
        C.dna1label = 'Oligo Model'
        C.dsDNA1 = 3431
        C.ssDNA1 = 4
        C.fitoffset = 1
        # fitrange = 200:260;
        C.fitrange = [10, 20]  # Force in pN
        C.fakeoffset = 0
    elif ConstructName == 'T4':
        C.hairpin = 0
        # dsDNA = 1.5e3;%first try
        # dsDNA1 = 3011+29+120;%second try, ds part + overhang + stuck part
        C.dsDNA1 = 3.4e3
        C.ssDNA1 = 0
        C.fitoffset = 1
        C.fitrange = [4, 6]
        # fitrange = [2 30]; % Force in pN
        C.fakeoffset = 0
        C.stdF = [2, 15]
    elif ConstructName == 'Oligo':
        C.hairpin = 0
        C.dna1label = 'Oligo Model'
        C.dsDNA1 = 3263
        C.ssDNA1 = 35
        C.fitoffset = 1
        # fitrange = 200:260;
        C.fitrange = [2, 10]  # Force in pN
        C.fakeoffset = 0
    elif ConstructName == 'fork':
        C.hairpin = 0
        C.dna1label = 'Oligo Model'
        C.dsDNA1 = 3105
        C.ssDNA1 = 4
        C.fitoffset = 1
        # fitrange = 200:260;
        C.fitrange = [10, 20]  # Force in pN
        C.fakeoffset = 0
    elif ConstructName == 'RecA':
        C.hairpin = 0
        C.dna1label = 'RecA Model'
        C.dsDNA1 = 3263
        C.ssDNA1 = 140
        C.fitoffset = 1
        C.fitrange = [4, 10]  # Force in pN
        C.fakeoffset = 0
    elif ConstructName == 'Hairpin-RH-dT-19':
        C.hairpin = 1
        C.dna1label = 'Hairpin Closed Model'
        C.ssDNA1 = 19
        C.dsDNA1 = 3.05e3
        C.dna2label = 'Hairpin Open Model'
        C.dsDNA2 = 3.05e3
        C.ssDNA2 = 19+89*2+4
        C.fitoffset = 1
        # fitrange = 100:150;
        C.fitrange = [2, 10]  # Force in pN
        C.fakeoffset = 0
        C.stdF = [3, 12]
        C.stdF2 = [17.5, 25]
    elif ConstructName == 'Hairpin-RH-dT-38':
        C.hairpin = 1
        C.dna1label = 'Hairpin Closed Model'
        C.ssDNA1 = 38
        C.dsDNA1 = 3.05e3
        C.dna2label = 'Hairpin Open Model'
        C.dsDNA2 = 3.05e3
        C.ssDNA2 = 2*89+4+38  # 220
        C.fitoffset = 1
        # fitrange = 100:150;
        C.fitrange = [2, 10]  # Force in pN
        C.fakeoffset = 0
        C.stdF = [2, 10]
        C.stdF2 = [15, 16]
    elif ConstructName == 'Hairpin-RH-dT-10':
        C.hairpin = 1
        C.dna1label = 'Hairpin Closed Model'
        C.ssDNA1 = 10  # 0;
        C.dsDNA1 = 3.055e3
        C.dna2label = 'Hairpin Open Model'
        C.dsDNA2 = 3.055e3
        C.ssDNA2 = 2*89+4+C.ssDNA1  # 2*50+4+10; %
        C.fitoffset = 1
        # fitrange = 100:150;
        C.fitrange = [2, 6]  # Force in pN
        C.fakeoffset = 0
        C.stdF = [2, 10]
        C.stdF2 = [15, 16]
        C.fkorhp = 'hairpin'
    elif ConstructName == 'Hairpin-LambdaRH-dT-10':
        C.hairpin = 1
        C.dna1label = 'Hairpin Closed Model'
        C.ssDNA1 = 10  # 0;
        C.dsDNA1 = 3.25e3
        C.dna2label = 'Hairpin Open Model'
        C.dsDNA2 = 3.25e3
        C.ssDNA2 = 2*89+4+C.ssDNA1  # 2*50+4+10; %
        C.fitoffset = 1
        # fitrange = 100:150;
        C.fitclosed = 1
        C.fitrange = [2, 10]  # Force in pN
        C.fakeoffset = 0
        C.stdF = [2, 10]
        C.stdF2 = [15, 16]
    elif ConstructName == 'nick_test':
        C.hairpin = 1
        C.dna1label = 'Hairpin Closed Model'
        C.ssDNA1 = 1010  # 0;
        C.dsDNA1 = 2.25e3
        C.dna2label = 'Hairpin Open Model'
        C.dsDNA2 = 2.25e3
        C.ssDNA2 = 2*86+4+C.ssDNA1  # 2*50+4+10; %
        C.fitoffset = 1
        # fitrange = 100:150;
        C.fitclosed = 1
        C.fitrange = [2, 5]  # Force in pN
        C.fakeoffset = 0
        C.stdF = [2, 10]
        C.stdF2 = [15, 16]
    elif ConstructName == 'GQ':
        C.hairpin = 1
        C.dna1label = 'Hairpin Closed Model'
        C.ssDNA1 = 34  # 0;
        C.dsDNA1 = 2.26e3
        C.dna2label = 'Hairpin Open Model'
        C.dsDNA2 = 2.26e3
        C.ssDNA2 = 22+C.ssDNA1  # 2*50+4+10; %
        C.fitoffset = 1
        # fitrange = 100:150;
        C.fitclosed = 1
        C.fitrange = [3, 6]  # Force in pN
        C.fakeoffset = 0
        C.stdF = [2, 10]
        C.stdF2 = [15, 16]
    elif ConstructName == 'MarcodsDNA':
        C.hairpin = 0
        # dsDNA = 1.5e3;%first try
        C.dsDNA1 = 1.565e3
        C.ssDNA1 = 0
        C.fitoffset = 1
        # fitrange = 200:350;
        C.fitrange = [2, 30]  # Force in pN
        C.fakeoffset = 0
        C.stdF = [2, 15]
    else:
        raise ValueError(f'Unknown construct: {ConstructName}')

    return C

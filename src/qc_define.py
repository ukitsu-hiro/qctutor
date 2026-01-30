# -*- coding: utf-8 -*-

# author: Ukitsu

import numpy as np

class Define():
    # Gates の U(2x2 Unitary matrix) を定義
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)
    H = np.array([[1, 1], [1, -1]], dtype=complex)/np.sqrt(2)
    I = np.eye(2, dtype=complex)
    T   = np.array([[1, 0], [0, np.exp(1j  * np.pi/4)]], dtype=complex) # Z周り45°
    Tdg = np.array([[1, 0], [0, np.exp(-1j * np.pi/4)]], dtype=complex)
    S   = np.array([[1, 0], [0, 1j]],  dtype=complex) # Z周り90°
    Sdg = np.array([[1, 0], [0, -1j]], dtype=complex)

    # 処理用 list の定義
    # 全ゲート名と、U を関連付け
    proc_gate_list = [["I", I], ["H", H], ["X", X], ["Y", Y], ["Z", Z], ["T", T], 
                    ["T+", Tdg], ["S", S], ["S+", Sdg], ["M", I], ["CNOT", X], ["CZ", Z]]
    proc_gate_list.extend([["R"+str(i+1), np.array([[1, 0], [0,  np.exp(1j *2* np.pi/(2**(i+1)))]])]
                        for i in range(7)])
    #print("define", proc_gate_list)
    # 全ゲート名のリスト
    gate_name_list = [s for s,_ in proc_gate_list] # name のみ
    # 全制御ゲート名のリスト
    ctrl_gate_list = ["CNOT", "CZ"]
    ctrl_gate_list.extend(["R"+str(n+1) for n in range(7)]) # CNOT以外を追加
    # 全Z軸回転ゲート名のリスト
    z_axis_list = ["Z", "T", "T+", "S", "S+", "CZ"] # Bloch 球のZ軸回転
    z_axis_list.extend(["R"+str(n+1) for n in range(7)])
    #print("define", ctrl_gate_list)

    def __init__(self):
        pass
        #print("len=", len(self.proc_gate_list))
        #print("Define", self.gate_name_list)

    def get_gates_group(self, name):
        """ 回路の定義 """
        msg = ""
        gm = ["measure", ["M", "all"]]
        if name == "glover": # glover 回路
            init = "0000"
            g0 = ["init", ["H", "all"]]
            g1 = ["oracle", ["X", 1], ["H", 0], ["CNOT", 0, [3,2,1]], ["H", 0], ["X", 1]]
            g2 = ["diffusion", ["H", "all"], ["X", "all"], ["H", 0], ["CNOT", 0, [3,2,1]], ["H", 0], ["X", "all"], ["H", "all"]]
            gates = [g0, g1, g2, g1, g2, g1, g2, gm]
            msg = "4 bits, 3 loops, oracle: |13)=|1101>"
            #print(gates)
        elif name == "toffoli": # CCNOT (2 制御ビットXゲート)  
            init = "0" * 3
            g0 = ["toffoli", ["H", 0], ["CNOT", 0, 1], ["T+",0], ["CNOT", 0, 2], ["T",0],
                    ["CNOT", 0, 1], ["T+",0], ["CNOT", 0, 2], ["T",[1,0]],
                    ["CNOT", 1, 2], ["T",2], ["T+",1], ["CNOT", 1, 2], ["H",0],] # OK CCNOT
            gates = [g0, gm]
            msg = "CCNOT: 2 bits(q2,q1) controlled X gate(q0)"
        elif name == "fourier_2": # 量子離散フーリエ変換 2bit
            init = "0" * 2
            g0 = ["fourier2", ["H", 1], ["R2", 1, 0], ["H", 0]]
            g1 = ["swap", ["CNOT", 1, 0], ["CNOT", 0, 1], ["CNOT", 1, 0]]
            gates = [g0, g1, gm]
            msg = "2 bits DFT(Discrete Fourier Transform)"
        elif name == "fourier_3":
            init = "0" * 3
            g0 = ["fourier4", ["H", 2], ["R2", 2, 1], ["R3", 2, 0], 
                            ["H", 1], ["R2", 1, 0],  ["H", 0]]
            g1 = ["swap", ["CNOT", 2, 0], ["CNOT", 0, 2], ["CNOT", 2, 0]]
            gates = [g0, g1, gm]
            msg = "3 bits DFT(Discrete Fourier Transform)"
        elif name == "fourier_4":
            init = "0" * 4
            g0 = ["fourier4", ["H", 3], ["R2", 3, 2], ["R3", 3, 1], ["R4", 3, 0],
                            ["H", 2], ["R2", 2, 1], ["R3", 2, 0], 
                            ["H", 1], ["R2", 1, 0], ["H", 0]]
            g1 = ["swap", ["CNOT", 3, 0], ["CNOT", 0, 3], ["CNOT", 3, 0],
                            ["CNOT", 2, 1], ["CNOT", 1, 2], ["CNOT", 2, 1]]
            gates = [g0, g1, gm]
            msg = "4 bits DFT(Discrete Fourier Transform)"
        elif name == "swap": # 2bit 交換
            init = "00"
            g0 = ["swap", ["CNOT", 1, 0], ["CNOT", 0, 1], ["CNOT", 1, 0]]
            gates = [g0, gm]
            msg = "swap q1 and q0"
        elif name == "fredkin": # 制御付き 2bit 交換
            init = "000"
            g0 = ["fredkin", ["CNOT", 0, [2,1]], ["CNOT", 1, [2,0]], ["CNOT", 0, [2,1]]]
            gates = [g0, gm]
            msg = "q2 control swap q1 and q0"
        elif name == "and": # AND 回路: q2 = q1 and q0
            init = "110"
            g0 = ["and", ["CNOT", 2, [0,1]]]
            gates = [g0, gm]
            msg = "q1 and q0 => q2 (valid input 0-3)\n(see by matrix logical mode)"
        elif name == "xor": # XOR 回路: q2 = q1 xor q0
            init = "111"
            g0 = ["xor", ["CNOT", 2, 0], ["CNOT", 2, 1]]
            gates = [g0, gm]
            msg = "q1 xor q0 => q2 (valid input 0-3)\n(see by matrix logical mode)"
        elif name == "or":  # OR 回路: q2 = q1 or q0
            init = "000"
            g0 = ["or", ["CNOT", 2, 0], ["CNOT", 2, 1], ["CNOT", 2, [0,1]]]
            gates = [g0, gm]
            msg = "q1 or q0 => q2 (valid input 0-3)\n(see by matrix logical mode)"
        elif name == "add1+1":  # 加算回路: q3q2 = q1 + q0 OK
            init = "1100"
            g0 = ["add1+1", ["CNOT", 2, 0], ["CNOT", 2, 1], ["CNOT", 3, [0,1]]]
            gates = [g0, gm]
            msg = "q1 + q0 => q3q2 (valid input 0-3)\n(see by matrix logical mode)"
        elif name == "add1+1a":  # 加算回路: q3q2 = q1 + q0
            init = "1100"
            g0 = ["add1+1a", ["CNOT", 3, [0,1]], ["CNOT", 3, [0,2]], ["CNOT", 3, [1,2]],
                    ["CNOT", 2, 0], ["CNOT", 2, 1]]
            gates = [g0, gm]
            msg = "q1 + q0 => q3q2 (valid input 0-3)\n(see by matrix logical mode)"
        elif name == "add2+2":
            init = "0" * 6
            # q5:c, q4:s2, q3:s1=b1, q2:s0=b0, q1:a1, a0:a0
            g0 = ["add2+2", ["CNOT", 3, 1], ["CNOT", 5, 1], ["CNOT", 5, [0,2]],
                    ["CNOT", 1, [5,3]], ["X", [3,4]], ["CNOT", 4, 1],
                    ["CNOT", 3, 5], ["CNOT", 1, [3,5]], 
                    ["CNOT", 5, [0,2]], ["X", [3,4]],
                    ["CNOT", 5, 1], ["CNOT", 3, 1], ["CNOT", 2, 0], ]
            gates = [g0, gm]
            msg = "q3q2 + q1q0 => q4_q2 (valid input 0-f)\n(see by matrix logical mode)" ##???
        elif name == "add3+3": # q5_q3 + q2_q0 => q6_q3 OK
            init = "0" * 8
            # q7:c, q6:s3, q5:s2=b2, q4:s1=b1, q3:s0=b0, q2:a2, q1:a1, q0:a0
            g0 = ["add3+3", ["CNOT", 5, 2], ["CNOT", 4, 1], ["CNOT", 7, 1],
                    ["CNOT", 7, [0,3]], ["CNOT", 1, 2], ["CNOT", 1, [7,4]],
                    ["CNOT", 2, [1,5]], ["X", [4,5,6]], ["CNOT", 4, 7], ["CNOT", 5, 1],
                    ["CNOT", 6, 2], ["CNOT", 2, [1,5]], ["CNOT", 1, [7,4]],
                    ["CNOT", 7, [0,3]], ["X", [4,5,6]], ["CNOT", 1, 2],
                    ["CNOT", 7, 1], ["CNOT", 3, 0], ["CNOT", 4, 1], ["CNOT", 5, 2], ]
            gates = [g0, gm]
            msg = "q5_q3 + q2_q0 => q6_q3 (valid input 0-3f)\n(see by matrix logical mode)"
        elif name == "add3+3a":
            init = "0" * 8
            g0 = ["add3+3",
                    # MAJ(a0,b0,c)  => CNOT c->b0,  CNOT c->a0, CNOT a0,b0->c
                    ["CNOT", 3, 6], ["CNOT", 0, 6], ["CNOT", 6, [0,3]],
                    # MAJ(a1,b1,c)  => CNOT c->b1,  CNOT c->a1, CNOT a1,b1->c
                    ["CNOT", 4, 6], ["CNOT", 1, 6], ["CNOT", 6, [1,4]], 
                    # MAJ(a2,b2,c)  => CNOT c->b2,  CNOT c->a2, CNOT a2,b2->c
                    ["CNOT", 5, 6], ["CNOT", 2, 6], ["CNOT", 6, [2,5]], 
                    # set carry => CNOT c->s3
                    ["CNOT", 7, 6],
                    # UMA(a2,b2,c) => CNOT a2,b2->c, CNOT c->a2, CNOT a2->b2(s2)
                    ["CNOT", 6, [2,5]], ["CNOT", 2, 6], ["CNOT", 5, 2], 
                    # UMA(a1,b1,c) => CNOT a1,b1->c, CNOT c->a1, CNOT a1->b1(s1)
                    ["CNOT", 6, [1,4]], ["CNOT", 1, 6], ["CNOT", 4, 1], 
                    # UMA(a0,b0,c)  => CNOT a0,b0->c, CNOT c->a0, CNOT a0->b0(s0)
                    ["CNOT", 6, [0,3]], ["CNOT", 0, 6], ["CNOT", 3, 0],]
            g1 = ["swap", ["CNOT", 7, 6], ["CNOT", 6, 7], ["CNOT", 7, 6]]
            gates = [g0, g1, gm]
            msg = "q5_q3 + q2_q0 => q6_q3 (valid input 0-3f)\n(see by matrix logical mode)"
        elif name == "add3+3_10bit":
            init = "0" * 10
            g0 = ["add3+3", ["CNOT", 7, [0,3]], ["CNOT", 6, 0], ["CNOT", 6, 3],
                    ["CNOT", 8, [1,4]], ["CNOT", 8, [1,7]], ["CNOT", 8, [4,7]],
                    ["CNOT", 7, 1],     ["CNOT", 7, 4], ["CNOT", 9, [2,5]], ["CNOT", 9, [2,8]],
                    ["CNOT", 9, [5,8]], ["CNOT", 8, 2], ["CNOT", 8, 5]]
            gates = [g0, gm]
            msg = "q5_q3 + q2_q0 => q9_q6 (valid input 0-3f)\n(see by matrix logical mode)"
        elif name == "exam": # ゲート種類実験
            init = "100"
            g0 = ["exam", ["X", [1,0]], ["Y", [1,0]], ["Z", [1,0]], ["H", [1,0]], 
                ["S", [1,0]], ["S+", [1,0]], ["T", [1,0]], ["T+", [1,0]], ["CNOT", 0, [2,1]], 
                ["R2", 0, 1], ["R3", 0, 1], ["R4", 0, 1], ["R5", 0, 1], ["R6", 0, 1], ["R7", 0, 1]]
            gates = [g0, gm]
            msg = "listed up all gate types"
        else:
            init = "0" * 2
            g0 = ["exam", ["CZ", 1, 0]]
            gates = [g0, gm]
            msg = "controlled Z"
        return init, gates, name + ": " + msg

    def get_circuit_names(self):
        """ tab menu 用に回路名リストを渡す(階層化) """
        circuit_names = ["glover", "toffoli",  
                    "fourier_2", "fourier_3", "fourier_4", 
                    "swap", "fredkin", 
                    ["logical", "and", "xor", "or"], # sub menu
                    ["add", "add1+1", "add1+1a", "add2+2",
                        "add3+3","add3+3a","add3+3_10bit"], # sub menu
                    "exam", "exam1"]
        return circuit_names

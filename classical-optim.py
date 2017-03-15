#!/usr/bin/python
import re
import fileinput
import sys
import os
#import numpy as np

####################
#some constants
####################

output_file = open("optimized.qasm", "w")

GateList = {     'CNOT'     : ( 2 , 1),
                 'c-z'      : ( 2 , 1),
                 'c-x'      : ( 2 , 1),
                 'measure'  : ( 1 , 0),
                 'MeasX'    : ( 1 , 0),
                 'MeasY'    : ( 1 , 0),
                 'MeasZ'    : ( 1 , 0),
                 'dmeter'   : ( 1 , 0),
                 'h'        : ( 1 , 0),
                 'H'        : ( 1 , 0),
                 'X'        : ( 1 , 0),
                 'Y'        : ( 1 , 0),
                 'Z'        : ( 1 , 0),
                 'S'        : ( 1 , 0),
                 'T'        : ( 1 , 0),
                 'Tdag'     : ( 1 , 0),
                 'U'        : ( 1 , 0),
                 'ZZ'       : ( 2 , 0),
                 'SS'       : ( 2 , 0),
                 'zero'     : ( 1 , 0),
                 'nop'      : ( 1 , 0),
                 'PrepZ'    : ( 1 , 0),
                 'PrepY'    : ( 1 , 0),
                 'PrepX'    : ( 1 , 0),
                 'discard'  : ( 1 , 0),
                 'slash'    : ( 1 , 0),
                 'space'    : ( 1 , 0),
                 'swap'     : ( 2 , 0),
                 'Tof'      : ( 3 , 2),
                 'Utwo'     : ( 2 , 0)
                 }
######################

class Qgate:

    def __init__(self, gt, ops, ln):
        
        self.gatetype = gt #gate name
        self.ops =[op.strip(' ') for op in ops.split(',')] #qubits operand 
        self.id = ln    #line number as ID to locate the gate

        # check if the gate exists
        if not GateList.has_key(self.gatetype):
            tmp_str = (self.id, self.gatetype, ops)
            err_print("On line %d of the input file: unidentified gate type%s on $%s" % tmp_str )

        #check if the argument has the right number of bits
        if(len(self.ops) != GateList[self.gatetype][0]):
            tmp_str = (self.id, self.gatetype)
            err_print("On line %d of the input file: wrong number of qubits in %s" %s )

        #check for duplicate operands
        if (len(set(self.ops)) < len(self.ops)):
            tmp_var = (self.id, self.gatetype)
            err_print("On line %d of the input file: duplicate bit operands in %s" %s )

class qasm_parser: 

    def __init__(self, input_file):

        self.bits = []
        self.bittypes = []
        self.gates = []

        linenum = 0 # line number counting,  for error messages

        for line in input_file: 
            
            linenum += 1

            if(line[0] == '#'):
                continue

            #qubit spec 
            m = re.compile('\S').search(line)
            if(not m):
                continue
            #qubit spec 
            m = re.compile('\s*qubit\s+(\S+)').search(line)
            if(m):
                self.bits.append(m.group(1))
                self.bittypes.append(0)
                continue

            # cbit spec - syntax: cbit name
            m = re.compile('\s*cbit\s+(\S+)').search(line)
            if(m):
                self.bits.append(m.group(1))	# add name
                self.bittypes.append(1)		# add as cbit
                # print "cbit: %s" % m.group(1)
                continue

            # gate acting on qubits
            m = re.compile('\s*(\S+)\s+(\S+)').search(line)
            if(m):
                op = m.group(1)
                args = m.group(2)
                self.gates.append(Qgate(op,args,linenum))
                continue

class qcircuit:		# quantum circuit class

    def __init__(self, bits, bittypes):
        self.bits = bits
        self.bittypes = bittypes
        self.gates = []
        self.classical_subcircuits = []
        self.optimized_gates = []

    def add_gate(self, qgate): self.gates.append(qgate)

    def generate_classical_subcircuits(self):
        gate_id = -1
        init_id = 0
        used_bits = dict([(bit, False) for bit in self.bits])
        forbidden_bits = dict([(bit, False) for bit in self.bits])

        for qgate in self.gates:
            gate_id += 1
            if qgate.gatetype == 'CNOT' or qgate.gatetype == 'Tof':
                for op in qgate.ops:
                    if forbidden_bits[op] == True:
                        self.classical_subcircuits.append((init_id, used_bits))
                        init_id = gate_id
                        used_bits = dict([(bit, False) for bit in self.bits])
                        forbidden_bits = dict([(bit, False) for bit in self.bits])
                    else:    
                        used_bits[op] = True
            else:
                for op in qgate.ops:
                    forbidden_bits[op] = True

        self.classical_subcircuits.append((init_id, used_bits))




#--------------------------------------
#main
qp = qasm_parser(fileinput.input())	# parse the qasm file
qc = qcircuit(qp.bits,qp.bittypes)	# initialize the circuit
for g in qp.gates:			# add each gate to the circuit
    qc.add_gate(g)

qc.generate_classical_subcircuits()
for cir in qc.classical_subcircuits:
    print "init:"+str(cir[0])
    print "used bits:"
    for bit in cir[1]:
        if cir[1][bit] == True: print bit


####################
#Helper functions
####################
def err_print(msg):
#error handling 
    sys.stderr.write('ERROR:' + msg +'\n')
    sys.exit(-1)

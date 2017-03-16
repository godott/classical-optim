#!/usr/bin/python
import re
import fileinput
import sys
import os
import numpy as np

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
        self.qgates = []
        self.classical_subcircuits = []
        self.optimized_subcircuits = []

    def add_gate(self, qgate): self.gates.append(qgate)

    def generate_classical_subcircuits(self):
        init_id = 0
        used_bits = dict([(bit, False) for bit in self.bits])
        forbidden_bits = dict([(bit, False) for bit in self.bits])
        cls_ln = 0
        for gate_id, qgate in enumerate(self.gates):

            new_circuit = False

            if qgate.gatetype == 'CNOT' or qgate.gatetype == 'Tof':
                for op in qgate.ops:
                    used_bits[op] = True
                    if forbidden_bits[op]: new_circuit = True

                if new_circuit:
                    self.classical_subcircuits.append((init_id, used_bits, gate_id-init_id, cls_ln))
                    init_id = gate_id
                    used_bits = dict([(bit, False) for bit in self.bits])
                    for op in qgate.ops:
                        used_bits[op] = True
                    forbidden_bits = dict([(bit, False) for bit in self.bits])
                    cls_ln = 1

                else:    
                    cls_ln += 1
            else:
                for op in qgate.ops:
                    forbidden_bits[op] = True

        self.classical_subcircuits.append((init_id, used_bits, gate_id - init_id, cls_ln))

    def optimize_circuits(self):
        for subcircuit in self.classical_subcircuits:
            self.optimized_subcircuits.append(self.__optimize_circuit(subcircuit))
        print self.optimized_subcircuits[2]
    
    def __optimize_circuit(self, para):
    
        return self.__convert_to_matrix(para)
        #return  self.__optimized_matrix(self.__convert_to_matrix(para))

    #########################
    #private helper methods 
    #########################
    def __convert_to_matrix(self, para):

        start, all_bits, ln, cls_ln = para
        bits = np.array([qb for qb in all_bits if all_bits[qb] == True])
        col_num = len(bits)
        row_num = cls_ln
        end = start + ln
        qgates = []
        cls_ind = -1

        bits_dict = dict([(bits[i], i) for i in range(col_num)])

        circuit = np.zeros((row_num, col_num), dtype = np.int8)

        for ind in range(start, end):

            if self.gates[ind].gatetype != ("CNOT" or "Tof"):
                qgates.append(self.gates[ind])
                continue

            elif self.gates[ind].gatetype == "CNOT":
                cls_ind += 1
                circuit[cls_ind][bits_dict[self.gates[ind].ops[0]]]= 1
                circuit[cls_ind][bits_dict[self.gates[ind].ops[1]]]= -1

            else:
                cls_ind += 1
                circuit[cls_ind][bits_dict[self.gates[ind].ops[0]]]=1
                circuit[cls_ind][bits_dict[self.gates[ind].ops[1]]]=1
                circuit[cls_ind][bits_dict[self.gates[ind].ops[2]]]=-1
        self.qgates.append(qgates) 
        return circuit

 
    #def __optimize_matrix(circuit):

        



        
#--------------------------------------
#main
#--------------------------------------
qp = qasm_parser(fileinput.input())	# parse the qasm file
qc = qcircuit(qp.bits,qp.bittypes)	# initialize the circuit
for g in qp.gates:			# add each gate to the circuit
    qc.add_gate(g)

qc.generate_classical_subcircuits()
qc.optimize_circuits()

####################
#Helper functions
####################
def err_print(msg):
#error handling 
    sys.stderr.write('ERROR:' + msg +'\n')
    sys.exit(-1)

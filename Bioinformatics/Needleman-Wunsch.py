#Needleman-Wunsch alogorithm
#BLOSUM62 substitution matrix
#gap: d = 11
BLOSUM62 = [[9], 
            [-1, 4], 
            [-1, 1, 5], 
            [-3, -1, -1, 7], 
            [0, 1, 0, -1, 4], 
            [-3, 0, -2, -2, 0, 6], 
            [-3, 1, 0, -2, -2, 0, 6], 
            [-3, 0, -1, -1, -2, -1, 1, 6],
            [-4, 0, -1, -1, -1, -2, 0, 2, 5],
            [-3, 0, -1, -1, -1, -2, 0, 0, 2, 5],
            [-3, -1, -2, -2, -2, -2, 1, -1, 0, 0, 8], 
            [-3, -1, -1, -2, -1, -2, 0, -2, 0, 1, 0, 5],
            [-3, 0, -1, -1, -1, -2, 0, -1, 1, 1, -1, 2, 5],
            [-1, -1, -1, -2, -1, -3, -2, -3, -2, 0, -2, -1, -1, 5],
            [-1, -2, -1, -3, -1, -4, -3, -3, -3, -3, -3, -3, -3, 1, 4],
            [-1, -2, -1, -3, -1, -4, -3, -4, -3, -2, -3, -2, -2, 2, 2, 4],
            [-1, -2, 0, -2, 0, -3, -3, -3, -2, -2, -3, -3, -2, 1, 3, 1, 4],
            [-2, -2, -2, -4, -2, -3, -3, -3, -3, -3, -1, -3, -3, 0, 0, 0, -1, 6],
            [-2, -2, -2, -3, -2, -3, -2, -3, -2, -1, 2, -2, -2, -1, -1, -1, -1, 3, 7],
            [-2, -3, -2, -4, -3, -2, -4, -4, -3, -2, -2, -3, -3, -1, -3, -2, -3, 1, 2, 11]]

amino = {"C":0,
         "S":1,
         "T":2,
         "P":3,
         "A":4,
         "G":5,
         "N":6,
         "D":7,
         "E":8,
         "Q":9,
         "H":10,
         "R":11,
         "K":12,
         "M":13,
         "I":14,
         "L":15,
         "V":16,
         "F":17,
         "Y":18,
         "W":19}

class Amino_elem:
    """docstring for amino_elem"""
    def __init__(self, cur_i, cur_j, i, j, value):
        self.cur_i = cur_i
        self.cur_j = cur_j
        self.pre_i = i
        self.pre_j = j
        self.value = value

def swap(a, b):
    a = a + b
    b = a - b
    a = a - b
    return [a, b]

def NeedlemanWunsch(seq_a, seq_b, dot_matrix, i, j):
    sigma = BLOSUM62[amino[seq_a[i-1]]][amino[seq_b[j-1]]] if amino[seq_a[i-1]] >= amino[seq_b[j-1]] else BLOSUM62[amino[seq_b[j-1]]][amino[seq_a[i-1]]]
    temp_dict = [0, 0, 0]
    temp_dict[2] = dot_matrix[i-1][j-1].value + sigma
    temp_dict[0] = i - 1
    temp_dict[1] = j - 1
    if temp_dict[2] < dot_matrix[i-1][j].value - 11:
        temp_dict[2] = dot_matrix[i-1][j].value - 11
        temp_dict[0] = i - 1
        temp_dict[1] = j

    if temp_dict[2] < dot_matrix[i][j-1].value - 11:
        temp_dict[2] = dot_matrix[i][j-1].value - 11
        temp_dict[0] = i
        temp_dict[1] = j -1 

    return Amino_elem(i, j, temp_dict[0], temp_dict[1], temp_dict[2])

def dot_matrix_init(dot_matrix):
    for i in range(len(dot_matrix)):
        if i == 0:
            elem = Amino_elem(0, 0, None, None, 0)
            dot_matrix[0][0] = elem
        else:
            elem = Amino_elem(i, 0, i-1 , 0, -11 * i)
            dot_matrix[i][0] = elem
    for i in range(len(dot_matrix[0])):
        if i == 0:
            elem = Amino_elem(0, 0, None, None, 0)
            dot_matrix[0][0] = elem
        else:
            elem = Amino_elem(0, i, 0 , i-1 , -11 * i)
            dot_matrix[0][i] = elem
    return dot_matrix

def amino_seq(seq_a, seq_b):
    dot_matrix = [[ 0 for i in range(len(seq_b) + 1)] for i in range(len(seq_a) + 1)]
    dot_matrix = dot_matrix_init(dot_matrix)
    for i in range(1, len(dot_matrix)):
        for j in range(1, len(dot_matrix[0])):
            dot_matrix[i][j] = NeedlemanWunsch(seq_a, seq_b, dot_matrix, i, j)
    return dot_matrix

def nw_display(dot_matrix):
    for i in dot_matrix:
        for j in i:
            print("{:<7d}".format(j.value), end = "")
        print("\n")

def nw_solve(dot_matrix, seq_a, seq_b):
    a_len = len(seq_a)
    b_len = len(seq_b)
    a_li = [i for i in seq_a]
    b_li = [i for i in seq_b]
    elem = dot_matrix[a_len][b_len]
    while elem.pre_i != None and elem.pre_j != None:
        if elem.cur_i-1 == elem.pre_i and elem.cur_j == elem.pre_j:
            b_li.insert(elem.cur_j, "-")
        elif elem.cur_i == elem.pre_i and elem.cur_j-1 == elem.pre_j:
            a_li.insert(elem.cur_i, "-")
        print("i:{0}, j:{1}, value:{2}".format(elem.cur_i, elem.cur_j, elem.value))
        elem = dot_matrix[elem.pre_i][elem.pre_j]
    print(a_li)
    print(b_li)
    return[a_li, b_li]
if __name__ == "__main__":
    seq_a = ["V", "E", "S", "L", "C", "Y"]
    seq_b = ["V", "D", "S", "C", "Y"]
    dot_matrix = amino_seq(seq_a, seq_b)
    nw_display(dot_matrix)
    nw_solve(dot_matrix, seq_a, seq_b)

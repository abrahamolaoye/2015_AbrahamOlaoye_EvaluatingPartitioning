__author__ = 'Lovecraft'

script_dir = "/SAN/biosciences/oma/research/Abraham_Olaoye/Pipeline/"
alignments_dir = script_dir + "Data/"

min_cluster_size, max_cluster_size = 4, 8
sample_no = 1

def get_no(pair):
    return pair[1]

x = [("c", 3), ("a", 1), ("b", 2), ("d", 2)]

numbers = range(61, -1, -1)

print numbers

x.sort(key=get_no)
print(x)

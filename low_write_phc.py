import argparse
import math

parser = argparse.ArgumentParser()
parser.add_argument('-s', type=str, dest="s", default = "")
parser.add_argument('-reps', type=int, dest="reps", default = 10)
parser.add_argument('-m', type=int, dest="m", default=5)
parser.add_argument('-a1', type=int, dest="a1", default=3)
parser.add_argument('-a2', type=int, dest="a2", default=12)

args = vars(parser.parse_args())
reps = args["reps"]
a1 = args["a1"]
a2 = args["a2"]
script_filename = args["s"]
m = args["m"]

f = open('phc_eqs/'+str(script_filename)+'.sh','w')
for a in range(a1,a2+1):
	for t in range(reps):
		for j in range(m):
			md = math.ceil(math.log(a))
			base_name = 'ln'+str(a)+'md'+str(md)
			f.write('../.././phc -b '+base_name+'_'+str(t)+'_eqs.txt '+base_name+'_'+str(t)+'_roots'+str(j)+'.txt\n')
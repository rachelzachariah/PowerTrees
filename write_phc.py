import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-s', type=str, dest="s", default = "")
parser.add_argument('-w', type=str, dest="w", default = "")
parser.add_argument('-reps', type=int, dest="reps", default = 10)


args = vars(parser.parse_args())
reps = args["reps"]
phc_filename = args["w"]
script_filename = args["s"]

f = open('phc_eqs/'+str(script_filename)+'.sh','w')
for t in range(reps):
	f.write('../.././phc -b '+str(phc_filename)+'_'+str(t)+'_eqs.txt '+str(phc_filename)+'_'+str(t)+'_roots.txt')







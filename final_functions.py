import math

def low_grep():
	f = open('phc_eqs/low_grep.sh','w')
	for a in range(3,13):
		md = int(math.ceil(math.log(a)))
		for t in range(10):
			base_name = 'ln'+str(a)+'md'+str(md)+'_'+str(t)
			for j in range(5):
				f.write('grep "User time in seconds" '+base_name+'_'+str(t)+'_roots'+str(j)+'.txt > '+str(base_name)+'_time'+str(j)+'.txt\n')
	f.close()

def med_grep():
	f = open('phc_eqs/med_grep.sh','w')
	for a in range(3,13):
		md = int(math.ceil(math.sqrt(a)))
		for t in range(10):
			base_name = 'ln'+str(a)+'md'+str(md)+'_'+str(t)
			for j in range(5):
				f.write('grep "User time in seconds" '+base_name+'_'+str(t)+'_roots'+str(j)+'.txt > '+str(base_name)+'_time'+str(j)+'.txt\n')
	f.close()

def high_grep():
	f = open('phc_eqs/high_grep.sh','w')
	for a in range(3,13):
		md = a-1
		for t in range(10):
			base_name = 'ln'+str(a)+'md'+str(md)+'_'+str(t)
			for j in range(5):
				f.write('grep "User time in seconds" '+base_name+'_'+str(t)+'_roots'+str(j)+'.txt > '+str(base_name)+'_time'+str(j)+'.txt\n')
	f.close()	



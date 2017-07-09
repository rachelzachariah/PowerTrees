import math

def low_grep():
	f = open('phc_eqs/low_grep.sh','w')
	for a in range(3,13):
		md = int(math.ceil(math.log(a)))
		for t in range(10):
			base_name = 'ln'+str(a)+'md'+str(md)+'_'+str(t)
			for j in range(5):
				f.write('grep "User time in seconds" '+base_name+'_roots'+str(j)+'.txt > '+str(base_name)+'_time'+str(j)+'.txt\n')
	f.close()


def low_scan():
	timing_info = {}
	for a in range(3,13):
		md = int(math.ceil(math.log(a)))

		a_times = []
		t_failures = []

		for t in range(10):
			base_name = 'ln'+str(a)+'md'+str(md)+'_'+str(t)

			t_times = []

			for j in range(5):
				f = open(str(base_name)+'_time'+str(j)+'.txt')
				ar = f.readlines()
				if len(ar) == 4:
					try:
						last_line = ar[-1]
						w = last_line.strip()
						w2 = w.split(' ')
						w3 = filter(lambda aa: aa!='',w2)
						equals_ind = w3.index('=')
						this_time = float(w3[equals_ind-1])
						t_times.append(this_time)
					except Exception:
						this_time = 0
				f.close()

			if len(t_times) > 0:
				average_t_time = float(sum(t_times))/(float(len(t_times)))
			else:
				t_failures.append(t)
			a_times.append(average_t_time)
		if len(a_times) == 0:
			print str(a)+' nodes : No successes'
		else:
			average_a = float(sum(a_times))/(float(len(a_times)))
			print str(a)+' nodes : '+str(average_a)+' seconds'
			print 'Failed instances : '+str(t_failures)
				


		

def med_grep():
	f = open('phc_eqs/med_grep.sh','w')
	for a in range(3,13):
		md = int(math.ceil(math.sqrt(a)))
		for t in range(10):
			base_name = 'ln'+str(a)+'md'+str(md)+'_'+str(t)
			for j in range(5):
				f.write('grep "User time in seconds" '+base_name+'_roots'+str(j)+'.txt > '+str(base_name)+'_time'+str(j)+'.txt\n')
	f.close()


def med_scan():
	timing_info = {}
	for a in range(3,13):
		md = int(math.ceil(math.sqrt(a)))

		a_times = []
		t_failures = []

		for t in range(10):
			base_name = 'mn'+str(a)+'md'+str(md)+'_'+str(t)

			t_times = []

			for j in range(5):
				f = open(str(base_name)+'_time'+str(j)+'.txt')
				ar = f.readlines()
				if len(ar) == 4:
					try:
						last_line = ar[-1]
						w = last_line.strip()
						w2 = w.split(' ')
						w3 = filter(lambda aa: aa!='',w2)
						equals_ind = w3.index('=')
						this_time = float(w3[equals_ind-1])
						t_times.append(this_time)
					except Exception:
						this_time = 0
				f.close()

			if len(t_times) > 0:
				average_t_time = float(sum(t_times))/(float(len(t_times)))
			else:
				t_failures.append(t)
			a_times.append(average_t_time)
		if len(a_times) == 0:
			print str(a)+' nodes : No successes'
		else:
			average_a = float(sum(a_times))/(float(len(a_times)))
			print str(a)+' nodes : '+str(average_a)+' seconds'
			print 'Failed instances : '+str(t_failures)	

def high_grep():
	f = open('phc_eqs/high_grep.sh','w')
	for a in range(3,13):
		md = a-1
		for t in range(10):
			base_name = 'hn'+str(a)+'md'+str(md)+'_'+str(t)
			for j in range(5):
				f.write('grep "User time in seconds" '+base_name+'_roots'+str(j)+'.txt > '+str(base_name)+'_time'+str(j)+'.txt\n')
	f.close()	



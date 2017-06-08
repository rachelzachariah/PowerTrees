def write_equations(eqs,filename):
	disp = str(len(eqs))+'\n'
	for coeff in eqs:
		str_coeff = str(coeff)
		disp = disp+str_coeff.replace(' ','') + ';\n'
	f = open(filename+'_eqs.txt','w')
	f.write(disp)
	f.close()

def display_equations(eqs):
	disp = str(len(eqs))+'\n'
	for coeff in eqs:
		str_coeff = str(coeff)
		disp = disp+str_coeff.replace(' ','') + ';\n'
	print disp

def write_eqs(n,filename):
	f = open(prefix+'phc_script.sh','w')
	for conf in confs:
		x_term = conf[0]
		y_term = conf[1]
		z_term = conf[2]
		eqs = self.find_equations(x_term,y_term,z_term)
		filename = self.make_filename(x_term,y_term,z_term)
		self.write_equations(eqs,prefix+filename)
		f.write("../../.././phc -b "+filename+"_eqs.txt " +filename+"_roots.txt\n")
	f.close()
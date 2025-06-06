import grins.reg_funcs as regfn

def odesys(t,y,args):
	(A, B) = y
	(Prod_A, Prod_B, Deg_A, Deg_B, InhFld_B_A, Thr_B_A, Hill_B_A, InhFld_A_B, Thr_A_B, Hill_A_B) = args
	d_A = Prod_A*regfn.nsH(B, InhFld_B_A, Thr_B_A, Hill_B_A) - Deg_A*A
	d_B = Prod_B*regfn.nsH(A, InhFld_A_B, Thr_A_B, Hill_A_B) - Deg_B*B
	d_y = (d_A, d_B)
	return d_y

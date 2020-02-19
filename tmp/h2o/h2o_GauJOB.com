%chk=h2o_MO.chk
%nproc=4
%mem=100MW
#p HF/6-31G guess=read scf=(xqc,maxcon=128,maxcyc=512,conver=8) force IOp(1/10=10) IOp(1/18=10) IOp(2/15=3) IOp(2/12=1)

force calculation: 1 0 0 0 0

0 1
O 	 0	  -0.743799822636	  -0.020302731545	  -0.397490659011
H 	 0	  -1.241138584781	   0.636471191669	   0.074844309857
H 	 0	  -1.313015367583	  -0.616168460123	  -0.869439650850



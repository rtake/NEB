%chk=h2o_MO.chk
%nproc=4
%mem=100MW
#p HF/6-31G guess=read scf=(xqc,maxcon=128,maxcyc=512,conver=8) freq=noraman IOp(7/8=20000) IOp(1/10=10) IOp(1/18=10) IOp(2/15=3) IOp(2/12=1)

frequency calculation: 1 0 0 0 0

0 1
O 	 0	  -1.158815261413	   0.003416203666	  -0.397360325119
H 	 0	  -0.589605035440	   0.599132739306	   0.074778645992
H 	 0	  -0.661480029059	  -0.653505682108	  -0.869504313118



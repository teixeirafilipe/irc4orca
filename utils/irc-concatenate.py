#! /usr/bin/env python3
# -*- coding: utf8 -*-

##########################################################################
#                                                                        #
# Programa: irc-concatenate.py              Data: 08/12/2013             #
#                                                                        #
# Uso: irc-concatenate.py forward.trj backward.trj out.trj               #
#                                                                        #
# Inverts backward.trj and concatenates it to forward.trj, resulting in  #
# a single trj file for the whole reaction.                              #
#                                                                        #
# (c) Filipe Teixeira 2013                                               #
#                                                                        #
# Revised 01/02/2018: Updated to Python 3.                               #
#                                                                        #
#                                                                        #
##########################################################################

import sys

def xyzInvert(fdata):
	nlines=int(fdata[0])+2
	tmp=[]
	for i in range(0,len(fdata),nlines):
		tmp.append(fdata[i:i+nlines])
	tmp.reverse()
	o=[]
	for cc in tmp:
		for l in cc:
			o.append(l)
	return o

if (__name__=='__main__'):
	try:
		fn1=sys.argv[1]
		fn2=sys.argv[2]
		fn3=sys.argv[3]
	except:
		print('Usage: irc-concatenate.py forward.trj backward.trj out.trj')
		sys.exit(1)
	ffile=open(fn1,'r')
	bfile=open(fn2,'r')
	ffdata=ffile.readlines()
	bbdata=bfile.readlines()
	ffile.close()
	bfile.close()
	nbdata=xyzInvert(bbdata)
	of=open(fn3,'w')
	for line in nbdata:
		of.write(line)
	for line in ffdata:
		of.write(line)
	of.close()


	





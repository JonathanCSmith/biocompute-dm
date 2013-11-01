#!/usr/bin/python


import os

files=os.listdir(".")

pages=[]
TF=open("template_top","r")
template_top=TF.readlines()
TF.close()

BF=open("template_bottom","r")
template_bottom=BF.readlines()
BF.close()



for m in range(0,len(files)):
	if files[m][-5:]=="_html":
		pages.append(files[m])     

for wr in pages:
	CT=open(wr,"r")
	content=CT.readlines()
	CT.close()

	pageName="../html/CoreInSys/"+wr[:-5]+".html"
	print pageName
	page=[]
	for top in template_top:
		page.append(top)

	for co in content:
		page.append(co)

	for bottom in template_bottom:
		page.append(bottom)

	PA=open(pageName,"w")
	for p in page:
		PA.write(p)
	PA.close()


	


	





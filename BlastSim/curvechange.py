# make new curves

import re



f = open("change.txt")
content = f.read()
f.close()

f = open("write.txt", "w")

content = "".join(content.split())
number = "((?:-|\\+)?\\d*\\.\\d+(?:e-?\\d+)?)"
regex = "a1="+number+""+number+"\\*Y;" + \
"a2=\\("+number+""+number+"\\*Y\\)\\*Z;" + \
"a3=\\("+number+""+number+"\\*Z"+number+"\\*Y\\)\\*Y\\*Y;" + \
"a4=\\("+number+""+number+"\\*Y"+number+"\\*Z\\)\\*Z\\*Z;(?:" + \
"a5="+number+""+number+"\\*Y;" + \
"a6=\\("+number+""+number+"\\*Y\\)\\*Z;"+ \
"a7=\\("+number+""+number+"\\*Z+"+number+"\\*Y\\)\\*Y\\*Y;" + \
"a8=\\("+number+""+number+"\\*Y"+number+"\\*Z\\)\\*Z\\*Z;" + \
"a9=exp\\("+number+"(?:"+number+"\\*Y)?(?:"+number+"\\*Z)?(?:"+number+"\\*Y\\*Z)?\\);)?"

cregex = re.compile(regex)




start = 0
for a in re.findall(regex, content):
	a = [0, *a]

	for i in range(1, len(a)):
		if a[i] and a[i][0] not in ['-', '+']:
			a[i] = '+' + a[i]


	"""
	= a2 + a4 Z + a8 Z2 + 2 a5 Y + 2 a7 YZ + 3 a9 Y2
	+ (a12 +a14 Z + 2 al7 YZ + al8 Z2 + 2 a15 Y + 3 al9 Y2) /
	[1 ± exp(a21 + a22 Y + a23 Z + a24 YZ)]
	±± (a11 + a12 Y + a13 Z + a14 YZ
	+ a15 Y2 + al7 Y2Z + a19 Y3 + a16 Z2 + a18 YZ2 + a20 Z3)(a22 + a24 Z)
	[exp(a21 + a22 Y + a23 Z + a24 YZ)/
	[1 ± exp(a21 + a22 Y + a23 Z + a24 YZ)]2 
	"""

	f.write("""
		a1 = %s + (%s %s * Z) * Z;
		a2 = (%s * 2 %s * 2 * Z %s * 3 * Y) * Y;
		a3 = %s + (%s %s * Y %s * Z) * Z + (%s * 2 %s * 2 * Y) * Y;
		a4 = exp(%s %s * Y %s * Z %s * Y * Z);
		a5 = %s %s * Y + (%s %s * Y) * Z;
		a6 = (%s %s * Z %s * Y) * Y * Y;
		a7 = (%s %s * Y %s * Z) * Z * Z;
		a8 = %s %s * Z;
		a9 = exp(%s %s * Y %s * Z %s * Y * Z);
		U = a1 + a2 + a3 / (1 + a4) - (a5 + a6 + a7) * a8 * a9 / ((1 + a9) * (1 + a9));\n""" %
		(a[2], a[4], a[8], a[5], a[7], a[9],
		a[12], a[14], a[17], a[18], a[15], a[19],
		a[21], a[22], a[23], a[24],
		a[11], a[12], a[13], a[14],
		a[15], a[17], a[19], a[16], a[18], a[20],
		a[22], a[24],
		a[21], a[22], a[23], a[24]))



	"""
	= a3 + 2 a6 Z + a10 Z2 + a4 Y + 2 a8 YZ + a7 Y2
	+ (a13 + 2 a16 Z + 2 a18 YZ + 3 a20 Z2 + a14 Y + 2 a17 Y2) /
	[1 ± exp(a21 + a22 Y + a23 Z + a24 YZ)]
	±± (a11 + a12 Y + a13 Z + a14 YZ
	+ a15 Y2 + al7 Y2Z + a19 Y3 + a16 Z2 + a18 YZ2 + a20 Z3)(a23 + a24 Y)
	[exp(a21 + a22 Y + a23 Z + a24 YZ)/
	[1 ± exp(a21 + a22 Y + a23 Z + a24 YZ)]2 
	"""

	f.write("""
		a1 = %s + (%s * 2 %s * Z) * Z;
		a2 = (%s %s * Z %s * Y) * Y;
		a3 = %s + (%s * 2 %s * 2 * Y %s * 3 * Z) * Z + (%s * 2 %s * Y) * Y;
		a4 = exp(%s %s * Y %s * Z %s * Y * Z);
		a5 = %s %s * Y + (%s %s * Y) * Z;
		a6 = (%s %s * Z %s * Y) * Y * Y;
		a7 = (%s %s * Y %s * Z) * Z * Z;
		a8 = %s %s * Z;
		a9 = exp(%s %s * Y %s * Z %s * Y * Z);
		V = a1 + a2 + a3 / (1 + a4) - (a5 + a6 + a7) * a8 * a9 / ((1 + a9) * (1 + a9));\n""" %
		(a[3], a[6], a[10], a[4], a[8], a[7],
		a[13], a[16], a[18], a[20], a[14], a[17],
		a[21], a[22], a[23], a[24],
		a[11], a[12], a[13], a[14],
		a[15], a[17], a[19], a[16], a[18], a[20],
		a[22], a[24],
		a[21], a[22], a[23], a[24]))



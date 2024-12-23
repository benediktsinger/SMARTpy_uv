import sys, os
import pandas as pd

lines = []
with open('chimera_log.txt', 'r') as file:
	lines = file.readlines()

o = {}
out = []
name = ''
for i in range(len(lines)):
	line = lines[i]
	if line.startswith('RhRh'):
		o = {}
		name = line
		print(line)
		o['Name'] = name
	elif 'molmap res 2.8:  volume' in line:
		o['Volume'] = line.split()[-1]
	elif 'molmap res 2.8:  area' in line:
		o['Area'] = line.split()[-1]
	elif 'Contact area on molmap res 2.8 within distance 1.1' in line:
		nextline = lines[i+1]
		o['CSA'] = nextline.split()[-1]
	elif 'molmap res 2.8 masked:  volume' in line:
		o['Prox5_Volume'] = line.split()[-1]
	elif 'molmap res 2.8 masked:  area' in line:
		o['Prox5_Area'] = line.split()[-1]
	elif 'Contact area on molmap res 2.8 masked within distance 1.1' in line:
		nextline = lines[i+1]
		o['Prox5_CSA'] = nextline.split()[-1]
		out.append(o)

df = pd.DataFrame(out)
print(df)
df.to_excel('chimera_parsed.xlsx')

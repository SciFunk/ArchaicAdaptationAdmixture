import gzip
import pandas as pd
import numpy as np

file_name = "chr1.csv"

vcf = gzip.open('chr1.vcf.gz', "rt")
line = vcf.readline()
spline = line.split()
for line in range(253):
    line = vcf.readline()
spline = line.split()
header = spline
one_zero=[0]*(len(header)-9)

zero_one=[0]*(len(header)-9)

one_one=[0]*(len(header)-9)

zeros=[0]*(len(header)-9)

multiallelic=[0]*(len(header)-9)

for line in vcf:
        spline = line.split()
        
        genos = spline[9:]
        for entry in range(len(genos)):
                indgenos = genos[entry]
                if indgenos == '1|0':  
                    one_zero[entry] +=1  
                elif indgenos == '0|1':
                    zero_one[entry] +=1
                elif indgenos == '1|1':
                    one_one[entry] +=1
                elif indgenos == '0|0':
                    zeros[entry] +=1
                else: 
                    multiallelic[entry] +=1

# Index from call_samples
df_samples = pd.read_csv('call_samples.ALL.panel.csv', delimiter='\t')
samples = df_samples.drop(columns=["Unnamed: 4", "Unnamed: 5"], axis = 1)
samples.head()
pop = samples['pop']

# Combined data and save to .csv
data = pd.DataFrame(
        {'0|1' : zero_one,
         '1|0' : one_zero,
         '1|1' : one_one,
        })
prelim = data.sum(axis = 1)
data_sum = np.array(prelim)

final = pd.DataFrame(
        {'0|1' : zero_one,
         '1|0' : one_zero,
         '1|1' : one_one,
         'sum_biallelic' : data_sum,
         'multiallelic' : multiallelic,  
        }, index=pop)
    
final.to_csv(file_name)    
#______________________________________________________________________________

file_name = "chr2.csv"

vcf = gzip.open('chr2.vcf.gz', "rt")
line = vcf.readline()
spline = line.split()
for line in range(253):
    line = vcf.readline()
spline = line.split()
header = spline
one_zero=[0]*(len(header)-9)

zero_one=[0]*(len(header)-9)

one_one=[0]*(len(header)-9)

zeros=[0]*(len(header)-9)

multiallelic=[0]*(len(header)-9)

for line in vcf:
        spline = line.split()
        
        genos = spline[9:]
        for entry in range(len(genos)):
                indgenos = genos[entry]
                if indgenos == '1|0':  
                    one_zero[entry] +=1  
                elif indgenos == '0|1':
                    zero_one[entry] +=1
                elif indgenos == '1|1':
                    one_one[entry] +=1
                elif indgenos == '0|0':
                    zeros[entry] +=1
                else: 
                    multiallelic[entry] +=1

# Index from call_samples
df_samples = pd.read_csv('call_samples.ALL.panel.csv', delimiter='\t')
samples = df_samples.drop(columns=["Unnamed: 4", "Unnamed: 5"], axis = 1)
samples.head()
pop = samples['pop']

# Combined data and save to .csv
data = pd.DataFrame(
        {'0|1' : zero_one,
         '1|0' : one_zero,
         '1|1' : one_one,
        })
prelim = data.sum(axis = 1)
data_sum = np.array(prelim)

final = pd.DataFrame(
        {'0|1' : zero_one,
         '1|0' : one_zero,
         '1|1' : one_one,
         'sum_biallelic' : data_sum,
         'multiallelic' : multiallelic,  
        }, index=pop)
    
final.to_csv(file_name) 

#______________________________________________________________________________

file_name = "chr3.csv"

vcf = gzip.open('chr3.vcf.gz', "rt")
line = vcf.readline()
spline = line.split()
for line in range(253):
    line = vcf.readline()
spline = line.split()
header = spline
one_zero=[0]*(len(header)-9)

zero_one=[0]*(len(header)-9)

one_one=[0]*(len(header)-9)

zeros=[0]*(len(header)-9)

multiallelic=[0]*(len(header)-9)

for line in vcf:
        spline = line.split()
        
        genos = spline[9:]
        for entry in range(len(genos)):
                indgenos = genos[entry]
                if indgenos == '1|0':  
                    one_zero[entry] +=1  
                elif indgenos == '0|1':
                    zero_one[entry] +=1
                elif indgenos == '1|1':
                    one_one[entry] +=1
                elif indgenos == '0|0':
                    zeros[entry] +=1
                else: 
                    multiallelic[entry] +=1

# Index from call_samples
df_samples = pd.read_csv('call_samples.ALL.panel.csv', delimiter='\t')
samples = df_samples.drop(columns=["Unnamed: 4", "Unnamed: 5"], axis = 1)
samples.head()
pop = samples['pop']

# Combined data and save to .csv
data = pd.DataFrame(
        {'0|1' : zero_one,
         '1|0' : one_zero,
         '1|1' : one_one,
        })
prelim = data.sum(axis = 1)
data_sum = np.array(prelim)

final = pd.DataFrame(
        {'0|1' : zero_one,
         '1|0' : one_zero,
         '1|1' : one_one,
         'sum_biallelic' : data_sum,
         'multiallelic' : multiallelic,  
        }, index=pop)
    
final.to_csv(file_name) 

#______________________________________________________________________________

file_name = "chr4.csv"

vcf = gzip.open('chr4.vcf.gz', "rt")
line = vcf.readline()
spline = line.split()
for line in range(253):
    line = vcf.readline()
spline = line.split()
header = spline
one_zero=[0]*(len(header)-9)

zero_one=[0]*(len(header)-9)

one_one=[0]*(len(header)-9)

zeros=[0]*(len(header)-9)

multiallelic=[0]*(len(header)-9)

for line in vcf:
        spline = line.split()
        
        genos = spline[9:]
        for entry in range(len(genos)):
                indgenos = genos[entry]
                if indgenos == '1|0':  
                    one_zero[entry] +=1  
                elif indgenos == '0|1':
                    zero_one[entry] +=1
                elif indgenos == '1|1':
                    one_one[entry] +=1
                elif indgenos == '0|0':
                    zeros[entry] +=1
                else: 
                    multiallelic[entry] +=1

# Index from call_samples
df_samples = pd.read_csv('call_samples.ALL.panel.csv', delimiter='\t')
samples = df_samples.drop(columns=["Unnamed: 4", "Unnamed: 5"], axis = 1)
samples.head()
pop = samples['pop']

# Combined data and save to .csv
data = pd.DataFrame(
        {'0|1' : zero_one,
         '1|0' : one_zero,
         '1|1' : one_one,
        })
prelim = data.sum(axis = 1)
data_sum = np.array(prelim)

final = pd.DataFrame(
        {'0|1' : zero_one,
         '1|0' : one_zero,
         '1|1' : one_one,
         'sum_biallelic' : data_sum,
         'multiallelic' : multiallelic,  
        }, index=pop)
    
final.to_csv(file_name) 

#______________________________________________________________________________

file_name = "chr5.csv"

vcf = gzip.open('chr5.vcf.gz', "rt")
line = vcf.readline()
spline = line.split()
for line in range(253):
    line = vcf.readline()
spline = line.split()
header = spline
one_zero=[0]*(len(header)-9)

zero_one=[0]*(len(header)-9)

one_one=[0]*(len(header)-9)

zeros=[0]*(len(header)-9)

multiallelic=[0]*(len(header)-9)

for line in vcf:
        spline = line.split()
        
        genos = spline[9:]
        for entry in range(len(genos)):
                indgenos = genos[entry]
                if indgenos == '1|0':  
                    one_zero[entry] +=1  
                elif indgenos == '0|1':
                    zero_one[entry] +=1
                elif indgenos == '1|1':
                    one_one[entry] +=1
                elif indgenos == '0|0':
                    zeros[entry] +=1
                else: 
                    multiallelic[entry] +=1

# Index from call_samples
df_samples = pd.read_csv('call_samples.ALL.panel.csv', delimiter='\t')
samples = df_samples.drop(columns=["Unnamed: 4", "Unnamed: 5"], axis = 1)
samples.head()
pop = samples['pop']

# Combined data and save to .csv
data = pd.DataFrame(
        {'0|1' : zero_one,
         '1|0' : one_zero,
         '1|1' : one_one,
        })
prelim = data.sum(axis = 1)
data_sum = np.array(prelim)

final = pd.DataFrame(
        {'0|1' : zero_one,
         '1|0' : one_zero,
         '1|1' : one_one,
         'sum_biallelic' : data_sum,
         'multiallelic' : multiallelic,  
        }, index=pop)
    
final.to_csv(file_name) 

#______________________________________________________________________________

file_name = "chr6.csv"

vcf = gzip.open('chr6.vcf.gz', "rt")
line = vcf.readline()
spline = line.split()
for line in range(253):
    line = vcf.readline()
spline = line.split()
header = spline
one_zero=[0]*(len(header)-9)

zero_one=[0]*(len(header)-9)

one_one=[0]*(len(header)-9)

zeros=[0]*(len(header)-9)

multiallelic=[0]*(len(header)-9)

for line in vcf:
        spline = line.split()
        
        genos = spline[9:]
        for entry in range(len(genos)):
                indgenos = genos[entry]
                if indgenos == '1|0':  
                    one_zero[entry] +=1  
                elif indgenos == '0|1':
                    zero_one[entry] +=1
                elif indgenos == '1|1':
                    one_one[entry] +=1
                elif indgenos == '0|0':
                    zeros[entry] +=1
                else: 
                    multiallelic[entry] +=1

# Index from call_samples
df_samples = pd.read_csv('call_samples.ALL.panel.csv', delimiter='\t')
samples = df_samples.drop(columns=["Unnamed: 4", "Unnamed: 5"], axis = 1)
samples.head()
pop = samples['pop']

# Combined data and save to .csv
data = pd.DataFrame(
        {'0|1' : zero_one,
         '1|0' : one_zero,
         '1|1' : one_one,
        })
prelim = data.sum(axis = 1)
data_sum = np.array(prelim)

final = pd.DataFrame(
        {'0|1' : zero_one,
         '1|0' : one_zero,
         '1|1' : one_one,
         'sum_biallelic' : data_sum,
         'multiallelic' : multiallelic,  
        }, index=pop)
    
final.to_csv(file_name) 

#______________________________________________________________________________

file_name = "chr7.csv"

vcf = gzip.open('chr7.vcf.gz', "rt")
line = vcf.readline()
spline = line.split()
for line in range(253):
    line = vcf.readline()
spline = line.split()
header = spline
one_zero=[0]*(len(header)-9)

zero_one=[0]*(len(header)-9)

one_one=[0]*(len(header)-9)

zeros=[0]*(len(header)-9)

multiallelic=[0]*(len(header)-9)

for line in vcf:
        spline = line.split()
        
        genos = spline[9:]
        for entry in range(len(genos)):
                indgenos = genos[entry]
                if indgenos == '1|0':  
                    one_zero[entry] +=1  
                elif indgenos == '0|1':
                    zero_one[entry] +=1
                elif indgenos == '1|1':
                    one_one[entry] +=1
                elif indgenos == '0|0':
                    zeros[entry] +=1
                else: 
                    multiallelic[entry] +=1

# Index from call_samples
df_samples = pd.read_csv('call_samples.ALL.panel.csv', delimiter='\t')
samples = df_samples.drop(columns=["Unnamed: 4", "Unnamed: 5"], axis = 1)
samples.head()
pop = samples['pop']

# Combined data and save to .csv
data = pd.DataFrame(
        {'0|1' : zero_one,
         '1|0' : one_zero,
         '1|1' : one_one,
        })
prelim = data.sum(axis = 1)
data_sum = np.array(prelim)

final = pd.DataFrame(
        {'0|1' : zero_one,
         '1|0' : one_zero,
         '1|1' : one_one,
         'sum_biallelic' : data_sum,
         'multiallelic' : multiallelic,  
        }, index=pop)
    
final.to_csv(file_name) 

#______________________________________________________________________________

file_name = "chr8.csv"

vcf = gzip.open('chr8.vcf.gz', "rt")
line = vcf.readline()
spline = line.split()
for line in range(253):
    line = vcf.readline()
spline = line.split()
header = spline
one_zero=[0]*(len(header)-9)

zero_one=[0]*(len(header)-9)

one_one=[0]*(len(header)-9)

zeros=[0]*(len(header)-9)

multiallelic=[0]*(len(header)-9)

for line in vcf:
        spline = line.split()
        
        genos = spline[9:]
        for entry in range(len(genos)):
                indgenos = genos[entry]
                if indgenos == '1|0':  
                    one_zero[entry] +=1  
                elif indgenos == '0|1':
                    zero_one[entry] +=1
                elif indgenos == '1|1':
                    one_one[entry] +=1
                elif indgenos == '0|0':
                    zeros[entry] +=1
                else: 
                    multiallelic[entry] +=1

# Index from call_samples
df_samples = pd.read_csv('call_samples.ALL.panel.csv', delimiter='\t')
samples = df_samples.drop(columns=["Unnamed: 4", "Unnamed: 5"], axis = 1)
samples.head()
pop = samples['pop']

# Combined data and save to .csv
data = pd.DataFrame(
        {'0|1' : zero_one,
         '1|0' : one_zero,
         '1|1' : one_one,
        })
prelim = data.sum(axis = 1)
data_sum = np.array(prelim)

final = pd.DataFrame(
        {'0|1' : zero_one,
         '1|0' : one_zero,
         '1|1' : one_one,
         'sum_biallelic' : data_sum,
         'multiallelic' : multiallelic,  
        }, index=pop)
    
final.to_csv(file_name) 

#______________________________________________________________________________

file_name = "chr9.csv"

vcf = gzip.open('chr9.vcf.gz', "rt")
line = vcf.readline()
spline = line.split()
for line in range(253):
    line = vcf.readline()
spline = line.split()
header = spline
one_zero=[0]*(len(header)-9)

zero_one=[0]*(len(header)-9)

one_one=[0]*(len(header)-9)

zeros=[0]*(len(header)-9)

multiallelic=[0]*(len(header)-9)

for line in vcf:
        spline = line.split()
        
        genos = spline[9:]
        for entry in range(len(genos)):
                indgenos = genos[entry]
                if indgenos == '1|0':  
                    one_zero[entry] +=1  
                elif indgenos == '0|1':
                    zero_one[entry] +=1
                elif indgenos == '1|1':
                    one_one[entry] +=1
                elif indgenos == '0|0':
                    zeros[entry] +=1
                else: 
                    multiallelic[entry] +=1

# Index from call_samples
df_samples = pd.read_csv('call_samples.ALL.panel.csv', delimiter='\t')
samples = df_samples.drop(columns=["Unnamed: 4", "Unnamed: 5"], axis = 1)
samples.head()
pop = samples['pop']

# Combined data and save to .csv
data = pd.DataFrame(
        {'0|1' : zero_one,
         '1|0' : one_zero,
         '1|1' : one_one,
        })
prelim = data.sum(axis = 1)
data_sum = np.array(prelim)

final = pd.DataFrame(
        {'0|1' : zero_one,
         '1|0' : one_zero,
         '1|1' : one_one,
         'sum_biallelic' : data_sum,
         'multiallelic' : multiallelic,  
        }, index=pop)
    
final.to_csv(file_name) 

#______________________________________________________________________________

file_name = "chr10.csv"

vcf = gzip.open('chr10.vcf.gz', "rt")
line = vcf.readline()
spline = line.split()
for line in range(253):
    line = vcf.readline()
spline = line.split()
header = spline
one_zero=[0]*(len(header)-9)

zero_one=[0]*(len(header)-9)

one_one=[0]*(len(header)-9)

zeros=[0]*(len(header)-9)

multiallelic=[0]*(len(header)-9)

for line in vcf:
        spline = line.split()
        
        genos = spline[9:]
        for entry in range(len(genos)):
                indgenos = genos[entry]
                if indgenos == '1|0':  
                    one_zero[entry] +=1  
                elif indgenos == '0|1':
                    zero_one[entry] +=1
                elif indgenos == '1|1':
                    one_one[entry] +=1
                elif indgenos == '0|0':
                    zeros[entry] +=1
                else: 
                    multiallelic[entry] +=1

# Index from call_samples
df_samples = pd.read_csv('call_samples.ALL.panel.csv', delimiter='\t')
samples = df_samples.drop(columns=["Unnamed: 4", "Unnamed: 5"], axis = 1)
samples.head()
pop = samples['pop']

# Combined data and save to .csv
data = pd.DataFrame(
        {'0|1' : zero_one,
         '1|0' : one_zero,
         '1|1' : one_one,
        })
prelim = data.sum(axis = 1)
data_sum = np.array(prelim)

final = pd.DataFrame(
        {'0|1' : zero_one,
         '1|0' : one_zero,
         '1|1' : one_one,
         'sum_biallelic' : data_sum,
         'multiallelic' : multiallelic,  
        }, index=pop)
    
final.to_csv(file_name) 

#______________________________________________________________________________

file_name = "chr11.csv"

vcf = gzip.open('chr11.vcf.gz', "rt")
line = vcf.readline()
spline = line.split()
for line in range(253):
    line = vcf.readline()
spline = line.split()
header = spline
one_zero=[0]*(len(header)-9)

zero_one=[0]*(len(header)-9)

one_one=[0]*(len(header)-9)

zeros=[0]*(len(header)-9)

multiallelic=[0]*(len(header)-9)

for line in vcf:
        spline = line.split()
        
        genos = spline[9:]
        for entry in range(len(genos)):
                indgenos = genos[entry]
                if indgenos == '1|0':  
                    one_zero[entry] +=1  
                elif indgenos == '0|1':
                    zero_one[entry] +=1
                elif indgenos == '1|1':
                    one_one[entry] +=1
                elif indgenos == '0|0':
                    zeros[entry] +=1
                else: 
                    multiallelic[entry] +=1

# Index from call_samples
df_samples = pd.read_csv('call_samples.ALL.panel.csv', delimiter='\t')
samples = df_samples.drop(columns=["Unnamed: 4", "Unnamed: 5"], axis = 1)
samples.head()
pop = samples['pop']

# Combined data and save to .csv
data = pd.DataFrame(
        {'0|1' : zero_one,
         '1|0' : one_zero,
         '1|1' : one_one,
        })
prelim = data.sum(axis = 1)
data_sum = np.array(prelim)

final = pd.DataFrame(
        {'0|1' : zero_one,
         '1|0' : one_zero,
         '1|1' : one_one,
         'sum_biallelic' : data_sum,
         'multiallelic' : multiallelic,  
        }, index=pop)
    
final.to_csv(file_name) 

#______________________________________________________________________________

file_name = "chr12.csv"

vcf = gzip.open('chr12.vcf.gz', "rt")
line = vcf.readline()
spline = line.split()
for line in range(253):
    line = vcf.readline()
spline = line.split()
header = spline
one_zero=[0]*(len(header)-9)

zero_one=[0]*(len(header)-9)

one_one=[0]*(len(header)-9)

zeros=[0]*(len(header)-9)

multiallelic=[0]*(len(header)-9)

for line in vcf:
        spline = line.split()
        
        genos = spline[9:]
        for entry in range(len(genos)):
                indgenos = genos[entry]
                if indgenos == '1|0':  
                    one_zero[entry] +=1  
                elif indgenos == '0|1':
                    zero_one[entry] +=1
                elif indgenos == '1|1':
                    one_one[entry] +=1
                elif indgenos == '0|0':
                    zeros[entry] +=1
                else: 
                    multiallelic[entry] +=1

# Index from call_samples
df_samples = pd.read_csv('call_samples.ALL.panel.csv', delimiter='\t')
samples = df_samples.drop(columns=["Unnamed: 4", "Unnamed: 5"], axis = 1)
samples.head()
pop = samples['pop']

# Combined data and save to .csv
data = pd.DataFrame(
        {'0|1' : zero_one,
         '1|0' : one_zero,
         '1|1' : one_one,
        })
prelim = data.sum(axis = 1)
data_sum = np.array(prelim)

final = pd.DataFrame(
        {'0|1' : zero_one,
         '1|0' : one_zero,
         '1|1' : one_one,
         'sum_biallelic' : data_sum,
         'multiallelic' : multiallelic,  
        }, index=pop)
    
final.to_csv(file_name) 

#______________________________________________________________________________

file_name = "chr13.csv"

vcf = gzip.open('chr13.vcf.gz', "rt")
line = vcf.readline()
spline = line.split()
for line in range(253):
    line = vcf.readline()
spline = line.split()
header = spline
one_zero=[0]*(len(header)-9)

zero_one=[0]*(len(header)-9)

one_one=[0]*(len(header)-9)

zeros=[0]*(len(header)-9)

multiallelic=[0]*(len(header)-9)

for line in vcf:
        spline = line.split()
        
        genos = spline[9:]
        for entry in range(len(genos)):
                indgenos = genos[entry]
                if indgenos == '1|0':  
                    one_zero[entry] +=1  
                elif indgenos == '0|1':
                    zero_one[entry] +=1
                elif indgenos == '1|1':
                    one_one[entry] +=1
                elif indgenos == '0|0':
                    zeros[entry] +=1
                else: 
                    multiallelic[entry] +=1

# Index from call_samples
df_samples = pd.read_csv('call_samples.ALL.panel.csv', delimiter='\t')
samples = df_samples.drop(columns=["Unnamed: 4", "Unnamed: 5"], axis = 1)
samples.head()
pop = samples['pop']

# Combined data and save to .csv
data = pd.DataFrame(
        {'0|1' : zero_one,
         '1|0' : one_zero,
         '1|1' : one_one,
        })
prelim = data.sum(axis = 1)
data_sum = np.array(prelim)

final = pd.DataFrame(
        {'0|1' : zero_one,
         '1|0' : one_zero,
         '1|1' : one_one,
         'sum_biallelic' : data_sum,
         'multiallelic' : multiallelic,  
        }, index=pop)
    
final.to_csv(file_name) 

#______________________________________________________________________________

file_name = "chr14.csv"

vcf = gzip.open('chr14.vcf.gz', "rt")
line = vcf.readline()
spline = line.split()
for line in range(253):
    line = vcf.readline()
spline = line.split()
header = spline
one_zero=[0]*(len(header)-9)

zero_one=[0]*(len(header)-9)

one_one=[0]*(len(header)-9)

zeros=[0]*(len(header)-9)

multiallelic=[0]*(len(header)-9)

for line in vcf:
        spline = line.split()
        
        genos = spline[9:]
        for entry in range(len(genos)):
                indgenos = genos[entry]
                if indgenos == '1|0':  
                    one_zero[entry] +=1  
                elif indgenos == '0|1':
                    zero_one[entry] +=1
                elif indgenos == '1|1':
                    one_one[entry] +=1
                elif indgenos == '0|0':
                    zeros[entry] +=1
                else: 
                    multiallelic[entry] +=1

# Index from call_samples
df_samples = pd.read_csv('call_samples.ALL.panel.csv', delimiter='\t')
samples = df_samples.drop(columns=["Unnamed: 4", "Unnamed: 5"], axis = 1)
samples.head()
pop = samples['pop']

# Combined data and save to .csv
data = pd.DataFrame(
        {'0|1' : zero_one,
         '1|0' : one_zero,
         '1|1' : one_one,
        })
prelim = data.sum(axis = 1)
data_sum = np.array(prelim)

final = pd.DataFrame(
        {'0|1' : zero_one,
         '1|0' : one_zero,
         '1|1' : one_one,
         'sum_biallelic' : data_sum,
         'multiallelic' : multiallelic,  
        }, index=pop)
    
final.to_csv(file_name) 

#______________________________________________________________________________

file_name = "chr15.csv"

vcf = gzip.open('chr15.vcf.gz', "rt")
line = vcf.readline()
spline = line.split()
for line in range(253):
    line = vcf.readline()
spline = line.split()
header = spline
one_zero=[0]*(len(header)-9)

zero_one=[0]*(len(header)-9)

one_one=[0]*(len(header)-9)

zeros=[0]*(len(header)-9)

multiallelic=[0]*(len(header)-9)

for line in vcf:
        spline = line.split()
        
        genos = spline[9:]
        for entry in range(len(genos)):
                indgenos = genos[entry]
                if indgenos == '1|0':  
                    one_zero[entry] +=1  
                elif indgenos == '0|1':
                    zero_one[entry] +=1
                elif indgenos == '1|1':
                    one_one[entry] +=1
                elif indgenos == '0|0':
                    zeros[entry] +=1
                else: 
                    multiallelic[entry] +=1

# Index from call_samples
df_samples = pd.read_csv('call_samples.ALL.panel.csv', delimiter='\t')
samples = df_samples.drop(columns=["Unnamed: 4", "Unnamed: 5"], axis = 1)
samples.head()
pop = samples['pop']

# Combined data and save to .csv
data = pd.DataFrame(
        {'0|1' : zero_one,
         '1|0' : one_zero,
         '1|1' : one_one,
        })
prelim = data.sum(axis = 1)
data_sum = np.array(prelim)

final = pd.DataFrame(
        {'0|1' : zero_one,
         '1|0' : one_zero,
         '1|1' : one_one,
         'sum_biallelic' : data_sum,
         'multiallelic' : multiallelic,  
        }, index=pop)
    
final.to_csv(file_name) 

#______________________________________________________________________________

file_name = "chr16.csv"

vcf = gzip.open('chr16.vcf.gz', "rt")
line = vcf.readline()
spline = line.split()
for line in range(253):
    line = vcf.readline()
spline = line.split()
header = spline
one_zero=[0]*(len(header)-9)

zero_one=[0]*(len(header)-9)

one_one=[0]*(len(header)-9)

zeros=[0]*(len(header)-9)

multiallelic=[0]*(len(header)-9)

for line in vcf:
        spline = line.split()
        
        genos = spline[9:]
        for entry in range(len(genos)):
                indgenos = genos[entry]
                if indgenos == '1|0':  
                    one_zero[entry] +=1  
                elif indgenos == '0|1':
                    zero_one[entry] +=1
                elif indgenos == '1|1':
                    one_one[entry] +=1
                elif indgenos == '0|0':
                    zeros[entry] +=1
                else: 
                    multiallelic[entry] +=1

# Index from call_samples
df_samples = pd.read_csv('call_samples.ALL.panel.csv', delimiter='\t')
samples = df_samples.drop(columns=["Unnamed: 4", "Unnamed: 5"], axis = 1)
samples.head()
pop = samples['pop']

# Combined data and save to .csv
data = pd.DataFrame(
        {'0|1' : zero_one,
         '1|0' : one_zero,
         '1|1' : one_one,
        })
prelim = data.sum(axis = 1)
data_sum = np.array(prelim)

final = pd.DataFrame(
        {'0|1' : zero_one,
         '1|0' : one_zero,
         '1|1' : one_one,
         'sum_biallelic' : data_sum,
         'multiallelic' : multiallelic,  
        }, index=pop)
    
final.to_csv(file_name) 

#______________________________________________________________________________

file_name = "chr17.csv"

vcf = gzip.open('chr17.vcf.gz', "rt")
line = vcf.readline()
spline = line.split()
for line in range(253):
    line = vcf.readline()
spline = line.split()
header = spline
one_zero=[0]*(len(header)-9)

zero_one=[0]*(len(header)-9)

one_one=[0]*(len(header)-9)

zeros=[0]*(len(header)-9)

multiallelic=[0]*(len(header)-9)

for line in vcf:
        spline = line.split()
        
        genos = spline[9:]
        for entry in range(len(genos)):
                indgenos = genos[entry]
                if indgenos == '1|0':  
                    one_zero[entry] +=1  
                elif indgenos == '0|1':
                    zero_one[entry] +=1
                elif indgenos == '1|1':
                    one_one[entry] +=1
                elif indgenos == '0|0':
                    zeros[entry] +=1
                else: 
                    multiallelic[entry] +=1

# Index from call_samples
df_samples = pd.read_csv('call_samples.ALL.panel.csv', delimiter='\t')
samples = df_samples.drop(columns=["Unnamed: 4", "Unnamed: 5"], axis = 1)
samples.head()
pop = samples['pop']

# Combined data and save to .csv
data = pd.DataFrame(
        {'0|1' : zero_one,
         '1|0' : one_zero,
         '1|1' : one_one,
        })
prelim = data.sum(axis = 1)
data_sum = np.array(prelim)

final = pd.DataFrame(
        {'0|1' : zero_one,
         '1|0' : one_zero,
         '1|1' : one_one,
         'sum_biallelic' : data_sum,
         'multiallelic' : multiallelic,  
        }, index=pop)
    
final.to_csv(file_name) 

#______________________________________________________________________________

file_name = "chr18.csv"

vcf = gzip.open('chr18.vcf.gz', "rt")
line = vcf.readline()
spline = line.split()
for line in range(253):
    line = vcf.readline()
spline = line.split()
header = spline
one_zero=[0]*(len(header)-9)

zero_one=[0]*(len(header)-9)

one_one=[0]*(len(header)-9)

zeros=[0]*(len(header)-9)

multiallelic=[0]*(len(header)-9)

for line in vcf:
        spline = line.split()
        
        genos = spline[9:]
        for entry in range(len(genos)):
                indgenos = genos[entry]
                if indgenos == '1|0':  
                    one_zero[entry] +=1  
                elif indgenos == '0|1':
                    zero_one[entry] +=1
                elif indgenos == '1|1':
                    one_one[entry] +=1
                elif indgenos == '0|0':
                    zeros[entry] +=1
                else: 
                    multiallelic[entry] +=1

# Index from call_samples
df_samples = pd.read_csv('call_samples.ALL.panel.csv', delimiter='\t')
samples = df_samples.drop(columns=["Unnamed: 4", "Unnamed: 5"], axis = 1)
samples.head()
pop = samples['pop']

# Combined data and save to .csv
data = pd.DataFrame(
        {'0|1' : zero_one,
         '1|0' : one_zero,
         '1|1' : one_one,
        })
prelim = data.sum(axis = 1)
data_sum = np.array(prelim)

final = pd.DataFrame(
        {'0|1' : zero_one,
         '1|0' : one_zero,
         '1|1' : one_one,
         'sum_biallelic' : data_sum,
         'multiallelic' : multiallelic,  
        }, index=pop)
    
final.to_csv(file_name) 

#______________________________________________________________________________

file_name = "chr19.csv"

vcf = gzip.open('chr19.vcf.gz', "rt")
line = vcf.readline()
spline = line.split()
for line in range(253):
    line = vcf.readline()
spline = line.split()
header = spline
one_zero=[0]*(len(header)-9)

zero_one=[0]*(len(header)-9)

one_one=[0]*(len(header)-9)

zeros=[0]*(len(header)-9)

multiallelic=[0]*(len(header)-9)

for line in vcf:
        spline = line.split()
        
        genos = spline[9:]
        for entry in range(len(genos)):
                indgenos = genos[entry]
                if indgenos == '1|0':  
                    one_zero[entry] +=1  
                elif indgenos == '0|1':
                    zero_one[entry] +=1
                elif indgenos == '1|1':
                    one_one[entry] +=1
                elif indgenos == '0|0':
                    zeros[entry] +=1
                else: 
                    multiallelic[entry] +=1

# Index from call_samples
df_samples = pd.read_csv('call_samples.ALL.panel.csv', delimiter='\t')
samples = df_samples.drop(columns=["Unnamed: 4", "Unnamed: 5"], axis = 1)
samples.head()
pop = samples['pop']

# Combined data and save to .csv
data = pd.DataFrame(
        {'0|1' : zero_one,
         '1|0' : one_zero,
         '1|1' : one_one,
        })
prelim = data.sum(axis = 1)
data_sum = np.array(prelim)

final = pd.DataFrame(
        {'0|1' : zero_one,
         '1|0' : one_zero,
         '1|1' : one_one,
         'sum_biallelic' : data_sum,
         'multiallelic' : multiallelic,  
        }, index=pop)
    
final.to_csv(file_name) 

#______________________________________________________________________________

file_name = "chr20.csv"

vcf = gzip.open('chr20.vcf.gz', "rt")
line = vcf.readline()
spline = line.split()
for line in range(253):
    line = vcf.readline()
spline = line.split()
header = spline
one_zero=[0]*(len(header)-9)

zero_one=[0]*(len(header)-9)

one_one=[0]*(len(header)-9)

zeros=[0]*(len(header)-9)

multiallelic=[0]*(len(header)-9)

for line in vcf:
        spline = line.split()
        
        genos = spline[9:]
        for entry in range(len(genos)):
                indgenos = genos[entry]
                if indgenos == '1|0':  
                    one_zero[entry] +=1  
                elif indgenos == '0|1':
                    zero_one[entry] +=1
                elif indgenos == '1|1':
                    one_one[entry] +=1
                elif indgenos == '0|0':
                    zeros[entry] +=1
                else: 
                    multiallelic[entry] +=1

# Index from call_samples
df_samples = pd.read_csv('call_samples.ALL.panel.csv', delimiter='\t')
samples = df_samples.drop(columns=["Unnamed: 4", "Unnamed: 5"], axis = 1)
samples.head()
pop = samples['pop']

# Combined data and save to .csv
data = pd.DataFrame(
        {'0|1' : zero_one,
         '1|0' : one_zero,
         '1|1' : one_one,
        })
prelim = data.sum(axis = 1)
data_sum = np.array(prelim)

final = pd.DataFrame(
        {'0|1' : zero_one,
         '1|0' : one_zero,
         '1|1' : one_one,
         'sum_biallelic' : data_sum,
         'multiallelic' : multiallelic,  
        }, index=pop)
    
final.to_csv(file_name) 

#______________________________________________________________________________

file_name = "chr21.csv"

vcf = gzip.open('chr21.vcf.gz', "rt")
line = vcf.readline()
spline = line.split()
for line in range(253):
    line = vcf.readline()
spline = line.split()
header = spline
one_zero=[0]*(len(header)-9)

zero_one=[0]*(len(header)-9)

one_one=[0]*(len(header)-9)

zeros=[0]*(len(header)-9)

multiallelic=[0]*(len(header)-9)

for line in vcf:
        spline = line.split()
        
        genos = spline[9:]
        for entry in range(len(genos)):
                indgenos = genos[entry]
                if indgenos == '1|0':  
                    one_zero[entry] +=1  
                elif indgenos == '0|1':
                    zero_one[entry] +=1
                elif indgenos == '1|1':
                    one_one[entry] +=1
                elif indgenos == '0|0':
                    zeros[entry] +=1
                else: 
                    multiallelic[entry] +=1

# Index from call_samples
df_samples = pd.read_csv('call_samples.ALL.panel.csv', delimiter='\t')
samples = df_samples.drop(columns=["Unnamed: 4", "Unnamed: 5"], axis = 1)
samples.head()
pop = samples['pop']

# Combined data and save to .csv
data = pd.DataFrame(
        {'0|1' : zero_one,
         '1|0' : one_zero,
         '1|1' : one_one,
        })
prelim = data.sum(axis = 1)
data_sum = np.array(prelim)

final = pd.DataFrame(
        {'0|1' : zero_one,
         '1|0' : one_zero,
         '1|1' : one_one,
         'sum_biallelic' : data_sum,
         'multiallelic' : multiallelic,  
        }, index=pop)
    
final.to_csv(file_name) 

#______________________________________________________________________________

file_name = "chr22.csv"

vcf = gzip.open('chr22.vcf.gz', "rt")
line = vcf.readline()
spline = line.split()
for line in range(253):
    line = vcf.readline()
spline = line.split()
header = spline
one_zero=[0]*(len(header)-9)

zero_one=[0]*(len(header)-9)

one_one=[0]*(len(header)-9)

zeros=[0]*(len(header)-9)

multiallelic=[0]*(len(header)-9)

for line in vcf:
        spline = line.split()
        
        genos = spline[9:]
        for entry in range(len(genos)):
                indgenos = genos[entry]
                if indgenos == '1|0':  
                    one_zero[entry] +=1  
                elif indgenos == '0|1':
                    zero_one[entry] +=1
                elif indgenos == '1|1':
                    one_one[entry] +=1
                elif indgenos == '0|0':
                    zeros[entry] +=1
                else: 
                    multiallelic[entry] +=1

# Index from call_samples
df_samples = pd.read_csv('call_samples.ALL.panel.csv', delimiter='\t')
samples = df_samples.drop(columns=["Unnamed: 4", "Unnamed: 5"], axis = 1)
samples.head()
pop = samples['pop']

# Combined data and save to .csv
data = pd.DataFrame(
        {'0|1' : zero_one,
         '1|0' : one_zero,
         '1|1' : one_one,
        })
prelim = data.sum(axis = 1)
data_sum = np.array(prelim)

final = pd.DataFrame(
        {'0|1' : zero_one,
         '1|0' : one_zero,
         '1|1' : one_one,
         'sum_biallelic' : data_sum,
         'multiallelic' : multiallelic,  
        }, index=pop)
    
final.to_csv(file_name) 

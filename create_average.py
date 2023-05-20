f = open("result.txt")
data = f.read()

m = {}
result = data.split('\n')
for line in result :
    k = line.split()[0]
    v = line.split()[1]
    print(k)
    print(v)
    if k not in m: 
        m[k] = 0
    m[k] += float(v)

print(m.keys())
m['1'] /= 100
m['2'] /= 100
m['3'] /= 100
m['4'] /= 100


w = open("average.txt","a")
for k in m : 
    w.write("Average when k = " + k + " : " + str(m[k]) + "\n")
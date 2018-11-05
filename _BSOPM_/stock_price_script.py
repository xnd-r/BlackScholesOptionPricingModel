from subprocess import Popen, PIPE
import sys
import os

now = datetime.datetime.now()
dt = now.strftime("%d-%m-%Y-%H-%M-%S")

res_array = []
file  = sys.argv[1]
step  = sys.argv[2]
count = sys.argv[3]
name = os.path.join(os.getcwd(), file)

def one_res(name):
    proc = Popen(name, shell=True, stdout=PIPE, stderr=PIPE)
    proc.wait()
    res = proc.communicate()
    if proc.returncode:
        return res[1]
    return res[0]
	
def fill_array(count, step):
    for i in range(count):
        res_array.append(float(one_res(name + " " + str(step * i))))
		
def write_log():
	with open(str(dt) + 'An_Sol_res_array' + '.log', 'a') as f:
		for i in count:
			f.write(str(i * step) + "	" + str(res_array[i]))

if __name__ == '__main__':
	fill_array(int(count), int(step))
	write_log()
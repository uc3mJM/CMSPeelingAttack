from cProfile import label
from math import e
from math import exp
from math import log as ln
import pandas
import matplotlib.pyplot as plt
import re


from LogFile import LogFile


# Parameters for the teoric simulations
WHITE_BOX = 0
BLACK_BOX_PAIRS = 1
BLACK_BOX_SINGLE = 2
VARIANTS = [WHITE_BOX,BLACK_BOX_SINGLE,BLACK_BOX_PAIRS]

COLORS = ['red','green','blue']
NAMES = ['White Box','Black Box Pairs','Black Box Single']

k = None
cs = None
variant = None


def setK(kk, n, w):

    global k
    global cs
    k = kk

    #cs = ln(2) / k
    #print(cs)
    #cs = k / ln(2) * n
    #print(cs)
    #cs = n / (k * w)
    cs = n / (w*k)
    #cs = 0.22
    #print(cs)
    

def setVarian(v):
    global variant
    variant = v

def round(ps,pf,cf):
	global cs
	global variant
	
	λ1 = k*cf*(1-pf)**(k-1)
	λ2 = k*cs*(1-ps)**(k-1)
	if (variant == WHITE_BOX):
		ps = e**(-λ1)
	elif (variant == BLACK_BOX_SINGLE):
		ps = e**(-λ1)*e**(-λ2)
	elif (variant == BLACK_BOX_PAIRS):
		ps = e**(-λ1)*e**(-λ2)*(1+λ2)
	pf = e**(-λ2)
	return (ps,pf)

def getCoreDensity(cf):
	(ps,pf) = (0,0)
	for i in range(1000):
		(ps,pf) = round(ps,pf,cf)
	return (1-ps)**k

def find_threshold():
	(low,high) = (0,500)
	while(high - low > 0.0001):
		mid = (low + high)/2
		if (getCoreDensity(mid) < 0.000001):
			low = mid
		else:
			high = mid
	return low


def plot_data_false_positives(log_df, title):

    fig = plt.figure()


    variant = [BLACK_BOX_SINGLE]

    for v in variant:
        setVarian(v)
        #results[v] = plt.Line2D([(x/100) for x in range(100,600)], [(1-getCoreDensity(y/100*2**k*cs)) for y in range(100,600)], color=COLORS[v])
        #plt.plot([(x/100) for x in range(100,600)], [(1-getCoreDensity(y/100*2**k*cs)) for y in range(100,600)], color=COLORS[v], label=NAMES[v])
        plt.plot([(x/100) for x in range(100,600)], [(1-getCoreDensity(y/100*((1/pow(1-exp(-cs*k),k))*cs))) for y in range(100,600)], color=COLORS[v], label="Theoretical Threshold")
        #ax.add_line(results[v])
		#print(results[v])

    plt.plot((log_df['false_pos']/(log_df['total_inserted']/log_df['iterations'])), log_df['success_iter_perc']/100, color = "purple", marker='o', label="SSP")
    plt.plot((log_df['false_pos']/(log_df['total_inserted']/log_df['iterations'])), log_df['average_retrieved']/(log_df['total_inserted']/log_df['iterations']),color = "green", marker='v', label="ESP Average")
    #plt.plot((log_df['false_pos']/(log_df['total_inserted']/log_df['iterations'])), log_df['worst_retrieved']/(log_df['total_inserted']/log_df['iterations']), color = "red", marker='x', label="ESP Worst Case")





    plt.legend()
    plt.xlabel ("ratio of false positives to true positives (f/n)")
    plt.ylabel ("Success rate")

    plt.grid()
    plt.title(title)
    plt.show()

def plot_data_insertions(log_df, title):

    fig = plt.figure()

    plt.plot(log_df['total_stored']/log_df['iterations'], log_df['success_iter_perc']/100, color = "purple", marker='o', label="SSP")
    #plt.plot(log_df['total_stored']/log_df['iterations'], log_df['average_retrieved']/(log_df['total_inserted']/log_df['iterations']), color = "green", marker='o', label="ESP Average")
    plt.plot(log_df['total_inserted']/log_df['iterations'], log_df['success_iter_perc']/100, color = "green", marker='o', label="SSP unique")




    plt.xlabel ("Insertions numbers")
    plt.ylabel ("Success rate")

    plt.legend()
    plt.grid()
    plt.title(title)
    plt.show()

    return

def plot_data_universe(log_df, kk):

    fig = plt.figure()

    
    plt.plot(log_df['universe'], log_df['success_iter_perc'])
    plt.xlabel ("Universe size")
    plt.ylabel ("Success rate")

    plt.grid()
    plt.show()

    return


# Main method
def main():

    # logs File to read
    # Directory to write the logs
    directory = './logs/'
    logs_file = "result_m4096_d5_u16777216_hmd5_i10_fB"
    #logs = LogFile(logs_file)

    # for the theoric simulation
    kk = 5
    n =  4096
    w = 4096
    u = 16

    title = 'd = ' + str(kk) + ', m = ' + str(w) + ', n = ' + str(n)
    title_i = 'd = ' + str(kk) + ', m = ' + str(w) + ', mask = ' + str(u)

    setK(kk, n, w)

    log_df = pandas.read_csv(directory + logs_file , delimiter =";")
    #plot_data_universe(log_df)
    plot_data_insertions(log_df, title_i)
    #plot_data_false_positives(log_df, title)

if __name__ == "__main__":
    main()

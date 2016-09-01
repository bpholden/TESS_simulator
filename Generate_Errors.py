import numpy as np
from consts import JD0
import scipy.constants as sc

# needs to take output from sim_night
# generate the skeleton of a vels file for each target
# a file will be MJD error I2 Counts exptime


def compute_precision(i2counts,cols):

    blue = (cols < 1.2)
    red = (cols >= 1.2)
    gerr = .07/np.log(10) # fractional err
    merr = .07/np.log(10)

    A_gk = 4.47
    B_gk = -1.58
    A_m = 4.14
    B_m = -1.73

    if blue:
        mean = (np.log10(i2counts) - A_gk)/B_gk
    else:
        mean = (np.log10(i2counts) - A_m)/B_m                        
    #    means = np.zeros_like(i2counts)
    #    means[blue] += (np.log10(i2counts[blue]) - A_gk)/B_gk
    #    means[red] += (np.log10(i2counts[red]) - A_m)/B_m

    try:
        devs = np.random.normal(size=len(i2counts))
    except:
        devs = np.random.normal(size=1)
    if blue:
        devs *= gerr
    else:
        devs *= merr
    devs += mean
    return 10**devs[0]

def jitter(cols):
    # just use median
    
    jitter = np.zeros_like(cols)
    jitter[((cols > 0.4) & (cols < 0.7))] += 2.3 + 17.4*0.02
    jitter[((cols > 0.7) & (cols < 1.0))] += 2.1 + 4.7 *0.02
    jitter[((cols > 1.0) & (cols < 1.3))] += 1.6 - 0.003 * 0.3 
    jitter[((cols > 1.3) & (cols < 1.6))] += 2.1 + 2.7 * 0.5


    # if cols < 0.7:
    #     jitter =  2.3 + 17.4*0.02
    # elif cols > 0.7 and cols < 1.0:
    #     jitter = 2.1 + 4.7 *0.02
    # elif cols > 1.0 and cols < 1.3:
    #     jitter =  1.6 - 0.003 * 0.3 
    # else:
    #     jitter =  2.1 + 2.7 * 0.5
        
    return jitter

def generate_jitters(data):
    jitters = jitter(data['b-v'])
    known = dict()
    known['185144'] = 2.491
    known['10700'] = 2.273
    known['9407']=2.820
    known['219134'] = 1.570
    for ky in known.keys():
        jitters[data['starname'] == ky] = known[ky]
    jitters /= 2.
    return jitters

def sinnoise(dates,totjitter,periods,phases):
    rms = np.sqrt(totjitter**2 / len(periods))
#    Ks = np.zeros(len(periods))
    Ks = np.abs(np.random.normal(scale=rms,size=len(periods)))
    if type(dates) == float:
        currentphases = np.empty([1,len(Ks)])
        ldate = 1
    else:
        try:
            currentphases = np.empty([len(dates),len(Ks)])
            ldate = len(dates)
        except:
            currentphases = np.empty([1,len(Ks)])
            ldate = 1
    for i in range(0,len(Ks)):
       currentphases[:,i]=np.mod((np.mod(dates-JD0, periods[i])/periods[i])+phases[i],1)
    fs = currentphases * 2*sc.pi
    noise = np.sum(Ks*np.cos(fs),axis=1) + np.random.normal(scale=rms,size=ldate)
    return noise

def compute_real_uncertainty(i2counts,cols,star_jitter):

    precision = compute_precision(i2counts,cols)
    #    true= np.zeros_like(precision)
    true = precision**2
    true += star_jitter**2
    true = np.sqrt(true)
    true += 1.0
    dev = true*np.random.normal(size=1)
    
    return precision, dev, true

def compute_real_uncertainty_sinnoise(i2counts,date,star_tab):

    star_jitter = jitter(star_tab['b-v'])
        
    precision = compute_precision(i2counts,star_tab['b-v'])
    periods = [star_tab['shortperiod'],star_tab['longperiod']]
    phases = [star_tab['shortphase'],star_tab['longphase']]
    
    dev = np.random.normal(scale=np.sqrt(precision**2+1),size=1)
    if star_jitter > 0.0:
        dev += sinnoise(date,star_jitter,periods,phases)

    tot_error = np.sqrt(star_jitter**2 + precision**2 + 1)
    return precision, dev, tot_error

if __name__ == "__main__":
    file='../Datafiles/newgoogledex.csv'
    data=ascii.read(file,format='csv')
    jitters = jitter(data['b-v'])
    known = dict()
    known['185144'] = 2.491
    known['10700'] = 2.273
    known['9407']=2.820
    known['219134'] = 1.570
    for ky in known.keys():
        jitters[data['starname'] == ky] = known[ky]
    out = Table([data['starname'],jitters])
#    ascii.write(out,delimiter=",")
                     

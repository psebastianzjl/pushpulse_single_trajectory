import os
import sys
import numpy as np
import re


current_state = int(sys.argv[1]) # currently populated state (S0 = 1)

def eV2Ha(eV):
    """
    Function that converts energy from electron Volt (eV) to atomic units (hartree).
    Parameters
    ----------
    eV : float
    Energy in eV

    Returns
    -------
    Energy in atomic units (hartree)
    """
    Ha = eV / 27.21138602
    return Ha


def fs2aut(fs):
    """
     Function that converts time from femto-seconds to atomic units.
     Parameters
     ----------
     fs : float
     Time in femto-seconds that is to be converted to atomic units of time

     Returns
     -------
     Time in atomic units
     """
    aut = fs * 41.341373336
    return aut

def ricc2_energy(fname, model):
    search_string = "Final " + re.escape(model.upper()) + " energy"
    gs_energy = search_file(fname, search_string)
    split_columns(gs_energy, col=5, convert=np.float64)
    ex_energy = search_file(fname, "Energy:")
    ex_energy = split_columns(ex_energy, 1, convert=np.float64)
    energy = np.repeat(gs_energy, len(ex_energy) + 1)
    energy[1:] = energy[1:] + ex_energy
    return energy

def ricc2_oscill(fname):
    """ Read oscillator strengths from STDOUT file. """
    oscill = search_file(fname, r"oscillator strength \(length gauge\)")
    oscill = split_columns(oscill, col=5)
    try:
        oscill_ex = search_file(fname,r"oscillator strength \(velocity gauge\)")
        oscill_ex = split_columns(oscill_ex, col=5)
        return np.array(oscill, dtype=np.float64), np.array(oscill_ex, dtype=np.float64)
    except:
        return np.array(oscill, dtype=np.float64)

def search_file(file, search, after=0, max_res=None, close=True, stop_at=None):
    """ Find all occurrences of a regex in file.

    Parameters
    ----------
    file : str
        Name of the file
    search : str
        Regular expression to be searched
    after : int, optional
        Number of lines to be skipped after regular expression
    max_res : int, optional
        Number of maximum occurrences until loop is exited
    close : bool, optional
    stop_at : str, optional
        Regular expression that triggers loop to be exited
    """
    cfile = open_if_needed(file)
    search_reg = re.compile(search)
    if stop_at is not None:
        stop_reg = re.compile(stop_at)
    values = []
    n_hit = 0
    for line in cfile:
        if stop_at is not None:
            if stop_reg.search(line):
                break
        if search_reg.search(line):
            n_hit = n_hit + 1
            if after == 0:
                values.append(line.rstrip())
            else:
                for _ in range(after):
                    line = next(cfile).rstrip()
                    values.append(line.rstrip())
            if max_res is not None:
                if n_hit >= max_res:
                    break
    if close:
        cfile.close()
    if not values:
        raise ValueError("No matches for {} in file {}".format(
            search_reg.pattern, cfile.name))
    return values
    
    def split_columns(split_list, col=None, split=str.split, convert=None):
    """ Split all sublists of split_list and take specified columns.

    Parameters
    ----------
    split_list : list or array,
    col : int or list, int, optional
        Columns to be read
    split : func, optional
        Function to split lines
    convert : type, optional
        Function to convert values in split_list
    """
    if col is None and convert is not None:
        split_list = [convert(val) for val in split_list]
    elif isinstance(col, int):
        for i, line in enumerate(split_list):
            if convert is None:
                split_list[i] = split(line)[col]
            else:
                split_list[i] = convert(split(line)[col])
    elif isinstance(col, list):
        for i, line in enumerate(split_list):
            line = split(line)
            if convert is None:
                split_list[i] = [line[ci] for ci in col]
            else:
                split_list[i] = [convert(line[ci]) for ci in col]
    return split_list



def envelope(w, tau):
    """
    Function to predict the intensity of the signal using a Gaussian envelope.
    Parameters
    ----------
    w : float
        Energy of the probe pulse in atomic units
    tau : float
        Width of the envelope in atomic units

    Returns
    -------
    Gaussian envelope according to the parameters
    """
    E = np.exp(-(w * tau) ** 2 / 4.) * tau
    return E

def signal_sum(eV, Ha, Ha_0, f, fs, delta=False):
    """

    Parameters
    ----------
    eV : float
        Energy of the probe pulse in eV
    Ha : float
        Energy of the higher lying states in atomic units
    Ha_0 : float
        Energy of the populated state in atomic units
    f :  float
        Oscillator strength for the transition between the populated state and the electronic ground state
    fs : float
        Width of the laser pulse envelope
    delta :  bool, optional
        Switch to pulse width of 0 fs

    Returns
    -------

    """
    laser = eV2Ha(eV)
    tau = fs2aut(fs)
    S = []
    for (i, val) in enumerate(Ha):
        w_in = val - Ha_0
        dw = laser - w_in
        v2 = f[i] / 2. * 3 / w_in
        if not delta:
            S += [envelope(dw, tau) ** 2. * v2]
        else:
            S += 1. * v2
    return S / np.sum(S) # Normalize sum to 1 to ensure hop



w_pu = 2.8 # Energy of the push pulse in eV
tau = 1 # Width of the envelope in fs
num_states = 31 #Total number of states; GS + excited states
print(str(w_pu) + ' eV as push pulse')


qm_log = './qmdir/qm.log' # Directory containing the output of a TURBOMOLE calculation (ADC(2))
eaf = np.zeros(num_states)  # t, osc str between current and higher lying states

exE = turbomole.ricc2_energy(qm_log, 'mp2')
exE = np.append(exE, current_state)
_, eaf_temp = turbomole.ricc2_oscill(qm_log)
start_index = -num_states + 1  ### Allow for adjustments in same loop
end_index = 0
for i in range(1, int(exE[-1])):  ### S1 starts at position 0
    start_index += (num_states - i)
    end_index += (num_states - i - 1)
eaf[int(exE[-1]):] = [float(i) for i in eaf_temp[start_index:end_index]]
#print('eaf: ', eaf)
Def = signal_sum(w_pu, exE[int(exE[-1]):-1], exE[int(exE[-1] - 1)], eaf[int(exE[-1]):], tau)   ### Doorway function

for i in range(30 - len(Def)): # add zeros in front of array until its length is 30 for easier handling
    Def = np.insert(Def, 0, 0)
Def2 = Def * np.random.rand(30)
index = np.argmax(Def2) #Find state with the highest probability; index = 0 corresponds to S1

print('Current state, final state: ', int(exE[-1]), index + 2)
print('Energy Diff: ', 27.2114 * (exE[index + 1] - exE[int(exE[-1] - 1)]))
print('Osc Str: ', eaf[int(index) + 1])


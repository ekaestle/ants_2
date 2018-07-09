import numpy as np
from copy import deepcopy
from obspy.signal.filter import envelope
#from obspy.signal.util import next_pow_2
#from scipy import fftpack
from scipy.signal import iirfilter, zpk2sos, sosfilt
#from ants_2.tools.windows import my_centered

def bandpass(freqmin, freqmax, df, corners=4):
    """
    From obspy with modification.

    Butterworth-Bandpass Filter.

    Filter data from ``freqmin`` to ``freqmax`` using ``corners``
    corners.
    The filter uses :func:`scipy.signal.iirfilter` (for design)
    and :func:`scipy.signal.sosfilt` (for applying the filter).

    :type data: numpy.ndarray
    :param data: Data to filter.
    :param freqmin: Pass band low corner frequency.
    :param freqmax: Pass band high corner frequency.
    :param df: Sampling rate in Hz.
    :param corners: Filter corners / order.
    :param zerophase: If True, apply filter once forwards and once backwards.
        This results in twice the filter order but zero phase shift in
        the resulting filtered trace.
    :return: Filtered data.
    """
    fe = 0.5 * df
    low = freqmin / fe
    high = freqmax / fe
    # raise for some bad scenarios
    if high > 1:
        high = 1.0
        msg = "Selected high corner frequency is above Nyquist. " + \
              "Setting Nyquist as high corner."
        warnings.warn(msg)
    if low > 1:
        msg = "Selected low corner frequency is above Nyquist."
        raise ValueError(msg)
    z, p, k = iirfilter(corners, [low, high], btype='band',
                        ftype='butter', output='zpk')
    sos = zpk2sos(z, p, k)
    return sos


def whiten_taper(ind_fw1,ind_fw2,npts,taper_samples):
    
    
    if ind_fw1 - taper_samples >= 0:
        i_l = ind_fw1 - taper_samples
    else:
        i_l = 0
        print('** Could not fully taper during whitening. Consider using a \
            smaller frequency range for whitening.')

    if ind_fw2 + taper_samples < npts:
        i_h = ind_fw2 + taper_samples
    else:
        i_h = npts - 1
        print('** Could not fully taper during whitening. Consider using a \
            smaller frequency range for whitening.')
    
    
    taper_left = np.linspace(0.,np.pi/2,ind_fw1-i_l)
    taper_left = np.square(np.sin(taper_left))
    
    taper_right = np.linspace(np.pi/2,np.pi,i_h-ind_fw2)
    taper_right = np.square(np.sin(taper_right))
    
    taper = np.zeros(npts)
    taper[ind_fw1:ind_fw2] += 1.
    taper[i_l:ind_fw1] = taper_left
    taper[ind_fw2:i_h] = taper_right

    return taper


def whiten(spec,sampling_rate,freq1,freq2,taper_samples,white_waterlevel,whitening_taper):
    
    # zeropadding should make things faster
#    n_pad = next_pow_2(tr.stats.npts)
#
#    data = my_centered(tr.data,n_pad)
    
    if whitening_taper:

        freqaxis=np.fft.rfftfreq((len(spec)-1)*2,sampling_rate)
        # the freq axis has a one-sample error if its length is odd
        # this does hardly influence the taper, so I ignore it here

        ind_fw = np.where( ( freqaxis > freq1 ) & ( freqaxis < freq2 ) )[0]

        if len(ind_fw) == 0:
            return(np.zeros((len(spec)-1)*2))

        ind_fw1 = ind_fw[0]
    
        ind_fw2 = ind_fw[-1]
    
        # Build a cosine taper for the frequency domain
        #df = 1/(tr.stats.npts*tr.stats.delta)
    
        # Taper 
        white_tape = whiten_taper(ind_fw1,ind_fw2,len(freqaxis),taper_samples)
    
    if white_waterlevel:
        # Don't divide by 0
        tol = np.max(np.abs(spec)) / 1e5
        if tol == 0.0:
            spec = np.exp(1j * np.angle(spec))
        else:
            spec /= (np.abs(spec)+tol)
    else:
        # whiten. This elegant solution is from MSNoise: (but the above is faster)
        spec =  np.exp(1j * np.angle(spec))
    
    if whitening_taper:
        spec *= white_tape
    
    return spec
    # Go back to time domain
    # Difficulty here: The time fdomain signal might no longer be real.
    # I don't think it actually makes a difference in the result (Emanuel)
    # Hence, irfft cannot be used.
#    spec_neg = np.conjugate(spec)[::-1]
#    spec = np.concatenate((spec,spec_neg[1:-1]))
#
#    tr.data = np.real(np.fft.ifft(spec))
    
    
def cap(data,cap_thresh):
    
    std = np.std(data*1.e6)
    gllow = cap_thresh * std * -1
    glupp = cap_thresh * std
    return np.clip(data*1.e6,gllow,glupp)/1.e6

    #return tr
    
def ram_norm(data,sampling_rate,winlen,prefilt=None):
    
    data_orig = deepcopy(data)
    hlen = int(winlen*sampling_rate/2.)

    if 2*hlen >= len(data):
        return np.zeros(len(data))


    weighttrace = np.zeros(len(data))
    
    if prefilt is not None:
        sos = bandpass(freqmin=prefilt[0],freqmax=prefilt[1],
                df=sampling_rate,corners=prefilt[2])
        temp = sosfilt(sos,data)
        data = sosfilt(sos,temp[::-1])[::-1]
        
    envlp = envelope(data)

    for n in range(hlen,len(data)-hlen):
        weighttrace[n] = np.sum(envlp[n-hlen:n+hlen+1]/(2.*hlen+1))
        
    weighttrace[0:hlen] = weighttrace[hlen]
    weighttrace[-hlen:] = weighttrace[-hlen-1]
    
    return data_orig.data / weighttrace
   

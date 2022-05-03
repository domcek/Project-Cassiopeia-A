"""
The script takes the data measurements in radio from Parley et al (2014) 
and estimates the 2009 radio image flux at 4.8 GHz in 2009 using 
Trotter et al (2017)
"""
import lmfit
import matplotlib.pyplot as plt
# import matplotlib.ticker as ticker
# import mpld3
import numpy as np
from lmfit import minimize, Parameters


def setup(ax):
    # ax.spines['right'].set_color('none')
    # ax.spines['left'].set_color('none')
    # ax.yaxis.set_major_locator(ticker.NullLocator())
    # ax.spines['top'].set_color('none')
    # ax.xaxis.set_ticks_position('bottom')
    # ax.set_xlim(0, 5)
    # ax.set_ylim(0, 1)
    # ax.tick_params(labeltop=True, labelright=True)
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.label.set_size(20)
    ax.yaxis.label.set_size(20)
    ax.tick_params(axis='both', which='major', width=1.00, length=8, direction='in', labelsize=18)
    ax.tick_params(axis='both', which='minor', width=0.75, length=4, direction='in', labelsize=12)
    # ax.tick_params(axis='y',which='major', width=1.00, length=10,direction='in')#, labelsize=22)
    # ax.tick_params(axis='y',which='minor', width=0.75, length=5,direction='in', labelsize=22)
    # ax.patch.set_alpha(0.0)
    # ax.xaxis.set_major_locator(ticker.MultipleLocator(1.00))
    # ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.25))
    # ax.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
    # ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
    plt.rcParams.update({'font.size': 18, 'axes.labelsize': 20, 'ytick.labelsize': 20})


# # Radio calibration
def linear_res(linear_pars, x, data=None, eps=None):
    # unpack parameters:
    # extract .value attribute for each parameter
    # alfa = parvals['alfa']
    # k = parvals['k']
    parvals = linear_pars.valuesdict()
    model = -parvals['alfa'] * x + parvals['k']

    if data is None:
        return model
    if eps is None:
        return (model - data)
    return (model - data) / eps


def linear_mo(x, alfa, k):
    return -alfa * x + k


outpath = '../data/intermediate_products/'

lin_model = lmfit.Model(linear_mo)

linear_pars = Parameters()
linear_pars.add('alfa', value=0.77)
linear_pars.add('k', value=1)


perley = np.genfromtxt('../data/original/perley.txt', unpack=True)

out = minimize(linear_res, linear_pars, args=(np.log10(perley[0][:21]), np.log10(perley[-1][:21])))

print(lmfit.fit_report(out.params))
# print(out.params['alfa'].value, out.params['k'].value)
# print(lin_model.independent_vars)
f1ghz = 1**(-1 * out.params['alfa'].value) * 10**out.params['k'].value
f1_4ghz = 1.4 ** (-1 * out.params['alfa'].value) * 10 ** out.params['k'].value
f4_8ghz = 4.8 ** (-1 * out.params['alfa'].value) * 10 ** out.params['k'].value
print('--------------------------------------')
print('Value at 1 GHz=', f1ghz)
print('Value at 1.4 GHz=', f1ghz)
print('Value at 4.8 GHz=', f4_8ghz)
print('--------------------------------------')

c = f4_8ghz
d = f1_4ghz
print('year 2014 , flux density ', c)
for i in range(2013, 2008, -1):
    c = c * 1.0067  # 0.67 % decrease per year (Trotter et al 2017)
    d = d * 1.0067
    print('4.8 GHz: year', i, ', flux density ', c)
    print('1.4 GHz: year', i, ', flux density ', d)


x_model = np.linspace(perley[0][0], perley[0][21], 5)
y_model = x_model**(-1 * out.params['alfa'].value) * 10**out.params['k'].value

# plt.clf()
fig, ax1 = plt.subplots(figsize=(12, 8), dpi=150)  # ,dpi=100)
setup(ax1)
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_xlabel('Frequency [GHz]')
ax1.set_ylabel('Flux density [Jy]')
ax1.set_xlim(5e-2, 5)
ax1.plot(x_model, y_model, label=r'$\alpha$ =' + str(round(-1 * out.params['alfa'].value, 3)) + r'$\pm$' + str(round(out.params['alfa'].stderr, 3)))
ax1.plot(perley[0], perley[-1], '+', label='Perley & Butler (2014)', markersize=10)
ax1.legend()

ax1.set_title('radio calibration perley \n flux density (2009) = ' + str(round(c, 1)) + ' [fading 0.67%]')
plt.savefig(outpath + '2_radio_calibration.pdf', bbox_inches='tight', dpi=100)
plt.show()

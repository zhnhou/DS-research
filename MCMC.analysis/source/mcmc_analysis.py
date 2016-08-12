import numpy as np
import pandas as pd
import os
from scipy.interpolate import *
from scipy.special import *
from scipy.ndimage.filters import gaussian_filter
from sklearn import mixture

import matplotlib.pyplot as plt



class mcmc_analysis(object):
    def __init__(self, chain_csv_file=None, chain_path=None, feature=None, chain_prefix=None, is_cosmomc=False, num_burnin=1000, num_chain=8):
        
        if (chain_csv_file is None and chain_path is None):
            print "Either 'chain_csv_file' or 'chain_path' should be specified."
            print "Stop and exit."
            os._exit(0)
        self.is_cosmomc = is_cosmomc
        self.chain = dict()
        self.num_burnin = num_burnin

        if (chain_csv_file is None):
            if (os.path.isdir(chain_path)):

                self.chain_path = chain_path

                if (self.is_cosmomc):
                    self.chain_prefix = chain_prefix
                    self.read_chain_cosmomc(parameter=feature, num_chain=num_chain)
                else:

                    # check if weight is included
                    # if not, we assume uniform weight
    
                    if ('weight' in os.listdir(chain_path)):
                        self.has_weight = True
                    else:
                        self.has_weight = False

                    if feature is None:
                        parameter = os.listdir(chain_path)
                    else:
                        if (self.has_weight):
                            feature.append('weight')
                            parameter = np.unique(feature).tolist()
                    try:
                        parameter.remove('paramnames')
                    except ValueError:
                        pass

                    self.parameter = parameter
                    self.read_chain_single_column()

        else:
            if (os.path.isfile(chain_csv_file)):
                df = pd.read_csv(chain_csv_file)
        
        self.num_sample = dict()
        for param in self.parameter:
            self.num_sample[param] = np.shape(self.chain[param])[0]

    def reproduce_chain_step(self):
        if self.has_weight:
            self.chain['step'] = np.cumsum(self.chain['weight'])

    def read_chain_cosmomc(self, parameter=None, num_chain=8):

        self.chain_original = dict()

        self.num_element = np.zeros(num_chain, dtype=np.int)
    
        tmp = np.loadtxt(self.chain_path+'/'+self.chain_prefix+'.paramnames', dtype='str')[:,0]
        pname = np.insert(tmp, 0, ['weight', 'loglike'])

        if parameter is None:
            self.parameter = pname
        else:
            try:
                parameter.remove('weight')
            except ValueError:
                pass

            self.parameter = np.insert(parameter, 0, 'weight')

        self.has_weight=True

        for i in np.arange(num_chain)+1:
            filename = self.chain_path+'/'+self.chain_prefix+'_'+str(i)+'.txt'
            print 'reading '+self.chain_prefix+'_'+str(i)+'.txt'
            tmp = np.loadtxt(filename)

            self.num_element[i-1] = tmp.shape[0] 

            for p in self.parameter:
                ip = np.where(pname == p)[0]
                
                if i == 1:
                    self.chain_original[p] = tmp[:,ip]
                    self.chain[p] = self.chain_original[p][self.num_burnin:]
                else:
                    self.chain_original[p] = np.append(self.chain_original[p], tmp[:,ip])
                    self.chain[p] = self.chain_original[p][self.num_burnin:]

        self.chain_original['step'] = np.cumsum(self.chain_original['weight'])
        self.chain['step'] = np.cumsum(self.chain['weight'])

    def read_chain_single_column(self):
        
        for param in self.parameter:
            print "reading "+param
            tmp = np.loadtxt(self.chain_path+'/'+param)
            self.chain[param] = tmp

        self.reproduce_chain_step()

    
    def create_2d_posterior(self, x_param, y_param, xra=None, yra=None, nbins_raw=50, frate=5.0, gaussfilter_sigma=0):
        
        xmin = np.amin(self.chain[x_param])
        xmax = np.amax(self.chain[x_param])
        ymin = np.amin(self.chain[y_param])
        ymax = np.amax(self.chain[y_param])

        dx = (xmax - xmin) / nbins_raw
        dy = (ymax - ymin) / nbins_raw

        dx_fine = dx / frate
        dy_fine = dy / frate

        nbins_fine = int(nbins_raw * frate)
        nbins_fine = int(nbins_raw * frate)

        # the half pixel below prevents the offset bias of surf_raw

        xcoord_raw = (self.chain[x_param] - xmin + 0.5*dx) / dx
        ycoord_raw = (self.chain[y_param] - ymin + 0.5*dy) / dy

        
        coord_x, coord_y = np.mgrid[xmin:xmax:(nbins_raw+1)*1j,  ymin:ymax:(nbins_raw+1)*1j]
        grid_x,  grid_y  = np.mgrid[xmin:xmax:(nbins_fine+1)*1j, ymin:ymax:(nbins_fine+1)*1j]
        
        surf_raw = np.zeros((nbins_raw+1, nbins_raw+1))

        for p in [x_param, y_param]:
            if not (p in self.parameter):
                self.parameter.append(p)
                self.num_sample[p] = np.shape(self.chain[p])[0]

        num_sample = min(self.num_sample[x_param], self.num_sample[y_param], self.num_sample['weight'])
        for i in np.arange(num_sample):
            surf_raw[int(xcoord_raw[i]), int(ycoord_raw[i])] += self.chain['weight'][i]

        points = np.zeros(((nbins_raw+1)*(nbins_raw+1),3))
        points[:,0] = coord_x.reshape((nbins_raw+1)*(nbins_raw+1))
        points[:,1] = coord_y.reshape((nbins_raw+1)*(nbins_raw+1))
        points[:,2] = np.reshape(surf_raw, (nbins_raw+1)*(nbins_raw+1))

        posterior2d = gaussian_filter( griddata(points[:,0:2], points[:,2], (grid_x, grid_y), method='cubic'), gaussfilter_sigma)
        posterior = posterior2d.reshape((nbins_fine+1)*(nbins_fine+1))

        nm = np.sum(posterior2d)
        posterior2d = posterior2d / nm
        posterior   = posterior / nm

        nm = np.sum(surf_raw)
        surf_raw = surf_raw / nm

        psort = np.sort(posterior)[::-1]
        pcum  = np.cumsum(psort)

        clevels = []
        sigma = [2.0, 1.0]
        ptes = erf(sigma / np.sqrt(2.0))
        for p in ptes:
            ind = np.where(pcum > p)[0][0]
            clevels.append(psort[ind])

        d = {'grid_x':grid_x, 'grid_y':grid_y, 'posterior2d':posterior2d, 'levels':clevels, 
             'xmin':xmin, 'xmax':xmax, 'ymin':ymin, 'ymax':ymax, 
             'dx_raw':dx, 'dy_raw':dy, 'posterior2d_raw':surf_raw,
             'dx_fine':dx_fine, 'dy_fine':dy_fine}

        return d

    def create_2d_gmm(self, x_param, y_param, num_select=None, do_weight=False):

        xmin = np.amin(self.chain[x_param])
        xmax = np.amax(self.chain[x_param])
        ymin = np.amin(self.chain[y_param])
        ymax = np.amax(self.chain[y_param])
        
        num_sample = self.num_sample['weight']

        if (num_select is None):
            num_select_gmm = num_sample

            x = self.chain[x_param]
            y = self.chain[y_param]
            w = np.int64(self.chain['weight'])
        else:
            num_select_gmm = num_select
            n_interval = num_sample / num_select_gmm

            x = self.chain[x_param][::n_interval][0:num_select_gmm]
            y = self.chain[y_param][::n_interval][0:num_select_gmm]
            w = np.int64(self.chain['weight'][::n_interval][0:num_select_gmm])

        if (do_weight):
            num_steps = np.sum(w)

            X_train = np.zeros((num_steps,2))

            ip = 0
            for i in np.arange(num_select_gmm):
                istart = ip
                iend = ip+w[i]

                X_train[istart:iend,0] = x[i]
                X_train[istart:iend,1] = y[i]
                ip = iend + 1
        else:
            X_train = np.zeros((num_select_gmm,2))
            X_train[:,0] = x
            X_train[:,1] = y

        g = mixture.GMM(n_components=4, covariance_type='full', n_iter=200, min_covar=1.0e-10)
        g.fit(X_train)

        grid_x = np.linspace(xmin, xmax, num=400)
        grid_y = np.linspace(ymin, ymax, num=400)

        X, Y = np.meshgrid(grid_x,grid_y)
        XX = np.array([X.ravel(), Y.ravel()]).T

        Z = -g.score_samples(XX)[0]
        Z = Z.reshape(X.shape)

        return X, Y, Z


    def create_2d_scatter(self, x_param, y_param, z_param, num_select=1000):

        num_sample = self.num_sample['weight']

        n_interval = num_sample / num_select

        x = self.chain[x_param][::n_interval]
        y = self.chain[y_param][::n_interval]
        z = self.chain[z_param][::n_interval]

        return x, y, z


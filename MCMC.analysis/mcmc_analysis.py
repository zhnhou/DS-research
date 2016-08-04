import numpy as np
import pandas as pd
import os
from scipy.interpolate import *
from scipy.special import *


class mcmc_analysis(object):
    def __init__(self, chain_csv_file=None, chain_path=None, feature=None):
        
        if (chain_csv_file is None and chain_path is None):
            print "Either 'chain_csv_file' or 'chain_path' should be specified."
            print "Stop and exit."
            os._exit(0)

        self.chain = dict()

        if (chain_csv_file is None):
            if (os.path.isdir(chain_path)):

                self.chain_path = chain_path

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

    def read_chain_single_column(self):
        
        for param in self.parameter:
            print "reading "+param
            tmp = np.loadtxt(self.chain_path+'/'+param)
            self.chain[param] = tmp

    
    def create_2d_posterior(self, x_param, y_param, xra=None, yra=None, nbins_raw=50, frate=5.0):
        
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

        points = np.zeros(((nbins_raw+1)*(nbins_raw+1),2))
        points[:,0] = coord_x.reshape((nbins_raw+1)*(nbins_raw+1))
        points[:,1] = coord_y.reshape((nbins_raw+1)*(nbins_raw+1))
        
        surf_raw = np.zeros((nbins_raw+1, nbins_raw+1))

        for p in [x_param, y_param]:
            if not (p in self.parameter):
                self.parameter.append(p)
                self.num_sample[p] = np.shape(self.chain[p])[0]

        num_sample = min(self.num_sample[x_param], self.num_sample[y_param], self.num_sample['weight'])
        for i in np.arange(num_sample):
            surf_raw[int(xcoord_raw[i]), int(ycoord_raw[i])] += self.chain['weight'][i]

        posterior2d = griddata(points, np.reshape(surf_raw, (nbins_raw+1)*(nbins_raw+1)), (grid_x, grid_y), method='cubic')
        posterior = posterior2d.reshape((nbins_fine+1)*(nbins_fine+1))

        nm = np.sum(posterior2d)
        posterior2d = posterior2d / nm
        posterior   = posterior / nm

        '''
        nm = np.sum(surf_raw)
        surf_raw = surf_raw / nm
        '''

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

            

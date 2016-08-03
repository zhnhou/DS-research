import numpy as np
import pandas as pd
import os

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
  
                print parameter
                
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
            print param
            self.num_sample[param] = np.shape(self.chain[param])[0]

    def read_chain_single_column(self):
        
        for param in self.parameter:
            print self.chain_path+'/'+param
            tmp = np.loadtxt(self.chain_path+'/'+param)
            self.chain[param] = tmp

    
    def create_2d_posterior(self, x_param, y_param, xra=None, yra=None):
        a =1



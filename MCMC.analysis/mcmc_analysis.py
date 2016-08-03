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

                if 'weight' in os.listdir(chain_path)
                    self.has_weight = True
                else:
                    self.has_weight = False

                if features is None:
                    parameter = os.listdir(chain_path)
                else:
                    if (self.has_weight):
                        parameter = np.unique(feature.append('weight')).tolist()
                
                self.parameter = parameter
                self.chain = self.read_chain_single_column()

        else:
            if (os.path.isfile(chain_csv_file)):
                df = pd.read_csv(chain_csv_file)

    def read_chain_single_column(self):
        
        for param in self.parameter:
            tmp = np.loadtxt(param)
            self.chain[param] = tmp

        


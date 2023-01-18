#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  7 23:24:32 2021

@author: tomas
"""
import numpy as np
from scipy import stats
from scipy.sparse import dok_matrix
from scipy.sparse import issparse


# def lognorm_cdf(mu, sigma):
#     shape  = sigma
#     loc    = 0
#     scale  = np.exp(mu)
#     return stats.lognorm( shape, loc, scale)

# mu     = 1.#2.0785
# sigma  = 2#np.log(2)#.744
# dist      = lognorm_cdf( mu, sigma)
# dist=stats.lognorm([stddev],loc=mean)
# frozen_lognorm = stats.lognorm(s=sigma, scale=np.exp(mu))
# mu, sigma = 1., 1. # mean and standard deviation

# s.mean()
# np.exp(mu + sigma/2)
# dist.stats(moments='mvsk')


from time import time

from _meta_util import timer_func

def is_pos_def(x):
    return np.all(np.linalg.eigvals(x) > -10**-7)

class NormalIG:
    def __init__(self,prior):
        self.name="Independent Gaussian"
        self.mu0=prior["mu"]
        self.kappa0=prior["kappa"]
        self.alpha0=prior["alpha"]
        self.beta0=prior["beta"]
        self.mu = np.copy(self.mu0)
        self.kappa=np.copy(self.kappa0)
        self.alpha=np.copy(self.alpha0)
        self.beta=np.copy(self.beta0)
        self.sumy=np.zeros(len(self.mu))
        self.sumy2=np.zeros(len(self.mu))
        self.N=np.zeros(len(self.mu))
    
    def update_posterior(self,sumy,sumy2,N):
        for key in sumy:
            self.mu[key]= (self.kappa0[key]*self.mu0[key] + sumy[key] ) / (
            self.kappa0[key] + N[key])
            self.kappa[key]=self.kappa0[key] + N[key]
            self.alpha[key]=self.alpha0[key] + N[key]/2
            self.beta[key]=(self.beta0[key] + .5 * (
                sumy2[key] - (
                     (sumy[key]**2 / N[key]) if N[key]>0 else 0
                         ) ) + 
                (self.kappa0[key] * self.N[key] * (self.mu0[key]-
                            ( (sumy[key]/N[key]) if N[key]>0 else 0 ))**2 )/(
                            2 * (self.kappa0[key] + N[key]) )
                )
            self.sumy[key]=sumy[key]
            self.sumy2[key]=sumy2[key]
            self.N[key]=N[key]
    def update_posterior_one_observation_stationary(self,OBS):
        self.update_posterior(
            {key : self.sumy[key]+val  for key,val in OBS.items()}, 
            {key : self.sumy2[key] + val**2  for key,val in OBS.items()}, 
            {key : self.N[key]+1 for key,val in OBS.items()})
    def sample_lambda_posterior(self):
        return 1 / np.random.gamma(shape=self.alpha, scale=self.beta) #(stats.invgamma(self.alpha, loc=0, scale=self.beta)).rvs()
    def sample_poserior(self):
        lambda_ = self.sample_lambda_posterior()#(stats.invgamma(self.alpha, loc=0, scale=self.beta)).rvs()
        # print(f"min,max,sum(lambda_)={min(lambda_)},{max(lambda_)},{sum(lambda_)}")
        return (
            stats.norm(self.mu,
                       (lambda_/self.kappa)**.5
                       )).rvs() + .5*lambda_ 

class naive_approach(NormalIG):
    def __init__(self,prior):
        super().__init__(prior)
        self.name="Naive"
        self.epsilon=prior["epsilon"]#.1
    def update_posterior(self,sumy,N):
        for key in sumy:
            self.mu[key]= sumy[key] / N[key]
            self.sumy[key]=sumy[key]
            self.N[key]=N[key]
    def update_posterior_one_observation_stationary(self,OBS):
        self.update_posterior(
            {key : self.sumy[key]+val  for key,val in OBS.items()}, 
            {key : self.N[key]+1 for key,val in OBS.items()})
    def sample_poserior(self):
        mu_ = self.mu.copy()
        # if np.random.random()<=self.epsilon:
        #     mu_[np.argmax(self.N)] = - abs(min(mu_)) * 6
        return mu_
# @timer_func
# def get_submatrix(matrix,indexes1,indexes2):
#     return matrix[indexes1,:][:,indexes2]
assert_required = True

class Spatial_Metric(NormalIG):
    def __init__(self,prior):
        super().__init__(prior)
        self.name="Spatial_"+prior["name"]
        # self.theta=prior["theta"]
        self.phi=prior["phi"]#dok_matrix((len(self.mu), len(self.mu)), dtype=np.float32)
        self.sparse_correlation_matrix_flag = issparse(self.phi)
        # for (key,val) in prior["rho"].items():
        #     self.phi[key[0],key[1]]= np.exp(-prior["theta"] * val**2 )
        if assert_required: assert(is_pos_def( #np.matmul
                          (prior["phi"]*np.eye(len(self.mu))) ))
        self.observed_paths=[]
        self.phi_Pinv=[]
        self.influence_zones=[]
        self.reversed__=prior["reversed"]
    def update_posterior_one_observation_stationary(self,OBS):
        self.update_posterior(
            {key : self.sumy[key]+OBS[key]  for key in OBS}, 
            {key : self.sumy2[key] + OBS[key]**2  for key in OBS}, 
            {key : self.N[key]+1 for key in OBS})
        self.observed_paths.append(OBS)
    # @timer_func
    def computesqrtLxphixsqrtL_P(self,sqrt_Lambda_,path_arc_indexes,inf_zone):
        # sqrt_Lambda_inf_zone = get_submatrix(sqrt_Lambda_,inf_zone,inf_zone)
        # phi_inf_zone = get_submatrix(self.phi,inf_zone,inf_zone)
        # sqrt_Lambda_inf_zone_path_arc_indexes = get_submatrix(sqrt_Lambda_,inf_zone,path_arc_indexes)
        # return np.matmul(sqrt_Lambda_inf_zone ,( phi_inf_zone * 
        #                   sqrt_Lambda_inf_zone_path_arc_indexes ) )
        if self.sparse_correlation_matrix_flag:
            return np.matmul(sqrt_Lambda_[inf_zone,:][:,inf_zone] ,
                          ( self.phi[inf_zone,:][:,path_arc_indexes] * 
                          sqrt_Lambda_[path_arc_indexes,:][:,path_arc_indexes] ) ) 
        else:
            return np.matmul(sqrt_Lambda_[inf_zone,:][:,inf_zone] ,
                          np.matmul( self.phi[inf_zone,:][:,path_arc_indexes] , 
                          sqrt_Lambda_[path_arc_indexes,:][:,path_arc_indexes] ) )
    # @timer_func
    def calculate_influence_zone(self,path_arc_indexes,AmP):
        inf_zone=[]
        for i_p in path_arc_indexes:
            for i_np in AmP:
                if self.phi[i_p,i_np]>  .1**2 and not i_np in inf_zone:
                    inf_zone.append(i_np)
        return inf_zone
    # @timer_func
    def sample_poserior(self):
        lambda_ = self.sample_lambda_posterior()#(stats.invgamma(self.alpha, loc=0, scale=self.beta)).rvs()
        mu_ = (stats.norm(self.mu,(lambda_/self.kappa)**.5 )).rvs()
        lambda_ = np.diag(lambda_)
        # d__=stats.invgamma(self.alpha, loc=0, scale=self.beta)
        # lambda_ = np.diag(d__.rvs())
        # mean, var, skew, kurt = d__.stats( moments='mvsk')
        # assert (all(np.abs(d__.stats( moments='m') - self.beta/(self.alpha-1) )<=.001) and all(np.abs(d__.stats( moments='v') - (self.beta**2)/( (self.alpha-1)**2 * (self.alpha-2) ) )<=.001))
        # mu_ = (stats.norm(self.mu,1/(self.kappa*lambda_.diagonal())**.5 )).rvs()
        number_updates=np.zeros(len(self.mu))
        
        # for _k,obs_k in enumerate(self.observed_paths):
        if self.reversed__:
            for _k,obs_k_ in enumerate(
                    self.observed_paths[len(self.observed_paths) - 1::-1]):
                # if _k<len(self.observed_paths)-5:continue
                sqrt_Lambda_ = lambda_**.5
                path_arc_indexes = list(obs_k_.keys())
                AmP = [i for i in range(len(self.mu)) if 
                       not i in path_arc_indexes]
                # self.reversed__=True
                if _k==0 and len(self.observed_paths)>len(self.influence_zones):
                    self.influence_zones.append(self.calculate_influence_zone(
                        path_arc_indexes,AmP))
                sigma_AinfxP = self.computesqrtLxphixsqrtL_P(
                    sqrt_Lambda_,path_arc_indexes,
                    self.influence_zones[len(self.observed_paths)-_k - 1])
                if _k==0 and len(self.observed_paths)>len(self.phi_Pinv):#len(self.observed_paths)-1:
                    if self.sparse_correlation_matrix_flag:
                        self.phi_Pinv.append(np.linalg.inv(
                            self.phi[:,path_arc_indexes][path_arc_indexes,:]*
                            np.eye(len(path_arc_indexes))))
                    else:
                        self.phi_Pinv.append(np.linalg.inv(np.matmul(
                            self.phi[:,path_arc_indexes][path_arc_indexes,:],
                            np.eye(len(path_arc_indexes)))))
                inv_sqrt_Lambda_P = np.linalg.inv(
                    sqrt_Lambda_[:,path_arc_indexes][path_arc_indexes,:])
                sigma_P_inv = np.matmul(np.matmul(
                    inv_sqrt_Lambda_P,self.phi_Pinv[
                        len(self.observed_paths)-_k - 1]),inv_sqrt_Lambda_P)
                # sigma_P_inv=np.linalg.inv(np.matmul(sqrt_Lambda_[:,path_arc_indexes][path_arc_indexes,:],
                #     self.phi[:,path_arc_indexes][path_arc_indexes,:]*sqrt_Lambda_[:,path_arc_indexes][path_arc_indexes,:]))
    #            mu_[AmP]+= np.matmul(sigma_AxP[AmP,:],
    #                                 np.matmul(sigma_P_inv,
    #                                           (np.array([obs_k[i] for i in path_arc_indexes]) - mu_[path_arc_indexes] )))
                for i_,i in enumerate(
                        self.influence_zones[len(self.observed_paths)-_k - 1]):
                    if self.N[i]>0: continue
                    # if number_updates[i]>0: continue
                    # if np.matmul(sigma_AxP[i,:],
                    #                   np.matmul(sigma_P_inv,
                    #                             (np.array([obs_k[i] for i in path_arc_indexes]) - mu_[path_arc_indexes] )))>0.00001 or np.matmul(
                    #                                 np.matmul(sigma_AxP[i,:] , sigma_P_inv),np.transpose(sigma_AxP[i,:])>0.00001 ):
                    #     assert(False)
                    number_updates[i] += 1
                    mu_[i] += np.matmul(sigma_AinfxP[i_,:],
                                np.matmul(sigma_P_inv,
                                 (np.array([obs_k_[j] for j in path_arc_indexes]
                                    ) - mu_[path_arc_indexes] ))
                                )#/len(self.observed_paths)
                    cte_ = np.matmul(np.matmul(
                            sigma_AinfxP[i_,:] , sigma_P_inv),
                        np.transpose(sigma_AinfxP[i_,:]) )
                    if cte_>lambda_[i,i]:
                        print("")
                    lambda_[i,i] -= np.matmul(np.matmul(
                            sigma_AinfxP[i_,:] , sigma_P_inv),
                        np.transpose(sigma_AinfxP[i_,:]) )#/len(self.observed_paths)
                # if _k==0:
                #     print(f"Max error reduction it {_k}",
                #           max([lambda_[i,i] for i in
                #         self.influence_zones[len(self.observed_paths)-_k - 1]]))
        else:
            for _k,obs_k in enumerate(self.observed_paths):
                # if _k<len(self.observed_paths)-5:continue
                sqrt_Lambda_ = lambda_**.5
                path_arc_indexes = list(obs_k.keys())
                AmP = [i for i in range(len(self.mu)) if not i in path_arc_indexes]
                # self.reversed__=True
                if _k==len(self.observed_paths)-1 and len(
                        self.observed_paths)>len(self.influence_zones):
                    self.influence_zones.append(
                        self.calculate_influence_zone(path_arc_indexes,AmP))
                sigma_AinfxP = self.computesqrtLxphixsqrtL_P(
                    sqrt_Lambda_,path_arc_indexes,self.influence_zones[_k])
                if _k==len(self.observed_paths)-1 and len(
                        self.observed_paths)>len(self.phi_Pinv):#len(self.observed_paths)-1:
                    if self.sparse_correlation_matrix_flag:
                        self.phi_Pinv.append(np.linalg.inv(
                            self.phi[:,path_arc_indexes][path_arc_indexes,:]*
                            np.eye(len(path_arc_indexes))))
                    else:
                        self.phi_Pinv.append(np.linalg.inv(np.matmul(
                            self.phi[:,path_arc_indexes][path_arc_indexes,:],
                            np.eye(len(path_arc_indexes)))))
                inv_sqrt_Lambda_P = np.linalg.inv(
                    sqrt_Lambda_[:,path_arc_indexes][path_arc_indexes,:])
                sigma_P_inv = np.matmul(
                    np.matmul(inv_sqrt_Lambda_P,
                              self.phi_Pinv[_k]),inv_sqrt_Lambda_P)
                # sigma_P_inv=np.linalg.inv(np.matmul(sqrt_Lambda_[:,path_arc_indexes][path_arc_indexes,:],
                #     self.phi[:,path_arc_indexes][path_arc_indexes,:]*sqrt_Lambda_[:,path_arc_indexes][path_arc_indexes,:]))
    #            mu_[AmP]+= np.matmul(sigma_AxP[AmP,:],
    #                                 np.matmul(sigma_P_inv,
    #                                           (np.array([obs_k[i] for i in path_arc_indexes]) - mu_[path_arc_indexes] )))
                for i_,i in enumerate(self.influence_zones[_k]):
                    if self.N[i]>0: continue
                    # if number_updates[i]>0: continue
                    # if np.matmul(sigma_AxP[i,:],
                    #                   np.matmul(sigma_P_inv,
                    #                             (np.array([obs_k[i] for i in path_arc_indexes]) - mu_[path_arc_indexes] )))>0.00001 or np.matmul(
                    #                                 np.matmul(sigma_AxP[i,:] , sigma_P_inv),np.transpose(sigma_AxP[i,:])>0.00001 ):
                    #     assert(False)
                    number_updates[i] += 1
                    mu_[i] += np.matmul(sigma_AinfxP[i_,:],
                                np.matmul(
                                    sigma_P_inv,(
                                        np.array([obs_k[j] for j in 
                                            path_arc_indexes]) - 
                                        mu_[path_arc_indexes] ))
                                )#/len(self.observed_paths)
                    cte_ = np.matmul(np.matmul(
                            sigma_AinfxP[i_,:] , sigma_P_inv),
                        np.transpose(sigma_AinfxP[i_,:]) )
                    if cte_>lambda_[i,i]:
                        print("")
                    lambda_[i,i] -= np.matmul(np.matmul(
                            sigma_AinfxP[i_,:] , sigma_P_inv),
                        np.transpose(sigma_AinfxP[i_,:]) )#/len(self.observed_paths)
                # if _k==len(self.observed_paths)-1:
                #     print(f"Max error reduction it {_k}",
                #           max([lambda_[i,i] for i in
                #         self.influence_zones[len(self.observed_paths)-_k - 1]]))
        return mu_ + (lambda_.diagonal()/2)
        

class Spatial_Metric_2(NormalIG):#PA
    def __init__(self,prior):
        super().__init__(prior)
        self.name = "Spatial_"+prior["name"]
        # self.theta=prior["theta"]
        self.generate_cov_matrix = prior["generate_cov_matrix"]
        if not self.generate_cov_matrix and assert_required: 
            assert(is_pos_def( #np.matmul
                          (prior["phi"]*np.eye(len(self.mu))) ))
        if not self.generate_cov_matrix:
            self.phi = prior["phi"]#dok_matrix((len(self.mu), len(self.mu)), dtype=np.float32)
            self.sparse_correlation_matrix_flag = issparse(self.phi)
        else:
            self.kernel = prior["kernel"]
        self.observed_paths = []
        self.phi_Pinv = []
        self.influence_zones = []
        # self.reversed__ = prior["reversed"]
    def update_posterior_one_observation_stationary(self,OBS):
        self.update_posterior(
            {key : self.sumy[key]+OBS[key]  for key in OBS}, 
            {key : self.sumy2[key] + OBS[key]**2  for key in OBS}, 
            {key : self.N[key]+1 for key in OBS})
        self.observed_paths.append(OBS)
    # @timer_func
    def computesqrtLxphixsqrtL_P(self,sqrt_Lambda_,path_arc_indexes,inf_zone):
        if self.generate_cov_matrix:
            return np.matmul(sqrt_Lambda_[inf_zone,:][:,inf_zone] ,
                          np.matmul( self.kernel(inf_zone,path_arc_indexes) , 
                          sqrt_Lambda_[path_arc_indexes,:][:,path_arc_indexes] 
                          ) )
        elif self.sparse_correlation_matrix_flag:
            return np.matmul(sqrt_Lambda_[inf_zone,:][:,inf_zone] ,
                          ( self.phi[inf_zone,:][:,path_arc_indexes] * 
                          sqrt_Lambda_[path_arc_indexes,:][:,path_arc_indexes] 
                          ) ) 
        else:
            return np.matmul(sqrt_Lambda_[inf_zone,:][:,inf_zone] ,
                          np.matmul( self.phi[inf_zone,:][:,path_arc_indexes] , 
                          sqrt_Lambda_[path_arc_indexes,:][:,path_arc_indexes] 
                          ) )
    # @timer_func
    def sample_poserior(self):
        lambda_ = self.sample_lambda_posterior()#(stats.invgamma(self.alpha, loc=0, scale=self.beta)).rvs()
        mu_ = (stats.norm(self.mu,(lambda_/self.kappa)**.5 )).rvs()
        lambda_ = np.diag(lambda_)
        # number_updates = np.zeros(len(self.mu))
        hat_mu = np.zeros(len(self.mu))
        for i in range(len(self.sumy)):
            hat_mu[i] = (self.sumy[i]/self.N[i]) if self.N[i]>0 else 0.
        # self.sumy2[key]=sumy2[key]
        sqrt_Lambda_ = lambda_**.5
        path_arc_indexes = [i for i in range(len(self.mu)) if self.N[i]>0]
        AmP = [i for i in range(len(self.mu)) if 
                   not i in path_arc_indexes]
        sigma_AinfxP = self.computesqrtLxphixsqrtL_P(
                sqrt_Lambda_,path_arc_indexes,AmP)
        if self.generate_cov_matrix:
            phi_Pinv = np.linalg.inv(np.matmul(
                        self.kernel(path_arc_indexes,path_arc_indexes),
                        np.eye(len(path_arc_indexes))))
        elif self.sparse_correlation_matrix_flag:
            phi_Pinv = np.linalg.inv(
                        self.phi[:,path_arc_indexes][path_arc_indexes,:]*
                        np.eye(len(path_arc_indexes)))
        else:
            phi_Pinv = np.linalg.inv(np.matmul(
                        self.phi[:,path_arc_indexes][path_arc_indexes,:],
                        np.eye(len(path_arc_indexes))))
        inv_sqrt_Lambda_P = np.linalg.inv(
            sqrt_Lambda_[:,path_arc_indexes][path_arc_indexes,:])
        sigma_P_inv = np.matmul(np.matmul(
            inv_sqrt_Lambda_P,phi_Pinv),inv_sqrt_Lambda_P)
        for i_,i in enumerate(AmP):
            if self.N[i]>0:continue
            mu_[i] += np.matmul(sigma_AinfxP[i_,:],
                            np.matmul(sigma_P_inv,
                             (np.array([hat_mu[j] for j in path_arc_indexes]
                                ) - mu_[path_arc_indexes] ))
                            )#/len(self.observed_paths)
            # cte_ = np.matmul(np.matmul(
            #         sigma_AinfxP[i_,:] , sigma_P_inv),
            #     np.transpose(sigma_AinfxP[i_,:]) )
            lambda_[i,i] -= np.matmul(np.matmul(
                    sigma_AinfxP[i_,:] , sigma_P_inv),
                np.transpose(sigma_AinfxP[i_,:]) )
        if False:
            mu_tilde = {}
            lambda_tilde = {}
            for i_,i in enumerate(path_arc_indexes):
                path_arc_indexes_ = [k for k in path_arc_indexes if k!=i]
                sigma_AinfxP = self.computesqrtLxphixsqrtL_P(
                        sqrt_Lambda_,path_arc_indexes_,[i])
                if self.generate_cov_matrix:
                    phi_Pinv = np.linalg.inv(np.matmul(
                        self.kernel(path_arc_indexes_,path_arc_indexes_),
                        np.eye(len(path_arc_indexes_))))
                elif self.sparse_correlation_matrix_flag:
                    phi_Pinv = np.linalg.inv(
                        self.phi[:,path_arc_indexes_][path_arc_indexes_,:]*
                        np.eye(len(path_arc_indexes_)))
                else:
                    phi_Pinv = np.linalg.inv(np.matmul(
                        self.phi[:,path_arc_indexes_][path_arc_indexes_,:],
                        np.eye(len(path_arc_indexes_))))
                inv_sqrt_Lambda_P = np.linalg.inv(
                    sqrt_Lambda_[:,path_arc_indexes_][path_arc_indexes_,:])
                sigma_P_inv = np.matmul(np.matmul(
                    inv_sqrt_Lambda_P,phi_Pinv),inv_sqrt_Lambda_P)
                mu_tilde[i] = mu_[i] + np.matmul(sigma_AinfxP[0,:],
                                np.matmul(sigma_P_inv,
                                 (np.array([hat_mu[j] for j in path_arc_indexes_]
                                    ) - mu_[path_arc_indexes_] ))
                                )
                cte_ = np.matmul(np.matmul(
                        sigma_AinfxP[0,:] , sigma_P_inv),
                    np.transpose(sigma_AinfxP[0,:]) )
                if cte_>lambda_[i,i]:
                    print("")
                lambda_tilde[i] = lambda_[i,i] - np.matmul(np.matmul(
                        sigma_AinfxP[0,:] , sigma_P_inv),
                    np.transpose(sigma_AinfxP[0,:]) )
            for i in mu_tilde.keys():
                mu_[i] = mu_tilde[i]
                lambda_[i,i] = lambda_tilde[i]
        return mu_ + (lambda_.diagonal()/2)
        

        
        
class Spatial_Metric3(NormalIG):#PB
    def __init__(self,prior):
        super().__init__(prior)
        self.name="Spatial_"+prior["name"]
        # self.theta=prior["theta"]
        # for (key,val) in prior["rho"].items():
        #     self.phi[key[0],key[1]]= np.exp(-prior["theta"] * val**2 )
        self.generate_cov_matrix = prior["generate_cov_matrix"]
        if not self.generate_cov_matrix and assert_required: 
            assert(is_pos_def( #np.matmul
                          (prior["phi"]*np.eye(len(self.mu))) ))
        if not self.generate_cov_matrix:
            self.phi = prior["phi"]#dok_matrix((len(self.mu), len(self.mu)), dtype=np.float32)
            self.sparse_correlation_matrix_flag = issparse(self.phi)
        else:
            self.kernel = prior["kernel"]
        self.observed_paths=[]
        self.phi_Pinv=[]
        self.influence_zones=[]
        self.reversed__=prior["reversed"]
    def update_posterior_one_observation_stationary(self,OBS):
        self.update_posterior(
            {key : self.sumy[key]+OBS[key]  for key in OBS}, 
            {key : self.sumy2[key] + OBS[key]**2  for key in OBS}, 
            {key : self.N[key]+1 for key in OBS})
        self.observed_paths.append(OBS)
    # @timer_func
    def computesqrtLxphixsqrtL_P(self,sqrt_Lambda_,path_arc_indexes,inf_zone):
        # sqrt_Lambda_inf_zone = get_submatrix(sqrt_Lambda_,inf_zone,inf_zone)
        # phi_inf_zone = get_submatrix(self.phi,inf_zone,inf_zone)
        # sqrt_Lambda_inf_zone_path_arc_indexes = get_submatrix(sqrt_Lambda_,inf_zone,path_arc_indexes)
        # return np.matmul(sqrt_Lambda_inf_zone ,( phi_inf_zone * 
        #                   sqrt_Lambda_inf_zone_path_arc_indexes ) )
        if self.generate_cov_matrix:
            return np.matmul(sqrt_Lambda_[inf_zone,:][:,inf_zone] ,
                          np.matmul( self.kernel(inf_zone,path_arc_indexes) , 
                          sqrt_Lambda_[path_arc_indexes,:][:,path_arc_indexes] ) )
        elif self.sparse_correlation_matrix_flag:
            return np.matmul(sqrt_Lambda_[inf_zone,:][:,inf_zone] ,
                          ( self.phi[inf_zone,:][:,path_arc_indexes] * 
                          sqrt_Lambda_[path_arc_indexes,:][:,path_arc_indexes] ) ) 
        else:
            return np.matmul(sqrt_Lambda_[inf_zone,:][:,inf_zone] ,
                          np.matmul( self.phi[inf_zone,:][:,path_arc_indexes] , 
                          sqrt_Lambda_[path_arc_indexes,:][:,path_arc_indexes] ) )
    # @timer_func
    def calculate_influence_zone(self,path_arc_indexes,AmP):
        inf_zone=[]
        for i_p in path_arc_indexes:
            for i_np in AmP:
                if (self.generate_cov_matrix and 
                    self.kernel([i_p],[i_np])>  .1**2 and 
                    not i_np in inf_zone
                    ) or (not self.generate_cov_matrix and 
                    self.phi[i_p,i_np]>  .1**2 and 
                    not i_np in inf_zone):
                    inf_zone.append(i_np)
        return inf_zone
    
    def get_dense_submatrix(self,path_arc_indexes):
        if self.generate_cov_matrix:
            return self.kernel(path_arc_indexes,path_arc_indexes)
        elif self.sparse_correlation_matrix_flag:
            return self.phi[:,path_arc_indexes][path_arc_indexes,:
                                    ]*np.eye(len(path_arc_indexes))
        else:
            return np.matmul(
                self.phi[:,path_arc_indexes][path_arc_indexes,:],
                np.eye(len(path_arc_indexes)))
    # @timer_func
    def get_hat_sigma_2(self,path_arc_indexes,mu_,sqrt_Lambda_,AmP,
                        obs_k_):
        hat_sigma_2 = 0.
        for i in path_arc_indexes:
            path_arc_indexes_i = [i__ for i__ in path_arc_indexes 
                                  if i__!=i]
            Sigma_ixPmi = self.computesqrtLxphixsqrtL_P(
                sqrt_Lambda_,
                path_arc_indexes_i,[i])
            inv_sqrt_Lambda_P_i = np.linalg.inv(sqrt_Lambda_[
                    :,path_arc_indexes_i][path_arc_indexes_i,:])
            Sigma_P_i_inv = np.matmul(np.matmul(
                    inv_sqrt_Lambda_P_i,
                    np.linalg.inv(self.get_dense_submatrix(
                        path_arc_indexes_i))
                    ),inv_sqrt_Lambda_P_i)
            hat_y_i = mu_[i] + np.matmul(Sigma_ixPmi[0,:],
                    np.matmul(Sigma_P_i_inv,
                     (np.array([obs_k_[j] for j in path_arc_indexes_i]
                        ) - mu_[path_arc_indexes_i] )))
            hat_sigma_2 += (obs_k_[i]-hat_y_i)**2
        return hat_sigma_2/len(path_arc_indexes)
    def sample_poserior(self):
        lambda_ = self.sample_lambda_posterior()#(stats.invgamma(self.alpha, loc=0, scale=self.beta)).rvs()
        mu_ = (stats.norm(self.mu,(lambda_/self.kappa)**.5 )).rvs()
        lambda_ = np.diag(lambda_)
        # d__=stats.invgamma(self.alpha, loc=0, scale=self.beta)
        # lambda_ = np.diag(d__.rvs())
        # mean, var, skew, kurt = d__.stats( moments='mvsk')
        # assert (all(np.abs(d__.stats( moments='m') - self.beta/(self.alpha-1) )<=.001) and all(np.abs(d__.stats( moments='v') - (self.beta**2)/( (self.alpha-1)**2 * (self.alpha-2) ) )<=.001))
        # mu_ = (stats.norm(self.mu,1/(self.kappa*lambda_.diagonal())**.5 )).rvs()
        number_updates = np.zeros(len(self.mu))
        # for _k,obs_k in enumerate(self.observed_paths):
        if self.reversed__:
            for _k,obs_k_ in enumerate(
                    self.observed_paths[len(self.observed_paths) - 1::-1]):
                # if _k<len(self.observed_paths)-5:continue
                sqrt_Lambda_ = lambda_**.5
                path_arc_indexes = list(obs_k_.keys())
                AmP = [i for i in range(len(self.mu)) if 
                       not i in path_arc_indexes]
                hat_sigma_2 = self.get_hat_sigma_2(
                    path_arc_indexes,mu_,sqrt_Lambda_,AmP,obs_k_)
                # self.reversed__=True
                if _k==0 and len(self.observed_paths)>len(self.influence_zones
                                                          ):
                    self.influence_zones.append(self.calculate_influence_zone(
                        path_arc_indexes,AmP))
                sigma_AinfxP = self.computesqrtLxphixsqrtL_P(
                    sqrt_Lambda_,path_arc_indexes,
                    self.influence_zones[len(self.observed_paths)-_k - 1])
                if _k==0 and len(self.observed_paths)>len(self.phi_Pinv):#len(self.observed_paths)-1:
                    self.phi_Pinv.append(np.linalg.inv(
                        self.get_dense_submatrix(path_arc_indexes)+
                    hat_sigma_2*np.linalg.inv(
                        lambda_[:,path_arc_indexes][path_arc_indexes,:]) ))
                inv_sqrt_Lambda_P = np.linalg.inv(
                    sqrt_Lambda_[:,path_arc_indexes][path_arc_indexes,:])
                sigma_P_inv = np.matmul(np.matmul(
                    inv_sqrt_Lambda_P,self.phi_Pinv[
                        len(self.observed_paths)-_k - 1]),inv_sqrt_Lambda_P)
                for i_,i in enumerate(
                        self.influence_zones[len(self.observed_paths)-_k - 1]):
                    if self.N[i]>0: continue
                    # if number_updates[i]>0: continue
                    # if np.matmul(sigma_AxP[i,:],
                    #                   np.matmul(sigma_P_inv,
                    #                             (np.array([obs_k[i] for i in path_arc_indexes]) - mu_[path_arc_indexes] )))>0.00001 or np.matmul(
                    #                                 np.matmul(sigma_AxP[i,:] , sigma_P_inv),np.transpose(sigma_AxP[i,:])>0.00001 ):
                    #     assert(False)
                    number_updates[i] += 1
                    mu_[i] += np.matmul(sigma_AinfxP[i_,:],
                                np.matmul(sigma_P_inv,
                                 (np.array([obs_k_[j] for j in path_arc_indexes]
                                    ) - mu_[path_arc_indexes] ))
                                )#/len(self.observed_paths)
                    cte_ = np.matmul(np.matmul(
                            sigma_AinfxP[i_,:] , sigma_P_inv),
                        np.transpose(sigma_AinfxP[i_,:]) )
                    if cte_>lambda_[i,i]:
                        print("")
                    lambda_[i,i] -= np.matmul(np.matmul(
                            sigma_AinfxP[i_,:] , sigma_P_inv),
                        np.transpose(sigma_AinfxP[i_,:]) )#/len(self.observed_paths)
                # if _k==0:
                #     print(f"Max error reduction it {_k}",
                #           max([lambda_[i,i] for i in
                #         self.influence_zones[len(self.observed_paths)-_k - 1]]))
        else:
            for _k,obs_k in enumerate(self.observed_paths):
                # if _k<len(self.observed_paths)-5:continue
                sqrt_Lambda_ = lambda_**.5
                path_arc_indexes = list(obs_k.keys())
                AmP = [i for i in range(len(self.mu)) 
                       if not i in path_arc_indexes]
                hat_sigma_2 = self.get_hat_sigma_2(path_arc_indexes,
                                                   mu_,sqrt_Lambda_,AmP,obs_k)
                # self.reversed__=True
                if _k==len(self.observed_paths)-1 and len(
                        self.observed_paths)>len(self.influence_zones):
                    self.influence_zones.append(
                        self.calculate_influence_zone(path_arc_indexes,AmP))
                sigma_AinfxP = self.computesqrtLxphixsqrtL_P(
                    sqrt_Lambda_,path_arc_indexes,self.influence_zones[_k])
                if _k==len(self.observed_paths)-1 and len(
                        self.observed_paths)>len(self.phi_Pinv):#len(self.observed_paths)-1:
                    self.phi_Pinv.append(np.linalg.inv(
                        self.get_dense_submatrix(path_arc_indexes)+
                        hat_sigma_2*np.linalg.inv(
                            lambda_[:,path_arc_indexes][path_arc_indexes,:]) ))
                inv_sqrt_Lambda_P = np.linalg.inv(
                    sqrt_Lambda_[:,path_arc_indexes][path_arc_indexes,:])
                sigma_P_inv = np.matmul(
                    np.matmul(inv_sqrt_Lambda_P,
                              self.phi_Pinv[_k]),inv_sqrt_Lambda_P)
                # sigma_P_inv=np.linalg.inv(np.matmul(sqrt_Lambda_[:,path_arc_indexes][path_arc_indexes,:],
                #     self.phi[:,path_arc_indexes][path_arc_indexes,:]*sqrt_Lambda_[:,path_arc_indexes][path_arc_indexes,:]))
    #            mu_[AmP]+= np.matmul(sigma_AxP[AmP,:],
    #                                 np.matmul(sigma_P_inv,
    #                                           (np.array([obs_k[i] for i in path_arc_indexes]) - mu_[path_arc_indexes] )))
                for i_,i in enumerate(self.influence_zones[_k]):
                    if self.N[i]>0: continue
                    # if number_updates[i]>0: continue
                    # if np.matmul(sigma_AxP[i,:],
                    #                   np.matmul(sigma_P_inv,
                    #                             (np.array([obs_k[i] for i in path_arc_indexes]) - mu_[path_arc_indexes] )))>0.00001 or np.matmul(
                    #                                 np.matmul(sigma_AxP[i,:] , sigma_P_inv),np.transpose(sigma_AxP[i,:])>0.00001 ):
                    #     assert(False)
                    number_updates[i] += 1
                    mu_[i] += np.matmul(sigma_AinfxP[i_,:],
                                np.matmul(
                                    sigma_P_inv,(
                                        np.array([obs_k[j] for j in 
                                            path_arc_indexes]) - 
                                        mu_[path_arc_indexes] ))
                                )#/len(self.observed_paths)
                    cte_ = np.matmul(np.matmul(
                            sigma_AinfxP[i_,:] , sigma_P_inv),
                        np.transpose(sigma_AinfxP[i_,:]) )
                    if cte_>lambda_[i,i]:
                        print("")
                    lambda_[i,i] -= np.matmul(np.matmul(
                            sigma_AinfxP[i_,:] , sigma_P_inv),
                        np.transpose(sigma_AinfxP[i_,:]) )#/len(self.observed_paths)
                # if _k==len(self.observed_paths)-1:
                #     print(f"Max error reduction it {_k}",
                #           max([lambda_[i,i] for i in
                #         self.influence_zones[len(self.observed_paths)-_k - 1]]))
        return mu_ + (lambda_.diagonal())/2
        
 
        
        
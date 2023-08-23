import tensorflow as tf
        
class IsotropicApproximation(tf.keras.layers.Layer):
    def call(self, C):
        """
        Compute the best isotropic approximation of an elasticity tensor.

        Parameters
        ----------
        C : tensorflow.Tensor
            The elasticity tensor in Voigt notation.

        Returns
        -------
        tensorflow.Tensor
            The best isotropic approximation.

        """
        if C.shape[1::] == (3,3):
            return C  # TODO
        
        elif C.shape[1::] == (6,6):
            c11 = C[:,0,0,tf.newaxis]
            c22 = C[:,1,1,tf.newaxis]
            c33 = C[:,2,2,tf.newaxis]
            c44 = C[:,3,3,tf.newaxis]
            c55 = C[:,4,4,tf.newaxis]
            c66 = C[:,5,5,tf.newaxis]
            c12 = C[:,0,1,tf.newaxis]
            c13 = C[:,0,2,tf.newaxis]
            c23 = C[:,1,2,tf.newaxis]
            
            # isotropic parameters
            lam0 = 1/15*(c11 + c22 + c33 - 2*(c44 + c55 + c66) + 4*(c12 + c13 
                                                                    + c23))
            mu0 = 1/15*(c11 + c22 + c33 + 3*(c44 + c55 + c66) - (c12 + c13 
                                                                 + c23))
            
            # build isotropic approximation
            zeros = tf.zeros_like(lam0)
            C0 = tf.concat([lam0 + 2*mu0, lam0, lam0, zeros, zeros, zeros,
                            lam0, lam0 + 2*mu0, lam0, zeros, zeros, zeros,
                            lam0, lam0, lam0 + 2*mu0, zeros, zeros, zeros,
                            zeros, zeros, zeros, mu0, zeros, zeros,
                            zeros, zeros, zeros, zeros, mu0, zeros,
                            zeros, zeros, zeros, zeros, zeros, mu0], axis=0)
            C0 = tf.reshape(C0, (-1,6,6))
            return C0
        
        else:
            raise ValueError("input data has wrong shape")
    
class _IsotropicCholeskyMatrix(tf.keras.layers.Layer):
    def call(self, G_params):
        """
        Assemble the Cholesky decomposition of an isotropic tensor.

        Parameters
        ----------
        G_params : tensorflow.Tensor
            The entries of the Cholesky matrix. Order is (11,22,33,21) in 2d
            and (11,22,33,66,21,32) in 3d case.

        Returns
        -------
        G : tensorflow.Tensor
            The assembled Cholesky decomposition matrix.

        """
        if G_params.shape[1] == 4:
            g11 = G_params[:,0,tf.newaxis]
            g22 = G_params[:,1,tf.newaxis]
            g66 = G_params[:,2,tf.newaxis]
            g21 = G_params[:,3,tf.newaxis]
            zeros = tf.zeros_like(g11)
            
            G = tf.concat([g11, zeros, zeros,
                           g21, g22, zeros,
                           zeros, zeros, g66], axis=1)
            G = tf.reshape(G, (-1,3,3))
            return G
        
        elif G_params.shape[1] == 6:
            g11 = G_params[:,0,tf.newaxis]
            g22 = G_params[:,1,tf.newaxis]
            g33 = G_params[:,2,tf.newaxis]
            g66 = G_params[:,3,tf.newaxis]
            g21 = G_params[:,4,tf.newaxis]
            g32 = G_params[:,5,tf.newaxis]
            zeros = tf.zeros_like(g11)
            
            G = tf.concat([g11, zeros, zeros, zeros, zeros, zeros,
                           g21, g22, zeros, zeros, zeros, zeros,
                           g21, g32, g33, zeros, zeros, zeros,
                           zeros, zeros, zeros, g66, zeros, zeros,
                           zeros, zeros, zeros, zeros, g66, zeros,
                           zeros, zeros, zeros, zeros, zeros, g66], axis=1)
            G = tf.reshape(G, (-1,6,6))
            return G
    
        else:
            raise ValueError("input data has wrong shape")
    
class CholeskyToElasticity(tf.keras.layers.Layer):
    def call(self, G):
        """
        Compute the elasticity tensor from its Cholesky decomposition.

        Parameters
        ----------
        G : tensorflow.Tensor
            Cholesky decomposition of the elasticity tensor in Voigt notation.

        Returns
        -------
        C : tensorflow.Tensor
            The elasticity tensor in Voigt notation.

        """
        if G.shape[1::] != (3,3) and G.shape[1::] != (6,6):
            raise ValueError("input data has wrong shape")
        C = tf.matmul(G, G, transpose_b=True)
        return C

class IsotropicElasticity2D(tf.keras.Model):
    def __init__(self, nparams=1, nlayers=1, units=32, activation="softplus"):
        super(IsotropicElasticity2D, self).__init__()
        self.ls = [tf.keras.layers.Dense(units, activation=activation, 
                                         input_shape=(nparams,))]
        for l in range(nlayers - 1):
            self.ls += [tf.keras.layers.Dense(units, activation=activation)]
        self.ls += [tf.keras.layers.Dense(4)]
        self.activation = activation
        
    def call(self, t):
        # evaluate hidden layers
        for l in self.ls:
            t = l(t)
            
        # get non-diagonal elements
        non_diag = tf.identity(t)
        non_diag = non_diag[:,3,tf.newaxis]
        
        # apply activation function again on diagonal elements
        t = tf.keras.layers.Activation(self.activation)(t)
        diag = t[:,0:3]
        
        # construct cholesky matrix
        G = tf.concat([diag, non_diag], axis=1)
        G = _IsotropicCholeskyMatrix()(G)
        
        # compute elasticity matrix
        C = CholeskyToElasticity()(G)

        return C
    
class IsotropicElasticity3D(tf.keras.Model):
    def __init__(self, nparams=1, nlayers=1, units=32, activation="softplus"):
        super(IsotropicElasticity3D, self).__init__()
        self.ls = [tf.keras.layers.Dense(units, activation=activation, 
                                         input_shape=(nparams,))]
        for l in range(nlayers - 1):
            self.ls += [tf.keras.layers.Dense(units, activation=activation)]
        self.ls += [tf.keras.layers.Dense(6)]
        self.activation = activation
        
    def call(self, t):
        # evaluate hidden layers
        for l in self.ls:
            t = l(t)
            
        # get non-diagonal elements
        non_diag = tf.identity(t)
        non_diag = non_diag[:,4::]
        
        # apply activation function again on diagonal elements
        t = tf.keras.layers.Activation(self.activation)(t)
        diag = t[:,0:4]
        
        # construct cholesky matrix
        G = tf.concat([diag, non_diag], axis=1)
        G = _IsotropicCholeskyMatrix()(G)
        
        # compute elasticity matrix
        C = CholeskyToElasticity()(G)
        
        return C
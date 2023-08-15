import numpy as np
import tensorflow as tf
import matplotlib.pyplot as plt

C = lambda E, nu: E/((1 + nu)*(1 - 2*nu))*np.array([[1 - nu, nu, 0.0],
                                                    [nu, 1 - nu, 0.0],
                                                    [0.0, 0.0, (1 - 2*nu)/2]])


class ConstitutiveTensor2D(tf.keras.Model):
    def __init__(self, nlayers=1, units=32, activation="softplus"):
        super(ConstitutiveTensor2D, self).__init__()
        self.ls = [tf.keras.layers.Dense(units, activation=activation, 
                                         input_shape=(1,))]
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
        non_diag = tf.convert_to_tensor(non_diag[:,3])
        non_diag = tf.reshape(non_diag, (-1,1))
        
        # apply activation function again on diagonal elements
        t = tf.keras.layers.Activation(self.activation)(t)
        diag = t[:,0:3]
        
        cholesky = tf.concat([diag, non_diag], axis=1)
        return cholesky



t = np.arange(100) + 1

data = []
for i in t:
    E = 100/i
    nu = 0.3
    
    constitutive = C(E, nu)
    G = np.linalg.cholesky(constitutive)
    data.append(np.array([G[0,0], G[1,1], G[2,2], G[1,0]]))
    
data = np.array(data)
# data = tf.convert_to_tensor(data)

training_in = t[0::2]
training_out = data[0::2]

test_data = t[1::2]
test_data = tf.convert_to_tensor(test_data)
test_data = tf.reshape(test_data, (-1,1))

# make tensors out of training_params
training_in = tf.convert_to_tensor(training_in)
training_in = tf.reshape(training_in, (-1,1))
training_out = tf.convert_to_tensor(training_out)

# set model parameters
kwargs = {"nlayers": 3, "units": 16, "activation": "softplus"}
epochs = 100

# train model
model = ConstitutiveTensor2D(**kwargs)
model.compile("adam", "mse")
tf.keras.backend.set_value(model.optimizer.learning_rate, 0.002)
h = model.fit(training_in, training_out, epochs=epochs, verbose=2)

# test model
GG = model.predict(test_data)

# Data for plotting
fig, ax = plt.subplots()

x_data = training_in.numpy()
y_data = training_out.numpy()[:,0]
ax.plot(x_data, y_data)

x_data = test_data.numpy()
y_data = GG[:,0]

ax.plot(x_data, y_data, lw=0, marker=".")


# ax.set(xlabel='time (s)', ylabel='voltage (mV)',
#        title='About as simple as it gets, folks')
ax.grid()

# fig.savefig("test.png")
plt.show()
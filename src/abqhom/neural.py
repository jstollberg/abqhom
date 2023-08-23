import numpy as np
import os
import tensorflow as tf
from abqhom.models import IsotropicApproximation, IsotropicElasticity2D
from abqhom.utils import find_all_files, read_csv_file
    
def read_data_files(path):
    csv_files = find_all_files(path, extension=".csv")
    data_in = []
    data_out = []
    for file in csv_files:
        data_in.append(float(file.split("_")[-1][0:-4]))
        data_out.append(read_csv_file(file, dtype=float))
    data_in = np.array(data_in)
    data_out = np.array(data_out)
        
    return data_in, data_out

def prepare_training_data(data_in, data_out, batch_fraction):
    # split into training and test data indices
    n_items = data_in.shape[0]
    n_training = int(np.floor(n_items*batch_fraction))
    all_indices = np.arange(n_items, dtype=int)
    training_indices = np.random.choice(all_indices, size=n_training, 
                                        replace=False)
    training_indices.sort()
    test_indices = np.setdiff1d(all_indices, training_indices)

    # training data
    training_in = tf.convert_to_tensor(data_in[training_indices])
    training_in = tf.reshape(training_in, (-1,1))
    training_out = tf.convert_to_tensor(data_out[training_indices])
    
    # test data
    test_in = tf.convert_to_tensor(data_in[test_indices])
    test_in = tf.reshape(test_in, (-1,1))
    test_out = tf.convert_to_tensor(data_out[test_indices])
    
    # make isotropic approximations
    training_out = IsotropicApproximation()(training_out)
    test_out = IsotropicApproximation()(test_out)

    return training_in, training_out, test_in, test_out
    
# ---------------------------------------------------------------
workdir = os.path.join("C:/Users/jonat/Documents/", 
                       "Institute for Mechanics/abqhom/examples",
                       "lattice_foam_2d")
os.chdir(workdir)
data_path = os.path.join(workdir, "results")

# read csv files
data_in, data_out = read_data_files(data_path)
(training_in, training_out, 
  test_in, test_out) = prepare_training_data(data_in, data_out, 0.5)
    
# set model parameters
kwargs = {"nlayers": 1, "units": 32, "activation": "softplus"}
epochs = 100

# train model
model = IsotropicElasticity2D(**kwargs)
model.compile("adam", "mse")
tf.keras.backend.set_value(model.optimizer.learning_rate, 0.002)
history = model.fit(training_in, training_out, epochs=epochs, verbose=2)

# export and load again
# model.save("lattice", save_format="tf")
# model = tf.keras.models.load_model("lattice")


# test model
C = model.predict(test_in)


#%% plotting
color_map = {0: "#ecbc00", 1: "#bf3c3c", 2: "#324379"}
label_map = {0: "$C_{11}$", 1: "$C_{66}$", 2: "$C_{12}$",
              3: "$G_{21}$"}
index_map = {0: [0,0], 1: [2,2], 2: [0,1]}

import matplotlib.pyplot as plt
plt.rcParams.update({"text.usetex": True,
                      "font.size": 14
                      })

# # plot results
fig1, ax1 = plt.subplots(dpi=600)
for i in range(3):
    inds = index_map[i]
    ax1.plot(test_in, C[:,inds[0],inds[1]], linewidth=2, linestyle="-",
            color=color_map[i], label=label_map[i])
    ax1.plot(test_in, test_out[:,inds[0],inds[1]], linewidth=0, markevery=3, 
            markersize=4.5, marker="o", markerfacecolor=color_map[i], 
            color=color_map[i])

plt.legend(ncol=2, handlelength=1.2, columnspacing=0.7, loc="upper left")
plt.xlabel("$d/l$")
plt.ylabel("$C_{ij}$")
fig1.tight_layout(pad=0.2)
# plt.savefig("cholesky.pdf")

fig2, ax2 = plt.subplots(dpi=600)
ax2.semilogy(history.history["loss"], label="training loss", color="black")
plt.grid(which="both")
plt.xlabel("calibration epoch")
plt.ylabel("log$_{10}$ MSE")
fig2.tight_layout(pad=0.2)
# plt.savefig("loss.pdf")